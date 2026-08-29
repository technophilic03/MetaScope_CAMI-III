#!/bin/bash -l
# Split the CAMI III RefSeq database into one bowtie2 library per taxonomic group.
#
# MetaScope aligns to every library named in `libs` and merges the results, so a
# group costs nothing extra as its own index. Keeping groups apart means a later
# run can add viruses without rebuilding bacteria, and a recall gain can be
# attributed to the group that produced it.
#
# Bacteria do not fit one index. bowtie2-build already handled 124.61 Gbp on
# these nodes, so a shard is capped at 100 Gbp here. Assemblies are dealt
# largest-first into whichever shard is currently emptiest, which holds the
# shards within a few percent of each other. The array's wall time is then set
# by the largest shard rather than by the total.
#
#   db/make_group_manifests.sh db/assembly_summary.txt refdb/genomes db/manifests
#
# Writes one <group>_<NN>.manifest per shard, plus libraries.txt naming them all.
# Pass those basenames to args[8] of code/run_metascope_metag.R, comma joined.

set -euo pipefail

if [[ $# -lt 3 || $# -gt 5 ]]; then
    echo "Usage: $0 <assembly_summary.txt> <genomes_dir> <out_dir> [bacteria_rule] [plasmids_dir]" >&2
    echo "  bacteria_rule: species (default) keeps the best assembly per species taxid;" >&2
    echo "                 reference keeps only refseq_category == 'reference genome'." >&2
    exit 1
fi

summary=$1
genomes_dir=$2
out_dir=$3
bacteria_rule=${4:-species}
plasmids_dir=${5:-}

# 100 Gbp per shard. The one index known to build and align on these nodes holds
# 124.61 Gbp, so this stays inside a proven envelope.
shard_gbp=100

case "$bacteria_rule" in
    species|reference) ;;
    *) echo "bacteria_rule must be 'species' or 'reference', got '$bacteria_rule'" >&2; exit 1 ;;
esac

for f in "$summary" "$genomes_dir"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

mkdir -p "$out_dir"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# The archive names each file after the basename of column 20 (ftp_path). That
# was checked against the zip central directory: 7461 of 7461 sampled entries
# matched. Rebuilding the name from columns 1 and 16 does not, because NCBI
# rewrites the punctuation of asm_name.
echo "Selecting assemblies (bacteria rule: $bacteria_rule)"
awk -F'\t' -v rule="$bacteria_rule" -v dir="$genomes_dir" '
{
    name = $20
    sub(/\/$/, "", name)
    sub(/.*\//, "", name)
    path = dir "/" name "_genomic.fna.gz"

    if ($25 != "bacteria") { print $25 "\t" $26 "\t" path; next }

    if (rule == "reference") {
        if ($5 == "reference genome") print $25 "\t" $26 "\t" path
        next
    }

    # One genome per species keeps every species taxid the database knows, which
    # is what the gold standard is scored on. Assembly level breaks the tie so a
    # closed genome wins over a pile of contigs.
    lvl = ($12 == "Complete Genome") ? 4 : ($12 == "Chromosome") ? 3 : ($12 == "Scaffold") ? 2 : 1
    if (!($7 in best_lvl) || lvl > best_lvl[$7]) {
        best_lvl[$7] = lvl; best_bp[$7] = $26; best_path[$7] = path
    }
}
END {
    for (t in best_path) print "bacteria\t" best_bp[t] "\t" best_path[t]
}' "$summary" | sort -k1,1 -k2,2nr > "$tmp/selected.tsv"

n_sel=$(wc -l < "$tmp/selected.tsv")
[[ $n_sel -gt 0 ]] || { echo "no assemblies selected from $summary" >&2; exit 1; }
echo "  $n_sel assemblies selected"

# A missing FASTA would silently shrink a library and quietly cost recall, so
# stop here instead of at bowtie2-build several hours in.
echo "Checking every selected file is present"
find "$genomes_dir" -name '*_genomic.fna.gz' | sed 's#.*/##' | sort -u > "$tmp/present.txt"
cut -f3 "$tmp/selected.tsv" | sed 's#.*/##' | sort -u > "$tmp/wanted.txt"
comm -23 "$tmp/wanted.txt" "$tmp/present.txt" > "$tmp/absent.txt"
n_absent=$(wc -l < "$tmp/absent.txt")
if [[ $n_absent -gt 0 ]]; then
    echo "$n_absent selected files are not in $genomes_dir; first three:" >&2
    head -3 "$tmp/absent.txt" >&2
    exit 1
fi
echo "  all present"

awk -F'\t' -v cap=$((shard_gbp * 1000000000)) '
    { total[$1] += $2 }
    END { for (g in total) print g "\t" int((total[g] - 1) / cap) + 1 "\t" total[g] }
' "$tmp/selected.tsv" > "$tmp/shards.tsv"

echo "Dealing assemblies into shards"
awk -F'\t' -v out="$out_dir" '
    NR == FNR { n[$1] = $2; next }
    {
        g = $1
        k = n[g]
        pick = 1
        for (i = 2; i <= k; i++) if (load[g, i] < load[g, pick]) pick = i
        load[g, pick] += $2
        count[g, pick]++
        print $3 > sprintf("%s/%s_%02d.manifest", out, g, pick)
    }
    END {
        printf "%-14s %10s %10s\n", "library", "assemblies", "Gbp"
        for (key in load) {
            split(key, p, SUBSEP)
            printf "%-14s %10d %10.2f\n", sprintf("%s_%02d", p[1], p[2]), count[key], load[key] / 1e9
            libs[sprintf("%s_%02d", p[1], p[2])] = 1
            grand_n += count[key]; grand_bp += load[key]
        }
        printf "%-14s %10d %10.2f\n", "TOTAL", grand_n, grand_bp / 1e9
    }
' "$tmp/shards.tsv" "$tmp/selected.tsv"

# Plasmids ship as one FASTA per replicon rather than per assembly, so they are
# listed straight from the directory. They are a separate library on purpose:
# a plasmid is not a cell, and its short length would turn a handful of reads
# into a large cell count if it sat in the same library as the chromosomes.
if [[ -n "$plasmids_dir" ]]; then
    [[ -d "$plasmids_dir" ]] || { echo "MISSING INPUT: $plasmids_dir" >&2; exit 1; }
    find "$plasmids_dir" -name '*.fna.gz' | sort > "$out_dir/plasmid_01.manifest"
    printf "%-14s %10d %10s\n" "plasmid_01" "$(wc -l < "$out_dir/plasmid_01.manifest")" "-"
fi

(cd "$out_dir" && ls -1 ./*.manifest | sed 's#^\./##; s/\.manifest$//' | sort) > "$out_dir/libraries.txt"
echo
echo "Wrote $(wc -l < "$out_dir/libraries.txt") manifests to $out_dir"
echo "Library basenames are in $out_dir/libraries.txt"
