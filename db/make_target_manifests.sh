#!/bin/bash
#
# Put the whole CAMI III RefSeq database into one library, target_reference.
#
# This is the alternative to db/make_group_manifests.sh, which splits the same
# database into bacteria_01..07, archaea_01, fungi_01 and viral_01. Both keep all
# 150209 assemblies. They differ only in how the assemblies are named and grouped.
#
#   db/make_target_manifests.sh db/assembly_summary.txt refdb/genomes db/manifests
#   db/make_target_manifests.sh db/assembly_summary.txt refdb/genomes db/manifests refdb/plasmids
#
# WHY THE TWO SCRIPTS DO NOT GIVE THE SAME ANSWER
# The bowtie2 option -k 9 reports up to 9 alignments per index, not per run. Ten
# group libraries can hand the EM step up to 90 alignments for one read. One
# library hands it 9. So the two layouts feed metascope_id() different input, and
# the profiles will differ. That is the reason to build both.
#
# WHY THIS STILL SHARDS
# 652.76 Gbp does not fit one bowtie2 index on these nodes. bowtie2 holds the
# whole index in memory while aligning, and 652.76 Gbp is roughly 0.8 TB. So the
# database is cut into target_reference_01..07 of 93.25 Gbp each, which is the
# same envelope as the 124.61 Gbp index that already runs here.
#
# To get one true index, raise SHARD_GBP past 653 and give bowtie2 a node with
# about 0.8 TB of memory. With one shard the library is named target_reference,
# with no number.
#
# Every assembly is kept. No row is dropped, and no row is tested against a
# column. In particular refseq_category is never read: NCBI marks no virus as a
# reference genome, so testing it would silently discard all 14990 viral rows.

set -euo pipefail

# assembly_summary.txt columns used. This is the 38-column NCBI layout, and the
# CAMI copy carries no header line and no comment lines.
#   20 ftp_path    26 genome_size

# One shard holds at most this many Gbp. 124.61 Gbp already builds and aligns on
# these nodes, so 100 stays inside a proven envelope.
SHARD_GBP=100

LIBRARY=target_reference


### ARGUMENTS ################################################################

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: $0 <assembly_summary.txt> <genomes_dir> <out_dir> [plasmids_dir]" >&2
    echo "  Builds ${LIBRARY}_01..07 from all 150209 assemblies, 652.76 Gbp." >&2
    echo "  Give plasmids_dir to add plasmid_01 as a further library." >&2
    exit 1
fi

summary=$1
genomes_dir=$2
out_dir=$3
plasmids_dir=${4:-}

for f in "$summary" "$genomes_dir"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

mkdir -p "$out_dir"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# These two tables stay next to the manifests, so a surprising library can be
# traced back to the rows that built it. Everything else is scratch.
selected="$out_dir/selected.tsv"     # bases, path
assigned="$out_dir/assigned.tsv"     # library, bases, path


### STEP 1  LIST EVERY ASSEMBLY ##############################################
#
# The archive names each file after the basename of column 20. That was checked
# against the zip central directory: 7461 of 7461 sampled entries matched.
# Rebuilding the name from columns 1 and 16 does not match, because NCBI rewrites
# the punctuation of asm_name.

echo "Step 1: listing every assembly in the database"

awk -F'\t' -v dir="$genomes_dir" '
    {
        name = $20; sub(/\/$/, "", name); sub(/.*\//, "", name)
        print $26 "\t" dir "/" name "_genomic.fna.gz"
    }
' "$summary" > "$selected"

n_selected=$(wc -l < "$selected")
n_rows=$(wc -l < "$summary")
[[ $n_selected -eq $n_rows ]] || {
    echo "read $n_rows rows but listed $n_selected assemblies" >&2; exit 1; }
echo "  $n_selected assemblies listed"


### STEP 2  CHECK EVERY FILE IS THERE ########################################
#
# A missing FASTA would quietly shrink the library and cost recall in every later
# run. Stop here, not several hours into bowtie2-build.

echo "Step 2: checking every listed file is present"

find "$genomes_dir" -name '*_genomic.fna.gz' | sed 's#.*/##' | sort > "$work/present.txt"

awk -F'\t' '
    NR == FNR { present[$0] = 1; next }
    { file = $2; sub(/.*\//, "", file); if (!(file in present)) print file }
' "$work/present.txt" "$selected" > "$work/absent.txt"

n_absent=$(wc -l < "$work/absent.txt")
if [[ $n_absent -gt 0 ]]; then
    echo "  $n_absent listed files are missing from $genomes_dir; first three:" >&2
    head -3 "$work/absent.txt" >&2
    exit 1
fi
echo "  all present"


### STEP 3  PACK INTO SHARDS #################################################

total_bases=$(awk -F'\t' '{ s += $1 } END { print s }' "$selected")
n_shards=$(( (total_bases - 1) / (SHARD_GBP * 1000000000) + 1 ))

echo "Step 3: packing $(awk -v b="$total_bases" 'BEGIN{printf "%.2f", b/1e9}') Gbp into $n_shards shard(s) of at most ${SHARD_GBP} Gbp"

# Largest genome first into whichever shard is currently emptiest. That holds the
# shards within a few percent of each other, so the array's wall time is set by
# one shard rather than by the total.
#
# A single shard takes the bare library name, so a one-index build is called
# target_reference and not target_reference_01.
sort -k1,1nr "$selected" | awk -F'\t' -v lib="$LIBRARY" -v n="$n_shards" '
    {
        pick = 1
        for (i = 2; i <= n; i++) if (load[i] < load[pick]) pick = i
        load[pick] += $1

        name = (n == 1) ? lib : sprintf("%s_%02d", lib, pick)
        print name "\t" $1 "\t" $2
    }
' > "$assigned"


### STEP 4  WRITE THE MANIFESTS ##############################################

echo "Step 4: writing manifests"

rm -f "$out_dir"/*.manifest
awk -F'\t' -v dir="$out_dir" '{ print $3 > (dir "/" $1 ".manifest") }' "$assigned"

# Plasmids ship as one FASTA per replicon, not per assembly, so they are listed
# straight from the directory. They stay a library of their own on purpose: a
# plasmid is not a cell, and its short length would turn a handful of reads into
# a large cell count if it sat beside the chromosomes.
if [[ -n "$plasmids_dir" ]]; then
    [[ -d "$plasmids_dir" ]] || { echo "MISSING INPUT: $plasmids_dir" >&2; exit 1; }
    find "$plasmids_dir" -name '*.fna.gz' | sort > "$out_dir/plasmid_01.manifest"
fi

(cd "$out_dir" && ls -1 ./*.manifest | sed 's#^\./##; s/\.manifest$//' | sort) > "$out_dir/libraries.txt"


### SUMMARY ##################################################################

echo
printf "%-22s %10s %10s\n" library assemblies Gbp

# Plasmids are counted as replicons, not assemblies, so they are printed after
# the total rather than inside it.
awk -F'\t' '
    { n[$1]++; bases[$1] += $2 }
    END { for (lib in n) printf "%-22s %10d %10.2f\n", lib, n[lib], bases[lib] / 1e9 }
' "$assigned" | sort

awk -F'\t' '
    { n++; bases += $2 }
    END { printf "%-22s %10d %10.2f\n", "TOTAL", n, bases / 1e9 }
' "$assigned"

if [[ -e "$out_dir/plasmid_01.manifest" ]]; then
    printf "%-22s %10d %10s   (replicons, not assemblies)\n" \
        plasmid_01 "$(wc -l < "$out_dir/plasmid_01.manifest")" -
fi

echo
echo "Wrote $(wc -l < "$out_dir/libraries.txt") manifests to $out_dir"
echo "Library basenames are in $out_dir/libraries.txt"
