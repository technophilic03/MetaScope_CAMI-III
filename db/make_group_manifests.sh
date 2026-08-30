#!/bin/bash
#
# Split the CAMI III RefSeq database into one bowtie2 library per group.
#
# Run this after unpacking refseq_microbial_genomes.zip and plasmids.zip:
#
#   db/make_group_manifests.sh db/assembly_summary.txt refdb/genomes db/manifests
#   db/make_group_manifests.sh db/assembly_summary.txt refdb/genomes db/manifests refdb/plasmids
#
# It writes one <group>_<NN>.manifest per library, listing the FASTA files that
# belong in it, plus libraries.txt naming them all. db/build_bowtie_indices.sbatch
# turns each manifest into an index. Pass the names to args[8] and args[9] of
# code/run_metascope_metag.R, comma joined.
#
# EVERY ASSEMBLY IS KEPT
# All 150209 rows of assembly_summary.txt end up in a library. CAMI publishes
# this database for its challenge datasets, so the reference is theirs and this
# script does not choose a subset of it. Sharding below is only about index size.
# No row is ever dropped, and no row is ever tested against a column.
#
# In particular, refseq_category is never read. NCBI marks no virus as a
# reference genome: all 15093 viral rows carry "na" in column 5. So this filter,
# which looks like it keeps viruses, keeps none:
#     $5 == "reference genome" && $25 ~ /bacteria|archaea|fungi|protozoa|viral/
# Column 5 rejects every viral row before column 25 is ever read. Nearly half of
# the CAMI toy reads are viral, so that filter costs about half the dataset.
#
# Why one library per group: align_target_bowtie() already loops over `libs`, so
# groups cost nothing kept apart. Apart, a group can be dropped from a run
# without rebuilding the others, and a change in recall points at the group that
# caused it.
#
# Why bacteria are sharded: they are 624.00 Gbp, and no single index holds that.
# The one index known to build and align on these nodes holds 124.61 Gbp.

set -euo pipefail

# assembly_summary.txt columns used. This is the 38-column NCBI layout, and the
# CAMI copy carries no header line and no comment lines.
#   20 ftp_path    25 group    26 genome_size

# One shard holds at most this many Gbp. 124.61 Gbp already builds and aligns on
# these nodes, so 100 stays inside a proven envelope.
SHARD_GBP=100


### ARGUMENTS ################################################################

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: $0 <assembly_summary.txt> <genomes_dir> <out_dir> [plasmids_dir]" >&2
    echo "  Builds 10 libraries from all 150209 assemblies: bacteria_01..07 (624.00 Gbp)," >&2
    echo "  archaea_01 (6.29), fungi_01 (21.92), viral_01 (0.55). Total 652.76 Gbp." >&2
    echo "  Give plasmids_dir to add plasmid_01 as an eleventh library." >&2
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
selected="$out_dir/selected.tsv"     # group, bases, path
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
        print $25 "\t" $26 "\t" dir "/" name "_genomic.fna.gz"
    }
' "$summary" > "$selected"

n_selected=$(wc -l < "$selected")
n_rows=$(wc -l < "$summary")
[[ $n_selected -eq $n_rows ]] || {
    echo "read $n_rows rows but listed $n_selected assemblies" >&2; exit 1; }
echo "  $n_selected assemblies listed"


### STEP 2  CHECK EVERY FILE IS THERE ########################################
#
# A missing FASTA would quietly shrink a library and cost recall in every later
# run. Stop here, not several hours into bowtie2-build.

echo "Step 2: checking every selected file is present"

find "$genomes_dir" -name '*_genomic.fna.gz' | sed 's#.*/##' | sort > "$work/present.txt"

awk -F'\t' '
    NR == FNR { present[$0] = 1; next }
    { file = $3; sub(/.*\//, "", file); if (!(file in present)) print file }
' "$work/present.txt" "$selected" > "$work/absent.txt"

n_absent=$(wc -l < "$work/absent.txt")
if [[ $n_absent -gt 0 ]]; then
    echo "  $n_absent selected files are missing from $genomes_dir; first three:" >&2
    head -3 "$work/absent.txt" >&2
    exit 1
fi
echo "  all present"


### STEP 3  PACK EACH GROUP INTO SHARDS ######################################

echo "Step 3: packing groups into shards of at most ${SHARD_GBP} Gbp"

awk -F'\t' -v cap=$((SHARD_GBP * 1000000000)) '
    { total[$1] += $2 }
    END { for (group in total) print group "\t" int((total[group] - 1) / cap) + 1 }
' "$selected" > "$work/shard_counts.tsv"

# Largest genome first into whichever shard is currently emptiest. That holds the
# shards within a few percent of each other, so the array's wall time is set by
# one shard rather than by the total.
sort -k1,1 -k2,2nr "$selected" | awk -F'\t' '
    NR == FNR { shards[$1] = $2; next }
    {
        group = $1; bases = $2; path = $3

        pick = 1
        for (i = 2; i <= shards[group]; i++)
            if (load[group, i] < load[group, pick]) pick = i
        load[group, pick] += bases

        printf "%s_%02d\t%d\t%s\n", group, pick, bases, path
    }
' "$work/shard_counts.tsv" - > "$assigned"


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
printf "%-14s %10s %10s\n" library assemblies Gbp

# Plasmids are counted as replicons, not assemblies, so they are printed after
# the total rather than inside it.
awk -F'\t' '
    { n[$1]++; bases[$1] += $2 }
    END { for (lib in n) printf "%-14s %10d %10.2f\n", lib, n[lib], bases[lib] / 1e9 }
' "$assigned" | sort

awk -F'\t' '
    { n++; bases += $2 }
    END { printf "%-14s %10d %10.2f\n", "TOTAL", n, bases / 1e9 }
' "$assigned"

if [[ -e "$out_dir/plasmid_01.manifest" ]]; then
    printf "%-14s %10d %10s   (replicons, not assemblies)\n" \
        plasmid_01 "$(wc -l < "$out_dir/plasmid_01.manifest")" -
fi

echo
echo "Wrote $(wc -l < "$out_dir/libraries.txt") manifests to $out_dir"
echo "Library basenames are in $out_dir/libraries.txt"
