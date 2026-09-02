#!/bin/bash -l

# Profile one paired-end sample with mOTUs 4.
#
# The call is the one printed in the mOTUs 4 tutorial, "Profiling one sample":
#
#   motus profile -f sampleA_1.fastq -r sampleA_2.fastq -n sampleA \
#       -o sampleA.mOTUs4 -t 4
#
# Nothing else is passed. -g defaults to 3 marker genes, -l to 75 bp, -y to
# INSERT_SCALED, and -db to the database "motus downloadMGDB" installed. Those
# are the values this benchmark wants, and version 3 ran at the same ones.
#
# Usage:
#   run_motus4.sh read1 read2 sampleName sampleId outDir threads
#
#   sampleName   file stem for the output, e.g. sample_0
#   sampleId     the @SampleID OPAL pairs on, e.g. 0. mOTUs writes it as the
#                abundance column header, and motus4_to_cami.R reads it back
#                from there, so there is one source of truth.
#
# Outputs, both in outDir. One -o makes two files:
#   <sampleName>.motus4        integer counts against GTDB lineages
#   <sampleName>.motus4.relab  the same rows as relative abundances
#
# The counts file is the default; line 1 of it says value_type=counts. mOTUs
# derives the second name itself, as str(motu_file) + '.relab' in mentities.py.
# code/motus4_to_cami.R reads the .relab file, because a CAMI PERCENTAGE is a
# proportion, not a count.
#
# mOTUs 4 has no CAMI writer. Its profile parser accepts only -f -r -s -n -o
# -g -l -t -y --skip-pair-check -db, with no -C, and it reports GTDB lineages
# instead of NCBI taxids. The CAMI profile is therefore built in a second step
# by code/motus4_to_cami.R.

set -euo pipefail

if [[ $# -ne 6 ]]; then
    echo "usage: $0 read1 read2 sampleName sampleId outDir threads" >&2
    exit 1
fi

readPath1="$1"
readPath2="$2"
sampleName="$3"
sampleId="$4"
outDir="$5"
threads="$6"

for f in "$readPath1" "$readPath2"; do
    [[ -s "$f" ]] || { echo "MISSING OR EMPTY INPUT: $f" >&2; exit 1; }
done

command -v motus >/dev/null || { echo "motus is not on PATH" >&2; exit 1; }

mkdir -p "$outDir"
profileOut="${outDir}/${sampleName}.motus4"
relabOut="${profileOut}.relab"

now=$SECONDS

echo "PROFILING ${sampleName} AS SAMPLE ${sampleId}"
motus --help

motus profile -f "$readPath1" -r "$readPath2" \
    -n "$sampleId" -t "$threads" -o "$profileOut" -db /home/yl2800/wejlab/reflib/motus_db_4.1 

# Line 1 records the run parameters, line 2 is "mOTU<TAB>Taxonomy<TAB><sample>".
# The column header is the sample id the converter will read back.
for f in "$profileOut" "$relabOut"; do
    [[ -s "$f" ]] || { echo "EMPTY OUTPUT: $f" >&2; exit 1; }
    got=$(awk -F'\t' 'NR == 2 { print $3; exit }' "$f")
    [[ "$got" == "$sampleId" ]] \
        || { echo "WRONG SAMPLE COLUMN IN $f: got '${got}', want '${sampleId}'" >&2; exit 1; }
done

# Only the .relab file sums to 1. The counts file sums to the assigned inserts,
# so its total is reported but not asserted.
awk -F'\t' '
    NR > 2 { s += $3; n++ }
    END {
        if (n == 0) { print "NO TAXON ROWS" > "/dev/stderr"; exit 1 }
        d = s - 1; if (d < 0) d = -d
        if (d > 1e-4) { printf "RELATIVE ABUNDANCES SUM TO %.6f, NOT 1\n", s > "/dev/stderr"; exit 1 }
        printf "%d mOTUs, relative abundances sum to %.6f\n", n, s
    }' "$relabOut"

awk -F'\t' 'NR > 2 { s += $3 } END { printf "counts total %d\n", s }' "$profileOut"

echo "PROFILE STEP COMPLETE after $((SECONDS - now))s"
