#!/bin/bash -l

# Profile one paired-end sample with mOTUs 4.
#
# Usage:
#   run_motus4.sh read1 read2 sampleName sampleId outDir threads markerGenes [motusDB]
#
#   sampleName   file stem for the output, e.g. sample_0
#   sampleId     the @SampleID OPAL pairs on, e.g. 0. mOTUs writes it as the
#                abundance column header, and motus4_to_cami.R reads it back
#                from there, so there is one source of truth.
#   markerGenes  -g cutoff: marker genes needed to call a mOTU present.
#                1 is highest recall, 3 is the mOTUs default, 10 is strictest.
#   motusDB      optional path to a mOTUs marker gene database; empty means the
#                one downloaded by "motus downloadMGDB"
#
# Output, in outDir:
#   <sampleName>.motus4    relative abundances against GTDB lineages
#
# mOTUs 4 dropped the CAMI writer that version 3 had, and it reports GTDB
# lineages instead of NCBI taxids. The CAMI profile is therefore built in a
# second step by code/motus4_to_cami.R.

set -euo pipefail

if [[ $# -lt 7 || $# -gt 8 ]]; then
    echo "usage: $0 read1 read2 sampleName sampleId outDir threads markerGenes [motusDB]" >&2
    exit 1
fi

readPath1="$1"
readPath2="$2"
sampleName="$3"
sampleId="$4"
outDir="$5"
threads="$6"
markerGenes="$7"
motusDB="${8:-}"

for f in "$readPath1" "$readPath2"; do
    [[ -s "$f" ]] || { echo "MISSING OR EMPTY INPUT: $f" >&2; exit 1; }
done
if [[ -n "$motusDB" ]]; then
    [[ -d "$motusDB" ]] || { echo "MISSING DATABASE: $motusDB" >&2; exit 1; }
fi

command -v motus >/dev/null || { echo "motus is not on PATH" >&2; exit 1; }

mkdir -p "$outDir"
profileOut="${outDir}/${sampleName}.motus4"

now=$SECONDS

echo "PROFILING ${sampleName} AS SAMPLE ${sampleId}"
motus --version

if [[ -n "$motusDB" ]]; then
    motus profile -db "$motusDB" -f "$readPath1" -r "$readPath2" \
        -n "$sampleId" -t "$threads" -g "$markerGenes" -o "$profileOut"
else
    motus profile -f "$readPath1" -r "$readPath2" \
        -n "$sampleId" -t "$threads" -g "$markerGenes" -o "$profileOut"
fi

[[ -s "$profileOut" ]] || { echo "EMPTY OUTPUT: $profileOut" >&2; exit 1; }

# Line 1 records the run parameters, line 2 is "mOTU<TAB>Taxonomy<TAB><sample>".
# The column header is the sample id the converter will read back.
got=$(awk -F'\t' 'NR == 2 { print $3; exit }' "$profileOut")
[[ "$got" == "$sampleId" ]] \
    || { echo "WRONG SAMPLE COLUMN IN $profileOut: got '${got}', want '${sampleId}'" >&2; exit 1; }

# mOTUs reports relative abundances, so the column must sum to 1.
awk -F'\t' '
    NR > 2 { s += $3; n++ }
    END {
        if (n == 0) { print "NO TAXON ROWS" > "/dev/stderr"; exit 1 }
        d = s - 1; if (d < 0) d = -d
        if (d > 1e-4) { printf "ABUNDANCES SUM TO %.6f, NOT 1\n", s > "/dev/stderr"; exit 1 }
        printf "%d mOTUs, abundances sum to %.6f\n", n, s
    }' "$profileOut"

echo "PROFILE STEP COMPLETE after $((SECONDS - now))s"
