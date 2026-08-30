#!/bin/bash -l

# Profile one paired-end sample with mOTUs and write it in CAMI format.
#
# Usage:
#   run_motus.sh read1 read2 sampleName sampleId outDir threads [motusDB]
#
#   sampleName   file stem for the outputs, e.g. sample_0
#   sampleId     the @SampleID OPAL pairs on, e.g. 0
#   motusDB      optional path to a mOTUs database; empty means the one that
#                ships with the installed package
#
# Outputs, all in outDir:
#   <sampleName>.motus                   native relative-abundance profile
#   <sampleName>.mgc.tsv                 marker gene cluster read counts
#   <sampleName>.cami_precision.profile  CAMI profile, unambiguous taxids only
#   <sampleName>.cami_recall.profile     CAMI profile, every mOTU mapped
#
# One mOTU can map to several NCBI taxids, so mOTUs offers two ways to write a
# CAMI profile. "precision" drops the ambiguous ones, "recall" keeps them all.
# The choice moves OPAL's precision and recall in opposite directions, so both
# are written here. Only the first call reads the FASTQ files; the other two
# re-profile from the saved counts and take seconds.

set -euo pipefail

if [[ $# -lt 6 || $# -gt 7 ]]; then
    echo "usage: $0 read1 read2 sampleName sampleId outDir threads [motusDB]" >&2
    exit 1
fi

readPath1="$1"
readPath2="$2"
sampleName="$3"
sampleId="$4"
outDir="$5"
threads="$6"
motusDB="${7:-}"

for f in "$readPath1" "$readPath2"; do
    [[ -s "$f" ]] || { echo "MISSING OR EMPTY INPUT: $f" >&2; exit 1; }
done
if [[ -n "$motusDB" ]]; then
    [[ -d "$motusDB" ]] || { echo "MISSING DATABASE: $motusDB" >&2; exit 1; }
fi

command -v motus >/dev/null || { echo "motus is not on PATH" >&2; exit 1; }

mkdir -p "$outDir"

nativeOut="${outDir}/${sampleName}.motus"
mgcOut="${outDir}/${sampleName}.mgc.tsv"

# The -db flag is optional, so wrap the call instead of building an array.
run_profile() {
    if [[ -n "$motusDB" ]]; then
        motus profile -db "$motusDB" "$@"
    else
        motus profile "$@"
    fi
}

# Every CAMI profile must carry the @SampleID the gold standard uses, or OPAL
# pairs the sample with nothing and scores it as absent.
#
# mOTUs writes three '#' comment lines and a blank line before the header
# block, so the sample id is not on line 1. It also writes "@SampleID: 1" with
# a space, while the gold standard writes "@SampleID:0" without one. OPAL skips
# comments and trims the value, so both spellings pair. Compare the trimmed id.
check_cami() {
    local f="$1"
    [[ -s "$f" ]] || { echo "EMPTY OUTPUT: $f" >&2; exit 1; }
    local got
    got=$(sed -n 's/^@SampleID:[[:space:]]*//p' "$f" | head -1 | tr -d '[:space:]')
    [[ "$got" == "$sampleId" ]] \
        || { echo "WRONG @SampleID IN $f: got '${got}', want '${sampleId}'" >&2; exit 1; }
    local rows
    rows=$(awk 'NF && $0 !~ /^[@#]/' "$f" | wc -l | tr -d ' ')
    [[ "$rows" -gt 0 ]] || { echo "NO TAXON ROWS IN $f" >&2; exit 1; }
    echo "$f  ${rows} taxon rows"
}

now=$SECONDS

echo "PROFILING ${sampleName} AS @SampleID:${sampleId}"
run_profile -f "$readPath1" -r "$readPath2" -n "$sampleId" -t "$threads" \
    -M "$mgcOut" -o "$nativeOut"

[[ -s "$nativeOut" ]] || { echo "EMPTY OUTPUT: $nativeOut" >&2; exit 1; }
[[ -s "$mgcOut" ]] || { echo "NO MARKER GENE COUNTS: $mgcOut" >&2; exit 1; }
echo "PROFILE STEP COMPLETE after $((SECONDS - now))s"

for mode in precision recall; do
    camiOut="${outDir}/${sampleName}.cami_${mode}.profile"
    run_profile -m "$mgcOut" -n "$sampleId" -C "$mode" -o "$camiOut"
    check_cami "$camiOut"
done

echo "CAMI STEP COMPLETE"
echo "TOTAL $((SECONDS - now))s"
