#!/bin/bash -l

# Profile one paired-end sample with mOTUs 3.1.0, then write the same sample
# again in CAMI/bioboxes format.
#
# Reference: https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.218
#
#   motus profile -f input/ERR479298s.1.fq.gz -r input/ERR479298s.2.fq.gz \
#       -n ERR479298s -o ERR479298s-default.motus
#   motus profile -f input/ERR479298s.1.fq.gz -r input/ERR479298s.2.fq.gz \
#       -o ERR479298s-C_precision.cami -C precision
#
# Usage:
#   run_motus.sh read1 read2 sampleName sampleId outDir threads camiMode [motusDB]
#
#   sampleName   file stem for the outputs, e.g. sample_0
#   sampleId     the @SampleID OPAL pairs on, e.g. 0
#   camiMode     precision, recall or parenthesis
#   motusDB      optional path to a mOTUs database; empty means the one that
#                ships with the installed package
#
# Outputs, both in outDir:
#   <sampleName>.motus         native relative-abundance profile
#   <sampleName>.cami.profile  the same sample in CAMI format
#
set -euo pipefail

if [[ $# -lt 7 || $# -gt 8 ]]; then
    echo "usage: $0 read1 read2 sampleName sampleId outDir threads camiMode [motusDB]" >&2
    exit 1
fi

readPath1="$1"
readPath2="$2"
sampleName="$3"
sampleId="$4"
outDir="$5"
threads="$6"
camiMode="$7"
motusDB="${8:-}"

for f in "$readPath1" "$readPath2"; do
    [[ -s "$f" ]] || { echo "MISSING OR EMPTY INPUT: $f" >&2; exit 1; }
done
case "$camiMode" in
    precision|recall|parenthesis) ;;
    *) echo "BAD CAMI MODE: '${camiMode}'; use precision, recall or parenthesis" >&2; exit 1;;
esac
if [[ -n "$motusDB" ]]; then
    [[ -d "$motusDB" ]] || { echo "MISSING DATABASE: $motusDB" >&2; exit 1; }
fi

command -v motus >/dev/null || { echo "motus is not on PATH" >&2; exit 1; }

mkdir -p "$outDir"

nativeOut="${outDir}/${sampleName}.motus"
camiOut="${outDir}/${sampleName}.cami.profile"

# The -db flag is optional, so wrap the call instead of building an array.
run_profile() {
    if [[ -n "$motusDB" ]]; then
        motus profile -db "$motusDB" "$@"
    else
        motus profile "$@"
    fi
}

now=$SECONDS

echo "PROFILING ${sampleName} AS @SampleID:${sampleId}"
run_profile -f "$readPath1" -r "$readPath2" -n "$sampleId" -t "$threads" \
    -o "$nativeOut"
[[ -s "$nativeOut" ]] || { echo "EMPTY OUTPUT: $nativeOut" >&2; exit 1; }
echo "NATIVE PROFILE DONE after $((SECONDS - now))s"

run_profile -f "$readPath1" -r "$readPath2" -n "$sampleId" -t "$threads" \
    -C "$camiMode" -o "$camiOut"

# mOTUs names the top rank "superkingdom"; the CAMI III gold standard and
# MetaScope name it "domain". OPAL groups taxa by rank name, so an unrenamed
# profile scores zero at that rank. Rename it here, so the file on disk is
# ready for OPAL and the combine step downstream is a plain cat.
awk 'BEGIN { FS = OFS = "\t" }
     /^@Ranks:superkingdom/ { sub(/^@Ranks:superkingdom/, "@Ranks:domain") }
     /^[@#]/ || NF < 2 { print; next }
     $2 == "superkingdom" { $2 = "domain" }
     { print }' "$camiOut" > "${camiOut}.tmp"
mv "${camiOut}.tmp" "$camiOut"

# Every CAMI profile must carry the @SampleID the gold standard uses, or OPAL
# pairs the sample with nothing and scores it as absent.
#
# mOTUs writes three '#' comment lines and a blank line before the header
# block, so the sample id is not on line 1. It also writes "@SampleID: 1" with
# a space, while the gold standard writes "@SampleID:0" without one. OPAL skips
# comments and trims the value, so both spellings pair. Compare the trimmed id.
[[ -s "$camiOut" ]] || { echo "EMPTY OUTPUT: $camiOut" >&2; exit 1; }
got=$(sed -n 's/^@SampleID:[[:space:]]*//p' "$camiOut" | head -1 | tr -d '[:space:]')
[[ "$got" == "$sampleId" ]] \
    || { echo "WRONG @SampleID IN $camiOut: got '${got}', want '${sampleId}'" >&2; exit 1; }
rows=$(awk 'NF && $0 !~ /^[@#]/' "$camiOut" | wc -l | tr -d ' ')
[[ "$rows" -gt 0 ]] || { echo "NO TAXON ROWS IN $camiOut" >&2; exit 1; }

echo "${camiOut}  ${camiMode}, ${rows} taxon rows"
echo "TOTAL $((SECONDS - now))s"
