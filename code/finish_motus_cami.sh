#!/bin/bash -l

# Finish the mOTUs run and build the two OPAL inputs. No reads are read, so
# nothing is realigned.
#
# The array job aborted at a header check that was wrong. mOTUs writes three
# '#' comment lines and a blank line before "@SampleID: N", and the check read
# only line 1. The profiles mOTUs had already written are correct. OPAL skips
# '#' comments and trims the value after the colon, so the header itself needs
# no repair to be read.
#
# Two things still need doing, and neither needs the aligner:
#
#   1. The recall profiles were never written, because the job died first.
#      They come from the saved marker gene counts in about a second each.
#
#   2. NO RANK RENAME, unlike the CAMI III version of this script. mOTUs 3.1.0
#      calls the top rank "superkingdom", and so does the CAMI II strain-madness
#      gold standard, so the two already agree and the profiles are concatenated
#      unmodified. The MetaScope converters still say "domain"; that is
#      reconciled by "opal.py --collapse" at scoring time, not here.
#
# Usage:
#   conda activate motus
#   bash code/finish_motus_cami.sh

set -euo pipefail

proj=/projectsp/f_wj183_1/work/Yaoan/CAMI-II

resultsDir="${proj}/output/strain_madness/motus/results"
camiDir="${proj}/output/strain_madness/cami"
samples=({0..99})

[[ -d "$resultsDir" ]] || { echo "MISSING RESULTS DIR: $resultsDir" >&2; exit 1; }
mkdir -p "$camiDir"

# Read the @SampleID mOTUs wrote, whitespace and comment lines ignored.
sample_id_of() {
    sed -n 's/^@SampleID:[[:space:]]*//p' "$1" | head -1 | tr -d '[:space:]'
}

### 1. Fill in any missing recall profile from the saved marker gene counts.
made=0
for id in "${samples[@]}"; do
    sampleName="strmgCAMI2_sample_${id}"
    camiId="strmgCAMI2_short_read_sample_${id}"
    recallOut="${resultsDir}/${sampleName}.cami_recall.profile"
    mgc="${resultsDir}/${sampleName}.mgc.tsv"

    [[ -s "$recallOut" ]] && continue
    [[ -s "$mgc" ]] || { echo "NO MARKER GENE COUNTS FOR ${sampleName}: $mgc" >&2; exit 1; }

    command -v motus >/dev/null || { echo "motus is not on PATH; run 'conda activate motus'" >&2; exit 1; }
    motus profile -m "$mgc" -n "$camiId" -C recall -o "$recallOut"
    made=$((made + 1))
done
echo "RECALL PROFILES WRITTEN THIS RUN: ${made}"

### 2. Check every sample, then build the combined OPAL inputs.
# A plain cat: mOTUs 3's "superkingdom" already matches the CAMI II gold standard.

for mode in precision recall; do
    outFile="${camiDir}/motus_strmg_${mode}.profile"
    : > "$outFile"

    for id in "${samples[@]}"; do
        f="${resultsDir}/strmgCAMI2_sample_${id}.cami_${mode}.profile"
        [[ -s "$f" ]] || { echo "MISSING OR EMPTY: $f" >&2; exit 1; }

        want="strmgCAMI2_short_read_sample_${id}"
        got=$(sample_id_of "$f")
        [[ "$got" == "$want" ]] \
            || { echo "WRONG @SampleID IN $f: got '${got}', want '${want}'" >&2; exit 1; }

        rows=$(awk 'NF && $0 !~ /^[@#]/' "$f" | wc -l | tr -d ' ')
        [[ "$rows" -gt 0 ]] || { echo "NO TAXON ROWS IN $f" >&2; exit 1; }

        cat "$f" >> "$outFile"
    done

    blocks=$(grep -c '^@SampleID' "$outFile")
    sk=$(awk -F'\t' 'NF > 1 && $2 == "superkingdom"' "$outFile" | wc -l | tr -d ' ')
    echo "${outFile}  ${blocks} sample blocks, ${sk} superkingdom rows"
    [[ "$sk" -gt 0 ]] \
        || { echo "NO superkingdom ROWS IN ${outFile}; the rank rename was not removed" >&2; exit 1; }
    [[ "$blocks" -eq ${#samples[@]} ]] \
        || { echo "EXPECTED ${#samples[@]} BLOCKS, GOT ${blocks}" >&2; exit 1; }
done

echo
echo "Now score it:"
echo "  opal.py --collapse -g ${camiDir}/gold_standard_strmg_short.profile \\"
echo "          -o ${camiDir}/opal_motus \\"
echo "          ${camiDir}/motus_strmg_precision.profile \\"
echo "          ${camiDir}/motus_strmg_recall.profile \\"
echo "          ${camiDir}/metascope_b_strmg.profile"
