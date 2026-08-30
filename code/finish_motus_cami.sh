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
#   2. mOTUs 3.1.0 calls the top rank "superkingdom". The CAMI III gold
#      standard and metascope_to_cami.R both call it "domain". OPAL groups
#      taxa by rank name. Left alone, the mOTUs domain row would be scored
#      against nothing, and the gold standard domain rows against nothing.
#      The rename is applied only to the combined OPAL inputs. The per-sample
#      files mOTUs wrote are left exactly as they are.
#
# Usage:
#   conda activate motus
#   bash code/finish_motus_cami.sh

set -euo pipefail

stem=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

resultsDir="${stem}/output/toy/motus/results"
camiDir="${stem}/output/toy/cami"
samples=({0..19})

[[ -d "$resultsDir" ]] || { echo "MISSING RESULTS DIR: $resultsDir" >&2; exit 1; }
mkdir -p "$camiDir"

# Read the @SampleID mOTUs wrote, whitespace and comment lines ignored.
sample_id_of() {
    sed -n 's/^@SampleID:[[:space:]]*//p' "$1" | head -1 | tr -d '[:space:]'
}

### 1. Fill in any missing recall profile from the saved marker gene counts.
made=0
for id in "${samples[@]}"; do
    recallOut="${resultsDir}/sample_${id}.cami_recall.profile"
    mgc="${resultsDir}/sample_${id}.mgc.tsv"

    [[ -s "$recallOut" ]] && continue
    [[ -s "$mgc" ]] || { echo "NO MARKER GENE COUNTS FOR sample_${id}: $mgc" >&2; exit 1; }

    command -v motus >/dev/null || { echo "motus is not on PATH; run 'conda activate motus'" >&2; exit 1; }
    motus profile -m "$mgc" -n "$id" -C recall -o "$recallOut"
    made=$((made + 1))
done
echo "RECALL PROFILES WRITTEN THIS RUN: ${made}"

### 2. Check every sample, then build the combined OPAL inputs.
# Renaming happens here so the per-sample mOTUs output stays untouched.
rename_top_rank() {
    awk 'BEGIN { FS = OFS = "\t" }
         /^@Ranks:superkingdom/ { sub(/^@Ranks:superkingdom/, "@Ranks:domain"); print; next }
         /^[@#]/ || NF < 2 { print; next }
         $2 == "superkingdom" { $2 = "domain" }
         { print }'
}

for mode in precision recall; do
    outFile="${camiDir}/motus_toy_${mode}.profile"
    : > "$outFile"

    for id in "${samples[@]}"; do
        f="${resultsDir}/sample_${id}.cami_${mode}.profile"
        [[ -s "$f" ]] || { echo "MISSING OR EMPTY: $f" >&2; exit 1; }

        got=$(sample_id_of "$f")
        [[ "$got" == "$id" ]] \
            || { echo "WRONG @SampleID IN $f: got '${got}', want '${id}'" >&2; exit 1; }

        rows=$(awk 'NF && $0 !~ /^[@#]/' "$f" | wc -l | tr -d ' ')
        [[ "$rows" -gt 0 ]] || { echo "NO TAXON ROWS IN $f" >&2; exit 1; }

        rename_top_rank < "$f" >> "$outFile"
    done

    blocks=$(grep -c '^@SampleID' "$outFile")
    domains=$(awk -F'\t' 'NF > 1 && $2 == "domain"' "$outFile" | wc -l | tr -d ' ')
    echo "${outFile}  ${blocks} sample blocks, ${domains} domain rows"
    [[ "$blocks" -eq ${#samples[@]} ]] \
        || { echo "EXPECTED ${#samples[@]} BLOCKS, GOT ${blocks}" >&2; exit 1; }
done

echo
echo "Now score it:"
echo "  opal.py -g ${camiDir}/gold_standard.profile \\"
echo "          -o ${camiDir}/opal_motus \\"
echo "          ${camiDir}/motus_toy_precision.profile \\"
echo "          ${camiDir}/motus_toy_recall.profile \\"
echo "          ${camiDir}/metascope_toy.profile"
