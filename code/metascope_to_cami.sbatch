#!/bin/bash -l
#SBATCH --job-name=metascope_to_cami
#SBATCH --partition=p_wj183_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --mail-user=yl2800@rutgers.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=./%x_%j.out
#SBATCH --error=./%x_%j.err

# Turn every MetaScope result table into one CAMI III taxonomic profile, ready
# for OPAL. All samples go into a single file, one @SampleID block each.
#
# This is one job, not an array: the twenty blocks belong in one file, and the
# taxonomy dump is read once instead of twenty times.
#
# Prerequisite, built once (see db/build_genome_lengths.R for how):
#   db/genome_lengths.tsv    taxid <TAB> genome length in bases

set -euo pipefail

stem=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

taxDump="${stem}/db"
genomeLengths="${stem}/db/genome_lengths_new.tsv"
sampleMap="${stem}/code/sample_map_toy.tsv"
inDir="${stem}/output/toy/metascope/results"
outFile="${stem}/output/toy/cami/metascope_toy.profile"

for f in "${taxDump}/nodes.dmp" "${taxDump}/names.dmp" "${taxDump}/merged.dmp" \
         "$genomeLengths" "$sampleMap"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

shopt -s nullglob
inputs=("${inDir}"/*.metascope_id.csv)
[[ ${#inputs[@]} -gt 0 ]] || { echo "NO INPUTS IN: $inDir" >&2; exit 1; }

echo "taxonomy dump  : $taxDump"
echo "genome lengths : $genomeLengths"
echo "sample map     : $sampleMap"
echo "input files    : ${#inputs[@]} from $inDir"
echo "output         : $outFile"

module purge
module load R/4.5.2

export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"

mkdir -p "$(dirname "$outFile")"
Rscript --vanilla "${stem}/code/metascope_to_cami.R" \
    --taxdump "$taxDump" \
    --genome-lengths "$genomeLengths" \
    --sample-map "$sampleMap" \
    --out "$outFile" \
    "${inputs[@]}"

# Score against the gold standard. OPAL pairs a prediction to a gold standard by
# @SampleID, which is why the sample map exists.
echo
echo "Next, score it:"
echo "  cat ${stem}/toy_db/taxonomic_profiles/taxonomic_profile_*.txt > ${stem}/output/toy/cami/gold_standard.profile"
echo "  opal.py -g ${stem}/output/toy/cami/gold_standard.profile -o ${stem}/output/toy/cami/opal ${outFile}"
