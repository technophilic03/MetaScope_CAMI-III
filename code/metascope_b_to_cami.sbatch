#!/bin/bash -l
#SBATCH --job-name=metascope_b_to_cami
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

# Turn every MetaScope BLAST-reassigned table into CAMI III taxonomic profiles,
# ready for OPAL. All samples go into one file, one @SampleID block each.
#
# This is one job, not an array: the sample blocks belong in one file, and the
# taxonomy dump is read once instead of once per sample.
#
# Two profiles come out of one run:
#   metascope_b_toy.profile             every row of every input
#   metascope_b_toy_validated.profile   the blast_validated rows only
#
# The second one is stricter, not better. metascope_blast() BLASTs only the top
# num_results taxa, so a taxon below that rank is unvalidated for want of
# evidence. Treat that profile as an abundance cutoff at that rank.
#
# Prerequisite, built once (see db/build_genome_lengths.R for how):
#   db/genome_lengths_new.tsv    taxid <TAB> genome length in bases

set -euo pipefail

stem=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

taxDump="${stem}/db"
genomeLengths="${stem}/db/genome_lengths_new.tsv"
sampleMap="${stem}/code/sample_map_toy_b.tsv"
inDir="${stem}/output/toy/metascope_blast/results"
outFile="${stem}/output/toy/cami/metascope_b_toy.profile"
outValidated="${stem}/output/toy/cami/metascope_b_toy_validated.profile"

for f in "${taxDump}/nodes.dmp" "${taxDump}/names.dmp" "${taxDump}/merged.dmp" \
         "$genomeLengths" "$sampleMap"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

shopt -s nullglob
inputs=("${inDir}"/*.metascope_blast_reassigned.csv)
[[ ${#inputs[@]} -gt 0 ]] || { echo "NO INPUTS IN: $inDir" >&2; exit 1; }

echo "taxonomy dump  : $taxDump"
echo "genome lengths : $genomeLengths"
echo "sample map     : $sampleMap"
echo "input files    : ${#inputs[@]} from $inDir"
echo "output         : $outFile"
echo "output (valid) : $outValidated"

module purge
module load R/4.5.2

export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"

mkdir -p "$(dirname "$outFile")"
Rscript --vanilla "${stem}/code/metascope_b_to_cami.R" \
    --taxdump "$taxDump" \
    --genome-lengths "$genomeLengths" \
    --sample-map "$sampleMap" \
    --out "$outFile" \
    --out-validated "$outValidated" \
    "${inputs[@]}"

# Score against the gold standard. OPAL pairs a prediction to a gold standard by
# @SampleID, which is why the sample map exists.
echo
echo "Next, score it:"
echo "  cat ${stem}/toy_db/taxonomic_profiles/taxonomic_profile_*.txt > ${stem}/output/toy/cami/gold_standard.profile"
echo "  opal.py -g ${stem}/output/toy/cami/gold_standard.profile -o ${stem}/output/toy/cami/opal_metascope_b ${outFile} ${outValidated}"
