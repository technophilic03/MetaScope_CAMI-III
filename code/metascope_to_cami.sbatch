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
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/output/strain_madness/logs/%x_%j.out
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/output/strain_madness/logs/%x_%j.err

# Turn every MetaScope result table into one CAMI II taxonomic profile, ready
# for OPAL. All samples go into a single file, one @SampleID block each.
#
# This is one job, not an array: the hundred blocks belong in one file, and the
# taxonomy dump is read once instead of a hundred times.
#
# Prerequisite, built once (see db/build_genome_lengths.R for how):
#   db/genome_lengths_new.tsv    taxid <TAB> genome length in bases

set -euo pipefail

proj=/projectsp/f_wj183_1/work/Yaoan/CAMI-II
# The reference database stays where CAMI-III built it.
refdb=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

# The official-submission taxonomy, not the toy one. CAMI II has no taxonomy
# matching the CAMI III reference, and its predictions come from that
# reference, so the whole CAMI II chain uses this dump -- including
# motus4_to_cami.sbatch, which already pointed here.
taxDump="${refdb}/db/CAMI-III_reference"
genomeLengths="${refdb}/db/genome_lengths_new.tsv"
sampleMap="${proj}/code/sample_map_strmg.tsv"
inDir="${proj}/output/strain_madness/metascope_blast/results/"
outFile="${proj}/output/strain_madness/cami/metascope_strmg.profile"

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
Rscript --vanilla "${proj}/code/metascope_to_cami.R" \
    --taxdump "$taxDump" \
    --genome-lengths "$genomeLengths" \
    --sample-map "$sampleMap" \
    --out "$outFile" \
    "${inputs[@]}"

# Score against the gold standard. OPAL pairs a prediction to a gold standard by
# @SampleID, which is why the sample map exists.
echo
echo "Next, score it:"
echo "  # the gold standard is built once in stage 1; see output/strain_madness/cami/gold_standard_strmg_short.profile"
echo "  opal.py --collapse -g ${proj}/output/strain_madness/cami/gold_standard_strmg_short.profile -o ${proj}/output/strain_madness/cami/opal_metascope ${outFile}"
