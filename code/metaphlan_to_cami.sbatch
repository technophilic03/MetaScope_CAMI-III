#!/bin/bash -l
#SBATCH --job-name=metaphlan_to_cami
#SBATCH --partition=main,p_wj183_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --mail-user=yl2800@rutgers.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/output/strain_madness/logs/%x_%j.out
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/output/strain_madness/logs/%x_%j.err

# Turn every MetaPhlAn profile into one CAMI II taxonomic profile, ready for
# OPAL. All samples go into a single file, one @SampleID block each.
#
# One job, not an array: the blocks belong in one file and the taxonomy dump is
# read once.
#
# Run this after code/run_metaphlan_salloc.sh has finished.

set -euo pipefail

proj=/projectsp/f_wj183_1/work/Yaoan/CAMI-II
# The reference database stays where CAMI-III built it.
refdb=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

# The official-submission taxonomy, like every other CAMI II converter.
taxDump="${refdb}/db/CAMI-III_reference"
inDir="${proj}/output/strain_madness/metaphlan/results"
outFile="${proj}/output/strain_madness/cami/metaphlan_strmg.profile"

for f in "${taxDump}/nodes.dmp" "${taxDump}/merged.dmp"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

shopt -s nullglob
inputs=("${inDir}"/*.metaphlan.tsv)
[[ ${#inputs[@]} -gt 0 ]] || { echo "NO INPUTS IN: $inDir" >&2; exit 1; }

echo "taxonomy dump : $taxDump"
echo "input files   : ${#inputs[@]} from $inDir"
echo "output        : $outFile"

module purge
module load R/4.5.2

export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"

mkdir -p "$(dirname "$outFile")"
Rscript --vanilla "${proj}/code/metaphlan_to_cami.R" \
    --taxdump "$taxDump" \
    --out "$outFile" \
    "${inputs[@]}"

blocks=$(grep -c '^@SampleID' "$outFile")
echo "${outFile}  ${blocks} sample blocks"
[[ "$blocks" -eq ${#inputs[@]} ]] \
    || { echo "EXPECTED ${#inputs[@]} BLOCKS, GOT ${blocks}" >&2; exit 1; }
