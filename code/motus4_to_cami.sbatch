#!/bin/bash -l
#SBATCH --job-name=motus4_to_cami
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

# Turn every mOTUs 4 profile into one CAMI II taxonomic profile, ready for
# OPAL. All samples go into a single file, one @SampleID block each.
#
# This is one job, not an array: the hundred blocks belong in one file, the
# taxonomy dump is read once instead of a hundred times, and a GTDB lineage seen
# in several samples is resolved once.
#
# Run this after run_motus4.sbatch has finished all hundred array tasks.

set -euo pipefail

proj=/projectsp/f_wj183_1/work/Yaoan/CAMI-II
# The reference database stays where CAMI-III built it.
refdb=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

taxDump="${refdb}/db/CAMI-III_reference"
inDir="${proj}/output/strain_madness/motus4/results"
outFile="${proj}/output/strain_madness/cami/motus4_strmg.profile"

for f in "${taxDump}/nodes.dmp" "${taxDump}/names.dmp"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

shopt -s nullglob
# The .relab files, not the .motus4 files. One "motus profile -o" writes both:
# the named file holds integer counts, and "<name>.relab" holds the same rows as
# relative abundances. A CAMI PERCENTAGE is a proportion, so .relab is the input.
inputs=("${inDir}"/*.motus4.relab)
[[ ${#inputs[@]} -gt 0 ]] || { echo "NO INPUTS IN: $inDir" >&2; exit 1; }

echo "taxonomy dump : $taxDump"
echo "input files   : ${#inputs[@]} from $inDir"
echo "output        : $outFile"

module purge
module load R/4.5.2

export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"

mkdir -p "$(dirname "$outFile")"
Rscript --vanilla "${proj}/code/motus4_to_cami.R" \
    --taxdump "$taxDump" \
    --out "$outFile" \
    "${inputs[@]}"

blocks=$(grep -c '^@SampleID' "$outFile")
echo "${outFile}  ${blocks} sample blocks"
[[ "$blocks" -eq ${#inputs[@]} ]] \
    || { echo "EXPECTED ${#inputs[@]} BLOCKS, GOT ${blocks}" >&2; exit 1; }

# Score it against the cellular gold standard, the same one MetaScope and
# mOTUs 3 are scored on. mOTUs profiles prokaryotes only, and both profilers
# rescale each rank to 100, so the full gold standard would be on a different
# scale: its cellular part is only about half of each rank.
echo
echo "Now score it:"
echo "  opal.py --collapse -g ${proj}/output/strain_madness/cami/gold_standard_strmg_short.profile \\"
echo "          -o ${proj}/output/strain_madness/cami/opal_all \\"
echo "          -l \"MetaScope,mOTUs3,mOTUs4\" \\"
echo "          -r superkingdom,species \\"
echo "          ${proj}/output/strain_madness/cami/metascope_b_strmg.profile \\"
echo "          ${proj}/output/strain_madness/cami/motus_strmg_recall.profile \\"
echo "          ${outFile}"
