#!/bin/bash -l
#SBATCH --job-name=normalise_taxonomy
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

# Build the OPAL inputs for one benchmark: every profile followed through
# merged.dmp and cut to the samples that every profiler finished.
#
# Usage:
#   sbatch code/normalise_taxonomy.sbatch cami2
#   sbatch code/normalise_taxonomy.sbatch cami3
#
# Each benchmark is pushed into the taxonomy of the gold standard it is scored
# against, and under its own dump neither gold standard moves a single taxid:
#   cami2  ${refdb}/db/CAMI-III_reference, the official-submission dump
#          (2026-07-01). The strain-madness gold standard is current under it,
#          and it recovers 51 synonyms of gold-standard taxa against the toy
#          dump's 48.
#   cami3  ${refdb}/db, the toy dataset's own dump (2025-08-01). The official
#          dump has since merged 3424663 into 3680502, which would move 1.68%
#          of every toy gold-standard sample.
#
# Score the result with "opal.py --collapse": the MetaScope and mOTUs 4
# profiles call the top rank "domain", the gold standards and the others
# "superkingdom", and --collapse is what pairs the two.

set -euo pipefail

proj=/projectsp/f_wj183_1/work/Yaoan/CAMI-II
# The reference database stays where CAMI-III built it.
refdb=/projectsp/f_wj183_1/work/Yaoan/CAMI-III
# The CAMI III project itself, whose own profiles are normalised here too.
cami3="$refdb"

case "${1:-}" in
    cami2)
        taxDump="${refdb}/db/CAMI-III_reference"
        inDir="${proj}/output/strain_madness/cami"
        outDir="${inDir}/normalised_21"
        samples=$(printf 'strmgCAMI2_short_read_sample_%s,' \
            0 1 2 5 6 7 10 12 13 18 20 23 24 25 26 27 29 30 31 34 35)
        files=(gold_standard_strmg_short metascope_strmg metascope_b_strmg
               metascope_b_strmg_validated motus_strmg_precision
               motus_strmg_recall motus4_strmg metaphlan_strmg)
        ;;
    cami3)
        # Sample 19 never completed in either MetaScope pipeline.
        taxDump="${cami3}/db"
        inDir="${cami3}/output/toy/cami"
        outDir="${inDir}/normalised_19"
        samples=$(printf '%s,' $(seq 0 18))
        files=(gold_standard metascope_toy.fixlen metascope_toy_cami_refseq.fixlen
               metascope_b_toy.fixlen metascope_b_toy_validated.fixlen
               motus_toy_precision motus_toy_recall motus4_toy)
        ;;
    *)
        echo "usage: $0 {cami2|cami3}" >&2
        exit 1
        ;;
esac
samples=${samples%,}

inputs=()
for f in "${files[@]}"; do
    inputs+=("${inDir}/${f}.profile")
done
for f in "${taxDump}"/{nodes,names,merged}.dmp "${inputs[@]}"; do
    [[ -s "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

echo "taxonomy dump : $taxDump"
echo "samples       : $samples"
echo "output        : $outDir"

module purge
module load R/4.5.2

export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"

Rscript --vanilla "${proj}/code/normalise_taxonomy.R" \
    --taxdump "$taxDump" \
    --out-dir "$outDir" \
    --samples "$samples" \
    "${inputs[@]}"

want=$(tr ',' '\n' <<< "$samples" | wc -l)
for f in "${files[@]}"; do
    got=$(grep -c '^@SampleID' "${outDir}/${f}.profile")
    [[ "$got" -eq "$want" ]] \
        || { echo "EXPECTED ${want} BLOCKS IN ${f}.profile, GOT ${got}" >&2; exit 1; }
done
echo "${#files[@]} profiles of ${want} blocks each in ${outDir}"
