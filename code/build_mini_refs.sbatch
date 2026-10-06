#!/bin/bash -l
#SBATCH --job-name=mini_refs
#SBATCH --partition=p_wj183_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-III/output/mini/logs/%x_%j.out
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-III/output/mini/logs/%x_%j.err

set -euo pipefail

stem=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

# conda restores the PATH it saved at activation, so it must go before the modules
conda deactivate 2>/dev/null || true

module purge
module load R/4.5.2
module load blast/2.17.0-ty308
module list

command -v makeblastdb || { echo "makeblastdb NOT on PATH" >&2; exit 1; }

export R_LIBS_USER=/home/yl2800/wejlab/work/Yaoan/R-4.5

Rscript --vanilla "$stem"/code/build_mini_refs.R
