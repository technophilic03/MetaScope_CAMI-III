#!/bin/bash -l
#SBATCH --job-name=taxonomizr_db_CAMI2
#SBATCH --partition=p_wj183_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/db/scripts/taxonomizr_db_%j.log
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/db/scripts/taxonomizr_db_%j.log
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --mail-type=END,FAIL

set -eo pipefail


stem=/home/yl2800/wejlab/work/Yaoan/CAMI-II/db/scripts
export TMPDIR=/projectsp/f_wj183_1/work/Yaoan/tmp_taxonomizr
mkdir -p "$TMPDIR"

# deactivate conda if something breaks
conda deactivate || true

export MODULEPATH=$MODULEPATH:/projects/community/modulefiles

module purge
module load R/4.5.2

export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"

echo "Running R script"

Rscript --vanilla --max-ppsize=500000 \
"$stem"/build_accession_tax_2019.R

echo "Done!"

