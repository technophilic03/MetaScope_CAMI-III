#!/bin/bash
#SBATCH --job-name=bt2_cami2
#SBATCH --partition=p_wj183_1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/db/scripts/bt2_%j.log
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/db/scripts/bt2_%j.err
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yl2800@gsbs.rutgers.edu
set -eo pipefail

export MODULEPATH=$MODULEPATH:/projects/community/modulefiles
module purge
module load R/4.5.2
module load bowtie2
export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"

command -v python    # Rbowtie2 needs python; the job stops here if it is missing

cd /projectsp/f_wj183_1/work/Yaoan/CAMI-II/db
mkdir -p ref_cat
out=ref_cat/cami2_refseq.fna.gz

# 1. Join all .fna.gz files into one .gz file (skip if already done)
if [ ! -s "$out" ]; then
  find ncbi_refseq -maxdepth 1 -name '*.fna.gz' -print0 | sort -z \
    | xargs -0 cat > cami2_refseq.fna.gz.tmp
  mv cami2_refseq.fna.gz.tmp "$out"
fi

# 2. Check: the joined size must equal the sum of the input sizes
want=$(find ncbi_refseq -maxdepth 1 -name '*.fna.gz' -printf '%s\n' \
       | awk '{s += $1} END {printf "%.0f\n", s}')
have=$(stat -c %s "$out")
echo "input bytes: $want  joined bytes: $have"
[ "$want" -eq "$have" ]

# 3. Build the index
Rscript --vanilla /home/yl2800/wejlab/work/Yaoan/CAMI-II/db/scripts/build_bt2_index.R
