#!/bin/bash -l
#SBATCH --job-name=metascope_b_blastonly
#SBATCH --partition=main,mem
# The seven tasks of job 61290578 that the OOM killer took, all inside
# metascope_blast(). Their alignment and EM output survives in scratch.
#SBATCH --array=2,4,6,9,10,12,15
#SBATCH --nodes=1
#SBATCH --ntasks=1
# 32 workers put five of these tasks between 303G and 311G against a 320G cap.
# Peak tracks the worker count, so 8 cuts the fork cost about fourfold.
#SBATCH --cpus-per-task=8
#SBATCH --mem=240G
# The failed tasks reached this stage in under an hour of a 25-hour run.
#SBATCH --time=8:00:00
#SBATCH --mail-user=yl2800@rutgers.edu
#SBATCH --requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/output/strain_madness/logs/%x_%A_%a.out
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/output/strain_madness/logs/%x_%A_%a.err
#SBATCH --open-mode=append

set -euo pipefail

nproc --all

proj=/projectsp/f_wj183_1/work/Yaoan/CAMI-II
refdb=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

sampleName="strmgCAMI2_sample_${SLURM_ARRAY_TASK_ID}"

# Reuse the scratch dir job 61290578 left behind. Never delete it: it holds the
# .updated.bam this rerun reads, and rebuilding it costs 25 hours of alignment.
priorJob=61290578
workingDir="${proj}/output/strain_madness/metascope_blast/scratch/${sampleName}_${priorJob}"
tmpDir="${workingDir}/"

outDir="${proj}/output/strain_madness/metascope_blast/results"
expTag=$sampleName
taxDB="${refdb}/db/accessionTaxa.sql"
blast_db="${refdb}/db/blast_db/cami_refseq"

bamPath="${workingDir}/${sampleName}.updated.bam"
idCsv="${outDir}/${sampleName}.metascope_id.csv"

for f in "$bamPath" "$idCsv" "$taxDB" "${blast_db}.nal"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Sample Name: $sampleName"
echo "BAM: $bamPath"
echo "MetaScope ID CSV: $idCsv"

conda deactivate || true

export MODULEPATH=$MODULEPATH:/projects/community/modulefiles
module purge
module load samtools/1.23-ez82
module load R/4.5.2
module load blast
module list

export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"
Rscript --vanilla --max-ppsize=500000 \
"$proj"/code/run_metascope_blast_only.R \
${idCsv} ${bamPath} ${expTag} ${outDir} ${tmpDir} ${SLURM_CPUS_PER_TASK} \
${taxDB} ${blast_db}
