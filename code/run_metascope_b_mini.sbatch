#!/bin/bash -l
#SBATCH --job-name=metascope_b_mini
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
mini="${stem}/mini"

sampleName=mini_01
workingDir="${stem}/output/mini/scratch/${sampleName}_${SLURM_JOB_ID}"
rm -rf "$workingDir"
mkdir -p "$workingDir"

outDir="${stem}/output/mini/results"
mkdir -p "$outDir"

indexDir="${mini}/bowtie_indices"
target="mini_target"
filter="mini_filter"
taxDB="${mini}/accessionTaxa_mini.sql"
blast_db="${mini}/blast/mini_nt"
priors_df="NONE"

readPath1="${mini}/reads/${sampleName}_R1.fq.gz"
readPath2="${mini}/reads/${sampleName}_R2.fq.gz"

for f in "$readPath1" "$readPath2" "$taxDB" "${blast_db}.ndb" \
         "${indexDir}/${target}.1.bt2" "${indexDir}/${filter}.1.bt2"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

# conda restores the PATH it saved at activation, so it must go before the modules
conda deactivate 2>/dev/null || true

module purge
module load samtools/1.23-ez82
module load R/4.5.2
module load blast/2.17.0-ty308
module list

command -v samtools || { echo "samtools NOT on PATH" >&2; exit 1; }
command -v blastn   || { echo "blastn NOT on PATH" >&2; exit 1; }

export R_LIBS_USER=/home/yl2800/wejlab/work/Yaoan/R-4.5

Rscript --vanilla -e 'cat("MetaScope", as.character(packageVersion("MetaScope")), "\n"); cat("R sees blastn at:", Sys.which("blastn"), "\n"); cat("R sees samtools at:", Sys.which("samtools"), "\n")'

Rscript --vanilla --max-ppsize=500000 \
"$stem"/code/run_metascope_b.R \
${readPath1} ${readPath2} ${indexDir} ${sampleName} ${outDir} ${workingDir}/ ${SLURM_CPUS_PER_TASK} \
${target} ${filter} \
${taxDB} ${priors_df} ${blast_db}

echo "=== expected profile ==="
cat "${mini}/expected_profile.csv"
echo "=== ${sampleName}.metascope_id.csv ==="
cut -d, -f1-4 "${outDir}/${sampleName}.metascope_id.csv"
echo "=== ${sampleName}.metascope_blast_reassigned.csv ==="
cat "${outDir}/${sampleName}.metascope_blast_reassigned.csv"
