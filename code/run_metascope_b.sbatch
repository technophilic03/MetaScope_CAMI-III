#!/bin/bash -l
#SBATCH --job-name=metascope_b
#SBATCH --partition=p_wj183_1
#SBATCH --array=0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=yl2800@rutgers.edu
#SBATCH --requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-III/output/toy/logs/%x_%A_%a.out
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-III/output/toy/logs/%x_%A_%a.err
#SBATCH --open-mode=append

set -euo pipefail

nproc --all

stem=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

# One array task per CAMI III toy sample: sample_0 .. sample_19
sampleName="sample_${SLURM_ARRAY_TASK_ID}"
readsTar="${stem}/toy_db/${sampleName}_reads.tar.gz"

### MAKE WORKING DIR
workingDir="${stem}/output/toy/metascope_blast/scratch/${sampleName}_${SLURM_ARRAY_JOB_ID}"
rm -rf "$workingDir"
mkdir -p "$workingDir"


### Run MetaScope
indexDir="${stem}/db/bowtie_indices"
expTag=$sampleName
outDir="${stem}/output/toy/metascope_blast/results"
mkdir -p "$outDir"
tmpDir="${workingDir}/"

taxDB="${stem}/db/accessionTaxa.sql"

### Choose filter and targets to use
target=$(paste -sd, "${stem}/db/manifests/libraries.txt")

# CAMI provides no host genome, so there is nothing to filter against. "none"
filter="none"

## No priors file exists for the CAMI toy samples; the EM step starts from uniform weights
priors_df="NONE"

# BLAST db
blast_db="/home/yl2800/wejlab/reflib/2024_blast_nt/nt"

for f in "$readsTar" "$taxDB" "${blast_db}.nal"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

for lib in ${target//,/ }; do
    for ext in 1 2 3 4 rev.1 rev.2; do
        [[ -e "${indexDir}/${lib}.${ext}.bt2l" ]] || {
            echo "MISSING INDEX FILE: ${lib}.${ext}.bt2l" >&2; exit 1; }
    done
done

### SPLIT READS
# CAMI ships one interleaved anonymous_reads.fq.gz per sample; the R driver
# needs two files. Split in a single streaming pass and check the mates line up.
readPath1="${workingDir}/${sampleName}_R1.fq.gz"
readPath2="${workingDir}/${sampleName}_R2.fq.gz"

tar xzOf "$readsTar" --wildcards '*/anonymous_reads.fq.gz' \
  | gzip -dc \
  | awk -v o1="$readPath1" -v o2="$readPath2" '
      BEGIN { c1 = "gzip -1 > " o1; c2 = "gzip -1 > " o2 }
      { mate = int((NR - 1) / 4) % 2 }
      NR % 4 == 1 {
          if (mate == 0) {
              if ($0 !~ /\/1$/) { print "NOT MATE 1 at line " NR ": " $0 > "/dev/stderr"; bad = 1; exit 1 }
              n1++
          } else {
              if ($0 !~ /\/2$/) { print "NOT MATE 2 at line " NR ": " $0 > "/dev/stderr"; bad = 1; exit 1 }
              n2++
          }
      }
      mate == 0 { print | c1; next }
      { print | c2 }
      END {
          close(c1); close(c2)
          if (bad) exit 1
          if (n1 != n2) { print "MATE COUNTS DIFFER: " n1 " vs " n2 > "/dev/stderr"; exit 1 }
          print "SPLIT " n1 " READ PAIRS" > "/dev/stderr"
      }'

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Sample Name: $sampleName"
echo "R1 Read Path: $readPath1"
echo "R2 Read Path: $readPath2"

# deactivate conda if something breaks
conda deactivate || true

# Run MetaScope
export MODULEPATH=$MODULEPATH:/projects/community/modulefiles
module purge
module load samtools/1.23-ez82
module load R/4.5.2
module load blast
module list

export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"
Rscript --vanilla --max-ppsize=500000 \
"$stem"/code/run_metascope_b.R \
${readPath1} ${readPath2} ${indexDir} ${expTag} ${outDir} ${tmpDir} ${SLURM_CPUS_PER_TASK} \
${target} ${filter} \
${taxDB} ${priors_df} ${blast_db}
