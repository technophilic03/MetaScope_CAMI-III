#!/bin/bash -l
#SBATCH --job-name=metascope_metag
#SBATCH --partition=p_wj183_1
#SBATCH --array=0-19
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-III/output/toy/logs/%x_%A_%a.out
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-III/output/toy/logs/%x_%A_%a.err

set -euo pipefail

nproc --all

stem=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

# One array task per CAMI III toy sample: sample_0 .. sample_19
sampleName="sample_${SLURM_ARRAY_TASK_ID}"
readsTar="${stem}/toy_db/${sampleName}_reads.tar.gz"

### MAKE WORKING DIR
workingDir="${stem}/output/toy/metascope/scratch/${sampleName}_${SLURM_ARRAY_JOB_ID}"
rm -rf "$workingDir"
mkdir -p "$workingDir"


### Run MetaScope
indexDir="/home/yl2800/wejlab/reflib/2025_reference/bowtie_indices"
expTag=$sampleName
outDir="${stem}/output/toy/metascope/results"
mkdir -p "$outDir"
tmpDir="${workingDir}/"

taxDB="${stem}/db/accessionTaxa.sql"

### Choose filter and targets to use
target="target_reference"
filter="filter_reference"

for f in "$readsTar" "$taxDB" "${indexDir}/${target}.1.bt2l" "${indexDir}/${filter}.1.bt2l"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
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

# Run MetaScope
module purge
module load samtools/1.3.1
module load R/4.5.1

export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"
Rscript --vanilla --max-ppsize=500000 \
"$stem"/code/run_metascope_metag.R \
${readPath1} ${readPath2} ${indexDir} ${expTag} ${outDir} ${tmpDir} ${SLURM_CPUS_PER_TASK} \
${target} ${filter} \
${taxDB}
