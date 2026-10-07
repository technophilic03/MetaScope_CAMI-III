#!/bin/bash -l
#SBATCH --job-name=metascope_metag
# p_wj183_1 first: it is the lab partition at PriorityTier=40, above mem (20)
# and main (10), so it is the only place this job is not competing with the
# whole cluster. mem and main are listed as overflow.
#SBATCH --partition=p_wj183_1,mem,main
# Submit through code/submit_metascope_batches.sh metag, which passes --array
# per batch and overrides this.
#SBATCH --array=1-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
# 240G, not the 450G the BLAST job uses. This is what lets the job onto the
# lab nodes at all -- they are 244-250 GB, so 320G or more locks them out. It
# is also the shape CAMI III proved for this exact driver: 113 metascope_metag
# tasks ran at 240G, reaching 17 concurrent. The EM step fits; it is
# metascope_blast that needs more, and that now runs separately.
#
# At 240G a 504 GB main node also takes two of these instead of one.
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=yl2800@rutgers.edu
#SBATCH --requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/output/strain_madness/logs/%x_%A_%a.out
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/output/strain_madness/logs/%x_%A_%a.err
#SBATCH --open-mode=append

set -euo pipefail

nproc --all

proj=/projectsp/f_wj183_1/work/Yaoan/CAMI-II
# The reference database stays where CAMI-III built it; nothing here writes to it.
refdb=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

# One array task per CAMI II strain-madness short-read sample, 0 .. 99.
idx="${SLURM_ARRAY_TASK_ID}"
sampleName="strmgCAMI2_sample_${idx}"
readsTar="${proj}/strmg_db/${sampleName}_reads.tar.gz"

### MAKE WORKING DIR
# NOT deleted at the end of the job. It holds <sampleName>.updated.bam, which
# the BLAST stage reads; rebuilding it costs the whole alignment again.
workingDir="${proj}/output/strain_madness/metascope/scratch/${sampleName}_${SLURM_ARRAY_JOB_ID}"
rm -rf "$workingDir"
mkdir -p "$workingDir"


### Run MetaScope
indexDir="${refdb}/db/bowtie_indices"
expTag=$sampleName
outDir="${proj}/output/strain_madness/metascope/results"
mkdir -p "$outDir"
tmpDir="${workingDir}/"

taxDB="${refdb}/db/accessionTaxa.sql"

### Choose filter and targets to use
target=$(paste -sd, "${refdb}/db/manifests/libraries.txt")

# CAMI provides no host genome, so there is nothing to filter against. "none"
filter="none"

for f in "$readsTar" "$taxDB"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

for lib in ${target//,/ }; do
    for ext in 1 2 3 4 rev.1 rev.2; do
        [[ -e "${indexDir}/${lib}.${ext}.bt2l" ]] || {
            echo "MISSING INDEX FILE: ${lib}.${ext}.bt2l" >&2; exit 1; }
    done
done

### SPLIT READS
# CAMI II ships one interleaved anonymous_reads.fq.gz per sample, same /1 and /2
# header convention as CAMI III, so this splitter is unchanged. The R driver
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

export MODULEPATH=$MODULEPATH:/projects/community/modulefiles

# Run MetaScope
module purge
module load samtools/1.23-ez82
module load R/4.5.2

export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"
Rscript --vanilla --max-ppsize=500000 \
"$proj"/code/run_metascope_metag.R \
${readPath1} ${readPath2} ${indexDir} ${expTag} ${outDir} ${tmpDir} ${SLURM_CPUS_PER_TASK} \
${target} ${filter} \
${taxDB}
