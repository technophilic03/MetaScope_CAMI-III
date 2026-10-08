#!/bin/bash -l
#SBATCH --job-name=metascope_metag
# mem ONLY. Measured, not assumed:
#
#   main (tier 10) is PreemptMode=REQUEUE and our priority is ~9, so tasks there
#   were preempted after 00:00:09, 00:00:28, 00:00:38, 00:19:49, 00:36:41,
#   01:11:02 and 02:33:00. A 35 h job cannot finish on main.
#
#   p_wj183_1 (tier 40) is never preempted, but its nodes are 244-250 GB and the
#   EM step needs 237+ GiB, so it cannot hold this job. Three of its six nodes
#   also went down mid-run, each taking ~34 h of alignment with it.
#
#   mem (tier 20) is preemptible only by tier 40, which is rare, and its nodes are
#   1-2 TB. The three metascope_b jobs parked there have run 46+ h untouched.
#SBATCH --partition=mem
# Submit through code/submit_metascope_batches.sh metag, which passes --array
# per batch and overrides this.
#SBATCH --array=1-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
# 320G, raised from 240G, because 240G is demonstrably too small here.
#
# CAMI III ran this driver at 240G for 113 tasks, which is why it was the
# starting point. Strain madness is denser: metascope_id's footprint scales with
# alignment records, and sample_0 produces 6.52 alignments per read pair against
# plasmid_01 where the CAMI III toy samples produce 0.75.
#
# The measurements: task 8 was OOM-killed at MaxRSS 234.1 GiB against the 240G
# cap, while the pilot running the same EM step with a 320G cap sat at 236.8 GiB
# and kept going. The EM step therefore wants more than 240G, and 320G is the
# figure a surviving run has actually demonstrated.
#SBATCH --mem=320G
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

### ALREADY DONE?
# This job is requeued often -- 7 preemptions, 3 node failures and 2 OOM kills in
# the first batch alone -- and SLURM_ARRAY_JOB_ID is unchanged across a requeue,
# so without this guard a restart would re-run up to 35 h of finished alignment
# and, worse, the rm -rf below would delete the .updated.bam the BLAST stage
# needs. Both outputs must be present to count as done: the csv alone gives a
# MetaScope profile, but only the bam lets metascope_blast run.
idCsvDone="${proj}/output/strain_madness/metascope/results/${sampleName}.metascope_id.csv"
for d in $(ls -dt "${proj}/output/strain_madness/metascope/scratch/${sampleName}_"* 2>/dev/null); do
    if [[ -s "$idCsvDone" && -s "${d}/${sampleName}.updated.bam" ]]; then
        echo "ALREADY COMPLETE, NOTHING TO DO: ${sampleName}"
        echo "  id csv      : $idCsvDone"
        echo "  updated bam : ${d}/${sampleName}.updated.bam"
        exit 0
    fi
done

### MAKE WORKING DIR
# NOT deleted at the end of the job. It holds <sampleName>.updated.bam, which
# the BLAST stage reads; rebuilding it costs the whole alignment again.
#
# The rm -rf only runs when the guard above decided this sample really does need
# redoing, so it clears a partial attempt rather than a finished one.
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
