#!/bin/bash -l
#SBATCH --job-name=metascope_b
#SBATCH --partition=main,mem
#SBATCH --array=0-99%8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
# 450G, not the 320G the CAMI III script settled on.
#
# Alignment holds one bowtie2 index at a time, 126.6 GiB for a 93.25 Gbp shard;
# the R stages add the rest, and they scale with the number of ALIGNMENT RECORDS,
# not the number of reads. CAMI III measured MaxRSS 278 GiB on a 480G request
# (job 61355906_17), and five of its nineteen tasks still hit 303-311 GiB against
# a 320G cap inside metascope_blast.
#
# Strain madness is far denser: sample_0 produced 6.52 alignments per read pair
# against plasmid_01, where the CAMI III toy samples produced 0.75 -- 8.7x. So the
# R stages here need more headroom than CAMI III did, not less.
#
# 450G is free in scheduling terms: a 504 GB node (117 of them in main) holds one
# of these jobs either way, and that pool alone is six times the batch width.
#SBATCH --mem=450G
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
# The reference database is NOT copied into this project. It stays where CAMI-III
# built it and is read from there. Nothing below writes under $refdb.
refdb=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

# One array task per CAMI II strain-madness short-read sample, 0 .. 99.
# Two strings are needed, where CAMI III conflated them:
#   idx         the array index
#   sampleName  the file stem, which is also the download name
# The third string, the @SampleID the gold standard uses, is not needed here:
# metascope_to_cami.R reads it from code/sample_map_strmg*.tsv, which is the
# single source of truth for it. Restating it here would only invite drift.
idx="${SLURM_ARRAY_TASK_ID}"
sampleName="strmgCAMI2_sample_${idx}"
readsTar="${proj}/strmg_db/${sampleName}_reads.tar.gz"

### MAKE WORKING DIR
workingDir="${proj}/output/strain_madness/metascope_blast/scratch/${sampleName}_${SLURM_ARRAY_JOB_ID}"
rm -rf "$workingDir"
mkdir -p "$workingDir"


### Run MetaScope
indexDir="${refdb}/db/bowtie_indices"
expTag=$sampleName
outDir="${proj}/output/strain_madness/metascope_blast/results"
mkdir -p "$outDir"
tmpDir="${workingDir}/"

taxDB="${refdb}/db/accessionTaxa.sql"

### Choose filter and targets to use
target=$(paste -sd, "${refdb}/db/manifests/libraries.txt")

# CAMI provides no host genome, so there is nothing to filter against. "none"
filter="none"

## No priors file exists for the CAMI II strain-madness samples; the EM step starts from uniform weights
priors_df="NONE"

# BLAST db
blast_db="${refdb}/db/blast_db/cami_refseq"


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
# CAMI II ships one interleaved anonymous_reads.fq.gz per sample, with the same
# /1 and /2 header convention as CAMI III, so this splitter is unchanged. The R driver
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
"$proj"/code/run_metascope_b.R \
${readPath1} ${readPath2} ${indexDir} ${expTag} ${outDir} ${tmpDir} ${SLURM_CPUS_PER_TASK} \
${target} ${filter} \
${taxDB} ${priors_df} ${blast_db}
