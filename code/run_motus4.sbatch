#!/bin/bash -l
#SBATCH --job-name=motus4_metag
#SBATCH --partition=p_wj183_1,main
#SBATCH --array=0-99%12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=yl2800@rutgers.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/output/strain_madness/logs/%x_%A_%a.out
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/output/strain_madness/logs/%x_%A_%a.err

# Run mOTUs 4 on the CAMI II strain-madness short-read samples, one array task
# per sample.
#
# This is the newer database: 124,295 mOTUs built from 3.90 million prokaryotic
# genomes, against 33k in version 3. It is still prokaryotes only.
#
# It runs alongside the version 3 job in run_motus.sbatch. Keep both. Version 3
# writes CAMI profiles itself; version 4 does not, and reports GTDB lineages
# instead of NCBI taxids, so code/motus4_to_cami.R converts afterwards.
#
# Prerequisite, built once on a compute node, not the login node. The channel
# and pin come from the project's own mOTUs4-conda.yaml, which installs the
# tool from PyPI and its aligners from conda:
#   conda create -y -n mOTUs4 -c conda-forge -c bioconda \
#       python==3.12 bwa==0.7.19 vsearch pysam pip
#   conda activate mOTUs4
#   pip install motus-tool==4.1.0
#   motus downloadMGDB              # add --toy for a small test database
#
# pysam is installed through conda first on purpose, so conda pulls its C
# libraries. If "motus --version" still fails with
# "ImportError: libbz2.so.1.0: cannot open shared object file", run:
#   conda install -y -c conda-forge bzip2 libstdcxx-ng
#
# downloadMGDB installs the database where motus profile looks by default, so
# no -db flag is needed. The same is true of -g, -l and -y: their defaults are
# 3 marker genes, 75 bp and INSERT_SCALED, which is what version 3 ran at.

set -euo pipefail

nproc --all

proj=/projectsp/f_wj183_1/work/Yaoan/CAMI-II

motusEnv="mOTUs4"

# One array task per CAMI II strain-madness short-read sample, 0 .. 99.
# Unlike CAMI III, the @SampleID and the file stem are DIFFERENT strings here:
# the gold standard calls its blocks strmgCAMI2_short_read_sample_<N> while the
# reads tarballs are strmgCAMI2_sample_<N>. Both derive from the index.
#
# sampleId reaches "motus profile -n", which writes it as the abundance-column
# header, and motus4_to_cami.R reads the @SampleID back from there. That is the
# only path by which the CAMI II block name reaches the profile, because
# motus4_to_cami.R has no --sample-map flag.
idx="${SLURM_ARRAY_TASK_ID}"
sampleId="strmgCAMI2_short_read_sample_${idx}"
sampleName="strmgCAMI2_sample_${idx}"
readsTar="${proj}/strmg_db/${sampleName}_reads.tar.gz"

[[ -e "$readsTar" ]] || { echo "MISSING INPUT: $readsTar" >&2; exit 1; }

### MAKE WORKING DIR
workingDir="${proj}/output/strain_madness/motus4/scratch/${sampleName}_${SLURM_ARRAY_JOB_ID}"
rm -rf "$workingDir"
mkdir -p "$workingDir"

outDir="${proj}/output/strain_madness/motus4/results"
mkdir -p "$outDir"

### SPLIT READS
# CAMI ships one interleaved anonymous_reads.fq.gz per sample; mOTUs needs the
# two mates in separate files. Split in a single streaming pass and check the
# mates line up.
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
          if (n1 == 0) { print "NO READS EXTRACTED FROM THE TARBALL" > "/dev/stderr"; exit 1 }
          if (n1 != n2) { print "MATE COUNTS DIFFER: " n1 " vs " n2 > "/dev/stderr"; exit 1 }
          print "SPLIT " n1 " READ PAIRS" > "/dev/stderr"
      }'

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Sample Name: $sampleName"
echo "Sample ID: $sampleId"
echo "R1 Read Path: $readPath1"
echo "R2 Read Path: $readPath2"
echo "Output Dir: $outDir"

### Run mOTUs 4
module purge

command -v conda >/dev/null || { echo "conda is not on PATH" >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$motusEnv"

command -v motus >/dev/null \
    || { echo "motus is not in conda env '${motusEnv}'; see the header of this script" >&2; exit 1; }

bash "${proj}/code/run_motus4.sh" \
    "$readPath1" "$readPath2" "$sampleName" "$sampleId" \
    "$outDir" "$SLURM_CPUS_PER_TASK"

# The split FASTQ files are large and fully reproducible from the tarball.
rm -rf "$workingDir"

echo
echo "Once every array task is done, convert and score:"
echo "  sbatch ${proj}/code/motus4_to_cami.sbatch"
