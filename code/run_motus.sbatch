#!/bin/bash -l
#SBATCH --job-name=motus_metag
#SBATCH --partition=p_wj183_1
#SBATCH --array=0-19
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-user=yl2800@rutgers.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-III/output/toy/logs/%x_%A_%a.out
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-III/output/toy/logs/%x_%A_%a.err

# Run mOTUs on the CAMI III human-gut toy samples, one array task per sample.
# mOTUs is the reference profiler MetaScope is measured against here.
#
# The mOTUs calls follow Ruscheweyh et al. 2021, Current Protocols 1:e218
# (doi:10.1002/cpz1.218). Install, database and merge are Basic Protocol 1;
# the CAMI output flag is Support Protocol 3.
#
# Prerequisite, built once on a compute node, not the login node:
#   conda create -y -n motus -c bioconda -c conda-forge motus=3.1.0
#   conda activate motus
#   motus downloadDB                     # ~3.5 GB into the env's package dir
##
# This job is version 3.1.0 only.

set -euo pipefail

nproc --all

stem=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

motusEnv="motus"
motusDB=""

# What -C does with an mOTU that maps to several NCBI taxids.
#
#   precision    the ambiguous mOTU is deleted
#   recall       its relative abundance is split across the candidate taxids
#
camiModes="precision,recall"

# One array task per CAMI III toy sample: sample_0 .. sample_19.
# The gold standard names its samples by bare index, so that index is the
# @SampleID and the padded name is only a file stem.
sampleId="${SLURM_ARRAY_TASK_ID}"
sampleName="sample_${sampleId}"
readsTar="${stem}/toy_db/${sampleName}_reads.tar.gz"

[[ -e "$readsTar" ]] || { echo "MISSING INPUT: $readsTar" >&2; exit 1; }

### MAKE WORKING DIR
workingDir="${stem}/output/toy/motus/scratch/${sampleName}_${SLURM_ARRAY_JOB_ID}"
rm -rf "$workingDir"
mkdir -p "$workingDir"

outDir="${stem}/output/toy/motus/results"
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

### Run mOTUs
module purge

command -v conda >/dev/null || { echo "conda is not on PATH" >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$motusEnv"

command -v motus >/dev/null \
    || { echo "motus is not in conda env '${motusEnv}'; see the header of this script" >&2; exit 1; }
motus --version

bash "${stem}/code/run_motus.sh" \
    "$readPath1" "$readPath2" "$sampleName" "$sampleId" \
    "$outDir" "$SLURM_CPUS_PER_TASK" "$camiModes" "$motusDB"

# The split FASTQ files are large and fully reproducible from the tarball.
rm -rf "$workingDir"

echo
echo "Once every array task is done, merge the native profiles:"
echo "  motus merge -i \$(ls ${outDir}/sample_*.motus | paste -sd, -) \\"
echo "              -o ${stem}/output/toy/motus/motus_toy_merged.motus"
echo
echo "Then build the CAMI files and score them. A CAMI file holds one @SampleID"
echo "block per sample, so the twenty profiles are concatenated:"
echo "  cat ${outDir}/sample_*.cami_precision.profile > ${stem}/output/toy/cami/motus_toy_precision.profile"
echo "  cat ${outDir}/sample_*.cami_recall.profile    > ${stem}/output/toy/cami/motus_toy_recall.profile"
echo "  opal.py -g ${stem}/output/toy/cami/gold_standard.profile \\"
echo "          -o ${stem}/output/toy/cami/opal_motus \\"
echo "          -l \"MetaScope,mOTUs-precision,mOTUs-recall\" \\"
echo "          ${stem}/output/toy/cami/metascope_toy.profile \\"
echo "          ${stem}/output/toy/cami/motus_toy_precision.profile \\"
echo "          ${stem}/output/toy/cami/motus_toy_recall.profile"
