#!/bin/bash -l
#SBATCH --job-name=metascope_blast
# main and mem, NOT p_wj183_1. This is the half of the pipeline that needs
# memory rather than priority, and the lab nodes are only 244-250 GB.
#SBATCH --partition=main,mem
# Submit through code/submit_metascope_batches.sh blast, which passes --array
# per batch and overrides this.
#SBATCH --array=1-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
# 32 workers, 96G. Both now measured rather than inferred.
#
# This stage was first written at 450G because CAMI III ran the *combined*
# metascope_b driver at 32 workers against a 320G cap and five of nineteen tasks
# reached 303-311 GiB. Splitting the pipeline shows that peak belonged to
# metascope_id, not to metascope_blast: the first fifteen tasks of this stage
# peaked at 14.5-34.5 GiB and all fifteen completed. 450G was thirteen times what
# the stage uses, and it confined the job to the 117 main nodes of 504 GB.
#
# 96G is 2.8x the observed peak and fits every main node, so the stage can use
# all 509 of them. metascope_blast() still gets all 32 BiocParallel workers.
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
# Measured at 16:53-41:04 over fifteen samples, against the 3.13 h estimated from
# CAMI III result mtimes. 3 h is 4.4x the longest observed and still short enough
# to schedule around a maintenance reservation; a timeout here throws away the
# BLAST work but not the alignment.
#SBATCH --time=3:00:00
#SBATCH --mail-user=yl2800@rutgers.edu
#SBATCH --requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/output/strain_madness/logs/%x_%A_%a.out
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/output/strain_madness/logs/%x_%A_%a.err
#SBATCH --open-mode=append

# Second half of the split MetaScope pipeline.
#
# run_metascope_metag.sbatch does alignment and metascope_id at 240G so it can
# reach the lab partition, and leaves two things behind:
#   output/strain_madness/metascope/results/<sample>.metascope_id.csv
#   output/strain_madness/metascope/scratch/<sample>_<jobid>/<sample>.updated.bam
# This job reads both and runs metascope_blast + blast_reassignment, writing into
# the metascope_blast/ results directory that metascope_b_to_cami.sbatch reads.
#
# Run it only after the matching metag batch has finished.

set -euo pipefail

nproc --all

proj=/projectsp/f_wj183_1/work/Yaoan/CAMI-II
refdb=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

idx="${SLURM_ARRAY_TASK_ID}"
sampleName="strmgCAMI2_sample_${idx}"

idCsv="${proj}/output/strain_madness/metascope/results/${sampleName}.metascope_id.csv"
outDir="${proj}/output/strain_madness/metascope_blast/results"
taxDB="${refdb}/db/accessionTaxa.sql"
blast_db="${refdb}/db/blast_db/cami_refseq"

### LOCATE THE METAG SCRATCH
# The directory name embeds the metag ARRAY JOB id, which differs per batch and
# again on any resubmission, so it cannot be hardcoded the way
# rerun_blast_only.sbatch hardcodes priorJob. Take the newest directory that
# actually holds the bam.
#
# The glob is unambiguous despite the shared prefix: "..._sample_1_*" requires an
# underscore right after the 1, so it cannot match "..._sample_10_<jobid>".
#
# A truncated bam is skipped, not taken: a metag task OOM-killed while writing it
# leaves a non-empty file with no BGZF end-of-file block, and blasting that would
# silently profile part of a sample. The last 28 bytes are that block; samtools
# quickcheck tests the same thing but samtools is not loaded this early.
bgzfEof=1f8b08040000000000ff0600424302001b0003000000000000000000
bamPath=""
tmpDir=""
for d in $(ls -dt "${proj}/output/strain_madness/metascope/scratch/${sampleName}_"* 2>/dev/null); do
    bam="${d}/${sampleName}.updated.bam"
    [[ -s "$bam" ]] || continue
    [[ "$(tail -c 28 "$bam" | od -An -tx1 | tr -d ' \n')" == "$bgzfEof" ]] || continue
    bamPath="$bam"
    tmpDir="${d}/"
    break
done

[[ -n "$bamPath" ]] || {
    echo "NO .updated.bam FOR ${sampleName} under ${proj}/output/strain_madness/metascope/scratch/" >&2
    echo "Has run_metascope_metag.sbatch finished this sample? Is update_bam TRUE in run_metascope_metag.R?" >&2
    exit 1; }

for f in "$idCsv" "$bamPath" "$taxDB" "${blast_db}.nal"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

mkdir -p "$outDir"

echo "Task ID:      $SLURM_ARRAY_TASK_ID"
echo "Sample Name:  $sampleName"
echo "ID csv:       $idCsv"
echo "Updated bam:  $bamPath  ($(du -h "$bamPath" | cut -f1))"
echo "Tmp dir:      $tmpDir"
echo "Output Dir:   $outDir"

# deactivate conda if something breaks
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
${idCsv} ${bamPath} ${sampleName} ${outDir} ${tmpDir} ${SLURM_CPUS_PER_TASK} \
${taxDB} ${blast_db}

### ASSERT THE OUTPUT
# blast_reassignment writes the file metascope_b_to_cami.sbatch globs for, so a
# silent miss here would surface much later as a short profile.
for f in "${outDir}/${sampleName}.metascope_blast.csv" \
         "${outDir}/${sampleName}.metascope_blast_reassigned.csv"; do
    [[ -s "$f" ]] || { echo "MISSING OR EMPTY OUTPUT: $f" >&2; exit 1; }
    echo "$(wc -l < "$f") rows  $f"
done

echo "DONE ${sampleName}"
