#!/bin/bash -l
#SBATCH --job-name=bt2_build
#SBATCH --partition=p_wj183_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=yl2800@rutgers.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-III/output/logs/%x_%A_%a.out
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-III/output/logs/%x_%A_%a.err

# Build one bowtie2 index per library listed in db/manifests/libraries.txt.
#
# Set the array range to the number of lines in that file. For the whole CAMI III
# database that is 10: bacteria_01..07, archaea_01, fungi_01, viral_01. Add
# plasmid_01 and any host library and it grows from there.
#   sbatch --array=1-$(wc -l < db/manifests/libraries.txt) db/build_bowtie_indices.sbatch
#
# Every index lands in one directory because align_target_bowtie() takes a single
# lib_dir and looks up each name inside it.
#
# Host libraries join the same array. Write the path of the GRCh38 FASTA into
# db/manifests/human_01.manifest and it builds like any other library.

set -euo pipefail

stem=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

manifest_dir="${stem}/db/manifests"
index_dir="${stem}/db/bowtie_indices"
work_dir="/scratch/yl2800/bt2_build"

# The lab module name for bowtie2 varies; override it rather than editing here.
bowtie2_build=${BOWTIE2_BUILD:-bowtie2-build}
bowtie2_inspect=${BOWTIE2_INSPECT:-bowtie2-inspect}

command -v "$bowtie2_build" >/dev/null || {
    echo "bowtie2-build not on PATH. Load a module or set BOWTIE2_BUILD." >&2; exit 1; }

lib=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${manifest_dir}/libraries.txt")
[[ -n "$lib" ]] || { echo "no library on line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
manifest="${manifest_dir}/${lib}.manifest"
[[ -e "$manifest" ]] || { echo "MISSING INPUT: $manifest" >&2; exit 1; }

mkdir -p "$index_dir" "$work_dir"
fasta="${work_dir}/${lib}.fna"
trap 'rm -f "$fasta"' EXIT

echo "Library: $lib"
echo "Files in manifest: $(wc -l < "$manifest")"

# bowtie2-build reads the reference twice, so it cannot take a stream, and a
# 19000-name command line would overrun ARG_MAX. Concatenating first avoids both.
# Gzip members concatenate legally, so one decompression pass handles the lot.
echo "Concatenating to $fasta"
while read -r f; do cat "$f"; done < "$manifest" | gzip -dc > "$fasta"
ls -l "$fasta"

n_in=$(grep -c '^>' "$fasta")
echo "Sequences in: $n_in"

echo "Building index"
"$bowtie2_build" --threads "$SLURM_CPUS_PER_TASK" --large-index "$fasta" "${index_dir}/${lib}"

# A silently truncated index would cost recall in every later run and would look
# like a MetaScope result, so count the sequences back out before trusting it.
n_out=$("$bowtie2_inspect" -s "${index_dir}/${lib}" | grep -c '^Sequence-')
echo "Sequences out: $n_out"
[[ "$n_in" -eq "$n_out" ]] || {
    echo "SEQUENCE COUNT CHANGED: $n_in in, $n_out out" >&2; exit 1; }

[[ -e "${index_dir}/${lib}.1.bt2l" ]] || {
    echo "no .1.bt2l written for $lib" >&2; exit 1; }

echo "DONE $lib"
