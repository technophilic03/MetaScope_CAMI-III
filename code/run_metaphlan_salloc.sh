#!/bin/bash -l

# Run MetaPhlAn over a list of samples as job steps inside one salloc
# allocation, four at a time, instead of as an sbatch array.
#
# Why an allocation. The mem QOS caps a user at 100 submitted jobs, which the
# MetaScope metag chain already fills, and every array task counts as one. An
# allocation is a single job however many samples run inside it, and samples do
# not queue again one by one once it is granted.
#
# Usage, from a login node:
#   salloc --no-shell -p p_wj183_1 -N1 -n1 -c 64 --mem=240G --time=08:00:00 -J metaphlan_salloc
#   bash code/run_metaphlan_salloc.sh <allocation job id> 0 1 2 5 ...
#
# Each sample is one step of 16 cpus and 55G, so four fit a 64-cpu, 240G
# allocation. Memory is what sets the step size: bowtie2 holds the vJan25 index,
# 35 GB on disk, in memory once per process.
#
# One log per sample: output/strain_madness/logs/metaphlan_<sampleName>.log

set -euo pipefail

[[ $# -ge 2 ]] || { echo "usage: $0 allocJobId idx [idx ...]" >&2; exit 1; }
allocJob="$1"
shift

proj=/projectsp/f_wj183_1/work/Yaoan/CAMI-II
dbDir="${proj}/db/metaphlan_db_vJan25"
index=mpa_vJan25_CHOCOPhlAnSGB_202503
outDir="${proj}/output/strain_madness/metaphlan/results"
scratchRoot="${proj}/output/strain_madness/metaphlan/scratch"
logDir="${proj}/output/strain_madness/logs"
threads=16
parallel=4

# The steps inherit this shell's PATH, which is how metaphlan and its bowtie2
# reach the compute node.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate metaphlan
command -v metaphlan >/dev/null || { echo "metaphlan is not in conda env 'metaphlan'" >&2; exit 1; }

mkdir -p "$outDir" "$scratchRoot" "$logDir"
export proj dbDir index outDir scratchRoot logDir threads allocJob

# xargs exits non-zero if any sample fails, and set -e passes that on, so a
# failure is never reported as a finished run.
printf '%s\n' "$@" | xargs -P "$parallel" -I{} bash -c '
    i={}
    s="strmgCAMI2_sample_${i}"
    if srun --jobid="$allocJob" --exact -N1 -n1 -c "$threads" --mem=55G \
           bash "${proj}/code/run_metaphlan.sh" \
               "${proj}/strmg_db/${s}_reads.tar.gz" "$s" "strmgCAMI2_short_read_sample_${i}" \
               "$outDir" "${scratchRoot}/${s}" "$threads" "$dbDir" "$index" \
           > "${logDir}/metaphlan_${s}.log" 2>&1; then
        echo "DONE   $s"
    else
        echo "FAILED $s  (${logDir}/metaphlan_${s}.log)"
        exit 1
    fi
'
