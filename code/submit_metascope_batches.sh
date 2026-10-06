#!/bin/bash
#
# Submit the MetaScope+BLAST array in batches of 20 rather than as one 1-99
# array with a %8 throttle.
#
# WHY BATCHES. A throttled array keeps a trickle of tasks running for weeks and
# gives no natural checkpoint. Batches of 20 finish together, so after roughly
# every 1.8 days there is a complete, inspectable set of 20 samples: failures can
# be recovered with rerun_blast_only.sbatch before the next batch starts, and the
# run can be stopped between batches without abandoning work in flight.
#
# Batch width is 20 because that is the scale of the CAMI III run this is compared
# against, and because capacity is not the limit: main holds 117 nodes of 504 GB /
# 64 cpu and mem another 13 of 1-2 TB, which is 145+ concurrent slots at 450G.
#
# WHY afterany, NOT afterok. Five of CAMI III's nineteen tasks hit the memory cap
# inside metascope_blast. With afterok a single such task would stall the whole
# chain; afterany keeps it moving and the failures are recovered separately.
#
# Sample 0 is not here: it ran as the pilot (job in output/strain_madness/.pilot_jobid).
#
# Usage:
#   bash code/submit_metascope_batches.sh            # submit the whole chain
#   bash code/submit_metascope_batches.sh 21 99      # resubmit from a given point
#
set -euo pipefail
cd "$(dirname "$0")/.."

FIRST=${1:-1}
LAST=${2:-99}
WIDTH=20

dep=""
lo=$FIRST
while [[ $lo -le $LAST ]]; do
    hi=$(( lo + WIDTH - 1 ))
    [[ $hi -gt $LAST ]] && hi=$LAST

    # Every sample in the batch needs its tarball before the batch is queued.
    # Later batches are checked at submit time too; the download finishes long
    # before batch 2 starts, so this is a guard against a silent gap.
    for i in $(seq "$lo" "$hi"); do
        f="strmg_db/strmgCAMI2_sample_${i}_reads.tar.gz"
        [[ -s "$f" ]] || { echo "MISSING READS FOR BATCH ${lo}-${hi}: $f" >&2; exit 1; }
    done

    jid=$(sbatch --parsable ${dep} --array="${lo}-${hi}" code/run_metascope_b.sbatch)
    echo "batch ${lo}-${hi}  job ${jid}${dep:+  (starts after ${dep##*:})}"
    echo "${lo}-${hi} ${jid}" >> output/strain_madness/.metascope_batches
    dep="--dependency=afterany:${jid}"
    lo=$(( hi + 1 ))
done
