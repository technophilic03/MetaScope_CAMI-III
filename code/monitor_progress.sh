#!/bin/bash
# Append a progress snapshot every 15 min. Detached, so it survives a session
# restart; read output/strain_madness/logs/progress.tsv at any time.
set -uo pipefail
cd "$(dirname "$0")/.."
LOG=output/strain_madness/logs/progress.tsv
PIDFILE=output/strain_madness/.monitor.pid
# A pidfile, not pgrep: a pgrep pattern also matches the shell that runs it.
if [[ -f "$PIDFILE" ]] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
    echo "monitor already running as PID $(cat "$PIDFILE")" >&2; exit 0
fi
echo $$ > "$PIDFILE"
trap 'rm -f "$PIDFILE"' EXIT
[[ -s "$LOG" ]] || printf 'time\tb_run\tmetag_run\tlab\tpend\tdone_id\tdone_reass\tpilot_elapsed\tpilot_rss_GiB\tpilot_libs\n' > "$LOG"
JID=$(cat output/strain_madness/.pilot_jobid 2>/dev/null)
while :; do
  brun=$(squeue -h -u "$USER" -n metascope_b -t RUNNING -r 2>/dev/null | wc -l)
  mrun=$(squeue -h -u "$USER" -n metascope_metag -t RUNNING -r 2>/dev/null | wc -l)
  lab=$(squeue -h -u "$USER" -p p_wj183_1 -t RUNNING -r 2>/dev/null | wc -l)
  pend=$(squeue -h -u "$USER" -t PENDING -r -n metascope_b,metascope_metag 2>/dev/null | wc -l)
  did=$(ls output/strain_madness/metascope/results/*.metascope_id.csv output/strain_madness/metascope_blast/results/*.metascope_id.csv 2>/dev/null | wc -l)
  dre=$(ls output/strain_madness/metascope_blast/results/*.metascope_blast_reassigned.csv 2>/dev/null | wc -l)
  pel=$(squeue -h -j "${JID}_0" -o '%M' 2>/dev/null)
  prss=$(sstat -a -j "${JID}_0.batch" --format=MaxRSS -n -P 2>/dev/null | head -1 | tr -d 'K')
  prss=$(awk -v k="${prss:-0}" 'BEGIN{printf "%.1f", k/1048576}')
  plibs=$(ls output/strain_madness/metascope_blast/scratch/strmgCAMI2_sample_0_${JID}/*.bam 2>/dev/null | wc -l)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date '+%m-%d %H:%M')" "$brun" "$mrun" "$lab" "$pend" "$did" "$dre" "${pel:--}" "$prss" "$plibs" >> "$LOG"
  sleep 900
done
