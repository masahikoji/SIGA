#!/usr/bin/env bash
# Progress from the shard CSV files (robust to restarts) plus throughput within the current scenario.
# Usage: bash workflow/progress.sh LOGFILE [PROJECT_DIR]
log="${1:?usage: progress.sh LOGFILE [PROJECT_DIR]}"
proj="${2:-${R3_PROJECT_DIR:-$HOME/r3_runs}}"
if ! grep -q "^Scenarios:" "$log"; then
  echo "STATUS: calibration or start-up in progress (no scenario started yet)"
  echo "last log line       : $(tail -1 "$log")"
  echo "R processes         : $(pgrep -f simulation_r3_weak_null_crt.R | wc -l | tr -d ' ')"
  exit 0
fi
outdir=$(grep "^Output directory:" "$log" | tail -1 | sed 's/^Output directory: //')
[ -d "$outdir" ] || outdir=$(ls -d "$proj"/simulation_output/profile_manuscript_* 2>/dev/null | head -1)
n_outer=$(grep -m1 "^Outer trials:" "$log" | sed -E 's/Outer trials: ([0-9]+);.*/\1/')
ids=$(grep "^Scenarios:" "$log" | tail -1 | sed 's/^Scenarios: //')
total=$(echo "$ids" | tr ',' '\n' | wc -l | tr -d ' ')
done_total=0; finished=0
printf "%-4s %-60s %s\n" "id" "scenario" "completed"
for id in $(echo "$ids" | tr ',' ' '); do
  dir=$(ls -d "$outdir"/scenario_shards/scenario_$(printf "%03d" "$id")_* 2>/dev/null | head -1)
  if [ -n "$dir" ]; then rows=$(cat "$dir"/shard_*.csv 2>/dev/null | grep -c -v '^"scenario_id"'); else rows=0; fi
  done_total=$((done_total+rows))
  [ "$rows" -ge "$n_outer" ] && finished=$((finished+1))
  [ "$rows" -gt 0 ] && printf "%-4s %-60s %s/%s\n" "$id" "$(basename "$dir" | sed 's/scenario_[0-9]*_//')" "$rows" "$n_outer"
done
echo "scenarios finished : $finished / $total"
echo "trials completed   : $done_total / $((total*n_outer))"
# throughput within the current scenario (lines after the last scenario header)
rate=$(awk '/^Scenario /{buf=""} /assigned trials/{buf=buf $0 "\n"} END{printf "%s", buf}' "$log" | tail -20 | awk '
  { gsub("/.*","",$2); t=$2; split($0,a,"elapsed "); split(a[2],b," ");
    e=b[1]; if (b[2] ~ /^hour/) e=e*60; if (b[2] ~ /^second/) e=e/60;
    if (NR==1){t0=t;e0=e} tn=t; en=e }
  END { if (NR>1 && en>e0) printf "%.0f", (tn-t0)/(en-e0); else print "NA" }')
echo "recent throughput  : $rate trials/min (current scenario)"
if grep -q "Shard .* completed" "$log"; then echo "STATUS: FINISHED"; else echo "STATUS: running ($(pgrep -f simulation_r3_weak_null_crt.R | wc -l | tr -d ' ') R processes)"; fi
