#!/usr/bin/env bash
# Polls a run log every 5 minutes and raises a macOS notification (and a spoken message)
# when the shard has completed.  Usage: nohup bash workflow/notify_when_done.sh LOGFILE &
# Waits up to 10 minutes for the R process to appear before treating its absence as an error.
log="${1:?usage: notify_when_done.sh LOGFILE}"
grace=0
until grep -q "Shard .* completed" "$log" 2>/dev/null; do
  if ! pgrep -f simulation_r3_weak_null_crt.R >/dev/null; then
    grace=$((grace+1))
    if [ "$grace" -gt 20 ]; then
      osascript -e 'display notification "R process is no longer running; check the log" with title "R3 simulation"' 2>/dev/null
      echo "R process not running at $(date); check $log"; exit 1
    fi
    sleep 30; continue
  fi
  grace=0
  sleep 300
done
osascript -e 'display notification "Shard completed; run the aggregation step" with title "R3 simulation"' 2>/dev/null
say "R3 simulation finished" 2>/dev/null
echo "finished at $(date)"
