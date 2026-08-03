#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

R_SCRIPT="${SWIFT_R_SCRIPT:-$SCRIPT_DIR/swift_direct_siga_r_case_study_extension_v2.R}"
PAIR_ENGINE="${SWIFT_PAIR_ENGINE:-$SCRIPT_DIR/siga_pair_path_engine.R}"
WORKERS="${SWIFT_R_WORKERS:-24}"
N_OUTER="${PWRT_N_OUTER:-100000}"
N_PAIR_CALIBRATION="${PWRT_N_PAIR_CALIBRATION:-100000}"
BASE_SEED="${PWRT_SEED:-20260724}"
PROJECT_DIR="${PWRT_PROJECT_DIR:-$HOME/Dropbox/Tex/03_Post_PhD/029_exact_minimization/02_program}"
OUTPUT_DIR="${SWIFT_R_OUTPUT_DIR:-$PROJECT_DIR/swift_direct_siga_r_case_study_output}"

for file in "$R_SCRIPT" "$PAIR_ENGINE"; do
  if [ ! -f "$file" ]; then
    echo "Required file not found: $file" >&2
    exit 1
  fi
done

case "$WORKERS" in
  ''|*[!0-9]*) echo "SWIFT_R_WORKERS must be a positive integer." >&2; exit 1 ;;
esac
if [ "$WORKERS" -lt 1 ]; then
  echo "SWIFT_R_WORKERS must be at least 1." >&2
  exit 1
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1
export SWIFT_PAIR_ENGINE="$PAIR_ENGINE"
export PWRT_PROJECT_DIR="$PROJECT_DIR"
export SWIFT_R_OUTPUT_DIR="$OUTPUT_DIR"
export PWRT_N_OUTER="$N_OUTER"
export PWRT_N_PAIR_CALIBRATION="$N_PAIR_CALIBRATION"
export PWRT_SEED="$BASE_SEED"

mkdir -p "$OUTPUT_DIR/logs"

cat <<SETTINGS
SWIFT DIRECT SIGA-R case-study extension
  script                    : $R_SCRIPT
  project directory         : $PROJECT_DIR
  output directory          : $OUTPUT_DIR
  worker processes          : $WORKERS
  outer trials              : $N_OUTER
  pair calibration paths    : $N_PAIR_CALIBRATION
  base seed                 : $BASE_SEED
SETTINGS

echo "Creating or validating the pair-path calibration..."
SWIFT_R_MODE=calibrate Rscript "$R_SCRIPT" \
  > "$OUTPUT_DIR/logs/calibration.log" 2>&1

echo "Calibration complete. Launching $WORKERS SIGA-R shards..."
pids=""
for shard in $(seq 1 "$WORKERS"); do
  log_file=$(printf "%s/logs/shard_%04d_of_%04d.log" "$OUTPUT_DIR" "$shard" "$WORKERS")
  (
    export SWIFT_R_MODE=run
    export SWIFT_R_N_SHARDS="$WORKERS"
    export SWIFT_R_SHARD_ID="$shard"
    Rscript "$R_SCRIPT"
  ) > "$log_file" 2>&1 &
  pid=$!
  pids="$pids $pid"
  printf "  shard %d/%d started: PID=%d\n" "$shard" "$WORKERS" "$pid"
done

failure=0
for pid in $pids; do
  if ! wait "$pid"; then
    echo "Worker PID $pid failed. Inspect shard logs." >&2
    failure=1
  fi
done
if [ "$failure" -ne 0 ]; then
  echo "At least one worker failed. Re-running this launcher resumes completed shards." >&2
  exit 1
fi

echo "All shards completed. Aggregating..."
SWIFT_R_MODE=aggregate Rscript "$R_SCRIPT" \
  > "$OUTPUT_DIR/logs/aggregation.log" 2>&1

echo "Completed successfully."
echo "Summary: $OUTPUT_DIR/publication/swift_direct_siga_r_summary.csv"
echo "Table:   $OUTPUT_DIR/publication/swift_direct_siga_r_table.tex"
