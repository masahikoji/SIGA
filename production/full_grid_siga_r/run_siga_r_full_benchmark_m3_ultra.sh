#!/bin/bash
set -euo pipefail

# Comprehensive manuscript run on one Apple M3 Ultra.
# Place these files in the same directory before running:
#   simulation_siga_r_full_benchmark.R
#   aggregate_siga_r_full_benchmark.R
#   siga_pair_path_engine.R
# Optional verification files:
#   check_pair_path_engine.R
#   00_verify_rt_kernel_equivalence.R

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

SIM_SCRIPT="${PWRT_SIM_SCRIPT:-$SCRIPT_DIR/simulation_siga_r_full_benchmark.R}"
AGG_SCRIPT="${PWRT_AGG_SCRIPT:-$SCRIPT_DIR/aggregate_siga_r_full_benchmark.R}"
PAIR_ENGINE="${PWRT_PAIR_ENGINE:-$SCRIPT_DIR/siga_pair_path_engine.R}"
OUTPUT_DIR="${PWRT_OUTPUT_DIR:-$SCRIPT_DIR/siga_r_full_benchmark_output}"
WORKERS="${PWRT_WORKERS:-24}"
N_OUTER="${PWRT_N_OUTER:-100000}"
N_RERAND="${PWRT_N_RERAND:-4999}"
N_CALIBRATION="${PWRT_N_CALIBRATION:-100000}"
OUTER_BATCH="${PWRT_OUTER_BATCH:-10}"
CALIBRATION_BATCH="${PWRT_CALIBRATION_BATCH:-1000}"
BASE_SEED="${PWRT_SEED:-20260801}"
P_BIASED_COIN="${PWRT_P_BIASED_COIN:-0.80}"
EPSILON_EXPONENT="${PWRT_EPSILON_EXPONENT:-1.0}"
PBC_TAG="$(printf '%.3f' "$P_BIASED_COIN" | tr '.' 'p')"
EPS_TAG="$(printf '%.2f' "$EPSILON_EXPONENT" | tr '.' 'p')"
RUN_TAG="M${N_OUTER}_B${N_RERAND}_Bpsi${N_CALIBRATION}_pbc${PBC_TAG}_eps${EPS_TAG}_seed${BASE_SEED}"
RUN_DIR="$OUTPUT_DIR/$RUN_TAG"

for file in "$SIM_SCRIPT" "$AGG_SCRIPT" "$PAIR_ENGINE"; do
  if [ ! -f "$file" ]; then
    echo "Required file not found: $file" >&2
    exit 1
  fi
done

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript was not found in PATH." >&2
  exit 1
fi

case "$WORKERS" in
  ''|*[!0-9]*) echo "PWRT_WORKERS must be a positive integer." >&2; exit 1 ;;
esac
if [ "$WORKERS" -lt 1 ]; then
  echo "PWRT_WORKERS must be at least 1." >&2
  exit 1
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1
export PWRT_PAIR_ENGINE="$PAIR_ENGINE"
export PWRT_OUTPUT_DIR="$OUTPUT_DIR"
export PWRT_N_OUTER="$N_OUTER"
export PWRT_N_RERAND="$N_RERAND"
export PWRT_N_CALIBRATION="$N_CALIBRATION"
export PWRT_OUTER_BATCH="$OUTER_BATCH"
export PWRT_CALIBRATION_BATCH="$CALIBRATION_BATCH"
export PWRT_SEED="$BASE_SEED"
export PWRT_P_BIASED_COIN="$P_BIASED_COIN"
export PWRT_EPSILON_EXPONENT="$EPSILON_EXPONENT"

mkdir -p "$RUN_DIR/logs"

cat <<SETTINGS
Comprehensive SIGA-R benchmark
  script directory          : $SCRIPT_DIR
  output base               : $OUTPUT_DIR
  run directory              : $RUN_DIR
  worker processes          : $WORKERS
  outer trials/scenario     : $N_OUTER
  RT paths/outer trial      : $N_RERAND
  three-copy calibrations   : $N_CALIBRATION
  biased-coin probability   : $P_BIASED_COIN
  base seed                 : $BASE_SEED
  epsilon_n                 : n^(-$EPSILON_EXPONENT)
SETTINGS

# Optional unit checks when the files are available.
if [ -f "$SCRIPT_DIR/check_pair_path_engine.R" ]; then
  echo "Running pair-path engine check..."
  Rscript "$SCRIPT_DIR/check_pair_path_engine.R"
fi
if [ -f "$SCRIPT_DIR/00_verify_rt_kernel_equivalence.R" ]; then
  echo "Running RT-kernel equivalence check..."
  Rscript "$SCRIPT_DIR/00_verify_rt_kernel_equivalence.R"
fi

# Create calibrations once before workers start. This prevents simultaneous
# workers from waiting on calibration locks during the main run.
echo "Creating or validating the four reusable calibrations..."
PWRT_MODE=calibrate Rscript "$SIM_SCRIPT" \
  > "$RUN_DIR/logs/calibration.log" 2>&1

echo "Calibration complete. Launching $WORKERS simulation shards..."

pids=""
for shard in $(seq 1 "$WORKERS"); do
  log_file=$(printf "%s/logs/shard_%04d_of_%04d.log" "$RUN_DIR" "$shard" "$WORKERS")
  (
    export PWRT_MODE=run
    export PWRT_N_SHARDS="$WORKERS"
    export PWRT_SHARD_ID="$shard"
    Rscript "$SIM_SCRIPT"
  ) > "$log_file" 2>&1 &
  pid=$!
  pids="$pids $pid"
  printf "  shard %d/%d started: PID=%d, log=%s\n" "$shard" "$WORKERS" "$pid" "$log_file"
done

failure=0
for pid in $pids; do
  if ! wait "$pid"; then
    echo "Worker PID $pid failed. Inspect the shard logs." >&2
    failure=1
  fi
done

if [ "$failure" -ne 0 ]; then
  echo "At least one worker failed. Re-running this launcher will resume completed shards." >&2
  exit 1
fi

echo "All workers completed. Aggregating and creating TeX tables..."
Rscript "$AGG_SCRIPT" > "$RUN_DIR/logs/aggregation.log" 2>&1

echo "Completed successfully."
echo "Main summary: $RUN_DIR/publication_tables/main_siga_r_comprehensive_summary.tex"
echo "Full summary: $RUN_DIR/publication_tables/siga_r_full_benchmark_summary.csv"
