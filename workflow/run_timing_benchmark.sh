#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript
Rscript "${workflow_dir}/record_session_info.R" "${project_root}/sessionInfo.txt"
prepare_selected_parameters

cat <<'MSG'
The timing values in the manuscript are machine-dependent. For numerical
comparison with the reported table, run this script on the stated Apple M3
Ultra workstation with one CPU core per process and minimal background load.
The script itself forces common BLAS/OpenMP thread counts to one.
MSG

env \
  PWRT_PROJECT_DIR="${project_root}" \
  PWRT_ENGINE_SCRIPT="${code_dir}/simulation_binary_SIGA_small_n_margin_tuning_v3.R" \
  PWRT_TIMING_OUTPUT_DIR="${project_root}/standardized_timing_rerun_output" \
  PWRT_TIMING_SCOPE=all \
  PWRT_TIMING_REPEATS=3 \
  PWRT_TIMING_TRIALS=30 \
  PWRT_TIMING_RERAND=4999 \
  PWRT_TIMING_CALIBRATION=100000 \
  PWRT_SIGA_MIN_SECONDS=1.0 \
  PWRT_SIGA_MIN_CALLS=1000 \
  PWRT_TIME_CALIBRATION=1 \
  PWRT_CALIBRATION_REPEATS=3 \
  PWRT_CALIBRATION_WARMUP_PATHS=1000 \
  Rscript "${code_dir}/rerun_standardized_timing_benchmark_v3_with_calibration_fixed.R"

printf '\nTiming benchmark completed.\n'
