#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 N_SHARDS" >&2
  exit 2
fi
n_shards="$1"
if ! [[ "${n_shards}" =~ ^[0-9]+$ ]] || (( n_shards < 1 )); then
  echo "N_SHARDS must be a positive integer." >&2
  exit 2
fi

prepare_selected_parameters
common_env=(
  "PWRT_PROJECT_DIR=${project_root}"
  "PWRT_N_OUTER=100000"
  "PWRT_N_RERAND=4999"
  "PWRT_N_CALIBRATION=100000"
  "PWRT_N_SHARDS=${n_shards}"
  "PWRT_MODE=aggregate"
  "PWRT_STRICT_FINAL=1"
  "PWRT_ALLOW_PARTIAL=0"
)

printf '\nAggregating continuous unadjusted and adjusted results...\n'
env "${common_env[@]}" Rscript "${code_dir}/simulation_unadjusted_SIGA_vs_fixed_score_RT_power85.R"
env "${common_env[@]}" Rscript "${code_dir}/simulation_covariate_adjusted_SIGA_vs_fixed_score_RT_power85.R"
env PWRT_PROJECT_DIR="${project_root}" Rscript "${code_dir}/aggregate_all_siga_fixed_score_power85.R"

printf '\nAggregating binary large-sample results...\n'
env "${common_env[@]}" PWRT_BINARY_ANALYSIS=joint \
  Rscript "${code_dir}/simulation_binary_SIGA_weak_null_pureR_100K_4999_timing_split.R"

printf '\nAggregating binary small-sample results...\n'
env "${common_env[@]}" PWRT_STAGE=final \
  PWRT_SMALL_N_ALTERNATIVES_CSV="${project_root}/config/binary_small_n_selected_parameters.csv" \
  Rscript "${code_dir}/simulation_binary_SIGA_small_n_margin_tuning_v3.R"

printf '\nAggregating SWIFT DIRECT-inspired results...\n'
env "${common_env[@]}" \
  Rscript "${code_dir}/swift_direct_siga_case_study_power_pureR_all_in_one_v2.R"

printf '\nSimulation aggregation completed.\n'
if [[ -s "${project_root}/standardized_timing_rerun_output/clean_timing_benchmark_summary_split_complete.csv" ]]; then
  "${workflow_dir}/finalize_manuscript_outputs.sh"
else
  printf 'Timing output is not present. Run workflow/run_timing_benchmark.sh, then workflow/finalize_manuscript_outputs.sh.\n'
fi
