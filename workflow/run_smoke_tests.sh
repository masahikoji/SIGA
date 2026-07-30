#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript
Rscript "${workflow_dir}/record_session_info.R" "${project_root}/sessionInfo.txt"

smoke_root="${project_root}/smoke_run"
rm -rf "${smoke_root}"
mkdir -p "${smoke_root}"
cp "${project_root}/config/binary_small_n_selected_parameters.csv" \
   "${smoke_root}/binary_small_n_selected_parameters.csv"

common_env=(
  "PWRT_PROJECT_DIR=${smoke_root}"
  "PWRT_N_OUTER=12"
  "PWRT_N_RERAND=19"
  "PWRT_N_CALIBRATION=200"
  "PWRT_N_SHARDS=1"
  "PWRT_SHARD_ID=1"
  "PWRT_OUTER_BATCH=4"
  "PWRT_CORES=1"
  "PWRT_STRICT_FINAL=0"
  "PWRT_ALLOW_PARTIAL=1"
)

printf '\n[1/6] Verifying the two fixed-score randomization-test kernels...\n'
Rscript "${code_dir}/00_verify_rt_kernel_equivalence.R"

printf '\n[2/6] Continuous, unadjusted, one scenario...\n'
env "${common_env[@]}" PWRT_SCENARIO_IDS=1 \
  Rscript "${code_dir}/simulation_unadjusted_SIGA_vs_fixed_score_RT_power85.R"

printf '\n[3/6] Continuous, adjusted, one scenario...\n'
env "${common_env[@]}" PWRT_SCENARIO_IDS=1 \
  Rscript "${code_dir}/simulation_covariate_adjusted_SIGA_vs_fixed_score_RT_power85.R"

printf '\n[4/6] Binary large-sample engine, one scenario...\n'
env "${common_env[@]}" PWRT_SCENARIO_IDS=1 PWRT_BINARY_ANALYSIS=joint \
  Rscript "${code_dir}/simulation_binary_SIGA_weak_null_pureR_100K_4999_timing_split.R"
env "${common_env[@]}" PWRT_SCENARIO_IDS=1 PWRT_BINARY_ANALYSIS=joint PWRT_MODE=aggregate \
  Rscript "${code_dir}/simulation_binary_SIGA_weak_null_pureR_100K_4999_timing_split.R"

printf '\n[5/6] Binary small-sample final engine, one scenario...\n'
env "${common_env[@]}" PWRT_STAGE=final PWRT_SCENARIO_IDS=15 \
  PWRT_SMALL_N_ALTERNATIVES_CSV="${smoke_root}/binary_small_n_selected_parameters.csv" \
  Rscript "${code_dir}/simulation_binary_SIGA_small_n_margin_tuning_v3.R"
env "${common_env[@]}" PWRT_STAGE=final PWRT_SCENARIO_IDS=15 PWRT_MODE=aggregate \
  PWRT_SMALL_N_ALTERNATIVES_CSV="${smoke_root}/binary_small_n_selected_parameters.csv" \
  Rscript "${code_dir}/simulation_binary_SIGA_small_n_margin_tuning_v3.R"

printf '\n[6/6] SWIFT DIRECT-inspired engine...\n'
env "${common_env[@]}" \
  Rscript "${code_dir}/swift_direct_siga_case_study_power_pureR_all_in_one_v2.R"
env "${common_env[@]}" PWRT_MODE=aggregate \
  Rscript "${code_dir}/swift_direct_siga_case_study_power_pureR_all_in_one_v2.R"

required_files=(
  "${smoke_root}/unadjusted_siga_vs_fixed_score_rt_power85_output/unadjusted_siga_vs_fixed_score_rt_power85_summary.csv"
  "${smoke_root}/covariate_adjusted_siga_vs_fixed_score_rt_power85_output/covariate_adjusted_siga_vs_fixed_score_rt_power85_summary.csv"
  "${smoke_root}/binary_siga_weak_null_pureR_100k_4999_output/binary_siga_weak_null_summary.csv"
  "${smoke_root}/binary_siga_weak_null_small_n_100K_4999_output/binary_siga_weak_null_small_n_summary.csv"
  "${smoke_root}/swift_direct_siga_case_study_pureR_output/swift_direct_siga_case_study_summary.csv"
)
for path in "${required_files[@]}"; do
  [[ -s "${path}" ]] || { echo "Smoke-test output missing or empty: ${path}" >&2; exit 1; }
done

printf '\nPASS: all production engines completed a reduced smoke test.\n'
printf 'Smoke-test outputs: %s\n' "${smoke_root}"
