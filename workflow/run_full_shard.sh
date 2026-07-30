#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript
Rscript "${workflow_dir}/record_session_info.R" "${project_root}/sessionInfo.txt"

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 N_SHARDS SHARD_ID" >&2
  exit 2
fi
n_shards="$1"
shard_id="$2"
if ! [[ "${n_shards}" =~ ^[0-9]+$ && "${shard_id}" =~ ^[0-9]+$ ]] || \
   (( n_shards < 1 || shard_id < 1 || shard_id > n_shards )); then
  echo "N_SHARDS and SHARD_ID must satisfy 1 <= SHARD_ID <= N_SHARDS." >&2
  exit 2
fi

prepare_selected_parameters

common_env=(
  "PWRT_PROJECT_DIR=${project_root}"
  "PWRT_N_OUTER=100000"
  "PWRT_N_RERAND=4999"
  "PWRT_N_CALIBRATION=100000"
  "PWRT_N_SHARDS=${n_shards}"
  "PWRT_SHARD_ID=${shard_id}"
  "PWRT_STRICT_FINAL=1"
  "PWRT_ALLOW_PARTIAL=0"
)

printf '\nRunning exact manuscript simulations: shard %s of %s.\n' "${shard_id}" "${n_shards}"
printf 'Existing checkpoints are reused by the production scripts.\n'

env "${common_env[@]}" \
  Rscript "${code_dir}/simulation_unadjusted_SIGA_vs_fixed_score_RT_power85.R"
env "${common_env[@]}" \
  Rscript "${code_dir}/simulation_covariate_adjusted_SIGA_vs_fixed_score_RT_power85.R"
env "${common_env[@]}" PWRT_BINARY_ANALYSIS=joint \
  Rscript "${code_dir}/simulation_binary_SIGA_weak_null_pureR_100K_4999_timing_split.R"
env "${common_env[@]}" PWRT_STAGE=final \
  PWRT_SMALL_N_ALTERNATIVES_CSV="${project_root}/config/binary_small_n_selected_parameters.csv" \
  Rscript "${code_dir}/simulation_binary_SIGA_small_n_margin_tuning_v3.R"
env "${common_env[@]}" \
  Rscript "${code_dir}/swift_direct_siga_case_study_power_pureR_all_in_one_v2.R"

printf '\nShard %s of %s completed. Do not aggregate until all shards have completed.\n' "${shard_id}" "${n_shards}"
