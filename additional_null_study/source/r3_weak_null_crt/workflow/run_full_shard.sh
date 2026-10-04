#!/usr/bin/env bash
# One manuscript-scale shard.
# Usage: bash workflow/run_full_shard.sh N_SHARDS SHARD_ID [SCENARIO_SET]
#   SCENARIO_SET: null (default) | power | interaction | manuscript | all
# Environment overrides (optional): R3_PROJECT_DIR, R3_CORES, R3_PBC, R3_CONT_SD,
#   R3_CONT_MAX_D, R3_CONT_ETA_SD, R3_BIN_MAX_D, R3_BINARY_CONTINUITY, R3_SCENARIO_IDS.
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript
if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 N_SHARDS SHARD_ID [SCENARIO_SET]" >&2
  exit 2
fi
n_shards="$1"; shard_id="$2"; scenario_set="${3:-null}"
if ! [[ "${n_shards}" =~ ^[0-9]+$ && "${shard_id}" =~ ^[0-9]+$ ]] || (( n_shards < 1 || shard_id < 1 || shard_id > n_shards )); then
  echo "N_SHARDS and SHARD_ID must satisfy 1 <= SHARD_ID <= N_SHARDS." >&2
  exit 2
fi
export R3_PROFILE=manuscript
export R3_PROJECT_DIR="${R3_PROJECT_DIR:-${project_root}/r3_output}"
export R3_N_OUTER="${R3_N_OUTER:-100000}"
export R3_N_RERAND="${R3_N_RERAND:-4999}"
export R3_N_CALIBRATION="${R3_N_CALIBRATION:-100000}"
export R3_N_SHARDS="${n_shards}"
export R3_SHARD_ID="${shard_id}"
export R3_SCENARIO_SET="${scenario_set}"
export R3_MODE=run
printf '\nR3 simulation: scenario set %s, shard %s of %s (existing checkpoints are reused).\n' "${scenario_set}" "${shard_id}" "${n_shards}"
Rscript "${code_dir}/simulation_r3_weak_null_crt.R"
printf '\nShard %s of %s completed. Aggregate only after all shards have finished.\n' "${shard_id}" "${n_shards}"
