#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 N_SHARDS SHARD_ID [core|sensitivity|all]" >&2
  exit 2
fi

n_shards="$1"
shard_id="$2"
scenario_set="${3:-core}"

env \
  PWRT_PROFILE=manuscript \
  PWRT_MODE=run \
  PWRT_SCENARIO_SET="${scenario_set}" \
  PWRT_N_SHARDS="${n_shards}" \
  PWRT_SHARD_ID="${shard_id}" \
  Rscript "${code_dir}/simulation_pair_path_theory_validation.R"
