#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 N_SHARDS [core|sensitivity|all]" >&2
  exit 2
fi

n_shards="$1"
scenario_set="${2:-core}"

env \
  PWRT_PROFILE=manuscript \
  PWRT_MODE=aggregate \
  PWRT_SCENARIO_SET="${scenario_set}" \
  PWRT_N_SHARDS="${n_shards}" \
  PWRT_SHARD_ID=1 \
  Rscript "${code_dir}/simulation_pair_path_theory_validation.R"
