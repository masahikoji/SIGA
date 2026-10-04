#!/usr/bin/env bash
# Aggregate all shards of a manuscript-scale run and write CSV summaries and LaTeX tables.
# Usage: bash workflow/aggregate_full_results.sh N_SHARDS [SCENARIO_SET]
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 N_SHARDS [SCENARIO_SET]" >&2
  exit 2
fi
export R3_PROFILE=manuscript
export R3_PROJECT_DIR="${R3_PROJECT_DIR:-${project_root}/r3_output}"
export R3_N_OUTER="${R3_N_OUTER:-100000}"
export R3_N_RERAND="${R3_N_RERAND:-4999}"
export R3_N_CALIBRATION="${R3_N_CALIBRATION:-100000}"
export R3_N_SHARDS="$1"
export R3_SCENARIO_SET="${2:-null}"
export R3_MODE=aggregate
Rscript "${code_dir}/simulation_r3_weak_null_crt.R"
