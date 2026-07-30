#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript

scenario_set="${1:-all}"
env \
  PWRT_PROFILE="${PWRT_PROFILE:-manuscript}" \
  PWRT_MODE=calibrate \
  PWRT_SCENARIO_SET="${scenario_set}" \
  Rscript "${code_dir}/simulation_pair_path_theory_validation.R"
