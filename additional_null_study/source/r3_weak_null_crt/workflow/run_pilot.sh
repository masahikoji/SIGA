#!/usr/bin/env bash
# Pilot run on one machine: 2,000 outer trials, 999 regenerated paths, 20,000 calibration paths.
# Usage: bash workflow/run_pilot.sh [SCENARIO_SET]      (default: null)
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript
export R3_PROFILE=pilot
export R3_PROJECT_DIR="${R3_PROJECT_DIR:-${project_root}/r3_output_pilot}"
export R3_SCENARIO_SET="${1:-null}"
R3_MODE=run Rscript "${code_dir}/simulation_r3_weak_null_crt.R"
R3_MODE=aggregate R3_SKIP_DETAILS=1 Rscript "${code_dir}/simulation_r3_weak_null_crt.R"
