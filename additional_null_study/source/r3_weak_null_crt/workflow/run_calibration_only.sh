#!/usr/bin/env bash
# Create (or reuse) the three-copy allocation calibration for all four designs
# and write the scenario definitions, without simulating outcomes.
# Usage: bash workflow/run_calibration_only.sh [profile]   (default: manuscript)
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript
export R3_PROFILE="${1:-manuscript}"
export R3_PROJECT_DIR="${R3_PROJECT_DIR:-${project_root}/r3_output}"
R3_MODE=calibrate R3_SCENARIO_SET=all Rscript "${code_dir}/simulation_r3_weak_null_crt.R"
