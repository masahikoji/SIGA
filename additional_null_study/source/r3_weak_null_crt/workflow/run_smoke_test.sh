#!/usr/bin/env bash
# End-to-end smoke test: 20 outer trials, 99 regenerated paths, 500 calibration paths,
# one scenario per objective/outcome.  Verifies execution only.
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript
export R3_PROFILE=smoke
export R3_PROJECT_DIR="${R3_PROJECT_DIR:-${project_root}/r3_output_smoke}"
export R3_CORES="${R3_CORES:-1}"
ids="2,8,14,17,23,29,38,41"
Rscript "${code_dir}/check_r3_engine.R"
R3_MODE=run R3_SCENARIO_IDS="${ids}" Rscript "${code_dir}/simulation_r3_weak_null_crt.R"
R3_MODE=aggregate R3_SCENARIO_IDS="${ids}" R3_SCENARIO_SET=null R3_SKIP_DETAILS=1 \
  Rscript "${code_dir}/simulation_r3_weak_null_crt.R"
printf '\nSmoke test completed. Output: %s\n' "${R3_PROJECT_DIR}"
