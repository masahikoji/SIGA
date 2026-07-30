#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript

smoke_root="${project_root}/smoke_output"
rm -rf "${smoke_root}"
mkdir -p "${smoke_root}"

Rscript "${code_dir}/check_pair_path_engine.R"

env \
  PWRT_PROFILE=smoke \
  PWRT_MODE=run \
  PWRT_SCENARIO_SET=core \
  PWRT_SCENARIO_IDS=1,2,5,6 \
  PWRT_PROJECT_DIR="${smoke_root}" \
  PWRT_CORES=1 \
  Rscript "${code_dir}/simulation_pair_path_theory_validation.R"

env \
  PWRT_PROFILE=smoke \
  PWRT_MODE=aggregate \
  PWRT_SCENARIO_SET=core \
  PWRT_SCENARIO_IDS=1,2,5,6 \
  PWRT_PROJECT_DIR="${smoke_root}" \
  PWRT_CORES=1 \
  Rscript "${code_dir}/simulation_pair_path_theory_validation.R"

summary_file="$(find "${smoke_root}" -name pair_path_theory_validation_summary.csv -type f | head -n 1)"
[[ -n "${summary_file}" && -s "${summary_file}" ]] || {
  echo "Smoke-test summary was not created." >&2
  exit 1
}

echo "PASS: additional pair-path simulation smoke test completed."
echo "Summary: ${summary_file}"
