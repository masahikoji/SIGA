#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript

SIGA_PROJECT_DIR="${project_root}" Rscript "${workflow_dir}/build_manuscript_outputs.R"
SIGA_PROJECT_DIR="${project_root}" Rscript "${workflow_dir}/create_manuscript_tables.R"
SIGA_PROJECT_DIR="${project_root}" SIGA_VERIFY_TIMING="${SIGA_VERIFY_TIMING:-0}" \
  Rscript "${workflow_dir}/verify_reported_results.R"
printf '\nFinal manuscript-level CSV files are in %s/manuscript_outputs.\n' "${project_root}"
