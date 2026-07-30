#!/usr/bin/env bash
set -euo pipefail

workflow_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${workflow_dir}/.." && pwd)"
if [[ -d "${project_root}/R/production" ]]; then
  code_dir="${project_root}/R/production"
elif [[ -d "${project_root}/code" ]]; then
  code_dir="${project_root}/code"
else
  echo "Could not locate the production-code directory." >&2
  exit 1
fi

require_rscript() {
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript was not found. Install R before running this workflow." >&2
    exit 127
  fi
}

prepare_selected_parameters() {
  cp "${project_root}/config/binary_small_n_selected_parameters.csv" \
     "${project_root}/binary_small_n_selected_parameters.csv"
}
