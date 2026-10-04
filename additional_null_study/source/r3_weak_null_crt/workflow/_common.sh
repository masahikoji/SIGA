#!/usr/bin/env bash
set -euo pipefail
workflow_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${workflow_dir}/.." && pwd)"
code_dir="${project_root}/code"
require_rscript() {
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript was not found. Install R (base R only is required)." >&2
    exit 127
  fi
}
