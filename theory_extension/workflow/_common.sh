#!/usr/bin/env bash
set -euo pipefail

workflow_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${workflow_dir}/.." && pwd)"
code_dir="${project_root}/code"

require_rscript() {
  command -v Rscript >/dev/null 2>&1 || {
    echo "Rscript was not found. Install R before running this workflow." >&2
    exit 127
  }
}
