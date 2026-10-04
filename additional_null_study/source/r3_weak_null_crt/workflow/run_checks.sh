#!/usr/bin/env bash
# Internal consistency checks of the R3 additions (a few seconds).
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript
Rscript "${code_dir}/check_r3_engine.R"
