#!/usr/bin/env bash
# Recompute the 16 binary null scenarios with the tie-inclusive corrected test (the continuous
# scenarios are unaffected because ties have probability zero there).  Deletes the shard files of
# those scenarios in the manuscript-profile output directory, reruns them, and re-aggregates.
# Usage: R3_CORES=24 bash workflow/rerun_binary_null.sh
source "$(cd "$(dirname "$0")" && pwd)/_common.sh"
require_rscript
export R3_PROFILE=manuscript
export R3_PROJECT_DIR="${R3_PROJECT_DIR:-$HOME/r3_runs}"
ids="23,29,35,38,65,71,77,80,107,113,119,122,149,155,161,164"
outdir=$(ls -d "$R3_PROJECT_DIR"/simulation_output/profile_manuscript_M100000_B4999_B0100000_S1_pbc0.8_sd1_dC0.5_dB0.15 2>/dev/null | head -1)
if [ -z "$outdir" ]; then echo "manuscript-profile output directory not found under $R3_PROJECT_DIR" >&2; exit 1; fi
for id in $(echo "$ids" | tr ',' ' '); do
  for d in "$outdir"/scenario_shards/scenario_$(printf "%03d" "$id")_*; do
    [ -d "$d" ] && rm -rf "$d" && echo "removed $(basename "$d")"
  done
done
R3_MODE=run R3_SCENARIO_IDS="$ids" R3_SCENARIO_SET=null Rscript "${code_dir}/simulation_r3_weak_null_crt.R"
R3_MODE=aggregate R3_SCENARIO_SET=null R3_SKIP_DETAILS=1 Rscript "${code_dir}/simulation_r3_weak_null_crt.R"
