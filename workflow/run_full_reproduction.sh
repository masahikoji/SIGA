#!/usr/bin/env bash
set -euo pipefail
workflow_dir="$(cd "$(dirname "$0")" && pwd)"
cat <<'MSG'
WARNING: this is the exact unsharded manuscript workflow:
  100,000 outer trials per scenario,
  4,999 regenerated allocations per outer trial,
  100,000 allocation-only calibration paths,
plus the complete timing benchmark.
It is computationally very expensive. A sharded run is normally preferable.
MSG
"${workflow_dir}/run_full_shard.sh" 1 1
"${workflow_dir}/aggregate_full_results.sh" 1
"${workflow_dir}/run_timing_benchmark.sh"
SIGA_VERIFY_TIMING="${SIGA_VERIFY_TIMING:-0}" "${workflow_dir}/finalize_manuscript_outputs.sh"
