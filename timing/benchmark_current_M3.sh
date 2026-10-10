#!/bin/bash
# Benchmark only: do not alter or rerun the operating-characteristic simulations.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE="${SIGA_CODE:-$HOME/Dropbox/Tex/03_Post_PhD/029_exact_minimization/02_program/SIGA_M3Ultra_d012_20261009}"
RUN="${SIGA_RUN:-$HOME/SIGA_runs/main_d012_20261009}"
PY="${SIGA_PYTHON:-$HOME/.venvs/siga_m3_d012/bin/python}"
if [ ! -x "$PY" ]; then
  printf 'Existing simulation Python not found: %s\n' "$PY" >&2
  exit 1
fi
if [ ! -f "$CODE/run.py" ] || [ ! -f "$RUN/plan.json" ]; then
  printf 'Existing code or frozen production plan not found.\nCode: %s\nRun: %s\n' "$CODE" "$RUN" >&2
  exit 1
fi
export NUMBA_CACHE_DIR="$HOME/.cache/siga_timing_numba"
export PYTHONPYCACHEPREFIX="$HOME/.cache/siga_timing_python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1
if command -v caffeinate >/dev/null 2>&1; then
  caffeinate -is "$PY" -u "$HERE/benchmark_methods.py" --code "$CODE" --run "$RUN" "$@"
else
  "$PY" -u "$HERE/benchmark_methods.py" --code "$CODE" --run "$RUN" "$@"
fi
