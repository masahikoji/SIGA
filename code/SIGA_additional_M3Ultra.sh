#!/usr/bin/env bash
# Auxiliary SIGA validation suites for the existing d012 code, on one Apple Silicon Mac.
# Does not modify the finished 32-scenario primary run or any source .py file.
set -euo pipefail
ACTION="${1:-help}"
CODE="$HOME/Dropbox/Tex/03_Post_PhD/029_exact_minimization/02_program/SIGA_M3Ultra_d012_20261009"
RUN="$HOME/SIGA_runs/additional_d012_20261009"
PY="$HOME/.venvs/siga_m3_d012/bin/python"
SUITES='controls,sample_size,pbc95'
export PYTHONNOUSERSITE=1
unset PYTHONHOME PYTHONPATH
export SIGA_BLAS_THREADS=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMBA_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
export NUMBA_CACHE_DIR="$HOME/.cache/siga_m3_d012/numba"
export PYTHONPYCACHEPREFIX="$HOME/.cache/siga_m3_d012/python"
if [ "$ACTION" = help ]; then
  echo 'Usage: bash SIGA_additional_M3Ultra.sh prepare|run|status|aggregate'
  echo 'Optional: SIGA_WORKERS=16 bash SIGA_additional_M3Ultra.sh run'
  echo "Separate results: $RUN"
  exit 0
fi
if [ ! -f "$CODE/run.py" ] || [ ! -x "$PY" ]; then
  echo 'Original d012 source or its Python environment was not found.' >&2
  exit 1
fi
mkdir -p "$NUMBA_CACHE_DIR" "$PYTHONPYCACHEPREFIX"
cd "$CODE"
case "$ACTION" in
  prepare)
    if [ ! -f "$RUN/plan.json" ]; then
      "$PY" -u run.py init --out "$RUN" --suites "$SUITES" --profile production
    fi
    "$PY" -u run.py check --out "$RUN"
    "$PY" -u run.py calibrate --out "$RUN" --workers 2
    echo 'Additional 20-scenario design is ready. Next: bash SIGA_additional_M3Ultra.sh run'
    ;;
  run)
    "$PY" -u run.py check --out "$RUN"
    WORKERS="${SIGA_WORKERS:-16}"
    if [ -x /usr/bin/caffeinate ]; then
      /usr/bin/caffeinate -is "$PY" -u run.py run --out "$RUN" --workers "$WORKERS"
    else
      "$PY" -u run.py run --out "$RUN" --workers "$WORKERS"
    fi
    "$PY" -u run.py aggregate --out "$RUN"
    echo "Additional simulations complete: $RUN/summary"
    ;;
  status)
    "$PY" -u run.py status --out "$RUN"
    ;;
  aggregate)
    "$PY" -u run.py aggregate --out "$RUN"
    ;;
  *) echo "Unknown action: $ACTION" >&2; exit 2;;
esac
