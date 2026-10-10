#!/bin/bash
# Run this with: bash m3.sh setup | smoke | prepare | benchmark | run | status | archive
set -euo pipefail
CODE="$(cd "$(dirname "$0")" && pwd -P)"
cd "$CODE"
ACTION="${1:-help}"
if [ "$ACTION" = "help" ]; then
    echo 'Usage: bash m3.sh setup|smoke|prepare|benchmark|run|resume|status|aggregate|archive|precision'
    echo 'Optional: bash m3.sh run --workers 16'
    echo 'Default results: $HOME/SIGA_runs/main_d012_20261009 (outside Dropbox)'
    exit 0
fi
# Avoid accidentally importing packages from a Conda base or an old project.
unset PYTHONHOME PYTHONPATH
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMBA_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 SIGA_BLAS_THREADS=1
export NUMBA_CACHE_DIR="$HOME/.cache/siga_m3_d012/numba"
export PYTHONPYCACHEPREFIX="$HOME/.cache/siga_m3_d012/python"
mkdir -p "$NUMBA_CACHE_DIR" "$PYTHONPYCACHEPREFIX"
VENV="${SIGA_VENV:-$HOME/.venvs/siga_m3_d012}"
PY="$VENV/bin/python"
if [ "$ACTION" = "setup" ]; then
    if [ -n "${SIGA_PYTHON:-}" ]; then
        BOOT="$SIGA_PYTHON"
    elif command -v python3.13 >/dev/null 2>&1; then
        BOOT="$(command -v python3.13)"
    elif [ -x /Library/Frameworks/Python.framework/Versions/3.13/bin/python3.13 ]; then
        BOOT=/Library/Frameworks/Python.framework/Versions/3.13/bin/python3.13
    else
        echo 'Python 3.13 was not found. Install the standard macOS universal2 installer from python.org, then retry.' >&2
        exit 1
    fi
    "$BOOT" -c 'import sys,platform,sysconfig; assert sys.version_info[:2]==(3,13), "Use Python 3.13"; assert sys.platform=="darwin" and platform.machine()=="arm64", "Use native arm64 Python, not Rosetta"; assert not sysconfig.get_config_var("Py_GIL_DISABLED"), "Use standard, not free-threaded, CPython"'
    if [ ! -x "$PY" ]; then
        mkdir -p "$(dirname "$VENV")"
        "$BOOT" -m venv "$VENV"
    fi
    "$PY" -m pip install --only-binary=:all: -r "$CODE/requirements.txt"
    "$PY" -m pip check
    "$PY" mac.py env --require-mac
    "$PY" -u run.py test
    echo 'Setup complete. Next: bash m3.sh smoke'
    exit 0
fi
if [ ! -x "$PY" ]; then
    echo 'The project environment is missing. Run: bash m3.sh setup' >&2
    exit 1
fi
# Prevent idle/system sleep while this foreground workflow is running on AC power.
if [ -x /usr/bin/caffeinate ] && [ "$ACTION" != "status" ]; then
    exec /usr/bin/caffeinate -is "$PY" -u mac.py "$@" --require-mac
else
    exec "$PY" -u mac.py "$@" --require-mac
fi
