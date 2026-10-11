#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
MODE="${1:-smoke}"
case "$MODE" in
  smoke) OUT="${SIGA_FOCUS_OUT:-$HOME/SIGA_runs/jrssb_focused_smoke_320k_20261011}"; WORKERS="${SIGA_WORKERS:-2}"; SEED=2026101191 ;;
  full) OUT="${SIGA_FOCUS_OUT:-$HOME/SIGA_runs/jrssb_focused_full_320k_20261011}"; WORKERS="${SIGA_WORKERS:-12}"; SEED=2026101107 ;;
  *) echo "Usage: bash m3_focus.sh smoke|full" >&2; exit 2 ;;
esac
PY="${SIGA_PYTHON:-$HOME/.venvs/siga_m3_d012/bin/python}"
if [ ! -x "$PY" ]; then
  echo "Python environment not found: $PY" >&2
  echo "Set SIGA_PYTHON to the Python executable in your existing SIGA environment." >&2
  exit 2
fi
"$PY" -c 'import numpy, scipy, numba; print("Python dependencies: OK")'
"$PY" run_focus.py test
if [ ! -f "$OUT/plan.json" ]; then
  "$PY" run_focus.py init --out "$OUT" --profile "$MODE" --suite focused --seed "$SEED"
fi
"$PY" - "$OUT" "$MODE" "$SEED" <<'PY'
import json,sys
from pathlib import Path
p=json.loads((Path(sys.argv[1])/'plan.json').read_text())
if p['profile'] != sys.argv[2] or p['suite'] != 'focused' or p['seed'] != int(sys.argv[3]):
    raise SystemExit('Existing plan differs from requested run. Choose a new SIGA_FOCUS_OUT.')
print('Output:',sys.argv[1]);print('Profile:',p['profile'],'; scenarios:',len(p['scenarios']))
PY
run_awake() {
  if command -v caffeinate >/dev/null 2>&1; then
    caffeinate -i "$@"
  else
    "$@"
  fi
}
run_awake "$PY" run_focus.py calibrate --out "$OUT" --workers "$WORKERS"
run_awake "$PY" run_focus.py run --out "$OUT" --workers "$WORKERS"
"$PY" run_focus.py report --out "$OUT"
"$PY" export_results.py --out "$OUT"
echo "Complete: $OUT/report"
echo "Upload: $OUT/results_for_review.zip"
