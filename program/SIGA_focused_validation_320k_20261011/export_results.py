#!/usr/bin/env python3
"""Export compact results without raw checkpoints or source manuscripts."""
import argparse
from pathlib import Path
import json,zipfile
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',required=True);a=p.parse_args()
r=Path(a.out).expanduser().resolve()
status=json.loads((r/'report/status.json').read_text())
if not status['all_scenarios_complete']:
    raise SystemExit('Run is incomplete. Inspect report/completion.csv; no final archive produced.')
files=[]
for name in ['plan.json','plan.csv','resolved_design.json','work_estimate.json','calibration_execution.json','last_execution.json']:
    if (r/name).is_file():files.append(r/name)
files+=list((r/'report').glob('*'))+list((r/'calibration').glob('*.npz'))
with zipfile.ZipFile(r/'results_for_review.zip','w',zipfile.ZIP_DEFLATED) as z:
    for f in sorted(files):
        if f.is_file():z.write(f,f.relative_to(r))
print(r/'results_for_review.zip')
