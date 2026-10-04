#!/usr/bin/env python3
"""Archive existing plans/environments/calibrations without rerunning trials."""
from __future__ import annotations
import argparse, hashlib, json, platform, zipfile
from pathlib import Path

def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,default=Path.home()/'SIGA_runs')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();root=args.runs.expanduser().resolve();out=args.output.expanduser().resolve()
    if not root.is_dir():p.error(f'Run directory does not exist: {root}')
    if out.exists():p.error(f'Refusing to overwrite existing archive: {out}')
    records=[];missing=[];selected=[]
    for suite in ('factorial','controls','confirm'):
        run=root/f'results_{suite}_pilot'
        if not run.is_dir():missing.append(str(run));continue
        for name in ('plan.json','plan.csv','environment.json'):
            q=run/name
            if q.is_file():selected.append(q)
            else:missing.append(str(q))
        for folder,pattern in [('calibration','*.npz'),('summary','*')]:
            found=[x for x in (run/folder).glob(pattern) if x.is_file() and not x.is_symlink()]
            if not found:missing.append(str(run/folder/pattern))
            selected.extend(found)
    if not selected:p.error('No existing metadata or summaries found; nothing was changed.')
    for q in selected:
        if q.is_symlink():p.error(f'Symlink not followed: {q}')
        records.append({'path':str(q.relative_to(root)), 'bytes':q.stat().st_size,
                        'sha256':hashlib.sha256(q.read_bytes()).hexdigest()})
    out.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(out,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for q in selected:z.write(q,q.relative_to(root).as_posix())
        z.writestr('COLLECTION.json',json.dumps({'files':records,'missing':missing,
          'purpose':'Existing execution metadata; no simulations rerun; large trial checkpoints excluded.',
          'collector_platform':platform.platform()},indent=2))
    print(f'Created: {out}\nFiles: {len(selected)}')
    for item in missing:print('MISSING:',item)
    print('Keep large original trial checkpoints separately for possible reaggregation.')
if __name__=='__main__':main()
