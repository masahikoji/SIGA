#!/usr/bin/env python3
"""Run a read-only-upstream, separately recorded SIGA revision study."""
from __future__ import annotations
import os
# Control numerical threads BEFORE importing NumPy/Numba. Workers are processes.
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'vendor'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    init=commands.add_parser('init',help='Freeze a new plan. Refuses to overwrite nonempty directories.')
    init.add_argument('--out',required=True)
    init.add_argument('--profile',choices=('smoke','full'),default='smoke')
    init.add_argument('--suite',choices=('focused','gap','a13'),default='focused')
    init.add_argument('--seed',type=int,default=2026101107)
    for name in ('outer-null','outer-power','diagnostic-trials','reference','calibration','direction-calibration',
                 'calibration-repeats','audit-trials','high-reference','chunk'):
        init.add_argument('--'+name,type=int)
    for name in ('calibrate','run','report','status'):
        sub=commands.add_parser(name)
        sub.add_argument('--out',required=True)
        if name in ('calibrate','run'): sub.add_argument('--workers',type=int,default=4)
        if name == 'run': sub.add_argument('--scenarios',default=None,help='Comma-separated literal ids or glob patterns; quote wildcards.')
        if name=='run': sub.add_argument('--max-chunks',type=int,default=None,help='Optional partial execution; report retains incomplete status.')
    commands.add_parser('test')
    audit=commands.add_parser('audit-upstream',help='Read-only comparison to a local GitHub checkout.')
    audit.add_argument('--repo',required=True)
    audit.add_argument('--save',default=None)
    args=parser.parse_args()
    if args.command=='test':
        import unittest
        suite=unittest.defaultTestLoader.discover(str(ROOT/'tests'))
        result=unittest.TextTestRunner(verbosity=2).run(suite)
        return 0 if result.wasSuccessful() else 1
    if args.command=='audit-upstream':
        from revision.utils import verify_vendor,file_hash,atomic_json
        upstream=verify_vendor();repo=Path(args.repo).expanduser().resolve();rows=[]
        for rel,v in upstream['files'].items():
            p=repo/v['upstream_path']
            rows.append(dict(path=str(p),exists=p.is_file(),
                             matches_reviewed_snapshot=p.is_file() and file_hash(p)==v['sha256']))
        report=dict(reviewed_commit=upstream['commit'],files=rows,read_only=True)
        if args.save: atomic_json(Path(args.save),report)
        print(json.dumps(report,indent=2));return 0 if all(r['matches_reviewed_snapshot'] for r in rows) else 2
    from revision.utils import load_plan,run_lock,atomic_json
    from revision.runner import initialise,execute,run_trials
    from revision.calibration import calibration_jobs,calibrate_job,resolve_and_seal
    root=Path(args.out).expanduser().resolve()
    if args.command=='init':
        keys=('outer_null','outer_power','diagnostic_trials','reference','calibration','direction_calibration',
              'calibration_repeats','audit_trials','high_reference','chunk')
        _,work=initialise(root,args.profile,args.suite,args.seed,{k:getattr(args,k) for k in keys})
        print(json.dumps(work,indent=2));print('Frozen plan:',root/'plan.json');return 0
    if args.command=='status':
        plan=load_plan(root)
        print(json.dumps(dict(profile=plan['profile'],suite=plan['suite'],plan_hash=plan['plan_hash'],
                              scenarios=[s['id'] for s in plan['scenarios']],configuration=plan['configuration']),indent=2))
        status=root/'report/status.json'
        if status.exists(): print(status.read_text())
        return 0
    with run_lock(root,args.command):
        if args.command=='calibrate':
            plan=load_plan(root);jobs=calibration_jobs(plan,root);rows=[]
            for i,r in enumerate(execute(jobs,calibrate_job,args.workers),1):
                rows.append(r);print(f'Calibration {i}/{len(jobs)}: {r}',flush=True)
            resolve_and_seal(root,plan);atomic_json(root/'calibration_execution.json',rows)
        elif args.command=='run':
            run_trials(root,args.workers,args.scenarios,args.max_chunks)
        elif args.command=='report':
            from revision.reporting import make_report
            print(json.dumps(make_report(root),indent=2))
    return 0


if __name__=='__main__':
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print('\nInterrupted. Completed checkpoints are retained; rerun the same command to resume.',file=sys.stderr)
        raise SystemExit(130)
