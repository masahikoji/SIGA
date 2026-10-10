#!/usr/bin/env python3
"""Single-Mac orchestration. Numerical results are produced only by run.py."""
from __future__ import annotations
import argparse
import contextlib
import datetime as dt
import fcntl
import hashlib
import json
import os
import platform
import shutil
import signal
import subprocess
import sys
import sysconfig
import time
import zipfile
from importlib.metadata import version
from pathlib import Path

CODE=Path(__file__).resolve().parent
DEFAULT_RUN=Path.home()/'SIGA_runs'/'main_d012_20261009'
PINNED={'numpy':'2.3.5','scipy':'1.17.0','numba':'0.65.1','llvmlite':'0.47.0'}


def environment(require_mac=False):
    if sys.version_info[:2]!=(3,13) or sysconfig.get_config_var('Py_GIL_DISABLED'):
        raise RuntimeError('Use standard CPython 3.13, not a free-threaded build.')
    found={k:version(k) for k in PINNED}
    if found!=PINNED: raise RuntimeError(f'Wrong dependencies: {found}; expected {PINNED}. Run setup.')
    if require_mac and (sys.platform!='darwin' or platform.machine()!='arm64'):
        raise RuntimeError('Use native arm64 Python on this Apple Silicon Mac; do not use Rosetta.')
    result=dict(python=platform.python_version(),system=platform.system(),machine=platform.machine(),
                cpu_count=os.cpu_count(),**found)
    print(json.dumps(result,indent=2),flush=True)
    return result


@contextlib.contextmanager
def locked(root):
    lock_dir=root.parent/'.siga_locks';lock_dir.mkdir(parents=True,exist_ok=True)
    name=hashlib.sha256(str(root).encode()).hexdigest()[:20]+'.lock'
    with (lock_dir/name).open('a+') as f:
        try: fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another M3 workflow is using this output. Use status, or wait for it to finish.')
        f.seek(0);f.truncate();f.write(str(os.getpid())+'\n');f.flush()
        try: yield
        finally: fcntl.flock(f,fcntl.LOCK_UN)


def command(root,*arguments):
    logs=root.parent/'logs'/root.name;logs.mkdir(parents=True,exist_ok=True)
    stamp=dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    log=logs/(stamp+'_'+str(arguments[0])+'.log')
    cmd=[sys.executable,'-u',str(CODE/'run.py'),*map(str,arguments)]
    print('\nCOMMAND: '+' '.join(cmd)+'\nLOG: '+str(log),flush=True)
    with log.open('w') as f:
        proc=subprocess.Popen(cmd,cwd=CODE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                              text=True,bufsize=1)
        try:
            for line in proc.stdout:
                print(line,end='',flush=True);f.write(line);f.flush()
            code=proc.wait()
        except KeyboardInterrupt:
            if proc.poll() is None:
                proc.send_signal(signal.SIGINT)
                try: proc.wait(timeout=60)
                except subprocess.TimeoutExpired:
                    proc.terminate();proc.wait()
            raise
    if code: raise RuntimeError(f'Command failed (exit {code}); read {log}. Completed chunks were not deleted.')


def existing_primary(root):
    from siga_validation.experiment import load_plan
    from siga_validation.preflight import audit_plan
    p=load_plan(root)
    if p['profile']!='production' or p['suite']!='main' or p['directions']!='both':
        raise RuntimeError('M3 workflow requires the fixed production main plan with both directions.')
    if (p['configuration']['outer_null'],p['configuration']['outer_power'],
        p['configuration']['reference'],p['configuration']['calibration'])!=(100000,10000,4999,100000):
        raise RuntimeError('Unexpected production budgets. Do not change the frozen primary plan.')
    if len(p['scenarios'])!=32: raise RuntimeError('Expected 32 primary scenarios.')
    report=audit_plan(root)
    print('Preflight passed: 32 primary scenarios; binary amplitude 0.12; no shrinking.',flush=True)
    return p


def prepare(root):
    if not (root/'plan.json').exists():
        command(root,'init','--out',root,'--suites','main','--profile','production')
    existing_primary(root)
    command(root,'calibrate','--out',root,'--workers','2')
    print('\nPrepared. Outcome trials have not started. Next: bash m3.sh benchmark',flush=True)


def worker_count(root,override=None):
    from siga_validation.experiment import load_plan
    p=load_plan(root)
    if override is not None:
        w=override
    elif (root/'recommended_workers.json').exists():
        r=json.loads((root/'recommended_workers.json').read_text())
        if r['plan_hash']!=p['plan_hash']: raise RuntimeError('Worker recommendation belongs to a different plan.')
        w=r['workers']
    else:
        w=min(12,os.cpu_count() or 1)
        print('No local benchmark found; using a conservative default of '+str(w)+' workers.',flush=True)
    if not 1<=w<=(os.cpu_count() or 1): raise ValueError('Workers must be between 1 and the available CPU count.')
    return w


def tune(root,trials,candidates):
    from siga_validation.experiment import atomic_json
    from siga_validation.report import write_csv
    p=existing_primary(root);cpus=os.cpu_count() or 1
    workers=sorted(set(int(x) for x in candidates.split(',') if int(x)<=cpus))
    if not workers or workers[0]<1: raise ValueError('Invalid benchmark worker candidates.')
    stamp=dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f');dest=root/'benchmarks'/stamp;dest.mkdir(parents=True)
    rows=[]
    for w in workers:
        command(root,'benchmark','--out',root,'--trials',trials,'--workers',w)
        data=json.loads((root/'benchmark.json').read_text())
        if data['scenarios_timed']!=32: raise RuntimeError('Incomplete timing coverage.')
        rows.append(dict(workers=w,projected_hours=data['ideal_parallel_hours'],
                         benchmark_seconds=data['actual_benchmark_seconds']))
        for ext in ('json','csv'): shutil.copy2(root/('benchmark.'+ext),dest/f'workers_{w}.{ext}')
    best=min(rows,key=lambda x:x['projected_hours'])
    write_csv(dest/'worker_comparison.csv',rows)
    atomic_json(root/'recommended_workers.json',dict(workers=best['workers'],plan_hash=p['plan_hash'],
         projected_hours=best['projected_hours'],benchmarks=str(dest),
         note='Timing estimate only; excludes full-run I/O, calibration, and future workload.'))
    print('\nMeasured worker comparison:',flush=True)
    for row in rows: print(f"  {row['workers']:2d} workers: projected {row['projected_hours']:.2f} hours",flush=True)
    print(f"Selected {best['workers']} workers from this local benchmark. This is not a runtime guarantee.",flush=True)


def archive(root,destination):
    existing_primary(root)
    # Revalidate actual checkpoint coverage instead of trusting a stale STATUS file.
    command(root,'aggregate','--out',root)
    status_path=root/'summary'/'STATUS.json'
    if not status_path.exists(): raise RuntimeError('Run and aggregate before archiving.')
    status=json.loads(status_path.read_text())
    if not status.get('complete') or not status.get('production_budget'):
        raise RuntimeError('Refusing to archive an incomplete or smoke run as a final result.')
    # A unique local archive is completed before the single file is copied to Dropbox.
    stamp=dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    local=root.parent/'exports';local.mkdir(parents=True,exist_ok=True)
    output=local/f'SIGA_main_d012_{stamp}.zip';temp=output.with_suffix('.zip.partial')
    with zipfile.ZipFile(temp,'w',compression=zipfile.ZIP_STORED,allowZip64=True) as z:
        for file in sorted(root.rglob('*')):
            if file.is_file(): z.write(file,'run/'+file.relative_to(root).as_posix())
        for file in sorted(CODE.rglob('*')):
            if not file.is_file() or any(x.startswith('.') or x=='__pycache__' for x in file.relative_to(CODE).parts): continue
            if file.suffix in ('.py','.sh','.md','.txt','.json','.csv','.patch') or file.name in ('LICENSE','SHA256SUMS'):
                z.write(file,'code/'+file.relative_to(CODE).as_posix())
        logs=root.parent/'logs'/root.name
        if logs.exists():
            for file in sorted(logs.glob('*.log')): z.write(file,'logs/'+file.name)
    os.replace(temp,output)
    destination.mkdir(parents=True,exist_ok=True)
    final=destination/output.name;copy=final.with_suffix('.zip.partial')
    shutil.copyfile(output,copy);os.replace(copy,final)
    with final.open('rb') as stream:
        digest=hashlib.file_digest(stream,'sha256').hexdigest()
    final.with_suffix('.zip.sha256').write_text(digest+'  '+final.name+'\n')
    print('\nFinal archive: '+str(final)+'\nSHA256: '+digest,flush=True)
    print('The local trial checkpoints are retained. No input or output was deleted.',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('env','smoke','prepare','benchmark','run','resume','status','aggregate','archive','precision'))
    p.add_argument('--out',type=Path,default=Path(os.environ.get('SIGA_RUN_DIR',DEFAULT_RUN)))
    p.add_argument('--workers',type=int)
    p.add_argument('--trials',type=int,default=None)
    p.add_argument('--candidates',default='8,12,16,20,24')
    p.add_argument('--destination',type=Path,default=CODE.parent/'SIGA_results')
    p.add_argument('--require-mac',action='store_true')
    a=p.parse_args();root=a.out.expanduser().resolve()
    environment(a.require_mac)
    if a.action=='env': return
    if a.action=='status':
        command(root,'status','--out',root);return
    with locked(root):
        if a.action=='smoke':
            smoke=root.with_name(root.name+'_smoke')
            command(smoke,'test')
            if not (smoke/'plan.json').exists():
                command(smoke,'init','--out',smoke,'--suites','all','--profile','smoke')
            command(smoke,'check','--out',smoke)
            command(smoke,'calibrate','--out',smoke,'--workers','2')
            command(smoke,'run','--out',smoke,'--workers','2')
            command(smoke,'aggregate','--out',smoke)
            print('\nSMOKE completed. These are NOT scientific results. Next: bash m3.sh prepare')
        elif a.action=='prepare': prepare(root)
        elif a.action=='benchmark': tune(root,a.trials or 10,a.candidates)
        elif a.action in ('run','resume'):
            existing_primary(root);w=worker_count(root,a.workers)
            command(root,'run','--out',root,'--workers',w)
            command(root,'aggregate','--out',root)
            print('\nAll requested trials completed and aggregated. Results: '+str(root/'summary'))
            print('To copy a complete archive to the project folder: bash m3.sh archive')
        elif a.action=='aggregate':
            existing_primary(root);command(root,'aggregate','--out',root)
        elif a.action=='precision':
            existing_primary(root);w=worker_count(root,a.workers)
            command(root,'precision','--out',root,'--scenarios','0,8,16,24','--trials',a.trials or 1000,
                    '--reference','9999','--calibration-repeats','3','--workers',w)
        elif a.action=='archive': archive(root,a.destination.expanduser().resolve())


if __name__=='__main__':
    try: main()
    except KeyboardInterrupt:
        print('\nInterrupted. Completed chunks remain; restart with bash m3.sh run.',file=sys.stderr)
        sys.exit(130)
    except Exception as exc:
        print(f'ERROR: {exc}',file=sys.stderr);sys.exit(1)
