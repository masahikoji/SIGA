#!/usr/bin/env python3
"""Reduced SIGA rerun. Start with `python run.py test` and a smoke plan."""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[_name]=os.environ.get('SIGA_BLAS_THREADS','1')
import argparse
import concurrent.futures as cf
import json
import multiprocessing as mp
import platform
import sys
import time
from pathlib import Path
import numpy as np
import scipy
import numba
from siga_validation.experiment import load_plan,calibrate_one,run_chunk,atomic_json
from siga_validation.report import aggregate,write_csv
from siga_validation.resimulation import (make_resimulation_plan,freeze,work_estimate,
    select_indices,benchmark_one,save_benchmark,prepare_precision,precision_one,SEED)


def parallel_results(fun,jobs,workers):
    if workers<1: raise ValueError('Workers must be positive.')
    if workers==1:
        for job in jobs: yield fun(job)
        return
    # Bound queued jobs: no millions of futures or reference paths in memory.
    it=iter(jobs)
    with cf.ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn')) as pool:
        pending=set()
        for _ in range(2*workers):
            try: pending.add(pool.submit(fun,next(it)))
            except StopIteration: break
        while pending:
            done,pending=cf.wait(pending,return_when=cf.FIRST_COMPLETED)
            for future in done:
                yield future.result()
                try: pending.add(pool.submit(fun,next(it)))
                except StopIteration: pass


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('init')
    q.add_argument('--out',required=True)
    q.add_argument('--suites',default='main',help='Comma-separated main,controls,reconstruction,sample_size,pbc95; or all.')
    q.add_argument('--profile',choices=('production','smoke'),default='production')
    q.add_argument('--directions',choices=('both','first_factor','interaction'),default='both')
    q.add_argument('--seed',type=int,default=SEED)
    q.add_argument('--calibration',type=int,default=None,help='Allocation triples per setting; default 100000 in production, 300 in smoke.')
    q.add_argument('--chunk',type=int,default=100)
    for command in ('calibrate','run','benchmark','precision'):
        q=sub.add_parser(command);q.add_argument('--out',required=True)
        q.add_argument('--workers',type=int,default=1)
        if command in ('run','benchmark','precision'): q.add_argument('--scenarios',default='')
        if command=='run':
            q.add_argument('--shard',default='1/1');q.add_argument('--max-chunks',type=int,default=None)
        if command=='benchmark': q.add_argument('--trials',type=int,default=5)
        if command=='precision':
            q.add_argument('--trials',type=int,default=200)
            q.add_argument('--reference',type=int,default=9999)
            q.add_argument('--calibration-repeats',type=int,default=3,help='Total batches including original; 1 skips additional calibration.')
    q=sub.add_parser('aggregate');q.add_argument('--out',required=True);q.add_argument('--allow-partial',action='store_true')
    q=sub.add_parser('status');q.add_argument('--out',required=True)
    q=sub.add_parser('check');q.add_argument('--out',required=True)
    sub.add_parser('test')
    a=p.parse_args()
    if a.command=='test':
        from siga_validation.tests import run_tests
        from siga_validation.resimulation_tests import run_resimulation_tests
        from siga_validation.d012_tests import run_d012_tests
        run_tests();run_resimulation_tests();run_d012_tests();return
    root=Path(a.out).expanduser().resolve()
    if a.command=='init':
        plan=make_resimulation_plan(a.suites,a.profile,a.directions,a.seed,a.calibration,a.chunk)
        plan=freeze(root,plan)
        atomic_json(root/'environment.json',dict(python=sys.version,platform=platform.platform(),
            numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__))
        print(json.dumps(work_estimate(plan),indent=2));print('Plan hash:',plan['plan_hash'])
        if a.profile=='smoke': print('SMOKE ONLY: these budgets must not be reported as a scientific study.')
        else: print('Production plan frozen; no simulation has run yet. Run calibrate and benchmark before run.')
        return
    plan=load_plan(root)
    if a.command=='check':
        from siga_validation.preflight import audit_plan
        print(json.dumps(audit_plan(root),indent=2));return
    if a.command=='aggregate':
        print(aggregate(root,a.allow_partial));return
    if a.command=='status':
        print(json.dumps(work_estimate(plan),indent=2))
        complete_chunks=0;total_chunks=0
        for sc in plan['scenarios']:
            for start in range(0,sc['outer'],plan['configuration']['chunk']):
                stop=min(start+plan['configuration']['chunk'],sc['outer']);total_chunks+=1
                complete_chunks+=int((root/'trials'/sc['id']/f'{start:07d}_{stop:07d}.npz').exists())
        print(f'Checkpoint files present: {complete_chunks}/{total_chunks}. Aggregate verifies actual coverage and hashes.')
        return
    if a.workers<1: raise ValueError('Workers must be positive.')
    if a.command=='calibrate':
        keys=sorted({s['calibration_key'] for s in plan['scenarios']})
        for result in parallel_results(calibrate_one,[(str(root),k) for k in keys],a.workers): print(result,flush=True)
        return
    requested=select_indices(plan,a.scenarios)
    for i in requested:
        key=plan['scenarios'][i]['calibration_key']
        if not (root/'calibration'/f'{key}.npz').exists():
            raise FileNotFoundError('Run calibrate first; a required calibration is missing.')
    if a.command=='benchmark':
        if a.trials<2: raise ValueError('Use at least 2 timed trials.')
        start=time.perf_counter()
        rows=list(parallel_results(benchmark_one,[(str(root),i,a.trials) for i in requested],a.workers))
        save_benchmark(root,rows,a.workers,time.perf_counter()-start)
        print((root/'benchmark.json').read_text());return
    if a.command=='precision':
        audit,jobs=prepare_precision(root,a.trials,a.reference,a.calibration_repeats,requested)
        rows=[]
        for values in parallel_results(precision_one,jobs,a.workers): rows.extend(values)
        rows.sort(key=lambda r:(r['index'],r['audit'],r['score'],r['method']))
        write_csv(audit/'precision_comparisons.csv',rows)
        print('Precision diagnostics:',audit)
        print('This audit is not the 100000-trial type I error study and has no automatic acceptance threshold.')
        return
    k,K=map(int,a.shard.split('/'))
    if not 1<=k<=K: raise ValueError('Shard must be k/K with 1 <= k <= K.')
    if a.max_chunks is not None and a.max_chunks<1: raise ValueError('max-chunks must be positive.')
    jobs=[];counter=0
    for sc in plan['scenarios']:
        for start in range(0,sc['outer'],plan['configuration']['chunk']):
            stop=min(start+plan['configuration']['chunk'],sc['outer'])
            if sc['index'] in requested and counter%K==k-1: jobs.append((str(root),sc['index'],start,stop))
            counter+=1
    if a.max_chunks is not None: jobs=jobs[:a.max_chunks]
    print(f'Running/verifying {len(jobs)} chunks with {a.workers} workers. Repeating this command resumes checkpoints.',flush=True)
    for result in parallel_results(run_chunk,jobs,a.workers): print(result,flush=True)

if __name__=='__main__':
    try: main()
    except KeyboardInterrupt:
        print('Interrupted. Completed checkpoints remain valid; rerun the same command.',file=sys.stderr)
        sys.exit(130)
