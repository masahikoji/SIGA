#!/usr/bin/env python3
"""See README.md before running manuscript-scale simulations."""
# Keep BLAS threads bounded before importing NumPy. Parallelise independent chunks.
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[_name]=os.environ.get('SIGA_BLAS_THREADS','1')
import argparse,csv,json,platform,sys
from pathlib import Path
import concurrent.futures as cf
import multiprocessing as mp
import numpy as np
import scipy,numba
from siga_validation.design import make_plan,stable_id,MASTER_SEED
from siga_validation.experiment import (source_signature,atomic_json,load_plan,calibrate_one,run_chunk)
from siga_validation.report import aggregate


def workers_map(fun,jobs,workers):
    if workers==1:
        for x in jobs: print(fun(x),flush=True)
    else:
        with cf.ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn')) as pool:
            for answer in pool.map(fun,jobs,chunksize=1): print(answer,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='mode',required=True)
    init=sub.add_parser('init');init.add_argument('--out',required=True)
    init.add_argument('--suite',choices=['smoke','confirm','factorial','controls','correlated','all'],default='confirm')
    init.add_argument('--profile',choices=['smoke','pilot','full'],default='pilot')
    init.add_argument('--seed',type=int,default=MASTER_SEED)
    init.add_argument('--profile-source',choices=['known_design','external_estimated'],default='known_design')
    for command in ('calibrate','run'):
        q=sub.add_parser(command);q.add_argument('--out',required=True);q.add_argument('--workers',type=int,default=1)
        if command=='run':
            q.add_argument('--shard',default='1/1',help='Shard k/K, e.g. 1/4. Run disjoint shards and combine after completion.')
            q.add_argument('--scenarios',default='',help='Optional comma-separated zero-based scenario indices; no silent scenario selection in final aggregation.')
    ag=sub.add_parser('aggregate');ag.add_argument('--out',required=True);ag.add_argument('--allow-partial',action='store_true')
    sub.add_parser('test')
    a=p.parse_args()
    if a.mode=='test':
        from siga_validation.tests import run_tests
        run_tests();return
    root=Path(a.out).expanduser().resolve()
    if a.mode=='init':
        if (root/'plan.json').exists(): raise SystemExit('Plan exists. Resume it or choose a new directory; do not overwrite it.')
        if root.exists() and any(root.iterdir()): raise SystemExit('Use an empty output directory.')
        plan=make_plan(a.suite,a.profile,a.seed,a.profile_source)
        plan['source_signature']=source_signature();plan['plan_hash']=stable_id(plan)
        atomic_json(root/'plan.json',plan)
        columns=['index','id','group','F','n','outcome','pbc','outcome_sd','individual_effect_sd',
                 'direction','requested_amplitude','actual_amplitude','feasibility_scale','role','delta','alpha','calibration_key']
        with (root/'plan.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows({k:s[k] for k in columns} for s in plan['scenarios'])
        atomic_json(root/'environment.json',dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__))
        cfg=plan['configuration'];S=len(plan['scenarios'])
        print(f'Frozen {S} scenarios, {cfg["outer"]} outer trials/scenario, {cfg["reference"]} reference draws/trial.')
        print(f'Up to {S*cfg["outer"]*cfg["reference"]:,} complete inner allocation sequences. Start with smoke testing.')
        print('Plan hash:',plan['plan_hash']);return
    plan=load_plan(root)
    if a.mode=='aggregate': print(aggregate(root,a.allow_partial));return
    if a.workers<1: raise SystemExit('--workers must be positive')
    if a.mode=='calibrate':
        keys=sorted({s['calibration_key'] for s in plan['scenarios']})
        workers_map(calibrate_one,[(str(root),k) for k in keys],a.workers);return
    k,K=map(int,a.shard.split('/'))
    if not 1<=k<=K: raise SystemExit('Shard must be k/K with 1 <= k <= K.')
    requested=set(map(int,a.scenarios.split(','))) if a.scenarios else set(range(len(plan['scenarios'])))
    if not requested.issubset(set(range(len(plan['scenarios'])))): raise SystemExit('Invalid scenario indices.')
    cfg=plan['configuration'];tasks=[];counter=0
    for sc in plan['scenarios']:
        for start in range(0,cfg['outer'],cfg['chunk']):
            stop=min(start+cfg['chunk'],cfg['outer'])
            if sc['index'] in requested and counter%K==k-1:
                tasks.append((str(root),sc['index'],start,stop))
            counter+=1
    workers_map(run_chunk,tasks,a.workers)

if __name__=='__main__':
    main()
