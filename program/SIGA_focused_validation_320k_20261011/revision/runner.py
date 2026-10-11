from __future__ import annotations
import copy
import fnmatch
import json
import multiprocessing as mp
import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from .utils import (ROOT, atomic_json, atomic_npz, environment, load_plan, stable_hash,
                    source_signature, verify_vendor, write_csv)
from .designs import PROFILES, make_scenarios
from .calibration import load_resolved, load_calibration
from .engine import one_trial, METHODS, METRICS, EFFECT_METRICS


def initialise(root,profile,suite,seed,overrides=None):
    root=Path(root).expanduser().resolve()
    if root==ROOT or ROOT in root.parents:
        raise ValueError('Put results outside the source package, e.g. ~/SIGA_runs/revision_20261011.')
    if root.exists() and any(root.iterdir()):
        raise FileExistsError('Use a new empty result directory. Existing data will not be overwritten.')
    upstream=verify_vendor()
    cfg=copy.deepcopy(PROFILES[profile])
    for k,v in (overrides or {}).items():
        if v is not None: cfg[k]=int(v)
    for k in ('outer_null','outer_power','diagnostic_trials','reference','calibration','direction_calibration','chunk'):
        if cfg[k]<2: raise ValueError(k+' must be at least 2.')
    if cfg['audit_trials']<0 or cfg['calibration_repeats']<1:
        raise ValueError('Invalid audit settings.')
    if cfg['high_reference']<cfg['reference']:
        raise ValueError('high_reference must be at least reference.')
    scenarios=make_scenarios(suite)
    for s in scenarios:
        s['outer']=(cfg['diagnostic_trials'] if s['family']=='diagnostic' else
                    cfg['outer_null'] if s['null_evaluation'] else cfg['outer_power'])
    plan=dict(format_version=1,extension_version='focused-1.0.0',seed=int(seed),profile=profile,
              suite=suite,configuration=cfg,scenarios=scenarios,
              upstream_commit=upstream['commit'],source_signature=source_signature(),
              environment=environment(),
              design_context='Follow-up designs developed after manuscript review. Not an original preregistration.',
              principal_effect_source='observed data only',
              oracle_policy='Truth is used only in explicit diagnostic columns and DGP generation.',
              inference_status='Only original RT, SIGA-S, SIGA-R and plug-in CRT are principal. Oracle and the existing noise subtraction are diagnostics only.',
              target='iid population average effect; nominal allocation 1:1; fixed binary factors')
    plan['plan_hash']=stable_hash(plan)
    atomic_json(root/'plan.json',plan)
    write_csv(root/'plan.csv',scenarios)
    from .calibration import calibration_jobs
    jobs=calibration_jobs(plan,root)
    main=sum(s['outer'] for s in scenarios)
    audit=sum(min(s['outer'],cfg['audit_trials']) for s in scenarios if s['audit'])
    work=dict(scenarios=len(scenarios),null_scenarios=sum(s['null_evaluation'] for s in scenarios),
              outcome_trials=main,additional_reference_audit_trials=audit,
              base_reference_paths=main*cfg['reference'],
              extra_reference_paths=audit*(cfg['high_reference']-cfg['reference']),
              calibration_jobs=len(jobs),calibration_triples=sum(j[-1] for j in jobs),
              statement='Scores, methods and calibration sensitivity reuse the same outcomes. They do not multiply the independent trial count.')
    atomic_json(root/'work_estimate.json',work)
    return plan,work


def selected_scenarios(resolved,text=None):
    if not text: return resolved['scenarios']
    patterns=[x.strip() for x in text.split(',') if x.strip()]
    chosen=[s for s in resolved['scenarios'] if any(fnmatch.fnmatchcase(s['id'],p) for p in patterns)]
    if not chosen: raise ValueError('No scenario matches '+text)
    unmatched=[p for p in patterns if not any(fnmatch.fnmatchcase(s['id'],p) for s in resolved['scenarios'])]
    if unmatched: raise ValueError('Unmatched scenario selector(s): '+','.join(unmatched))
    return chosen


def chunk_path(root,id,start,end):
    return Path(root)/'trials'/id/f'{start:08d}_{end:08d}.npz'


def warm_compile():
    from .calibration import simulate_calibration
    simulate_calibration(2,8,.8,2,0,['warmup'])
    from siga_validation.allocation import reference_batch
    reference_batch(np.zeros(8,dtype=np.int64),np.array([[0,0],[1,0],[0,1],[1,1]],dtype=np.int64),
                    .8,np.zeros((8,2)),np.full((2,8),.25))


def _run_job(job):
    root,id,start,end=job
    plan,resolved=load_resolved(root)
    sc=next(s for s in resolved['scenarios'] if s['id']==id)
    path=chunk_path(root,id,start,end)
    meta=dict(plan_hash=plan['plan_hash'],resolved_hash=sc['resolved_hash'],
              id=id,start=start,end=end,methods=METHODS,metrics=METRICS,effect_metrics=EFFECT_METRICS)
    meta=json.loads(json.dumps(meta))
    if path.exists():
        with np.load(path,allow_pickle=False) as f:
            if json.loads(str(f['metadata']))!=meta or not np.array_equal(f['trial'],np.arange(start,end)):
                raise RuntimeError('Incompatible checkpoint '+str(path))
        return dict(id=id,start=start,end=end,cached=True)
    cfg=plan['configuration']
    rep_count=cfg['calibration_repeats'] if sc['audit'] and start<cfg['audit_trials'] else 1
    cals=[load_calibration(root,sc['calibration_key'],i) for i in range(rep_count)]
    rows=[]; audits=[]; audit_indices=[]
    t=time.perf_counter()
    for trial in range(start,end):
        try:
            r=one_trial(sc,plan,trial,cals)
            if not np.isfinite(r['p'][:,:,:4]).all():
                raise FloatingPointError('A feasible SIGA/RT p-value was non-finite.')
            if not np.all((r['p']>=0)&(r['p']<=1)):
                raise FloatingPointError('A p-value was outside [0,1].')
            if r['audit']:
                audits.append(r['audit']); audit_indices.append(trial)
            rows.append(r)
        except Exception:
            atomic_json(path.parent/f'error_{trial:08d}.json',
                        dict(id=id,trial=trial,traceback=traceback.format_exc(),plan_hash=plan['plan_hash']))
            raise
    arrays={key:np.array([r[key] for r in rows]) for key in ('p','metrics','effect','statistic','fallback','seconds')}
    arrays.update(trial=np.arange(start,end),audit_trial=np.array(audit_indices,dtype=np.int64),
                  metadata=np.array(json.dumps(meta,sort_keys=True)),elapsed=np.array(time.perf_counter()-t))
    if audits:
        for key in ('high_p','high_metrics','alt_p','alt_metrics','alt_effect'):
            arrays['audit_'+key]=np.array([a[key] for a in audits])
    atomic_npz(path,**arrays)
    return dict(id=id,start=start,end=end,cached=False,seconds=time.perf_counter()-t)


def execute(jobs,worker_function,workers):
    if workers<1: raise ValueError('workers must be positive.')
    warm_compile()
    if workers==1:
        yield from map(worker_function,jobs)
    else:
        # spawn matches macOS and avoids forked BLAS state.
        with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn')) as pool:
            # Yield completed jobs promptly, including when earlier jobs take longer.
            from concurrent.futures import as_completed
            pending=[pool.submit(worker_function,job) for job in jobs]
            for future in as_completed(pending):
                yield future.result()


def run_trials(root,workers=4,selectors=None,max_chunks=None):
    plan,resolved=load_resolved(root,verify_calibrations=True)
    selected=selected_scenarios(resolved,selectors)
    chunk=plan['configuration']['chunk']
    jobs=[(str(root),s['id'],a,min(a+chunk,s['outer'])) for s in selected for a in range(0,s['outer'],chunk)]
    # A partial run is permitted, but report.py always shows actual completion.
    if max_chunks is not None: jobs=jobs[:max_chunks]
    t=time.perf_counter(); completed=0; cached=0
    requested_trials=sum(j[3]-j[2] for j in jobs); done_trials=0
    print(f'Running {requested_trials:,} requested trials in {len(jobs)} chunks with {workers} workers.',flush=True)
    for result in execute(jobs,_run_job,workers):
        completed+=1; cached+=int(result['cached'])
        done_trials+=result['end']-result['start']
        if completed%5==0 or completed==len(jobs):
            print(f'Trials {done_trials:,}/{requested_trials:,} ({100*done_trials/requested_trials:.1f}%); chunks {completed}/{len(jobs)}; cached={cached}; elapsed={(time.perf_counter()-t)/60:.1f} min; last={result["id"]}',flush=True)
    record=dict(workers=workers,chunks_requested=len(jobs),chunks_completed=completed,
                cached=cached,seconds=time.perf_counter()-t,selectors=selectors,
                source_signature=plan['source_signature'])
    atomic_json(Path(root)/'last_execution.json',record)
    return record
