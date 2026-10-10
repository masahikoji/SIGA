"""Reduced, prospectively frozen rerun plan and execution diagnostics.

Numerical inference and data-generation kernels are vendored from SIGA commit
1685daefb13c8a86602526e8828c321076fd6993. See UPSTREAM.json and docs/CHANGES.md.
"""
from __future__ import annotations
import copy
import csv
import json
import math
import time
from pathlib import Path
import numpy as np
from .design import make_scenario, stable_id, patterns_and_prob
from .experiment import (numerical_environment, source_signature, atomic_json, atomic_npz, load_plan,
                         calibrate_one, calibration_for, trial_evaluation, METHODS)
from .report import write_csv, se

UPSTREAM_COMMIT = '1685daefb13c8a86602526e8828c321076fd6993'
SEED = 2026100901
BINARY_AMPLITUDE = 0.12
SUITES = ('main','controls','reconstruction','sample_size','pbc95')
ROLES = ('superiority_null','superiority_power',
         'noninferiority_null','noninferiority_power')


def make_resimulation_plan(suites=('main',), profile='production', directions='both',
                           seed=SEED, calibration=None, chunk=100):
    if isinstance(suites,str): suites=tuple(suites.split(','))
    if 'all' in suites: suites=SUITES
    if not suites or not set(suites).issubset(SUITES):
        raise ValueError('Suites must be main, controls, reconstruction, sample_size, pbc95, or all.')
    suites=tuple(s for s in SUITES if s in suites)
    if profile not in ('production','smoke'): raise ValueError('Invalid profile.')
    if directions not in ('both','first_factor','interaction'): raise ValueError('Invalid direction.')
    if chunk<1: raise ValueError('Chunk must be positive.')
    if profile=='production':
        cfg=dict(outer_null=100000,outer_power=10000,reference=4999,
                 calibration=100000 if calibration is None else int(calibration),chunk=int(chunk))
    else:
        cfg=dict(outer_null=6,outer_power=4,reference=49,
                 calibration=300 if calibration is None else int(calibration),chunk=4)
    if cfg['calibration']<2: raise ValueError('At least two calibration replicates are required.')
    drs=('first_factor','interaction') if directions=='both' else (directions,)
    candidates=[]
    def add(F,n,outcome,pbc,sd,eta,amp,dr,role,group,sharp=False):
        candidates.append(make_scenario(F,n,outcome,pbc,sd,eta,amp,dr,role,group,sharp=sharp,binary_policy="error"))
    if 'main' in suites:
        for F,n in ((2,200),(5,400)):
            for outcome in ('continuous','binary'):
                for dr in drs:
                    for role in ROLES:
                        add(F,n,outcome,.8,1.,.25,.5 if outcome=='continuous' else BINARY_AMPLITUDE,dr,role,'main')
    if 'controls' in suites:
        for F,n in ((2,200),(5,400)):
            for outcome in ('continuous','binary'):
                add(F,n,outcome,.8,1.,0.,0.,'zero','superiority_null','controls',True)
    if 'reconstruction' in suites:
        # Three stratum-effect patterns, not three individual-effect laws.
        # eta=.25 is retained even when the stratum-effect vector is zero.
        for sd in (.25,1.):
            for amp,dr in ((0.,'zero'),(.5,'first_factor'),(.5,'interaction')):
                for role in ('superiority_null','superiority_power'):
                    add(5,400,'continuous',.8,sd,.25,amp,dr,role,'reconstruction')
    if 'sample_size' in suites:
        for F,n in ((2,1000),(5,2000)):
            for outcome in ('continuous','binary'):
                for role in ('superiority_null','superiority_power'):
                    add(F,n,outcome,.8,1.,.25,.5 if outcome=='continuous' else BINARY_AMPLITUDE,
                        'first_factor',role,'sample_size')
    if 'pbc95' in suites:
        # Same outcome law as the .8 reconstruction study; change only pbc.
        for sd in (.25,1.):
            for amp,dr in ((0.,'zero'),(.5,'first_factor')):
                for role in ('superiority_null','superiority_power'):
                    add(5,400,'continuous',.95,sd,.25,amp,dr,role,'pbc95')
    unique={}
    for s in candidates:
        # Presentation group and replication budget must not create duplicate DGPs.
        physics={k:v for k,v in s.items() if k not in ('id','group','calibration_key')}
        uid=stable_id(physics)[:20]
        if uid in unique:
            unique[uid]['memberships'].append(s['group'])
            continue
        s['id']=uid
        s['memberships']=[s['group']]
        s['evaluation']='type1' if s['null_evaluation'] else 'power'
        s['outer']=cfg['outer_null'] if s['null_evaluation'] else cfg['outer_power']
        unique[uid]=s
    scenarios=list(unique.values())
    for i,s in enumerate(scenarios):
        s['index']=i
        pi=np.array(s['profile_prob']);d=np.array(s['dgp_deviations'])
        if abs(pi@d)>1e-12: raise ValueError('Nonzero mean effect deviation.')
        if s['null_evaluation'] and abs(s['delta']-s['boundaries'][0])>1e-12:
            raise ValueError('Null scenario does not lie on its tested boundary.')
        if s['outcome']=='binary':
            p0=np.array(s['control_risks']);p1=np.array(s['treatment_risks'])
            if not (np.all((p0>0)&(p0<1)) and np.all((p1>.02)&(p1<.98))):
                raise ValueError('Invalid Bernoulli probabilities.')
            if abs(pi@(p1-p0)-s['delta'])>1e-12: raise ValueError('Binary marginal effect changed.')
            expected=0.0 if s['sharp_binary'] else BINARY_AMPLITUDE
            if s['feasibility_scale']!=1.0 or abs(s['actual_amplitude']-expected)>1e-12:
                raise ValueError('Binary amplitude was changed from the frozen d012 specification.')
    return dict(format_version=3,binary_amplitude=BINARY_AMPLITUDE,binary_feasibility="error",suite=','.join(suites),profile=profile,seed=int(seed),
                calibration_seed=int(seed),profile_source='known_design',directions=directions,
                configuration=cfg,upstream_commit=UPSTREAM_COMMIT,
                inference_effect_source='observed_data_only',
                plan_context='Revised design after previous simulations; new random streams. Not represented as an original preregistration.',
                scenarios=scenarios)


def work_estimate(plan):
    cfg=plan['configuration']; ss=plan['scenarios']; B=cfg['reference']
    keys={s['calibration_key']:s for s in ss}
    return dict(scenarios=len(ss),null_scenarios=sum(s['null_evaluation'] for s in ss),
                power_scenarios=sum(not s['null_evaluation'] for s in ss),
                outer_trials=sum(s['outer'] for s in ss),
                reference_paths=sum(s['outer']*B for s in ss),
                reference_participant_assignments=sum(s['outer']*B*s['n'] for s in ss),
                calibration_settings=len(keys),
                calibration_triples=len(keys)*cfg['calibration'],
                calibration_paths=3*len(keys)*cfg['calibration'],
                paired_score_analyses=2,methods=len(METHODS),
                note='Methods and adjusted/unadjusted scores share reference paths; do not multiply path counts by method or score counts.')


def freeze(root,plan):
    root=Path(root)
    if root.exists() and any(root.iterdir()): raise FileExistsError('Use a new empty output directory.')
    plan=copy.deepcopy(plan)
    plan['numerical_environment']=numerical_environment()
    plan['source_signature']=source_signature();plan['plan_hash']=stable_id(plan)
    atomic_json(root/'plan.json',plan)
    columns=['index','id','group','memberships','evaluation','outer','F','n','outcome','pbc',
             'outcome_sd','individual_effect_sd','direction','requested_amplitude',
             'actual_amplitude','feasibility_scale','role','delta','boundaries','tails','alpha','calibration_key']
    write_csv(root/'plan.csv',[{k:s[k] for k in columns} for s in plan['scenarios']])
    atomic_json(root/'work_estimate.json',work_estimate(plan))
    return plan


def select_indices(plan,text=''):
    if not text: return list(range(len(plan['scenarios'])))
    ids=sorted(set(int(x) for x in text.split(',')))
    if not ids or ids[0]<0 or ids[-1]>=len(plan['scenarios']): raise ValueError('Invalid zero-based scenario index.')
    return ids


def benchmark_one(args):
    root,index,trials=args;plan=load_plan(root);sc=plan['scenarios'][index]
    cal=calibration_for(root,sc['calibration_key'])
    if str(cal['plan_hash'])!=plan['plan_hash']: raise RuntimeError('Calibration mismatch.')
    # Independent, nonproduction trial indices; warm compilation in each worker.
    trial_evaluation(sc,plan,1_000_000_000,cal)
    elapsed=[];parts=[]
    for i in range(trials):
        t=time.perf_counter();value=trial_evaluation(sc,plan,1_000_000_001+i,cal)
        elapsed.append(time.perf_counter()-t);parts.append(value['seconds'])
    mean=float(np.mean(elapsed));p=np.mean(parts,axis=0)
    return dict(index=index,group=sc['group'],F=sc['F'],n=sc['n'],outcome=sc['outcome'],
                role=sc['role'],direction=sc['direction'],trials_timed=trials,
                reference_draws=plan['configuration']['reference'],
                mean_seconds_per_trial=mean,median_seconds_per_trial=float(np.median(elapsed)),
                construction_seconds=float(p[0]),reference_seconds=float(p[1]),tail_seconds=float(p[2]),
                outer_trials_planned=sc['outer'],estimated_service_hours=mean*sc['outer']/3600)


def save_benchmark(root,rows,workers,wall_seconds):
    root=Path(root);rows=sorted(rows,key=lambda x:x['index'])
    write_csv(root/'benchmark.csv',rows)
    total=sum(r['estimated_service_hours'] for r in rows)
    atomic_json(root/'benchmark.json',dict(workers=workers,actual_benchmark_seconds=wall_seconds,
        estimated_sum_of_worker_service_hours=total,
        ideal_parallel_hours=total/workers,
        scenarios_timed=len(rows),
        note='Measured on this computer at the requested worker count after warm-up. Ideal parallel time excludes calibration, scheduling, checkpoint I/O and other workload. This is a projection, not a completed full run.'))


def prepare_precision(root,trials,extra_reference,calibration_repeats,scenario_indices):
    root=Path(root);plan=load_plan(root)
    if extra_reference<=plan['configuration']['reference']: raise ValueError('Extra reference budget must exceed the base budget.')
    if trials<2 or calibration_repeats<1: raise ValueError('Use >=2 trials and >=1 calibration repeat.')
    settings=dict(trials=trials,extra_reference=extra_reference,calibration_repeats=calibration_repeats,
                  scenario_indices=scenario_indices,parent_plan=plan['plan_hash'])
    audit=root/'precision'/stable_id(settings)[:16]
    atomic_json(audit/'audit_settings.json',settings)
    alt_roots=[]
    keys=sorted({plan['scenarios'][i]['calibration_key'] for i in scenario_indices})
    for k in range(1,calibration_repeats):
        ar=audit/f'calibration_batch_{k}'
        alt=copy.deepcopy(plan);alt.pop('plan_hash')
        alt['profile']='precision_calibration'
        alt['calibration_seed']=plan['calibration_seed']+1_000_003*k
        alt['plan_hash']=stable_id(alt)
        if (ar/'plan.json').exists():
            if load_plan(ar)['plan_hash']!=alt['plan_hash']: raise RuntimeError('Audit plan mismatch.')
        else: atomic_json(ar/'plan.json',alt)
        for key in keys: print(calibrate_one((str(ar),key)),flush=True)
        alt_roots.append(str(ar))
    return audit,[(str(root),str(audit),i,trials,extra_reference,alt_roots) for i in scenario_indices]


def precision_one(args):
    root,audit,index,trials,extra_reference,alt_roots=args
    plan=load_plan(root);sc=plan['scenarios'][index]
    cal=calibration_for(root,sc['calibration_key'])
    higher=copy.deepcopy(plan);higher['configuration']['reference']=extra_reference
    alt_cals=[calibration_for(ar,sc['calibration_key']) for ar in alt_roots]
    path=Path(audit)/f'scenario_{index:03d}.npz'
    # Deterministic audit prefix, separate output; never substituted for the main run.
    if not path.exists():
        values=[]
        for t in range(trials):
            base=trial_evaluation(sc,plan,t,cal)
            high=trial_evaluation(sc,higher,t,cal)
            alternates=[trial_evaluation(sc,plan,t,c) for c in alt_cals]
            values.append([base,high,*alternates])
        atomic_npz(path,p=np.array([[x['p'] for x in row] for row in values]),
                   reject=np.array([[x['reject'] for x in row] for row in values]),
                   diag=np.array([[x['diagnostic'] for x in row] for row in values]),
                   trial_index=np.arange(trials),parent_plan_hash=np.array(plan['plan_hash']))
    with np.load(path,allow_pickle=False) as f:
        if str(f['parent_plan_hash'])!=plan['plan_hash'] or len(f['trial_index'])!=trials:
            raise RuntimeError('Audit checkpoint mismatch.')
        p=f['p'];r=f['reject'];d=f['diag']
    labels=[f'reference_{extra_reference}',*[f'independent_calibration_{k+1}' for k in range(len(alt_roots))]]
    rows=[]
    for v,label in enumerate(labels,1):
        for sj,score in enumerate(('unadjusted','adjusted')):
            for m,method in enumerate(METHODS):
                diff=r[:,v,sj,m].astype(float)-r[:,0,sj,m].astype(float)
                pdiff=p[:,v,sj,m]-p[:,0,sj,m]
                rows.append(dict(index=index,F=sc['F'],n=sc['n'],outcome=sc['outcome'],role=sc['role'],
                    direction=sc['direction'],audit=label,score=score,method=method,trials=trials,
                    base_rejection_rate=float(r[:,0,sj,m].mean()),audit_rejection_rate=float(r[:,v,sj,m].mean()),
                    paired_difference=float(diff.mean()),paired_mcse=se(diff),
                    discordant_trials=int(np.count_nonzero(diff)),mean_p_difference=float(pdiff.mean()),
                    mean_absolute_p_difference=float(np.abs(pdiff).mean())))
    return rows
