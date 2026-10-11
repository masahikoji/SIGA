from __future__ import annotations
import json
import time
from pathlib import Path
import numpy as np
from .utils import rng_for, atomic_json, atomic_npz, load_plan, file_hash
from .designs import patterns_prob, maximum_ratio_direction, resolve_scenario
from siga_validation.allocation import calibration_batch
from siga_validation.core import contrast_matrix, psd_roundoff


class CovAccumulator:
    """Stable mergeable sample covariance; mean subtraction is NOT omitted."""
    def __init__(self, dimension):
        self.n=0
        self.mean=np.zeros(dimension)
        self.m2=np.zeros((dimension,dimension))

    def add(self, x):
        if len(x)==0: return
        mu=x.mean(axis=0)
        c=x-mu
        n2=len(x)
        delta=mu-self.mean
        self.m2 += c.T@c + np.outer(delta,delta)*(self.n*n2/(self.n+n2))
        self.mean += delta*n2/(self.n+n2)
        self.n += n2

    def covariance(self):
        if self.n<2: raise ValueError('At least two independent triples are needed.')
        return self.m2/(self.n-1)


def simulate_calibration(F,n,pbc,replicates,seed,stream_key,batch=128):
    pat,pi=patterns_prob(F)
    J=len(pi)
    acc=[CovAccumulator(J) for _ in range(4)]
    rp=rng_for(seed,'calibration',stream_key,'profiles')
    ra=[rng_for(seed,'calibration',stream_key,'signs',a) for a in range(3)]
    treated=np.zeros(n+1,dtype=np.int64)
    for first in range(0,replicates,batch):
        m=min(batch,replicates-first)
        sids=np.searchsorted(np.cumsum(pi),rp.random((m,n))).astype(np.int64)
        u=np.stack([r.random((m,n)) for r in ra],axis=1)
        U,D,G1,G2,nt=calibration_batch(sids,pat,pbc,u)
        for a,x in zip(acc,(U,D,G1,G2)): a.add(x)
        treated += np.bincount(nt,minlength=n+1)
    L=contrast_matrix(pat)
    cd=acc[1].covariance()
    return dict(gamma=psd_roundoff(acc[0].covariance()),
                xi=psd_roundoff(L@cd@L.T),
                psi=psd_roundoff((acc[2].covariance()+acc[3].covariance())/2),
                sigma_D_finite=psd_roundoff(cd/n),
                treated_weights=treated/replicates,replicates=np.array(replicates),
                profile_prob=pi)


def calibration_jobs(plan,root):
    cfg=plan['configuration']
    specs={s['calibration_key']:s for s in plan['scenarios']}
    audit_keys={s['calibration_key'] for s in plan['scenarios'] if s['audit'] and cfg['audit_trials']>0}
    direction_keys={s['calibration_key'] for s in plan['scenarios'] if s['direction']=='maximum_ratio'}
    jobs=[]
    for key,s in sorted(specs.items()):
        for r in range(cfg['calibration_repeats'] if key in audit_keys else 1):
            jobs.append((str(root),key,s['F'],s['n'],s['pbc'],'analysis',r,cfg['calibration']))
        if key in direction_keys:
            jobs.append((str(root),key,s['F'],s['n'],s['pbc'],'direction',0,cfg['direction_calibration']))
    return jobs


def cal_path(root,key,stage='analysis',rep=0):
    return Path(root)/'calibration'/f'{key}_{stage}_{rep:02d}.npz'


def calibrate_job(job):
    root,key,F,n,pbc,stage,rep,R=job
    plan=load_plan(root)
    path=cal_path(root,key,stage,rep)
    meta=dict(plan_hash=plan['plan_hash'],stage=stage,rep=rep,key=key,replicates=R)
    if path.exists():
        with np.load(path,allow_pickle=False) as f:
            if json.loads(str(f['metadata']))!=meta:
                raise RuntimeError(f'Calibration checkpoint mismatch: {path}')
        return dict(key=key,stage=stage,rep=rep,cached=True)
    t=time.perf_counter()
    cal=simulate_calibration(F,n,pbc,R,plan['seed'],[key,stage,rep])
    elapsed=time.perf_counter()-t
    atomic_npz(path,**cal,metadata=np.array(json.dumps(meta,sort_keys=True)),seconds=np.array(elapsed))
    return dict(key=key,stage=stage,rep=rep,cached=False,seconds=elapsed)


def load_calibration(root,key,rep=0,stage='analysis'):
    with np.load(cal_path(root,key,stage,rep),allow_pickle=False) as f:
        return {k:f[k].copy() for k in ('gamma','xi','psi','sigma_D_finite','treated_weights','profile_prob')}


def resolve_and_seal(root,plan):
    directions={}
    for s in plan['scenarios']:
        key=s['calibration_key']
        if s['direction']=='maximum_ratio' and key not in directions:
            c=load_calibration(root,key,stage='direction')
            v,diagnostic=maximum_ratio_direction(c['psi'],c['profile_prob'])
            directions[key]=dict(direction=v.tolist(),**diagnostic)
    resolved=[resolve_scenario(s,directions) for s in plan['scenarios']]
    files={str(cal_path(root,j[1],j[5],j[6]).relative_to(root)):
           file_hash(cal_path(root,j[1],j[5],j[6])) for j in calibration_jobs(plan,root)}
    payload=dict(plan_hash=plan['plan_hash'],directions=directions,scenarios=resolved,calibration_files=files)
    path=Path(root)/'resolved_design.json'
    if path.exists():
        if json.loads(path.read_text())!=payload:
            raise RuntimeError('Calibrations or resolved designs changed after sealing.')
    else:
        atomic_json(path,payload)
    return payload


def load_resolved(root,verify_calibrations=False):
    plan=load_plan(root)
    path=Path(root)/'resolved_design.json'
    resolved=json.loads(path.read_text())
    if resolved['plan_hash']!=plan['plan_hash']:
        raise RuntimeError('Resolved design does not match frozen plan.')
    if verify_calibrations:
        for rel,h in resolved['calibration_files'].items():
            if file_hash(Path(root)/rel)!=h:
                raise RuntimeError('Calibration file changed: '+rel)
    return plan,resolved
