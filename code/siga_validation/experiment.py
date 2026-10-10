"""Frozen-plan calibration, shared-trial evaluation, and checkpointed output."""
from __future__ import annotations
import json, os, time, hashlib, tempfile, platform
from importlib.metadata import version
from pathlib import Path
import numpy as np
from .allocation import calibration_batch, reference_batch
from .core import (contrast_matrix, psd_roundoff, scores_and_variances,
                   mc_pvalue, normal_pvalue, lattice_pvalue)
from .design import rng_for, patterns_and_prob, generate_observation, stable_id

METHODS = ['RT', 'original_S_normal', 'refined_S_normal',
           'original_R_normal', 'refined_R_normal', 'original_CRT_data', 'refined_CRT_data',
           'original_R_lattice', 'refined_R_lattice',
           'original_CRT_lattice', 'refined_CRT_lattice',
           'original_CRT_noise', 'refined_CRT_noise']
DIAGNOSTICS = ['T', 'ref_variance', 'S_original', 'S_refined', 'R_original', 'R_refined',
               'ref_between', 'ref_within', 'ref_twice_cross',
               'pred_between_original', 'pred_between_refined',
               'pred_within_original', 'pred_within_refined',
               'lambda_original', 'lambda_refined', 'both_arms', 'represented',
               'kappa_floor_original', 'kappa_floor_refined', 'R_floor_original', 'R_floor_refined']
_CAL_CACHE = {}


def source_signature():
    root=Path(__file__).resolve().parent.parent
    h=hashlib.sha256()
    for p in sorted(root.glob('siga_validation/*.py'))+[root/'run.py',root/'mac.py']:
        h.update(p.relative_to(root).as_posix().encode());h.update(p.read_bytes())
    return h.hexdigest()


def atomic_json(path, data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix=path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as f: json.dump(data,f,indent=2,sort_keys=True,allow_nan=False)
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)


def atomic_npz(path, **arrays):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix=path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'wb') as f: np.savez_compressed(f,**arrays)
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)


def numerical_environment():
    """Freeze numerical software and interpreter/architecture for a single run."""
    return dict(python=platform.python_version(),system=platform.system(),
                machine=platform.machine(),
                **{name:version(name) for name in ('numpy','scipy','numba','llvmlite')})


def load_plan(root):
    root=Path(root)
    plan=json.loads((root/'plan.json').read_text())
    check=dict(plan);given=check.pop('plan_hash')
    if stable_id(check)!=given: raise RuntimeError('Plan was edited after freezing. Use a new output directory.')
    if plan['source_signature'] != source_signature():
        raise RuntimeError('Analysis source changed after freezing. Use a new output directory.')
    expected=plan.get('numerical_environment')
    if expected is not None and expected != numerical_environment():
        raise RuntimeError('Numerical environment changed after init. Restore the recorded '
                           'Python/library versions and architecture; do not mix checkpoints.')
    return plan


def covariance(sum_x, sum_xx, B):
    return psd_roundoff((sum_xx-np.outer(sum_x,sum_x)/B)/(B-1))


def calibrate_one(args):
    root,key=args;root=Path(root);plan=load_plan(root)
    path=root/'calibration'/f'{key}.npz'
    if path.exists():
        with np.load(path,allow_pickle=False) as old:
            if str(old['plan_hash']) != plan['plan_hash']: raise RuntimeError('Incompatible calibration.')
        return f'{key}: existing calibration verified'
    sc=next(s for s in plan['scenarios'] if s['calibration_key']==key)
    n,F=sc['n'],sc['F'];pat,pi=patterns_and_prob(F,sc['profile_law']);J=len(pi)
    true_pi=pi.copy()
    if plan['profile_source']=='external_estimated':
        size=max(5000,10*n)
        cnt=rng_for(plan['seed'],'external_profile_sample',key).multinomial(size,pi)
        pi=(cnt+1)/(size+J)
    elif plan['profile_source']!='known_design': raise ValueError('Unknown profile source.')
    L=contrast_matrix(pat);q=F+1;B=plan['configuration']['calibration']
    cal_seed=plan.get('calibration_seed',plan['seed'])
    rp=rng_for(cal_seed,'calibration',key,'profiles')
    ra=rng_for(cal_seed,'calibration',key,'assignments')
    sums=[np.zeros(J),np.zeros(q),np.zeros(J),np.zeros(J)]
    seconds=[np.zeros((J,J)),np.zeros((q,q)),np.zeros((J,J)),np.zeros((J,J))]
    weights=np.zeros(n+1)
    started=time.perf_counter()
    for begin in range(0,B,256):
        m=min(256,B-begin)
        sid=np.searchsorted(np.cumsum(pi),rp.random((m,n))).astype(np.int64)
        u,D,g1,g2,treated=calibration_batch(sid,pat,sc['pbc'],ra.random((m,3,n)))
        # X is obtained from the unnormalised first-copy D, not from Gamma.
        arrays=[u,D@L.T,g1,g2]
        for k,x in enumerate(arrays):
            sums[k]+=x.sum(axis=0);seconds[k]+=x.T@x
        weights+=np.bincount(treated,minlength=n+1)
    covs=[covariance(s,x,B) for s,x in zip(sums,seconds)]
    elapsed=time.perf_counter()-started
    atomic_npz(path,gamma=covs[0],xi=covs[1],psi=.5*(covs[2]+covs[3]),
               weights=weights/B,planning_pi=pi,generating_pi=true_pi,
               B0=np.array(B),n=np.array(n),pbc=np.array(sc['pbc']),
               plan_hash=np.array(plan['plan_hash']),elapsed_seconds=np.array(elapsed))
    return f'{key}: {B} allocation-only replicates, {elapsed:.1f} s'


def calibration_for(root,key):
    path=Path(root)/'calibration'/f'{key}.npz'
    ckey=str(path.resolve())
    if ckey not in _CAL_CACHE:
        with np.load(path,allow_pickle=False) as z: _CAL_CACHE[ckey]={k:z[k] for k in z.files}
    return _CAL_CACHE[ckey]


def trial_evaluation(sc, plan, trial, cal):
    obs=generate_observation(sc,plan['seed'],trial)
    clock=time.perf_counter()
    entries,geo=scores_and_variances(obs,sc['boundaries'],cal,noise_adjustment=True)
    construction_seconds=time.perf_counter()-clock
    scores=np.column_stack([e['score'] for e in entries])
    e_scores=np.column_stack([e['original']['S']['residual'] for e in entries])
    both=np.column_stack((scores,e_scores))
    L=len(entries);B=plan['configuration']['reference']
    ref=np.empty((B,2*L))
    rr=rng_for(plan['seed'],'reference',sc['id'],int(trial))
    clock=time.perf_counter()
    for begin in range(0,B,256):
        m=min(256,B-begin)
        ref[begin:begin+m]=reference_batch(obs.stratum,obs.patterns,sc['pbc'],both,rr.random((m,sc['n'])))
    reference_seconds=time.perf_counter()-clock
    conditional=ref[:,:L];residual_ref=ref[:,L:];stratum_ref=conditional-residual_ref
    p_boundary=np.ones((L,len(METHODS)))
    diagnostic=np.full((L,len(DIAGNOSTICS)),np.nan)
    discrete_used=np.zeros(2,dtype=np.int8)
    clock=time.perf_counter()
    for l,e in enumerate(entries):
        tail=sc['tails'][l//2];tr=conditional[:,l];t=e['statistic']
        p_boundary[l,0]=mc_pvalue(tr,t,tail)
        diagnostic[l,0]=t;diagnostic[l,1]=tr.var(ddof=1)
        vb=stratum_ref[:,l].var(ddof=1);ve=residual_ref[:,l].var(ddof=1)
        diagnostic[l,6:9]=[vb,ve,diagnostic[l,1]-vb-ve]
        diagnostic[l,15:17]=[e['both_arms'],int((e['counts']>0).sum())]
        applicable=sc['outcome']=='binary' and e['boundary']==0 and e['score_name']=='unadjusted' and tail=='two'
        if applicable: discrete_used[l%2]=1
        for version,tag in enumerate(('original','refined')):
            sv=e[tag]['S'];rv=e[tag]['R'];vS=sv['variance'];vR=rv['variance']
            lam=np.sqrt(vR/vS) if vS>0 and vR>0 else 1.0
            pnS=normal_pvalue(t,vS,tail);pnR=normal_pvalue(t,vR,tail)
            pcrt=mc_pvalue(tr,t,tail,lam,True)
            # Lattice-labelled columns mean lattice where applicable, normal
            # elsewhere. This rule is frozen identically for both constructions.
            if applicable:
                yplus=int(obs.y.sum())
                prd=lattice_pvalue(t,vR,sc['n'],yplus,cal['weights'])
                pcd=lattice_pvalue(t,vR,sc['n'],yplus,cal['weights'],lam,True)
            else: prd,pcd=pnR,pnS
            noise=e[tag]['R_noise']['variance']
            nl=np.sqrt(noise/vS) if vS>0 and noise>0 else 1.0
            pnoise=mc_pvalue(tr,t,tail,nl,True)
            for offset,value in [(1,pnS),(3,pnR),(5,pcrt),(7,prd),(9,pcd),(11,pnoise)]:
                p_boundary[l,offset+version]=value
            diagnostic[l,2+version]=vS;diagnostic[l,4+version]=vR
            diagnostic[l,9+version]=sv['between'];diagnostic[l,11+version]=sv['within']+rv['correction']
            diagnostic[l,13+version]=lam
            diagnostic[l,17+version]=int(sv['kappa_active']);diagnostic[l,19+version]=int(rv['floor'])
    tail_seconds=time.perf_counter()-clock
    # Entries alternate unadjusted/adjusted within a boundary.
    # TOST uses max of both one-sided p-values, not its boundary component.
    p=np.stack([p_boundary[j::2,:].max(axis=0) for j in (0,1)])
    reject=(p <= sc['alpha']).astype(np.int8)
    arm0=np.bincount(obs.stratum[obs.assignment==0],minlength=len(obs.patterns))
    arm1=np.bincount(obs.stratum[obs.assignment==1],minlength=len(obs.patterns))
    missing=(arm0==0)|(arm1==0)
    missing_mass=float(np.asarray(cal['planning_pi'])[missing].sum())
    return dict(p=p,reject=reject,diagnostic=diagnostic,
                fallback=int(geo['fallback']),balance_error=geo['balance_error'],
                discrete=discrete_used,missing_strata=int(missing.sum()),missing_profile_mass=missing_mass,
                seconds=np.array([construction_seconds,reference_seconds,tail_seconds]))


def run_chunk(args):
    root,index,start,stop=args;root=Path(root);plan=load_plan(root)
    sc=plan['scenarios'][index]
    path=root/'trials'/sc['id']/f'{start:07d}_{stop:07d}.npz'
    if path.exists():
        with np.load(path,allow_pickle=False) as saved:
            if str(saved['plan_hash'])!=plan['plan_hash'] or not np.array_equal(saved['trial_index'],np.arange(start,stop)):
                raise RuntimeError('Incompatible existing chunk.')
        return f'{index}:{start}-{stop}: existing chunk verified'
    cal=calibration_for(root,sc['calibration_key'])
    if str(cal['plan_hash'])!=plan['plan_hash']: raise RuntimeError('Calibration/plan mismatch.')
    # Exceptions stop this chunk; no failed trials are silently removed.
    values=[trial_evaluation(sc,plan,j,cal) for j in range(start,stop)]
    out={k:np.stack([v[k] for v in values]) for k in values[0]}
    out.update(trial_index=np.arange(start,stop),plan_hash=np.array(plan['plan_hash']),scenario_id=np.array(sc['id']))
    atomic_npz(path,**out)
    return f'{index}:{start}-{stop}: saved {len(values)} paired trials'
