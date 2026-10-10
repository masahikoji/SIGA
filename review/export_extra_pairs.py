#!/usr/bin/env python3
"""Read-only extra paired comparisons for completed SIGA d012 checkpoints.

Does not simulate, import the simulation package, or modify its plan/checkpoints.
Use the actual Mac plan, not reconstructed scenario IDs.
"""
from __future__ import annotations
import argparse
import ast
import csv
import hashlib
import json
import math
from pathlib import Path
import tempfile
import numpy as np

METHODS=['RT','original_S_normal','refined_S_normal','original_R_normal',
         'refined_R_normal','original_CRT_data','refined_CRT_data',
         'original_R_lattice','refined_R_lattice','original_CRT_lattice',
         'refined_CRT_lattice','original_CRT_noise','refined_CRT_noise']
PAIRS=[(f'{v}_{a}',b if b=='RT' else f'{v}_{b}')
       for v in ('original','refined')
       for a,b in [('CRT_data','RT'),('S_normal','RT'),('CRT_noise','CRT_data')]]

def stable_id(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def validate_source(code:Path,expected:str):
    files=sorted(code.glob('siga_validation/*.py'))+[code/'run.py',code/'mac.py']
    h=hashlib.sha256()
    for p in files:
        h.update(p.relative_to(code).as_posix().encode());h.update(p.read_bytes())
    if h.hexdigest()!=expected:raise ValueError('Simulation source differs from the frozen plan.')
    tree=ast.parse((code/'siga_validation/experiment.py').read_text())
    got=None
    for node in tree.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='METHODS' for t in node.targets):
            got=ast.literal_eval(node.value)
    if got!=METHODS:raise ValueError('Checkpoint method order is not the expected 13-method order.')

def collect(run:Path,code:Path):
    plan=json.loads((run/'plan.json').read_text())
    raw=dict(plan);given=raw.pop('plan_hash')
    if stable_id(raw)!=given:raise ValueError('Frozen plan hash mismatch.')
    validate_source(code,plan['source_signature'])
    if plan.get('format_version')!=3:raise ValueError('Only the verified d012 format_version=3 is supported.')
    rows=[];coverage=[]
    for sc in plan['scenarios']:
        n=int(sc['outer']);seen=np.zeros(n,dtype=bool)
        stat={(s,a,b):[0,0,0,0,0.,0.] for s in (0,1) for a,b in PAIRS}
        paths=sorted((run/'trials'/sc['id']).glob('*.npz'))
        for path in paths:
            with np.load(path,allow_pickle=False) as z:
                if str(z['plan_hash'])!=given or str(z['scenario_id'])!=sc['id']:
                    raise ValueError(f'Checkpoint hash/id mismatch: {path}')
                ind=np.asarray(z['trial_index'])
                if ind.ndim!=1 or ind.dtype.kind not in 'iu' or np.any(ind<0) or np.any(ind>=n):
                    raise ValueError(f'Invalid trial indices: {path}')
                if len(np.unique(ind))!=len(ind) or seen[ind].any():raise ValueError(f'Duplicate trial indices: {path}')
                r=np.asarray(z['reject']);p=np.asarray(z['p'])
                if r.shape!=(len(ind),2,13) or p.shape!=r.shape:
                    raise ValueError(f'Wrong array shape: {path}')
                if not np.isfinite(p).all() or np.any((p<0)|(p>1)) or not np.isin(r,[0,1]).all():
                    raise ValueError(f'Invalid p-values/decisions: {path}')
                if not np.array_equal(r,p<=sc['alpha']):raise ValueError(f'p-value/decision mismatch: {path}')
                seen[ind]=True
                for (s,a,b),v in stat.items():
                    ra=r[:,s,METHODS.index(a)].astype(bool);rb=r[:,s,METHODS.index(b)].astype(bool)
                    v[0]+=int(ra.sum());v[1]+=int(rb.sum())
                    v[2]+=int((ra&~rb).sum());v[3]+=int((~ra&rb).sum())
                    d=p[:,s,METHODS.index(a)]-p[:,s,METHODS.index(b)]
                    v[4]+=float(d.sum());v[5]+=float(np.abs(d).sum())
        if not seen.all():raise ValueError(f'Incomplete scenario {sc["index"]}: {seen.sum()}/{n}; nothing is published.')
        coverage.append({'index':sc['index'],'id':sc['id'],'n_trials':n})
        meta={k:sc[k] for k in ('index','id','group','F','n','outcome','pbc','direction','role','delta','alpha','evaluation')}
        for (s,a,b),(ka,kb,n10,n01,dp,adp) in stat.items():
            delta=(n10-n01)/n
            mcse=math.sqrt(max(0.,(n10+n01-n*delta*delta)/(n*(n-1)))) if n>1 else float('nan')
            rows.append(dict(**meta,score=('unadjusted','adjusted')[s],method=a,reference=b,
                n_trials=n,rejections_a=ka,rejections_b=kb,rate_a=ka/n,rate_b=kb/n,
                paired_rejection_difference=delta,paired_mcse=mcse,a_only=n10,b_only=n01,
                mean_pvalue_difference=dp/n,mean_absolute_pvalue_error=adp/n,plan_hash=given))
        print(f'Scenario {sc["index"]}: checked {n} saved trials',flush=True)
    return rows,coverage,plan

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',required=True,type=Path)
    parser.add_argument('--code',required=True,type=Path,help='Unchanged simulation-code folder for source signature check')
    parser.add_argument('--output',required=True,type=Path,help='New output directory outside the simulation code and run')
    args=parser.parse_args();run=args.run.expanduser().resolve();code=args.code.expanduser().resolve();out=args.output.expanduser().resolve()
    if out.exists():parser.error('Use a new output folder; existing output is never overwritten.')
    if out.is_relative_to(run) or out.is_relative_to(code):parser.error('Output must be outside the run and code folders.')
    rows,coverage,plan=collect(run,code)
    out.mkdir(parents=True,exist_ok=False)
    with (out/'extra_paired_comparisons.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'STATUS.json').write_text(json.dumps(dict(complete=True,reran_simulations=False,
        source_plan_hash=plan['plan_hash'],scenarios=coverage),indent=2)+'\n')
    (out/'README.txt').write_text('Rates and differences are probabilities; multiply by 100 for percentage points.\n'
        'Paired Monte Carlo standard errors condition on the original shared covariance estimates.\n'
        'All original checkpoint files were read only. No new simulated trial was generated.\n'
        'Zero empirical discordance is not a proof of identical procedures.\n')
    print(out)
if __name__=='__main__':main()
