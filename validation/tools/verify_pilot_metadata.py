#!/usr/bin/env python3
"""Verify archived original execution plans against the distributed code/results."""
from __future__ import annotations
import argparse, csv, io, json, sys, zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'code'))
from siga_validation.design import stable_id
from siga_validation.experiment import source_signature

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive',type=Path)
    a=p.parse_args();failures=[];records=[]
    with zipfile.ZipFile(a.archive) as z:
        names=set(z.namelist())
        for suite in ('factorial','controls','confirm'):
            root=f'results_{suite}_pilot/'
            name=root+'plan.json'
            if name not in names:failures.append(name+' missing');continue
            plan=json.loads(z.read(name));data=dict(plan);declared=data.pop('plan_hash')
            if stable_id(data)!=declared:failures.append(name+' plan hash mismatch')
            if plan['source_signature']!=source_signature():failures.append(name+' source version mismatch')
            coverage=root+'summary/coverage.csv'
            if coverage not in names:failures.append(coverage+' missing');continue
            rows=list(csv.DictReader(io.StringIO(z.read(coverage).decode())))
            ids={int(r['scenario']):r['id'] for r in rows}
            for sc in plan['scenarios']:
                if ids.get(sc['index'])!=sc['id']:failures.append(name+f" scenario {sc['index']} ID mismatch")
            cfg=plan['configuration']
            if any(int(r['n_trials'])!=cfg['outer'] for r in rows):failures.append(name+' completed count mismatch')
            for key in {s['calibration_key'] for s in plan['scenarios']}:
                fn=root+f'calibration/{key}.npz'
                if fn not in names:failures.append(fn+' missing');continue
                with np.load(io.BytesIO(z.read(fn)),allow_pickle=False) as obj:
                    if str(obj['plan_hash'])!=declared:failures.append(fn+' incompatible plan')
                    if int(obj['B0'])!=cfg['calibration']:failures.append(fn+' calibration count mismatch')
            records.append({'suite':suite,'configuration':cfg,'seed':plan['seed'],
                'plan_hash':declared,'source_signature':plan['source_signature']})
    print(json.dumps({'runs':records,'failures':failures},indent=2))
    if failures:raise SystemExit(1)
    print('PASS: original plan IDs, source version, declared counts and calibration headers.')
    print('This does not rerun simulations or reaggregate raw trial checkpoints.')
if __name__=='__main__':main()
