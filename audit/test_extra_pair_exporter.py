from pathlib import Path
import sys, tempfile, importlib.util, json, hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent.parent
spec=importlib.util.spec_from_file_location('exporter',ROOT/'review/export_extra_pairs.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
# Synthetic input tests only, not production simulation results.
with tempfile.TemporaryDirectory() as td:
    root=Path(td);code=root/'code';(code/'siga_validation').mkdir(parents=True)
    (code/'siga_validation/experiment.py').write_text('METHODS='+repr(m.METHODS)+'\n')
    for f in ['run.py','mac.py']:(code/f).write_text('# fixture\n')
    h=hashlib.sha256()
    for p in sorted(code.glob('siga_validation/*.py'))+[code/'run.py',code/'mac.py']:
        h.update(p.relative_to(code).as_posix().encode());h.update(p.read_bytes())
    sc=dict(index=0,id='fixture',outer=6,group='test',F=2,n=200,outcome='continuous',pbc=.8,direction='first_factor',role='superiority_null',delta=0.,alpha=.05,evaluation='type1')
    plan=dict(format_version=3,source_signature=h.hexdigest(),scenarios=[sc]);plan['plan_hash']=m.stable_id(plan)
    run=root/'run';run.mkdir();(run/'plan.json').write_text(json.dumps(plan));tr=run/'trials/fixture';tr.mkdir(parents=True)
    rng=np.random.default_rng(164);p=rng.uniform(0,.1,(6,2,13));r=(p<=.05).astype('int8')
    def save(path,ind,pr,re):np.savez(path,p=pr,reject=re,trial_index=ind,plan_hash=plan['plan_hash'],scenario_id='fixture')
    save(tr/'a.npz',np.arange(3),p[:3],r[:3]);save(tr/'b.npz',np.arange(3,6),p[3:],r[3:])
    rows,_,_=m.collect(run,code);assert len(rows)==12
    for row in rows:
        si=0 if row['score']=='unadjusted' else 1
        a=m.METHODS.index(row['method']);b=m.METHODS.index(row['reference']);d=r[:,si,a].astype(float)-r[:,si,b]
        assert np.isclose(row['paired_rejection_difference'],d.mean())
        assert np.isclose(row['paired_mcse'],d.std(ddof=1)/np.sqrt(6))
    save(tr/'duplicate.npz',np.arange(3),p[:3],r[:3])
    try:m.collect(run,code);raise AssertionError('Duplicate not detected')
    except ValueError as e:assert 'Duplicate' in str(e)
    (tr/'duplicate.npz').unlink();(tr/'b.npz').unlink()
    try:m.collect(run,code);raise AssertionError('Incomplete not detected')
    except ValueError as e:assert 'Incomplete' in str(e)
    bad=r[:3].copy();bad[0,0,0]=1-bad[0,0,0];save(tr/'a.npz',np.arange(3),p[:3],bad)
    try:m.collect(run,code);raise AssertionError('Mismatch not detected')
    except ValueError as e:assert 'mismatch' in str(e)
print('Synthetic exporter checks passed: paired statistics; duplicate, incomplete, and inconsistent decisions rejected.')
