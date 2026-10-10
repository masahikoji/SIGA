"""Tests specific to the reduced plan; no simulation performance claims."""
from __future__ import annotations
import copy
import itertools
import tempfile
from pathlib import Path
import numpy as np
from .allocation import reference_batch
from .core import Observation, scores_and_variances, effect_deviations, mc_pvalue
from .design import patterns_and_prob, rng_for, generate_observation
from .resimulation import make_resimulation_plan,work_estimate,freeze
from .experiment import load_plan


def test_reduced_plan():
    p=make_resimulation_plan()
    w=work_estimate(p)
    assert (w['scenarios'],w['null_scenarios'],w['power_scenarios'])==(32,16,16)
    assert w['outer_trials']==1760000
    assert w['reference_paths']==8798240000
    assert w['reference_participant_assignments']==2639472000000
    assert len({s['id'] for s in p['scenarios']})==32
    for s in p['scenarios']:
        assert (s['F'],s['n']) in ((2,200),(5,400)) and s['pbc']==.8
        assert not s['role'].startswith('equivalence')
        assert s['outer']==(100000 if s['null_evaluation'] else 10000)
        pi=np.array(s['profile_prob']);d=np.array(s['dgp_deviations'])
        assert abs(pi@d)<1e-12
        if s['outcome']=='binary':
            assert abs(pi@(np.array(s['treatment_risks'])-s['control_risks'])-s['delta'])<1e-12
    allp=make_resimulation_plan('all')
    assert len(allp['scenarios'])==60
    assert len(make_resimulation_plan('main,reconstruction')['scenarios'])==40
    assert len(make_resimulation_plan('main,controls')['scenarios'])==36
    assert len(make_resimulation_plan(directions='first_factor')['scenarios'])==16
    main_ids={s['id'] for s in p['scenarios']}
    assert main_ids.issubset({s['id'] for s in allp['scenarios']})
    # The stress study changes p only, not sd/eta/heterogeneity simultaneously.
    pair=make_resimulation_plan('reconstruction,pbc95')['scenarios']
    for s in pair:
        if s['pbc']==.95:
            matches=[t for t in pair if t['pbc']==.8 and all(t[k]==s[k] for k in
                ('F','n','outcome','outcome_sd','individual_effect_sd','direction','actual_amplitude','role','delta'))]
            assert len(matches)==1


def test_preserved_dgp_choices():
    p=make_resimulation_plan()
    for s in p['scenarios']:
        if s['outcome']=='continuous':
            assert s['outcome_sd']==1. and s['individual_effect_sd']==.25
            assert s['requested_amplitude']==.5
            expected={'superiority_null':0.,'noninferiority_null':-.2,
                'superiority_power':.43 if s['F']==2 else .3,
                'noninferiority_power':.23 if s['F']==2 else .1}
        else:
            assert s['requested_amplitude']==.12
            assert abs(s['actual_amplitude']-.12)<1e-12
            assert s['feasibility_scale']==1.
            expected={'superiority_null':0.,'noninferiority_null':-.1,
                'superiority_power':.2 if s['F']==2 else .14,
                'noninferiority_power':.1 if s['F']==2 else .05}
        assert s['delta']==expected[s['role']]
    for s in make_resimulation_plan('controls')['scenarios']:
        assert s['delta']==0 and s['actual_amplitude']==0
        assert s['individual_effect_sd']==0
        if s['outcome']=='binary': assert s['sharp_binary']


def test_plan_integrity():
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)/'run'
        p=freeze(root,make_resimulation_plan(profile='smoke'))
        assert load_plan(root)['plan_hash']==p['plan_hash']
        text=(root/'plan.json').read_text()
        (root/'plan.json').write_text(text.replace('"outer_null": 6','"outer_null": 7'))
        try: load_plan(root)
        except RuntimeError: pass
        else: raise AssertionError('Edited plan was not detected.')


def test_reference_prefix_and_plus_one():
    sc=make_resimulation_plan(profile='smoke')['scenarios'][0]
    obs=generate_observation(sc,12345,1)
    scores=np.column_stack((obs.y-obs.y.mean(),obs.y*0))
    r1=rng_for(4321,'reference',sc['id'],1).random((49,len(obs.y)))
    r2=rng_for(4321,'reference',sc['id'],1).random((99,len(obs.y)))
    a=reference_batch(obs.stratum,obs.patterns,.8,scores,r1)
    b=reference_batch(obs.stratum,obs.patterns,.8,scores,r2)
    np.testing.assert_array_equal(a,b[:49])
    assert mc_pvalue(np.zeros(4999),1.,'upper')==1/5000
    assert mc_pvalue(np.zeros(4999),0.,'upper')==1.


def test_squared_covariance_three_participants():
    # Independent exact enumeration: verify the finite-n ratio used in Appendix E.
    F=3;n=3;pat,pi=patterns_and_prob(F,'uniform');v=np.prod(2*pat-1,axis=1)
    for p in (.8,.95):
        psi=np.zeros((8,8))
        for sidtuple in itertools.product(range(8),repeat=n):
            sid=np.array(sidtuple);C=np.zeros((n,n));total=0.
            for ztuple in itertools.product((-1,1),repeat=n):
                marginal=np.zeros((F,2),int);overall=0;prob=1.
                for i,z in enumerate(ztuple):
                    plus=abs(overall+1);minus=abs(overall-1)
                    for f,x in enumerate(pat[sid[i]]):
                        plus+=abs(marginal[f,x]+1);minus+=abs(marginal[f,x]-1)
                    pr=p if plus<minus else (1-p if minus<plus else .5)
                    prob*=pr if z==1 else 1-pr
                    overall+=z
                    for f,x in enumerate(pat[sid[i]]): marginal[f,x]+=z
                z=np.array(ztuple);C+=prob*np.outer(z,z);total+=prob
            assert abs(total-1)<1e-12
            H=np.eye(8)[sid]
            psi+=H.T@(C*C)@H/(n*8**n)
        ratio=float(v@psi@v/(v@(np.diag(pi))@v))
        np.testing.assert_allclose(ratio,1-p*(1-p)*(2*p-1)**2/6,atol=1e-12)


def run_resimulation_tests():
    for f in (test_reduced_plan,test_preserved_dgp_choices,test_plan_integrity,
              test_reference_prefix_and_plus_one,test_squared_covariance_three_participants):
        f();print('PASS',f.__name__,flush=True)
