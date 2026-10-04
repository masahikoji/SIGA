"""Algebra, independent rule enumeration, and data-only interface tests."""
from __future__ import annotations
import inspect,itertools
import numpy as np
from .allocation import path_from_uniforms,calibration_batch,reference_batch
from .core import (covariance_constructions,contrast_matrix,scores_and_variances,
                   effect_deviations,lattice_pvalue,mc_pvalue,Observation)
from .design import make_plan,patterns_and_prob,generate_observation,rng_for


def explicit_range_path(sid,patterns,p,uniforms):
    counts=np.zeros((patterns.shape[1],2),int);overall=0;z=[]
    for s,u in zip(sid,uniforms):
        plus=abs(overall+1);minus=abs(overall-1)
        for f,x in enumerate(patterns[s]):
            plus+=abs(counts[f,x]+1);minus+=abs(counts[f,x]-1)
        pp=p if plus<minus else (1-p if minus<plus else .5)
        zi=1 if u<pp else -1;z.append(zi);overall+=zi
        for f,x in enumerate(patterns[s]): counts[f,x]+=zi
    return np.array(z)


def test_allocation():
    rng=np.random.default_rng(71401)
    for F in (1,2,3,5):
        pat,_=patterns_and_prob(F)
        for p in (.5,.8,.95,1.):
            sid=rng.integers(len(pat),size=61);u=rng.random(61)
            np.testing.assert_array_equal(path_from_uniforms(sid,pat,p,u),explicit_range_path(sid,pat,p,u))
    pat,_=patterns_and_prob(2);sids=rng.integers(4,size=(10,7));u=rng.random((10,3,7))
    U,D,g1,g2,k=calibration_batch(sids,pat,.8,u)
    for i in range(10):
        z=[explicit_range_path(sids[i],pat,.8,u[i,j]) for j in range(3)]
        d=np.bincount(sids[i],weights=z[0],minlength=4);cnt=np.bincount(sids[i],minlength=4)
        np.testing.assert_allclose(D[i],d)
        np.testing.assert_allclose(U[i],d/np.sqrt(np.maximum(cnt,1)))
        np.testing.assert_allclose(g1[i],np.bincount(sids[i],weights=z[0]*z[1],minlength=4)/np.sqrt(7))
    scores=rng.normal(size=(7,3));unif=rng.random((5,7))
    tr=reference_batch(sids[0],pat,.8,scores,unif)
    for i in range(5):
        np.testing.assert_allclose(tr[i],explicit_range_path(sids[0],pat,.8,unif[i])@scores/2)


def test_geometry():
    rng=np.random.default_rng(4102)
    for F in (2,5):
        pat,pi=patterns_and_prob(F);J=len(pi);L=contrast_matrix(pat)
        x=rng.normal(size=(J,J));gamma=x@x.T/J
        x=rng.normal(size=(F+1,F+1));xi=x@x.T
        for empty in (False,True):
            counts=rng.integers(3,30,size=J).astype(float)
            if empty: counts[0]=0
            old,new,d=covariance_constructions(counts,gamma,xi,L)
            assert not d['fallback']
            np.testing.assert_allclose(L@new@L.T,xi,atol=5e-10)
            assert np.linalg.eigvalsh(new).min()>-1e-10
            recode=np.eye(F+1)+0.05*rng.normal(size=(F+1,F+1))
            _, recoded, _=covariance_constructions(counts,gamma,recode@xi@recode.T,recode@L)
            np.testing.assert_allclose(new,recoded,atol=5e-10)
            active=counts>0;B=np.sqrt(counts)[:,None]*L.T;ki=np.linalg.inv(B.T@B)
            trace=np.trace(np.diag(active.astype(float))@gamma)-np.trace(ki@B.T@gamma@B)+np.trace(xi@ki)
            np.testing.assert_allclose(trace,np.sum(np.diag(new)[active]/counts[active]),atol=1e-10)
        counts=np.zeros(J);counts[0]=100
        old,new,d=covariance_constructions(counts,gamma,xi,L)
        assert d['fallback'];np.testing.assert_array_equal(old,new)
        sigma=np.eye(J)-L.T@np.linalg.inv(L@L.T)@L
        gamma0=sigma/np.sqrt(pi)[:,None]/np.sqrt(pi)[None,:]
        errors=[]
        for n in (1000,10000,100000):
            counts=n*pi
            _,new,_=covariance_constructions(counts,gamma0+np.eye(J)/np.sqrt(n),xi,L)
            errors.append(np.linalg.norm(new/n-sigma))
        assert errors[2]<errors[1]<errors[0]


def test_no_oracle():
    sc=make_plan('smoke','smoke')['scenarios'][0]
    obs=generate_observation(sc,7140,3);J=len(obs.patterns);L=contrast_matrix(obs.patterns)
    cal={'gamma':np.eye(J),'xi':np.eye(L.shape[0]),'psi':np.eye(J)/J}
    e1,_=scores_and_variances(obs,[0.],cal)
    a=.713
    translated=Observation(obs.y+a*obs.assignment,obs.assignment,obs.stratum,obs.patterns)
    e2,_=scores_and_variances(translated,[a],cal)
    for x,y in zip(e1,e2):
        np.testing.assert_allclose(x['score'],y['score'],atol=1e-12)
        np.testing.assert_allclose(x['dhat'],y['dhat'],atol=1e-12)
        for k in ('original','refined'):
            np.testing.assert_allclose(x[k]['S']['variance'],y[k]['S']['variance'])
            np.testing.assert_allclose(x[k]['R']['variance'],y[k]['R']['variance'])
    assert list(inspect.signature(scores_and_variances).parameters)==['obs','boundaries','calibration','noise_adjustment']
    # Observed data schema has no potential outcome, true effect or true variance.
    assert set(Observation.__dataclass_fields__)=={'y','assignment','stratum','patterns'}
    again=generate_observation(sc,7140,3)
    np.testing.assert_array_equal(obs.y,again.y)
    assert not np.array_equal(rng_for(11,'a',1).random(8),rng_for(11,'b',1).random(8))


def test_lattice():
    w=np.zeros(201);w[100]=1
    assert abs(lattice_pvalue(0.,8.,200,101,w)-1)<1e-12
    for t in (0.,.5,12.5):
        p=lattice_pvalue(t,8.,200,101,w)
        q=lattice_pvalue(t,8.,200,101,w,1.,True)
        np.testing.assert_allclose(p,q,atol=1e-14)
    samples=np.array([-12.5,12.5,12.51,13.,-1.])
    t=12.5;lam=1.001
    no=mc_pvalue(samples,t,'two',lam,False)
    yes=mc_pvalue(samples,t,'two',lam,True)
    assert yes>no
    # A non-tied point can change classification even when lambda is close to 1.
    assert mc_pvalue(np.array([12.51]),t,'two') != mc_pvalue(np.array([12.51]),t,'two',lam,True)


def test_plans():
    counts={}
    for suite in ('confirm','factorial','controls','correlated','all'):
        p=make_plan(suite,'smoke');counts[suite]=len(p['scenarios'])
        for sc in p['scenarios']:
            pi=np.array(sc['profile_prob']);d=np.array(sc['dgp_deviations'])
            assert abs(pi@d)<1e-12
            assert abs(pi@(d+sc['delta'])-sc['delta'])<1e-12
            if sc['outcome']=='binary':
                p0=np.array(sc['control_risks']);p1=np.array(sc['treatment_risks'])
                assert np.min(p1)>=.02-1e-12 and np.max(p1)<=.98+1e-12
                assert abs(pi@(p1-p0)-sc['delta'])<1e-12
            if sc['role'].startswith('equivalence'):
                assert len(sc['boundaries'])==2 and sc['tails']==['upper','lower']
            if sc['null_evaluation']:
                assert min(abs(sc['delta']-b) for b in sc['boundaries'])<1e-12
    assert counts['confirm']==224
    assert counts['all']==sum(counts[x] for x in ('confirm','factorial','controls','correlated'))
    print('Scenario counts:',counts)


def run_tests():
    for fun in [test_allocation,test_geometry,test_no_oracle,test_lattice,test_plans]:
        fun();print('PASS',fun.__name__,flush=True)
