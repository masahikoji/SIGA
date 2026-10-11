"""New follow-up designs; never relabel them as the archived original trials.

Only the data generator and explicitly labelled oracle diagnostics see a true
outcome parameter. Feasible analysis functions accept Observation, b, and an
allocation-only calibration, not a Scenario.
"""
from __future__ import annotations
import copy
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import brentq
from scipy.special import expit
from .utils import rng_for, stable_hash
from siga_validation.core import Observation
from siga_validation.allocation import path_from_uniforms

PROFILES = {
    'smoke': dict(outer_null=8, outer_power=6, diagnostic_trials=8, reference=49,
                  calibration=500, direction_calibration=1000, calibration_repeats=3,
                  audit_trials=8, high_reference=99, chunk=4),
    'full': dict(outer_null=100000, outer_power=10000, diagnostic_trials=100000,
                 reference=4999, calibration=100000, direction_calibration=200000,
                 calibration_repeats=3, audit_trials=2000, high_reference=9999, chunk=100),
}


def patterns_prob(F):
    if F not in (2, 5):
        raise ValueError('The frozen follow-up covers F=2 and F=5 only.')
    patterns = ((np.arange(2**F)[:, None] >> np.arange(F)) & 1).astype(np.int64)
    marg = np.array([.5, .4, .3, .2, .1])[:F]
    pi = np.prod(np.where(patterns, marg, 1-marg), axis=1)
    return patterns, pi / pi.sum()


def baseline_means(patterns, pi):
    val = patterns @ np.array([.5, .4, .3, .2, .1])[:patterns.shape[1]]
    val += .25 * patterns[:, 0] * patterns[:, 1]
    return val - float(pi @ val)


def baseline_risks(patterns, pi):
    eta = patterns @ np.array([.35, -.25, .20, -.15, .10])[:patterns.shape[1]]
    eta += .20 * patterns[:, 0] * patterns[:, 1]
    a = brentq(lambda t: float(pi @ expit(t+eta))-.60, -30, 30, xtol=1e-13)
    return expit(a+eta)


def first_direction(patterns, pi):
    v = patterns[:, 0].astype(float)
    v -= float(pi @ v)
    return v / np.max(np.abs(v))


def maximum_ratio_direction(psi, pi):
    """Population-mean-zero constraint, with a deterministic sign convention."""
    basis = null_space(np.sqrt(pi)[None, :])
    white = psi / np.sqrt(pi[:, None] * pi[None, :])
    R = basis.T @ white @ basis
    eig, vectors = np.linalg.eigh((R+R.T)/2)
    v = (basis @ vectors[:, -1]) / np.sqrt(pi)
    v -= float(pi @ v)
    if v[np.argmax(np.abs(v))] < 0:
        v = -v
    v /= np.max(np.abs(v))
    return v, dict(minimum_ratio=float(eig[0]), maximum_ratio=float(eig[-1]),
                   eigengap=float(eig[-1]-eig[-2]) if len(eig)>1 else None,
                   note='Eigenvalues of an independently estimated finite-n covariance, not known limiting eigenvalues.')


def design_key(F, n, pbc):
    return f'F{F}_n{n}_p{int(round(100*pbc))}'


def make_scenarios(suite='focused'):
    """Four continuous large-gap settings and one 100,000-trial A13 evaluation.

    New follow-up simulations, not a reproduction of the archived random streams.
    No external comparator and no new CRT variant are evaluated.
    The inherited, fixed local alternative h=1.5 is not an 80% power target.
    """
    if suite not in ('focused', 'gap', 'a13'):
        raise ValueError(suite)
    rows=[]
    if suite in ('focused','gap'):
        for n,hist in ((200,'S09'),(1000,'S11')):
            for role in ('null','power'):
                b=-0.20
                delta=b if role=='null' else b+1.5/np.sqrt(n)
                rows.append(dict(id=f'G_C_n{n}_{role}',family='gap',F=2,n=n,pbc=0.95,
                    outcome='continuous',model='symmetric',boundary=b,delta=float(delta),
                    direction='maximum_ratio',amplitude=1.0,sigma=0.05,eta=0.0,
                    tail='upper',alpha=0.025,null_evaluation=(role=='null'),
                    historical_counterpart=hist if role=='null' else None,audit=False,
                    calibration_key=design_key(2,n,0.95)))
    if suite in ('focused','a13'):
        rows.append(dict(id='A13_diagnostic',family='diagnostic',F=5,n=400,pbc=0.95,
            outcome='continuous',model='prognostic',boundary=0.0,delta=0.0,direction='zero',
            amplitude=0.0,sigma=0.25,eta=0.25,tail='two',alpha=0.05,
            null_evaluation=True,historical_counterpart='A13',audit=True,
            calibration_key=design_key(5,400,0.95)))
    return rows


def resolve_scenario(sc, directions):
    s=copy.deepcopy(sc)
    pat,pi=patterns_prob(s['F'])
    if s['direction']=='maximum_ratio':
        v=np.asarray(directions[s['calibration_key']]['direction'])
    elif s['direction']=='first_factor':
        v=first_direction(pat,pi)
    elif s['direction']=='zero':
        v=np.zeros(len(pi))
    else:
        raise ValueError(s['direction'])
    amplitude=s['amplitude']
    if s['outcome']=='binary' and s['model']=='symmetric':
        # Original H.1 feasibility recipe is based on b, not on estimated data.
        b=s['boundary']
        upper=np.where(v>0,(.90-b)/np.maximum(v,1e-300),np.inf)
        lower=np.where(v<0,(.90+b)/np.maximum(-v,1e-300),np.inf)
        tmax=float(min(np.min(upper),np.min(lower)))
        amplitude=min(.80,.98*tmax)
    dev=amplitude*v
    if abs(pi@dev)>1e-12:
        raise AssertionError('Effect vector not probability centred.')
    tau=s['delta']+dev
    if s['model']=='symmetric':
        if s['outcome']=='continuous':
            mu0=-tau/2; mu1=tau/2
            var0=np.full(len(pi),s['sigma']**2); var1=var0.copy()
        else:
            mu0=(1-tau)/2; mu1=(1+tau)/2
            var0=mu0*(1-mu0); var1=mu1*(1-mu1)
    elif s['model']=='prognostic':
        if s['outcome']=='continuous':
            mu0=baseline_means(pat,pi); mu1=mu0+tau
            var0=np.full(len(pi),s['sigma']**2)
            var1=np.full(len(pi),s['sigma']**2+s['eta']**2)
        else:
            mu0=baseline_risks(pat,pi); mu1=mu0+tau
            var0=mu0*(1-mu0); var1=mu1*(1-mu1)
    else:
        raise ValueError(s['model'])
    if s['outcome']=='binary' and not (np.all((mu0>.02)&(mu0<.98)) and np.all((mu1>.02)&(mu1<.98))):
        raise ValueError(f'Invalid binary probabilities in {s["id"]}; no clipping or outcome-driven shrinking is allowed.')
    if abs(pi@(mu1-mu0)-s['delta'])>1e-12:
        raise AssertionError('Generating marginal effect changed.')
    s.update(profile_prob=pi.tolist(),dgp_deviations=dev.tolist(),actual_amplitude=float(amplitude),
             mu0=mu0.tolist(),mu1=mu1.tolist(),var0=var0.tolist(),var1=var1.tolist())
    s['resolved_hash']=stable_hash(s)
    return s


def generate_observation(sc, seed, trial):
    pat,pi=patterns_prob(sc['F']); n=sc['n']; key=sc['id']
    # Historical counterparts share distributions, NOT old random-number streams.
    sid=np.searchsorted(np.cumsum(pi),rng_for(seed,'outer',key,int(trial),'profiles').random(n)).astype(np.int64)
    u=rng_for(seed,'outer',key,int(trial),'assignment').random(n)
    z=path_from_uniforms(sid,pat,sc['pbc'],u)
    A=(z.astype(np.int64)+1)//2
    ry=rng_for(seed,'outer',key,int(trial),'outcome')
    m0=np.asarray(sc['mu0']); m1=np.asarray(sc['mu1'])
    if sc['outcome']=='binary':
        y=(ry.random(n)<np.where(A==1,m1[sid],m0[sid])).astype(float)
    elif sc['model']=='prognostic':
        base=m0[sid]+sc['sigma']*ry.normal(size=n)
        noise=sc['eta']*ry.normal(size=n)
        y=base+A*(m1[sid]-m0[sid]+noise)
    else:
        y=np.where(A==1,m1[sid],m0[sid])+sc['sigma']*ry.normal(size=n)
    return Observation(y.astype(float),A,sid,pat)
