"""Outcome-generating models and frozen scenario families, not analysis code."""
from __future__ import annotations
import hashlib
import json
import numpy as np
from scipy.special import expit
from scipy.optimize import brentq
from .core import Observation
from .allocation import path_from_uniforms

MASTER_SEED = 2026100319
DESIGNS = [(2,200), (2,1000), (5,400), (5,2000)]
PROFILES = dict(smoke=dict(outer=8, reference=49, calibration=300, chunk=8),
                pilot=dict(outer=2000, reference=999, calibration=20000, chunk=100),
                full=dict(outer=100000, reference=4999, calibration=100000, chunk=250))


def stable_id(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def rng_for(master_seed, *keys):
    """Schedule/shard-independent PCG64DXSM stream; keys include stage/trial."""
    digest = hashlib.sha256(json.dumps([int(master_seed), *keys], separators=(',', ':')).encode()).digest()
    entropy = np.frombuffer(digest, dtype='<u4').tolist()
    return np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence(entropy)))


def patterns_and_prob(F, law='independent'):
    pat = ((np.arange(2**F)[:,None] >> np.arange(F)) & 1).astype(np.int64)
    marginal = np.array([.5,.4,.3,.2,.1])[:F]
    pi = np.prod(np.where(pat, marginal, 1-marginal), axis=1)
    if law == 'correlated':
        # The correlated family already considered in the supplement.
        u = 2*pat-1
        mu = 2*marginal-1
        pi *= 1+.1*(u[:,0]-mu[0])*(u[:,1]-mu[1])
    elif law == 'uniform':
        pi[:] = 1/len(pi)
    elif law != 'independent': raise ValueError(law)
    pi /= pi.sum()
    return pat, pi


def baseline_means(patterns, pi):
    val = patterns@np.array([.5,.4,.3,.2,.1])[:patterns.shape[1]]
    val = val+.25*patterns[:,0]*patterns[:,1]
    return val-float(pi@val)


def baseline_risks(patterns, pi):
    eta = patterns@np.array([.35,-.25,.20,-.15,.10])[:patterns.shape[1]]
    eta += .20*patterns[:,0]*patterns[:,1]
    intercept = brentq(lambda a: float(pi@expit(a+eta))-.60, -30,30, xtol=1e-13)
    return expit(intercept+eta)


def direction(patterns, pi, name):
    if name == 'zero': return np.zeros(len(pi))
    if name == 'first_factor': val = patterns[:,0].astype(float)
    elif name == 'interaction': val = np.prod(2*patterns-1,axis=1).astype(float)
    else: raise ValueError("Only fixed, prespecified effect directions are allowed.")
    val -= float(pi@val)
    return val/np.max(np.abs(val))


def make_scenario(F,n,outcome,pbc,sd,eta,amplitude,dir_name,role,group,
                  profile_law='independent', sharp=False):
    idx = DESIGNS.index((F,n))
    if outcome == 'continuous':
        margin = .45; ni = -.20
        super_alt = [.430,.190,.300,.135][idx]
        ni_alt = [.230,-.010,.100,-.065][idx]
        eq_alt = [.000,.280,.180,.330][idx]
    else:
        margin = [.21,.10,.15,.10][idx]; ni = -.10
        super_alt = [.20,.10,.14,.07][idx]
        ni_alt = [.10,.00,.05,-.03][idx]
        eq_alt = [.00,.00,.00,.04][idx]
    roles = dict(superiority_null=(0., [0.], ['two'], .05),
                 superiority_power=(super_alt, [0.], ['two'], .05),
                 noninferiority_null=(ni, [ni], ['upper'], .025),
                 noninferiority_power=(ni_alt, [ni], ['upper'], .025),
                 equivalence_lower=(-margin, [-margin,margin], ['upper','lower'], .05),
                 equivalence_upper=(margin, [-margin,margin], ['upper','lower'], .05),
                 equivalence_power=(eq_alt, [-margin,margin], ['upper','lower'], .05))
    delta,boundaries,tails,alpha = roles[role]
    pat,pi = patterns_and_prob(F, profile_law)
    d = amplitude*direction(pat,pi,dir_name)
    scale = 1.0
    if outcome == 'binary':
        p0 = baseline_risks(pat,pi)
        # Feasibility is imposed before any trial is simulated. It does not use
        # a rejection rate, a covariance eigenvector, or the candidate method.
        if np.min(p0+delta) <= .02 or np.max(p0+delta) >= .98:
            raise ValueError("Invalid binary mean effect even without heterogeneity.")
        # Preserve the author's R3 feasibility convention: multiply the whole
        # centred direction by 0.9 until all treatment risks are in (.02,.98).
        # This is done when freezing the plan, before outcomes or test results.
        while not np.all((p0+delta+scale*d > .02) & (p0+delta+scale*d < .98)):
            scale *= .9
            if scale < 1e-3:
                raise ValueError("Binary heterogeneity cannot be made feasible.")
        d = d*scale
        p1 = p0+delta+d
    else:
        p0 = p1 = None
    actual = float(pi@(delta+d))
    if abs(actual-delta)>1e-12: raise AssertionError("Marginal effect not preserved.")
    sc = dict(F=F,n=n,outcome=outcome,pbc=pbc,outcome_sd=sd,
              individual_effect_sd=eta,requested_amplitude=amplitude,
              actual_amplitude=float(np.max(np.abs(d))),direction=dir_name,
              profile_law=profile_law,group=group,role=role,delta=delta,
              boundaries=boundaries,tails=tails,alpha=alpha,
              null_evaluation=role.endswith('_null') or role.startswith('equivalence_') and role!='equivalence_power',
              sharp_binary=bool(sharp and outcome=='binary'),
              feasibility_scale=scale,profile_prob=pi.tolist(),dgp_deviations=d.tolist(),
              control_risks=p0.tolist() if p0 is not None else None,
              treatment_risks=p1.tolist() if p1 is not None else None)
    sc['calibration_key'] = stable_id(dict(F=F,n=n,pbc=pbc,law=profile_law))[:20]
    sc['id'] = stable_id(sc)[:20]
    return sc


def make_plan(suite, profile, seed=MASTER_SEED, profile_source='known_design'):
    if profile not in PROFILES: raise ValueError(profile)
    roles = ['superiority_null','superiority_power','noninferiority_null',
             'noninferiority_power','equivalence_lower','equivalence_upper','equivalence_power']
    scenarios = []
    if suite in ('confirm','all'):
        for F,n in DESIGNS:
            for setting,p,sd,eta in [('moderate',.80,1.,.25),('strong',.95,.25,.10)]:
                for outcome in ('continuous','binary'):
                    amp = .15 if outcome=='binary' else (.5 if setting=='moderate' else 1.)
                    for dr in ('first_factor','interaction'):
                        for role in roles:
                            scenarios.append(make_scenario(F,n,outcome,p,sd,eta,amp,dr,role,'confirm_'+setting))
    if suite in ('factorial','all'):
        # Decouple the three changes that were combined in the earlier stress run.
        for p in (.80,.95):
            for sd in (.25,1.):
                for amp,dr in [(0.,'zero'),(.5,'first_factor'),(1.,'first_factor'),
                               (.5,'interaction'),(1.,'interaction')]:
                    for role in ('superiority_null','superiority_power'):
                        scenarios.append(make_scenario(5,400,'continuous',p,sd,.25,amp,dr,role,'factorial'))
    if suite in ('controls','all'):
        for F,n in DESIGNS:
            for p in (.80,.95):
                for outcome in ('continuous','binary'):
                    scenarios.append(make_scenario(F,n,outcome,p,1.,0.,0.,'zero','superiority_null','sharp_control',sharp=True))
    if suite in ('correlated','all'):
        for F,n in DESIGNS:
            for p in (.80,.95):
                for outcome in ('continuous','binary'):
                    for dr in ('first_factor','interaction'):
                        scenarios.append(make_scenario(F,n,outcome,p,1.,.25,.5 if outcome=='continuous' else .15,
                            dr,'superiority_null','correlated',profile_law='correlated'))
    if suite == 'smoke':
        # No smoke output is to be treated as a performance result.
        scenarios = [make_scenario(2,200,'continuous',.8,1.,.25,.5,'first_factor','superiority_null','smoke'),
                     make_scenario(2,200,'binary',.8,1.,.25,.15,'first_factor','superiority_null','smoke'),
                     make_scenario(5,400,'continuous',.95,.25,.1,1.,'first_factor','noninferiority_null','smoke'),
                     make_scenario(5,400,'binary',.95,1.,.1,.15,'interaction','equivalence_lower','smoke')]
    if not scenarios: raise ValueError("Unknown or empty suite.")
    if len({s['id'] for s in scenarios}) != len(scenarios): raise AssertionError("Duplicate scenarios.")
    for i,s in enumerate(scenarios): s['index']=i
    return dict(format_version=1,suite=suite,profile=profile,seed=int(seed),
                profile_source=profile_source,configuration=PROFILES[profile],
                inference_effect_source='observed_data_only',
                tuning='none; design defined after exploratory diagnostics; independent validation streams',
                scenarios=scenarios)


def generate_observation(sc, master_seed, trial):
    pat,pi = patterns_and_prob(sc['F'],sc['profile_law'])
    n=sc['n']; key=sc['id']
    rp = rng_for(master_seed,'outer',key,int(trial),'profiles')
    sid = np.searchsorted(np.cumsum(pi),rp.random(n)).astype(np.int64)
    ua = rng_for(master_seed,'outer',key,int(trial),'assignment').random(n)
    z = path_from_uniforms(sid,pat,sc['pbc'],ua)
    A = (z.astype(np.int64)+1)//2
    ry = rng_for(master_seed,'outer',key,int(trial),'outcome')
    if sc['outcome']=='continuous':
        base = baseline_means(pat,pi)[sid]+sc['outcome_sd']*ry.normal(size=n)
        effect_noise = sc['individual_effect_sd']*ry.normal(size=n)
        y = base+A*(sc['delta']+np.array(sc['dgp_deviations'])[sid]+effect_noise)
    elif sc['sharp_binary']:
        y = (ry.random(n)<np.array(sc['control_risks'])[sid]).astype(float)
    else:
        prob = np.where(A==1,np.array(sc['treatment_risks'])[sid],np.array(sc['control_risks'])[sid])
        y = (ry.random(n)<prob).astype(float)
    # No true effect field is passed into the data-only analysis API.
    return Observation(y.astype(float),A,sid,pat)
