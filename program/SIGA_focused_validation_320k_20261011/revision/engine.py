"""Feasible analyses are kept separate from explicitly labelled oracle diagnostics."""
from __future__ import annotations
import time
import numpy as np
from .utils import rng_for
from siga_validation.core import (scores_and_variances, normal_pvalue, mc_pvalue,
                                   effect_deviations, reference_variance, TIE_TOL)
from siga_validation.allocation import reference_batch

SCORES = ('unadjusted', 'adjusted')
CONSTRUCTIONS = ('base', 'reconstructed')
UPSTREAM_TAGS = ('original', 'refined')
METHODS = ('RT','SIGA_S','SIGA_R','CRT_plugin','SIGA_R_oracle_d','CRT_oracle_d','CRT_noise')
PRINCIPAL = METHODS[:4]
DIAGNOSTIC_METHODS = METHODS[4:]
METRICS = ('VS','VR_plugin_raw','VR_plugin','VR_noise_raw','VR_MC',
           'reference_mean','raw_tie_mass','kappa','between','within','floor_plugin',
           'floor_noise','VR_oracle_raw','VR_masked_raw')
EFFECT_METRICS = ('oracle_correction','masked_correction','missing_component',
                 'linear_noise_component','quadratic_noise_component',
                 'expected_quadratic_noise','estimated_noise_subtraction',
                 'singleton_expected_noise','missing_strata','empty_strata',
                 'missing_mass','weighted_squared_effect_error','weighted_dhat_mean')


def reference_statistics(obs,pbc,entries,seed,scenario_id,trial,draws,batch=256):
    scores=np.column_stack([e['score'] for e in entries])
    rng=rng_for(seed,'reference',scenario_id,int(trial))
    result=np.empty((draws,len(entries)))
    for start in range(0,draws,batch):
        m=min(batch,draws-start)
        u=rng.random((m,len(obs.y)))
        result[start:start+m]=reference_batch(obs.stratum,obs.patterns,pbc,scores,u)
    return result


def _crt(samples,t,vs,vr,tail):
    if not vs>0 or not vr>0:
        # Match upstream's degenerate denominator convention, and flag in metrics.
        return 1.0
    return mc_pvalue(samples,t,tail,np.sqrt(vr/vs),retain_raw_ties=True)


def feasible_from_entries(entries,samples,tail):
    """Original four procedures. No true DGP parameter is accepted here."""
    if samples.shape[0]<2 or samples.shape[1]!=2:
        raise ValueError('Expected >=2 reference draws for two score columns.')
    p=np.full((2,2,len(METHODS)),np.nan)
    metrics=np.full((2,2,len(METRICS)),np.nan)
    for si,e in enumerate(entries):
        t=e['statistic']; refs=samples[:,si]
        vmc=float(np.var(refs,ddof=1))  # Diagnostic ONLY, never a CRT denominator.
        tie=np.abs(np.abs(refs)-abs(t))<=TIE_TOL if tail=='two' else np.abs(refs-t)<=TIE_TOL
        rt=mc_pvalue(refs,t,tail)
        for ci,tag in enumerate(UPSTREAM_TAGS):
            sv=e[tag]['S']; rv=e[tag]['R']; rn=e[tag]['R_noise']; vs=sv['variance']
            p[si,ci,:4]=[rt,normal_pvalue(t,vs,tail),normal_pvalue(t,rv['variance'],tail),
                          _crt(refs,t,vs,rv['variance'],tail)]
            p[si,ci,METHODS.index('CRT_noise')]=_crt(refs,t,vs,rn['variance'],tail)
            values=dict(VS=vs,VR_plugin_raw=rv['raw'],VR_plugin=rv['variance'],
                        VR_noise_raw=rn['raw'],VR_MC=vmc,reference_mean=float(refs.mean()),
                        raw_tie_mass=float(tie.mean()),kappa=sv['kappa'],between=sv['between'],
                        within=sv['within'],floor_plugin=float(rv['floor']),floor_noise=float(rn['floor']))
            for k,v in values.items(): metrics[si,ci,METRICS.index(k)]=v
    return p,metrics


def analyse_feasible(obs,boundary,calibration,samples,tail):
    entries,geometry=scores_and_variances(obs,[boundary],calibration,noise_adjustment=True)
    p,metrics=feasible_from_entries(entries,samples,tail)
    return dict(p=p,metrics=metrics,entries=entries,geometry=geometry)


def oracle_diagnostics(obs,entries,samples,calibration,true_d,var0,var1,pi,tail):
    """SIMULATION ONLY. Never used by the feasible analysis functions above.

    Exact raw correction decomposition:
      Q(dhat)-Q(d) = [Q(I*d)-Q(d)] + 2(I*d)'M*error + error'M*error.
    All correction terms here are on the variance scale n/16, not p-value scale.
    """
    n=len(obs.y); J=len(obs.patterns)
    counts=np.bincount(obs.stratum,minlength=J)
    sizes=[np.bincount(obs.stratum[obs.assignment==a],minlength=J) for a in (0,1)]
    b=entries[0]['boundary']
    dhat,vd,both=effect_deviations(obs,b)
    d=np.asarray(true_d,dtype=float)
    masked=d*both
    error=dhat-masked
    M=calibration['psi']-np.diag(counts/n)
    true_vd=np.where(both,np.asarray(var0)/np.maximum(sizes[0],1)+np.asarray(var1)/np.maximum(sizes[1],1),0.)
    singleton_vd=np.where(both,np.asarray(var0)*(sizes[0]==1)+np.asarray(var1)*(sizes[1]==1),0.)
    fac=n/16.
    qo=fac*float(d@M@d); qm=fac*float(masked@M@masked)
    linear=fac*2*float(masked@M@error); quadratic=fac*float(error@M@error)
    qhat=fac*float(dhat@M@dhat)
    if not np.isclose(qhat-qo,(qm-qo)+linear+quadratic,atol=1e-9,rtol=1e-10):
        raise AssertionError('Raw correction decomposition failed.')
    dvals=dict(oracle_correction=qo,masked_correction=qm,missing_component=qm-qo,
               linear_noise_component=linear,quadratic_noise_component=quadratic,
               expected_quadratic_noise=fac*float(np.diag(M)@true_vd),
               estimated_noise_subtraction=fac*float(np.diag(M)@vd),
               singleton_expected_noise=fac*float(np.diag(M)@singleton_vd),
               missing_strata=float((~both).sum()),empty_strata=float((counts==0).sum()),
               missing_mass=float(np.asarray(pi)[~both].sum()),
               weighted_squared_effect_error=float(np.asarray(pi)@((dhat-d)**2)),
               weighted_dhat_mean=float(np.asarray(pi)@dhat))
    effect=np.array([dvals[k] for k in EFFECT_METRICS])
    po=np.empty((2,2,2)); vo=np.empty((2,2,2))
    for si,e in enumerate(entries):
        for ci,tag in enumerate(UPSTREAM_TAGS):
            vs=e[tag]['S']['variance']; t=e['statistic']
            full=reference_variance(vs,counts,d,calibration['psi'])
            mask=reference_variance(vs,counts,masked,calibration['psi'])
            po[si,ci]=[normal_pvalue(t,full['variance'],tail),
                       _crt(samples[:,si],t,vs,full['variance'],tail)]
            vo[si,ci]=[full['raw'],mask['raw']]
    return po,vo,effect


def evaluate_with_truth(obs,entries,samples,calibration,sc):
    # A wrapper for simulation only; principal values are computed BEFORE truth is read.
    p,metrics=feasible_from_entries(entries,samples,sc['tail'])
    true_d=np.asarray(sc['mu1'])-np.asarray(sc['mu0'])-sc['boundary']
    po,vo,effect=oracle_diagnostics(obs,entries,samples,calibration,true_d,
                                  sc['var0'],sc['var1'],sc['profile_prob'],sc['tail'])
    p[:,:,4:6]=po
    if np.all(true_d==0):
        if not np.array_equal(p[:,:,METHODS.index('CRT_oracle_d')],p[:,:,METHODS.index('RT')]):
            raise AssertionError('Zero-effect oracle CRT must agree with RT on every trial.')
    metrics[:,:,METRICS.index('VR_oracle_raw')]=vo[:,:,0]
    metrics[:,:,METRICS.index('VR_masked_raw')]=vo[:,:,1]
    return p,metrics,effect


def one_trial(sc,plan,trial,calibrations):
    from .designs import generate_observation
    cfg=plan['configuration']
    t0=time.perf_counter()
    obs=generate_observation(sc,plan['seed'],trial)
    entries,geometry=scores_and_variances(obs,[sc['boundary']],calibrations[0],noise_adjustment=True)
    t1=time.perf_counter()
    audit=bool(sc['audit'] and trial<cfg['audit_trials'])
    draws=cfg['high_reference'] if audit else cfg['reference']
    refs=reference_statistics(obs,sc['pbc'],entries,plan['seed'],sc['id'],trial,draws)
    t2=time.perf_counter()
    base_refs=refs[:cfg['reference']]
    p,metrics,effect=evaluate_with_truth(obs,entries,base_refs,calibrations[0],sc)
    extra={}
    if audit:
        ph,mh,eh=evaluate_with_truth(obs,entries,refs,calibrations[0],sc)
        alt_p=[]; alt_metrics=[]; alt_effect=[]
        for cal in calibrations[1:]:
            en,_=scores_and_variances(obs,[sc['boundary']],cal,noise_adjustment=True)
            for original,other in zip(entries,en):
                if not np.array_equal(original['score'],other['score']):
                    raise AssertionError('Calibration changed a supposedly fixed outcome score.')
            pa,ma,ea=evaluate_with_truth(obs,en,base_refs,cal,sc)
            alt_p.append(pa); alt_metrics.append(ma); alt_effect.append(ea)
        extra=dict(high_p=ph,high_metrics=mh,alt_p=np.array(alt_p),
                   alt_metrics=np.array(alt_metrics),alt_effect=np.array(alt_effect))
    t3=time.perf_counter()
    return dict(p=p,metrics=metrics,effect=effect,statistic=np.array([e['statistic'] for e in entries]),
                fallback=np.array(float(geometry['fallback'])),
                seconds=np.array([t1-t0,t2-t1,t3-t2]),audit=extra)
