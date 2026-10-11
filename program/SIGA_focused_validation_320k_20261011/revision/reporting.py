"""Trial-level paired summaries. Calibration repeats are NEVER pooled as iid trials."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from .utils import write_csv, atomic_json, se
from .calibration import load_resolved
from .engine import METHODS, METRICS, EFFECT_METRICS, SCORES, CONSTRUCTIONS, PRINCIPAL


def binomial_summary(p,alpha):
    p=np.asarray(p,dtype=float)
    finite=np.isfinite(p); n=int(finite.sum()); total=len(p)
    if n==0:
        return dict(trials=total,available=0,availability=0.,rejections=0,rate=np.nan,mcse=np.nan,
                    wilson_lower=np.nan,wilson_upper=np.nan)
    k=int((p[finite]<=alpha).sum()); r=k/n
    z=1.959963984540054
    den=1+z*z/n; mid=(r+z*z/(2*n))/den
    half=z*np.sqrt(r*(1-r)/n+z*z/(4*n*n))/den
    return dict(trials=total,available=n,availability=n/total,rejections=k,rate=r,
                mcse=np.sqrt(r*(1-r)/n),wilson_lower=mid-half,wilson_upper=mid+half)


def paired_summary(pa,pb,alpha):
    pa=np.asarray(pa); pb=np.asarray(pb)
    good=np.isfinite(pa)&np.isfinite(pb); n=int(good.sum())
    if n==0:
        return dict(paired_trials=0,rate_A=np.nan,rate_B=np.nan,difference=np.nan,
                    paired_mcse=np.nan,A_only=0,B_only=0,mean_p_difference=np.nan,mean_absolute_p_difference=np.nan)
    a=pa[good]<=alpha; b=pb[good]<=alpha
    delta=a.astype(float)-b.astype(float)
    pd=pa[good]-pb[good]
    return dict(paired_trials=n,rate_A=float(a.mean()),rate_B=float(b.mean()),difference=float(delta.mean()),
                paired_mcse=se(delta),A_only=int((a&~b).sum()),B_only=int((b&~a).sum()),
                mean_p_difference=float(pd.mean()),mean_absolute_p_difference=float(np.abs(pd).mean()))


def ratio_summary(x,y):
    x=np.asarray(x); y=np.asarray(y)
    good=np.isfinite(x)&np.isfinite(y)
    x=x[good]; y=y[good]
    if len(y)==0 or not y.mean()>0: return np.nan,np.nan
    ratio=float(x.mean()/y.mean())
    return ratio,se(x-ratio*y)/float(y.mean())


def load_scenario_results(root,sc,plan):
    folder=Path(root)/'trials'/sc['id']
    paths=sorted(folder.glob('*.npz'))
    if not paths: return None
    base={k:[] for k in ('trial','p','metrics','effect','statistic','fallback','seconds')}
    audit={k:[] for k in ('trial','high_p','high_metrics','alt_p','alt_metrics','alt_effect')}
    for path in paths:
        with np.load(path,allow_pickle=False) as f:
            m=json.loads(str(f['metadata']))
            if m['plan_hash']!=plan['plan_hash'] or m['resolved_hash']!=sc['resolved_hash']:
                raise RuntimeError('Checkpoint provenance mismatch '+str(path))
            if m['methods']!=list(METHODS) or m['metrics']!=list(METRICS):
                raise RuntimeError('Checkpoint schema mismatch.')
            for k in base: base[k].append(f[k].copy())
            if len(f['audit_trial']):
                audit['trial'].append(f['audit_trial'].copy())
                for k in ('high_p','high_metrics','alt_p','alt_metrics','alt_effect'):
                    audit[k].append(f['audit_'+k].copy())
    base={k:np.concatenate(v,axis=0) for k,v in base.items()}
    idx=np.argsort(base['trial'])
    base={k:v[idx] for k,v in base.items()}
    if len(np.unique(base['trial']))!=len(base['trial']):
        raise RuntimeError('Overlapping chunks; duplicate trials must not be pooled.')
    if np.any(base['trial']<0) or np.any(base['trial']>=sc['outer']):
        raise RuntimeError('Trial ids outside frozen plan.')
    if audit['trial']:
        audit={k:np.concatenate(v,axis=0) for k,v in audit.items()}
        ai=np.argsort(audit['trial']); audit={k:v[ai] for k,v in audit.items()}
    else:
        audit={}
    return base,audit


def make_report(root):
    plan,resolved=load_resolved(root)
    out=Path(root)/'report'; out.mkdir(exist_ok=True)
    rates=[]; diagnostic_rates=[]; pairs=[]; variance=[]; effects=[]; completion=[]
    audit_rows=[]; audit_pairs=[]; cal_spread=[]
    comparisons=[('SIGA_R','RT'),('SIGA_S','RT'),('CRT_plugin','RT'),
                 ('CRT_plugin','SIGA_S'),('CRT_noise','CRT_plugin'),
                 ('CRT_oracle_d','CRT_plugin'),('CRT_oracle_d','RT'),
                 ('SIGA_R_oracle_d','SIGA_R')]
    for sc in resolved['scenarios']:
        data=load_scenario_results(root,sc,plan)
        common={k:sc[k] for k in ('id','family','F','n','pbc','outcome','model','boundary','delta','alpha','null_evaluation','historical_counterpart')}
        if data is None:
            completion.append(dict(id=sc['id'],planned=sc['outer'],completed=0,complete=False)); continue
        base,audit=data; actual=len(base['trial'])
        complete=actual==sc['outer'] and np.array_equal(base['trial'],np.arange(sc['outer']))
        completion.append(dict(id=sc['id'],planned=sc['outer'],completed=actual,complete=complete))
        common.update(planned=sc['outer'],complete=complete)
        for si,score in enumerate(SCORES):
            for ci,construction in enumerate(CONSTRUCTIONS):
                group=dict(**common,score=score,construction=construction)
                for mi,method in enumerate(METHODS):
                    if method=='RT' and ci==1: continue
                    row=dict(**group,method=method,**binomial_summary(base['p'][:,si,ci,mi],sc['alpha']))
                    if method=='RT': row['construction']='none'
                    (rates if method in PRINCIPAL else diagnostic_rates).append(row)
                for a,b in comparisons:
                    pairs.append(dict(**group,method_A=a,method_B=b,subset='all',
                         **paired_summary(base['p'][:,si,ci,METHODS.index(a)],base['p'][:,si,ci,METHODS.index(b)],sc['alpha'])))
                mat=base['metrics'][:,si,ci,:]
                vs=mat[:,METRICS.index('VS')]
                vrmc=mat[:,METRICS.index('VR_MC')]
                sampling_emp=float(np.var(base['statistic'][:,si],ddof=1)) if actual>1 else np.nan
                vr=dict(**group,trials=actual,mean_VS=float(vs.mean()),
                        mean_VR_MC=float(vrmc.mean()),empirical_sampling_variance=sampling_emp,
                        mean_kappa=float(mat[:,METRICS.index('kappa')].mean()),
                        reconstruction_fallback_rate=float(base['fallback'].mean()),
                        mean_raw_tie_mass=float(mat[:,METRICS.index('raw_tie_mass')].mean()))
                for name in ('VR_plugin','VR_plugin_raw','VR_noise_raw','VR_oracle_raw','VR_masked_raw','VR_MC'):
                    v=mat[:,METRICS.index(name)]
                    ratio,rse=ratio_summary(v,vs)
                    vr[name+'_over_VS']=ratio; vr[name+'_over_VS_mcse']=rse
                    lam=np.sqrt(np.maximum(v,vs/sc['n'])/np.maximum(vs,1e-300))
                    for q in (.05,.5,.95,.99): vr[name+'_lambda_q'+str(int(q*100))]=float(np.quantile(lam,q))
                    if name!='VR_MC': vr[name+'_relative_MC_variance_discrepancy']=float(v.mean()/vrmc.mean()-1) if vrmc.mean()>0 else np.nan
                vr['VS_over_sampling_variance_at_null']=float(vs.mean()/sampling_emp) if sc['null_evaluation'] and sampling_emp>0 else np.nan
                vr['MC_over_sampling_variance_at_null']=float(vrmc.mean()/sampling_emp) if sc['null_evaluation'] and sampling_emp>0 else np.nan
                for name in ('floor_plugin','floor_noise'):
                    vr[name+'_rate']=float(mat[:,METRICS.index(name)].mean())
                variance.append(vr)
                for di,name in enumerate(EFFECT_METRICS):
                    x=base['effect'][:,di]
                    rr,rse=ratio_summary(x,vs)
                    effects.append(dict(**group,component=name,mean=float(x.mean()),mcse=se(x),
                                        ratio_to_mean_VS=rr,ratio_mcse=rse))
        if audit:
            positions=np.searchsorted(base['trial'],audit['trial'])
            if not np.array_equal(base['trial'][positions],audit['trial']): raise RuntimeError('Audit/main trial mismatch.')
            bpa=base['p'][positions]
            versions=[('base_audit_subset',bpa),('high_reference',audit['high_p'])]
            if audit['alt_p'].ndim>1:
                for ri in range(audit['alt_p'].shape[1]):
                    versions.append((f'independent_calibration_{ri+1:02d}',audit['alt_p'][:,ri]))
            for si,score in enumerate(SCORES):
                for ci,construction in enumerate(CONSTRUCTIONS):
                    group=dict(**common,score=score,construction=construction)
                    for mi,method in enumerate(METHODS):
                        if method not in PRINCIPAL: continue
                        cr=[]
                        for label,pa in versions:
                            summary=binomial_summary(pa[:,si,ci,mi],sc['alpha'])
                            audit_rows.append(dict(**group,method=method,audit=label,**summary))
                            if label!='high_reference': cr.append(summary['rate'])
                            if label!='base_audit_subset':
                                audit_pairs.append(dict(**group,method=method,audit=label,
                                   **paired_summary(pa[:,si,ci,mi],bpa[:,si,ci,mi],sc['alpha'])))
                        cal_spread.append(dict(**group,method=method,calibration_batches=len(cr),
                            same_outcome_trials=len(bpa),minimum_rate=min(cr),maximum_rate=max(cr),
                            sd_of_batch_rates=float(np.std(cr,ddof=1)) if len(cr)>1 else np.nan,
                            note='Descriptive spread on identical outcomes/reference paths; batches are not independent new trials.'))
    for name,rows in [('rates',rates),('diagnostic_rates',diagnostic_rates),('paired_comparisons',pairs),
                      ('variance_diagnostics',variance),('effect_decomposition',effects),('completion',completion),
                      ('audit_rates',audit_rows),('audit_paired',audit_pairs),('calibration_sensitivity',cal_spread)]:
        write_csv(out/(name+'.csv'),rows)
    total=sum(r['completed'] for r in completion)
    status=dict(profile=plan['profile'],plan_hash=plan['plan_hash'],
                planned_trials=sum(r['planned'] for r in completion),completed_trials=total,
                all_scenarios_complete=all(r['complete'] for r in completion),
                warning='SMOKE TEST ONLY; not performance evidence.' if plan['profile']=='smoke' else
                        'Pointwise conditional Monte Carlo uncertainty; inspect independent-calibration and high-reference audits separately.')
    atomic_json(out/'status.json',status)
    _write_readme(out,status)
    from .prediction import null_predictions
    write_csv(out/'first_order_diagnostics.csv',null_predictions(root,plan,resolved))
    return status


def _write_readme(out,status):
    text=f"""# Focused follow-up results

Profile: {status['profile']}; completed {status['completed_trials']:,} / {status['planned_trials']:,} unique trials.
{status['warning']}

- rates.csv: the four principal procedures only; RT plus base/reconstructed SIGA-S, SIGA-R and the original plug-in CRT.
- diagnostic_rates.csv: oracle-effect analyses and the existing exploratory noise subtraction, NOT new principal procedures.
- paired_comparisons.csv: matched trial comparisons, including CRT-RT and SIGA-S-RT. Differences are probabilities; multiply by 100 for percentage points.
- variance_diagnostics.csv: ratios of means and quantiles of scale factors. VR_MC is a Monte Carlo diagnostic of the reference variance, not a new CRT implementation.
- effect_decomposition.csv: signed decomposition of the RAW variance correction, not a decomposition of rejection rates.
- audit_rates.csv and audit_paired.csv: A13 on the SAME 2,000 observations across nested reference budgets and independent calibration batches. No pooling as independent trials.
- calibration_sensitivity.csv: descriptive spread across three calibration batches. It is not a confidence interval for calibration uncertainty.
- first_order_diagnostics.csv: null-only variance formulas using DGP quantities and independent finite-n calibration estimates. They are not exact limiting variances or a proof of asymptotic sign.
- completion.csv: all planned scenarios, including unrun ones.

A13_diagnostic is a fresh diagnostic sample, not a reanalysis of the original A13 trial checkpoints, and is not a high-precision size study. The four large-gap scenarios use independent new simulations and retain h=1.5 for the power alternatives; no 80% power target was fitted. Calibration and direction selection are independent of the outcomes. No external comparator and no CRT_MC are run.
"""
    (out/'README.md').write_text(text,encoding='utf-8')
