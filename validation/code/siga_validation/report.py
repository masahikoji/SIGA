"""Aggregate all predeclared scenarios, retaining unfavorable results."""
from __future__ import annotations
import csv,json, math
from pathlib import Path
import numpy as np
from scipy.stats import norm
from .experiment import load_plan, METHODS, DIAGNOSTICS


def write_csv(path,rows):
    if not rows: raise RuntimeError('No report rows.')
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def se(x):
    return float(np.std(x,ddof=1)/np.sqrt(len(x))) if len(x)>1 else float('nan')


def wilson(k,n):
    z=norm.ppf(.975);p=k/n;den=1+z*z/n
    cen=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return cen-half,cen+half


def aggregate(root,allow_partial=False):
    root=Path(root);plan=load_plan(root)
    cfg=plan['configuration'];expected=cfg['outer']
    rates=[];pairs=[];diags=[];coverage=[]
    for sc in plan['scenarios']:
        parts=sorted((root/'trials'/sc['id']).glob('*.npz'))
        blocks=[];indices=[]
        for p in parts:
            with np.load(p,allow_pickle=False) as z:
                if str(z['plan_hash'])!=plan['plan_hash'] or str(z['scenario_id'])!=sc['id']:
                    raise RuntimeError(f'Incompatible checkpoint: {p}')
                blocks.append({k:z[k] for k in ['p','reject','diagnostic','fallback','balance_error','seconds']})
                indices.extend(z['trial_index'].tolist())
        if len(indices)!=len(set(indices)): raise RuntimeError(f'Duplicate outer trials: {sc["id"]}')
        complete=sorted(indices)==list(range(expected))
        coverage.append(dict(scenario=sc['index'],id=sc['id'],n_trials=len(indices),expected=expected,complete=complete))
        if not complete and not allow_partial:
            raise RuntimeError(f'Scenario {sc["index"]}: {len(indices)}/{expected} trials. Finish all scenarios or explicitly use --allow-partial.')
        if not blocks: continue
        a={k:np.concatenate([x[k] for x in blocks],axis=0) for k in blocks[0]}
        num=len(indices)
        meta={k:sc[k] for k in ['index','id','group','F','n','outcome','pbc','direction','role','delta','alpha','actual_amplitude','feasibility_scale']}
        meta.update(n_trials=num,complete=complete)
        for score_idx,score in enumerate(['unadjusted','adjusted']):
            for m,name in enumerate(METHODS):
                k=int(a['reject'][:,score_idx,m].sum());ph=k/num;low,hi=wilson(k,num)
                rates.append(dict(**meta,score=score,method=name,rejections=k,rate=ph,
                        mcse=math.sqrt(ph*(1-ph)/num),ci_lower=low,ci_upper=hi,
                        difference_from_nominal=ph-sc['alpha'],
                        note='null rejection probability' if sc['null_evaluation'] else 'power under fixed alternative'))
            contrasts=[]
            for stem in ['S_normal','R_normal','CRT_data','R_lattice','CRT_lattice','CRT_noise']:
                contrasts.append((f'refined_{stem}',f'original_{stem}'))
            for tag in ('original','refined'):
                contrasts.extend([(f'{tag}_R_normal','RT'),(f'{tag}_R_lattice','RT'),
                                  (f'{tag}_CRT_lattice',f'{tag}_CRT_data'),
                                  (f'{tag}_S_normal',f'{tag}_CRT_data')])
            for x,y in contrasts:
                i,j=METHODS.index(x),METHODS.index(y)
                delta=a['reject'][:,score_idx,i].astype(float)-a['reject'][:,score_idx,j]
                pd=a['p'][:,score_idx,i]-a['p'][:,score_idx,j]
                pairs.append(dict(**meta,score=score,method=x,reference=y,
                       paired_rejection_difference=float(delta.mean()),paired_mcse=se(delta),
                       a_only=int((delta==1).sum()),b_only=int((delta==-1).sum()),
                       mean_pvalue_difference=float(pd.mean()),mean_absolute_pvalue_error=float(np.abs(pd).mean())))
        for bi,b in enumerate(sc['boundaries']):
            for sj,score in enumerate(['unadjusted','adjusted']):
                d=a['diagnostic'][:,2*bi+sj,:]
                empirical=d[:,1]
                for tag,vScol,vRcol in [('original',2,4),('refined',3,5)]:
                    vS,vR=d[:,vScol],d[:,vRcol]
                    good=empirical>0
                    ratio=np.divide(empirical,vR,out=np.full(num,np.nan),where=vR>0)
                    # Conditional moments are estimated from B reference sequences,
                    # NOT truth values supplied to the tested procedure.
                    row=dict(**meta,boundary=b,score=score,construction=tag,
                        mean_reference_variance=float(empirical.mean()),mean_estimated_S=float(vS.mean()),
                        mean_estimated_R=float(vR.mean()),
                        R_bias_ratio_of_means=float(vR.mean()/empirical.mean()-1) if empirical.mean()>0 else np.nan,
                        mean_reference_over_estimate=float(np.nanmean(ratio)),mcse_mean_ratio=se(ratio[np.isfinite(ratio)]),
                        empirical_sampling_variance=float(d[:,0].var(ddof=1)),
                        S_bias_ratio_of_means=(float(vS.mean()/d[:,0].var(ddof=1)-1) if abs(b-sc['delta'])<1e-12 and d[:,0].var(ddof=1)>0 else np.nan),
                        S_variance_diagnostic_valid=abs(b-sc['delta'])<1e-12,
                        mean_ref_between=float(d[:,6].mean()),mean_ref_within=float(d[:,7].mean()),
                        mean_ref_twice_cross=float(d[:,8].mean()),
                        mean_both_arms=float(d[:,15].mean()),mean_represented=float(d[:,16].mean()),
                        rank_fallback_count=int(a['fallback'].sum()),
                        zero_reference_variance_count=int((~good).sum()),
                        kappa_floor_count=int(d[:,17+(tag=='refined')].sum()),
                        R_floor_count=int(d[:,19+(tag=='refined')].sum()))
                    diags.append(row)
    out=root/('partial_summary' if allow_partial else 'summary')
    write_csv(out/'coverage.csv',coverage)
    write_csv(out/'rejection_rates.csv',rates);write_csv(out/'paired_comparisons.csv',pairs)
    write_csv(out/'variance_diagnostics.csv',diags)
    (out/'README.txt').write_text(
        'All requested scenarios and methods are retained. A partial summary is not a completed validation.\n'
        'Monte Carlo SEs condition on the shared allocation-covariance estimates.\n'
        'The primary methods use estimated stratum effects, never DGP effect values.\n'
        'R_lattice and CRT_lattice use the lattice only for unadjusted binary b=0 superiority.\n'
        'Elsewhere they use Gaussian R and S tails respectively. Normal-only columns are ablations.\n'
        'CRT_noise is exploratory and not the primary proposed procedure.\n'
        'Only null-boundary rows interpret mean S as the variance of the observed statistic.\n'
        'Reference variances are Monte Carlo estimates from the common inner draws, not exact values.\n'
        'Scenario-wise intervals are descriptive, not simultaneous multiplicity-adjusted guarantees.\n')
    return str(out)
