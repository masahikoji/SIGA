from pathlib import Path
import hashlib, json, math, sys
import numpy as np
import pandas as pd
from scipy.special import ndtr, ndtri

ROOT=Path(__file__).resolve().parent
R=ROOT/'results'
P=ROOT/'program/SIGA_focused_validation_320k_20261011'
checks=[]
def check(name,condition):
    if not bool(condition): raise AssertionError(name)
    checks.append(name)
def close(name,a,b,atol=1e-11):
    check(name,np.allclose(a,b,rtol=1e-10,atol=atol,equal_nan=True))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def stable(obj): return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
plan=json.loads((R/'plan.json').read_text()); resolved=json.loads((R/'resolved_design.json').read_text()); status=json.loads((R/'report/status.json').read_text())
check('plan hash',stable({k:v for k,v in plan.items() if k!='plan_hash'})==plan['plan_hash'])
check('resolved and report plan hashes',resolved['plan_hash']==status['plan_hash']==plan['plan_hash'])
paths=[P/'run_focus.py',P/'UPSTREAM.json']+sorted((P/'revision').glob('*.py'))+sorted((P/'vendor/siga_validation').glob('*.py'))
check('executable source signature',stable({str(p.relative_to(P)):sha(p) for p in paths})==plan['source_signature'])
for rel,d in json.loads((P/'UPSTREAM.json').read_text())['files'].items(): check('vendor '+rel,sha(P/rel)==d['sha256'])
for name,hashval in resolved['calibration_files'].items():
    check('calibration file '+name,sha(R/name)==hashval)
    with np.load(R/name,allow_pickle=False) as f:
        meta=json.loads(str(f['metadata']))
        check('calibration plan '+name,meta['plan_hash']==plan['plan_hash'])
        expected=200000 if meta['stage']=='direction' else 100000
        check('calibration budget '+name,int(f['replicates'])==meta['replicates']==expected)
        for field in ['gamma','psi','sigma_D_finite','xi']:
            a=f[field];close('symmetry '+name+field,a,a.T);check('PSD '+name+field,np.linalg.eigvalsh(a).min()>-1e-8)
        close('profile sum '+name,np.sum(f['profile_prob']),1)
check('distinct calibration objects',len(set(resolved['calibration_files'].values()))==7)
completion=pd.read_csv(R/'report/completion.csv');check('complete scenario counts',len(completion)==5 and completion.complete.all() and (completion.completed==completion.planned).all() and completion.completed.sum()==320000)
check('report completion flag',status['all_scenarios_complete'] and status['completed_trials']==320000)
for s in plan['scenarios']:
    check('planned trial budget '+s['id'],s['outer']==(100000 if s['null_evaluation'] else 10000))
check('reference budget',plan['configuration']['reference']==4999)
check('audit budget',plan['configuration']['audit_trials']==2000 and plan['configuration']['high_reference']==9999 and plan['configuration']['calibration_repeats']==3)
rates=pd.read_csv(R/'report/rates.csv');diag=pd.read_csv(R/'report/diagnostic_rates.csv');allrates=pd.concat([rates,diag],ignore_index=True)
check('all principal methods present',len(rates)==70 and set(rates.method)=={'RT','SIGA_S','SIGA_R','CRT_plugin'})
check('all diagnostic rows present',len(diag)==60)
for name,df in [('rates',allrates),('audit rates',pd.read_csv(R/'report/audit_rates.csv'))]:
    check(name+' complete and available',df.complete.all() and (df.trials==df.available).all())
    n=df.available.to_numpy();k=df.rejections.to_numpy();p=k/n
    close(name+' rates',df.rate,p);close(name+' MCSE',df.mcse,np.sqrt(p*(1-p)/n))
    z=ndtri(.975);den=1+z*z/n;mid=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    close(name+' Wilson lower',df.wilson_lower,mid-half);close(name+' Wilson upper',df.wilson_upper,mid+half)
for name,df in [('paired',pd.read_csv(R/'report/paired_comparisons.csv')),('audit paired',pd.read_csv(R/'report/audit_paired.csv'))]:
    n=df.paired_trials.to_numpy();a=df.A_only.to_numpy();b=df.B_only.to_numpy();delta=(a-b)/n
    close(name+' discordance difference',df.difference,delta)
    close(name+' rate difference',df.difference,df.rate_A-df.rate_B)
    close(name+' paired MCSE',df.paired_mcse,np.sqrt(np.maximum(0,(a+b)/n-delta**2)/(n-1)))
    check(name+' absolute p difference',np.all(np.abs(df.mean_p_difference)<=df.mean_absolute_p_difference+1e-13))
    if name=='paired':
        for _,row in df.iterrows():
            for side in ['A','B']:
                method=row['method_'+side];cons='none' if method=='RT' else row.construction
                found=allrates[(allrates.id==row.id)&(allrates.score==row.score)&(allrates.construction==cons)&(allrates.method==method)]
                check('paired marginal lookup',len(found)==1)
                close('paired marginal value',row['rate_'+side],found.iloc[0].rate)
e=pd.read_csv(R/'report/effect_decomposition.csv');v=pd.read_csv(R/'report/variance_diagnostics.csv')
check('variance rows',len(v)==20)
for _,row in v.iterrows():
    close('variance ratio reconciliation',row.VR_plugin_over_VS*row.mean_VS/row.mean_VR_MC-1,row.VR_plugin_relative_MC_variance_discrepancy)
    if row.null_evaluation:
        close('sampling ratio',row.mean_VS/row.empirical_sampling_variance,row.VS_over_sampling_variance_at_null)
        close('reference sampling ratio',row.mean_VR_MC/row.empirical_sampling_variance,row.MC_over_sampling_variance_at_null)
    ed=e[(e.id==row.id)&(e.score==row.score)&(e.construction==row.construction)].set_index('component')
    close('decomposition ratios',ed['mean']/row.mean_VS,ed.ratio_to_mean_VS)
    correction=ed.loc[['oracle_correction','missing_component','linear_noise_component','quadratic_noise_component'],'mean'].sum()
    close('raw correction decomposition',row.mean_VS+correction,row.VR_plugin_raw_over_VS*row.mean_VS)
    close('noise subtraction',row.VR_noise_raw_over_VS*row.mean_VS,row.VR_plugin_raw_over_VS*row.mean_VS-ed.loc['estimated_noise_subtraction','mean'])
    close('oracle correction',row.VR_oracle_raw_over_VS*row.mean_VS,row.mean_VS+ed.loc['oracle_correction','mean'])
    if row.id=='A13_diagnostic':
        close('zero true A13 correction',ed.loc[['oracle_correction','masked_correction','missing_component','linear_noise_component'],'mean'],0)
check('no recorded floor or fallback',np.max(v[['floor_plugin_rate','floor_noise_rate','reconstruction_fallback_rate']].to_numpy())==0)
first=pd.read_csv(R/'report/first_order_diagnostics.csv')
for sc in resolved['scenarios']:
    pi=np.asarray(sc['profile_prob']);d=np.asarray(sc['mu1'])-np.asarray(sc['mu0'])-sc['boundary']
    close('DGP effect mean '+sc['id'],pi@d,sc['delta']-sc['boundary'])
    if not sc['null_evaluation']:continue
    with np.load(R/'calibration'/f"{sc['calibration_key']}_analysis_00.npz") as cal:
        q=(np.array(sc['mu1'])+np.array(sc['mu0'])-sc['boundary'])/2
        F=sc['F'];pat=((np.arange(2**F)[:,None]>>np.arange(F))&1).astype(float)
        for score in ['unadjusted','adjusted']:
            C=np.ones((len(pi),1)) if score=='unadjusted' else np.column_stack([np.ones(len(pi)),pat])
            coef=np.linalg.lstsq(C*np.sqrt(pi[:,None]),q*np.sqrt(pi),rcond=None)[0]
            m=q-C@coef
            baseline=m@cal['sigma_D_finite']@m+pi@((np.array(sc['var0'])+np.array(sc['var1']))/2)
            vs=(baseline+pi@(d*d)/4)/4;vr=(baseline+d@cal['psi']@d/4)/4
            fact=2 if sc['tail']=='two' else 1
            pred=fact*ndtr(-ndtri(1-sc['alpha']/fact)*np.sqrt(vr/vs))
            r=first[(first.id==sc['id'])&(first.score==score)].iloc[0]
            close('first-order diagnostic '+sc['id']+score,[vs,vr,vr/vs,pred],[r.sampling_variance_per_n,r.reference_variance_per_n,r.reference_to_sampling_ratio,r.predicted_RT_null_rate])
pa=pd.read_csv(R/'report/paired_comparisons.csv')
a13or=pa[(pa.id=='A13_diagnostic')&(pa.method_A=='CRT_oracle_d')&(pa.method_B=='RT')]
check('A13 true zero oracle agrees with RT in supplied paired summaries',len(a13or)==4 and (a13or[['A_only','B_only','mean_absolute_p_difference']]==0).all().all())
notes=['All checks are on supplied calibration arrays, metadata and aggregate summaries. Trial-level checkpoints were not supplied, so no independent trial-level reconstruction is claimed.','The generated report/README.md contains an obsolete sentence saying A13 is not a high-precision size study. The sealed plan and counts establish 100,000 A13 trials; only its calibration/reference sensitivity subset has 2,000 trials. Raw input files have not been changed.','No simulation, theorem, or clinical dataset was newly generated for this audit.']
out={'checks_passed':len(checks),'unique_check_names':len(set(checks)),'main_trials':320000,'gap_trials':220000,'A13_trials':100000,'A13_paired_sensitivity_subset':2000,'base_reference_budget':4999,'source_signature':plan['source_signature'],'plan_hash':plan['plan_hash'],'results_zip_sha256':sha(ROOT/'results_for_review.zip'),'scope_notes':notes,'checks':checks}
(ROOT/'audit_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
