#!/usr/bin/env python3
"""Recreate every new manuscript table from the supplied result CSVs.
Run: python audit/rebuild_tables.py
Requires pandas and numpy. Does not rerun any outcome simulation.
"""
from pathlib import Path
import pandas as pd
import numpy as np
from decimal import Decimal, ROUND_HALF_UP
from fractions import Fraction
import statistics
import sys
import re
ROOT=Path(__file__).resolve().parent.parent
SC=pd.read_json(ROOT/'audit'/'reconstructed_scenario_specifications.json')
SC['paper_id']=SC.apply(lambda s: ('P' if s['batch']=='primary' else 'A')+f'{int(s["index"])+1:02}',axis=1)
frames={}
for b in ['primary','additional']:
    for n in ['rejection_rates','paired_comparisons','variance_diagnostics']:
        d=pd.read_csv(ROOT/'data'/b/(n+'.csv'));d['batch']=b;frames[b,n]=d
RR=pd.concat([frames[b,'rejection_rates'] for b in ['primary','additional']],ignore_index=True)
PC=pd.concat([frames[b,'paired_comparisons'] for b in ['primary','additional']],ignore_index=True)
VD=pd.concat([frames[b,'variance_diagnostics'] for b in ['primary','additional']],ignore_index=True)
MAP=SC.set_index(['batch','index'])['paper_id'].to_dict()
for d in [RR,PC,VD]:d['paper_id']=[MAP[b,i] for b,i in zip(d['batch'],d['index'])]
method_order=['RT','original_S_normal','refined_S_normal','original_R_normal','refined_R_normal','original_CRT_data','refined_CRT_data']
head_methods=r'RT & $S_0$ & $S_1$ & $R_0$ & $R_1$ & $C_0$ & $C_1$'
method_note=r'$S$, $R$ and $C$ denote ordinary Gaussian SIGA-S, ordinary Gaussian SIGA-R and simulated CRT, respectively; subscripts 0 and 1 denote the base and reconstructed covariance constructions. U and A denote unadjusted and baseline-adjusted scores.'

def decimal_round(value,k):
    q=Decimal(1).scaleb(-k)
    return format(value.quantize(q,rounding=ROUND_HALF_UP),f'.{k}f')
def pct(v,k=2):
    if isinstance(v,Fraction):
        x=Decimal(v.numerator)*100/Decimal(v.denominator)
    else: x=Decimal(str(v))*100
    return decimal_round(x,k)
def fmt(v,k=3):return f'{float(v):.{k}f}'
def pm(v,k=3):
 val=100*float(v)
 if abs(val)<.5*10**(-k):val=0.
 return f'{val:+.{k}f}'
def row(vals):return ' & '.join(map(str,vals))+r' \\'
def rate_values(sub,dec):
 return [pct(sub.loc[sub.method.eq(m),'rate'].iloc[0],dec) for m in method_order]
def table(caption,label,cols,headers,rows,notes='',font=r'\small',tabsep='4pt'):
 return ('\\begin{table}[!htbp]\n\\centering\n\\caption{'+caption+'}\n\\label{'+label+'}\n\\begin{threeparttable}\n\\begin{singlespace}\n'+font+'\n\\setlength{\\tabcolsep}{'+tabsep+'}\n\\renewcommand{\\arraystretch}{1.05}\n\\begin{tabular}{@{}'+cols+'@{}}\n\\toprule\n'+headers+r' \\'+'\n\\midrule\n'+'\n'.join(rows)+'\n\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\\footnotesize\n\\item[] '+notes+'\n\\end{tablenotes}\n\\end{singlespace}\n\\end{threeparttable}\n\\end{table}\n')
def longtable(caption,label,cols,headers,rows,notes='',font=r'\small',tabsep='4pt'):
 return ('\\begingroup\n\\begin{singlespace}\n'+font+'\n\\setlength{\\tabcolsep}{'+tabsep+'}\n\\renewcommand{\\arraystretch}{1.04}\n\\setlength{\\LTleft}{\\fill}\\setlength{\\LTright}{\\fill}\n\\begin{longtable}{@{}'+cols+'@{}}\n\\caption{'+caption+'}\\label{'+label+r'}\\'+'\n\\toprule\n'+headers+r' \\'+'\n\\midrule\n\\endfirsthead\n\\multicolumn{'+str(len(cols))+'}{c}{\\tablename~\\thetable{} -- continued}'+r'\\'+'\n\\toprule\n'+headers+r' \\'+'\n\\midrule\n\\endhead\n\\midrule\n\\multicolumn{'+str(len(cols))+'}{r}{Continued on next page}'+r'\\'+'\n\\endfoot\n\\bottomrule\n\\endlastfoot\n'+'\n'.join(rows)+'\n\\end{longtable}\n\\noindent\\footnotesize '+notes+'\n\\end{singlespace}\n\\endgroup\n')

# Main tables: all primary settings, no outcome-dependent selection of examples.
main_tables=[]
for evaluation,title,label,dec in [('type1','Primary null rejection probabilities (\\%)','tab:main-null-results',3),('power','Primary power estimates (\\%)','tab:main-power',2)]:
 rows=[]
 for test,testlabel in [('superiority','Two-sided superiority'),('noninferiority','Upper-tailed non-inferiority')]:
  rows.append(r'\multicolumn{11}{l}{\textit{'+testlabel+(r' ($\alpha=.05$)' if test=='superiority' else r' ($\alpha=.025$)')+r'}}\\')
  for s in SC[(SC.batch=='primary')&(SC.evaluation==evaluation)&SC.role.str.startswith(test)].sort_values('index').itertuples():
   for score in ['unadjusted','adjusted']:
    sub=RR[(RR.batch=='primary')&(RR['index']==s.index)&(RR.score==score)]
    rows.append(row([s.F,'Cont.' if s.outcome=='continuous' else 'Binary','First' if s.direction=='first_factor' else 'Inter.','U' if score=='unadjusted' else 'A']+rate_values(sub,dec)))
  if test=='superiority':rows.append(r'\addlinespace[3pt]')
 notes=method_note+r' All settings use $p_{\rm bc}=.80$; $n=200$ for $F=2$ and $n=400$ for $F=5$. First and Inter. denote the first-factor and highest-order interaction effect directions.'
 notes+=r' Each row is based on 100,000 trials. For a rate $r$, its Monte Carlo standard error is $\{r(1-r)/100000\}^{1/2}$ on the probability scale.' if evaluation=='type1' else r' Each row is based on 10,000 trials. The generating effects are specified in Supplementary Table~\ref*{S-tab:supp-primary-design}; for a rate $r$, its Monte Carlo standard error is $\{r(1-r)/10000\}^{1/2}$ on the probability scale.'
 main_tables.append(table(title,label,'rlllrrrrrrr',r'$F$ & Outcome & Direction & Score & '+head_methods,rows,notes,font=r'\footnotesize',tabsep='3.1pt'))
rows=[]
for F in [2,5]:
 for outcome in ['continuous','binary']:
  for score in ['unadjusted','adjusted']:
   for ev in ['type1','power']:
    x=PC[(PC.batch=='primary')&(PC.F==F)&(PC.outcome==outcome)&(PC.score==score)&(PC.evaluation==ev)&(PC.reference=='RT')]
    vals=[]
    for m in ['original_R_normal','refined_R_normal']:
     a=[]
     for z in x[x.method==m].itertuples():
      r=RR[(RR.batch==z.batch)&(RR['index']==z.index)&(RR.score==z.score)].set_index('method')
      a.append(Fraction(abs(int(r.loc[m,'rejections'])-int(r.loc['RT','rejections'])),int(z.n_trials)))
     assert len(a)==4
     vals += [pct(statistics.median(a),3),pct(max(a),3)]
    rows.append(row([F,'Cont.' if outcome=='continuous' else 'Binary','U' if score=='unadjusted' else 'A','Null' if ev=='type1' else 'Power']+vals))
main_tables.append(table(r'Absolute SIGA-R--RT rejection-rate differences in the primary evaluation (percentage points).','tab:main-paired','rlllrrrr',r'$F$ & Outcome & Score & Target & \multicolumn{2}{c}{Base} & \multicolumn{2}{c}{Reconstructed}\\\cmidrule(lr){5-6}\cmidrule(lr){7-8} & & & & Median & Maximum & Median & Maximum',rows,r'Each row summarises four scenario--score comparisons: two effect directions and two testing objectives. Null and power results are summarised separately. SIGA-R uses ordinary Gaussian tails here; the lattice approximation is reported separately. Medians and maxima are descriptive summaries, not confidence bounds. Differences are reconstructed from integer rejection counts; medians are computed before rounding, with exact decimal ties rounded away from zero. Scenario-specific paired Monte Carlo standard errors are in Supplementary Table~\ref*{S-tab:supp-primary-paired}. U and A denote unadjusted and baseline-adjusted scores.'))

# Supplementary design and result tables.
g_tables=[]
rows=[]
for F,n in [(2,200),(5,400)]:
 for outcome in ['continuous','binary']:
  s=SC[(SC.batch=='primary')&(SC.F==F)&(SC.outcome==outcome)&(SC.direction=='first_factor')]
  sup=s[s.role=='superiority_power'].iloc[0];ni=s[s.role=='noninferiority_power'].iloc[0]
  ranges=('P01--P08' if outcome=='continuous' else 'P09--P16') if F==2 else ('P17--P24' if outcome=='continuous' else 'P25--P32')
  rows.append(row([ranges,F,n,'Cont.' if outcome=='continuous' else 'Binary',fmt(sup.actual_amplitude,2),fmt(sup.delta,3),fmt(ni.boundaries[0],2),fmt(ni.delta,3)]))
g_tables.append(table(r'Primary scenario settings.','tab:supp-primary-design','lrrlrrrr',r'IDs & $F$ & $n$ & Outcome & $a_d$ & $\Delta_{\rm sup}$ & $b_{\rm NI}$ & $\Delta_{\rm NI}$',rows,r'All primary scenarios use $p_{\rm bc}=.80$. Within each eight-scenario range, the first four use the first-factor direction and the next four use the interaction direction; the order is superiority null, superiority alternative, non-inferiority null and non-inferiority alternative. Superiority tests $b=0$. Non-inferiority nulls have $\Delta=b_{\rm NI}$. The alternative effects are $\Delta_{\rm sup}$ and $\Delta_{\rm NI}$. The amplitude $a_d=\max_s|d_s^\circ|$ is centred around the generating average effect, not the tested boundary.'))
rows=[]
for s in SC[SC.batch=='additional'].sort_values('index').itertuples():
 rows.append(row([s.paper_id,s.F,s.n,'Cont.' if s.outcome=='continuous' else 'Binary',fmt(s.pbc,2),fmt(s.outcome_sd,2) if s.outcome=='continuous' else '--',fmt(s.individual_effect_sd,2) if s.outcome=='continuous' else '--',fmt(s.actual_amplitude,2),fmt(s.delta,3),'Null' if s.evaluation=='type1' else 'Power']))
g_tables.append(table(r'Additional scenario settings.','tab:supp-additional-design','lrrlrrrrrl',r'ID & $F$ & $n$ & Outcome & $p_{\rm bc}$ & $\sigma$ & $\sigma_\delta$ & $a_d$ & $\Delta$ & Target',rows,r'All additional scenarios use two-sided superiority at $b=0$, $\alpha=.05$. A01--A04 are sharp-null controls; A05--A12 increase sample size; A13--A20 use $p_{\rm bc}=.95$. Nonzero deviations use the first-factor direction. Zero stratum deviations in A13, A14, A17 and A18 do not impose constant individual effects because $\sigma_\delta=.25$. The power effects in A06, A08, A10 and A12 differ from their smaller-sample counterparts. The primary evaluation supplies a matched $p_{\rm bc}=.80$ comparison only for A19 and A20 (P17 and P18, respectively).',font=r'\footnotesize',tabsep='4pt'))
rows=[]
for s in SC[SC.batch=='additional'].sort_values('index').itertuples():
 for score in ['unadjusted','adjusted']:
  sub=RR[(RR.batch=='additional')&(RR['index']==s.index)&(RR.score==score)]
  rows.append(row([s.paper_id,'U' if score=='unadjusted' else 'A']+rate_values(sub,3 if s.evaluation=='type1' else 2)))
g_tables.append(table(r'Rejection probabilities in the 20 additional scenarios (\%).','tab:supp-additional-results','llrrrrrrr',r'ID & Score & '+head_methods,rows,method_note+r' Scenario IDs are defined in Supplementary Table~\ref{tab:supp-additional-design}. Null scenarios use 100,000 trials and alternatives 10,000 trials. Entries of 100.00 mean that all 10,000 simulated trials rejected; they are not claims of population power exactly one.',font=r'\footnotesize',tabsep='4pt'))
# Pairwise differences table: 64 rows, exact per scenario rather than selected successes.
rows=[]
for s in SC[SC.batch=='primary'].sort_values('index').itertuples():
 for score in ['unadjusted','adjusted']:
  x=PC[(PC.batch=='primary')&(PC['index']==s.index)&(PC.score==score)]
  v=[]
  for m,ref in [('original_R_normal','RT'),('refined_R_normal','RT'),('refined_CRT_data','original_CRT_data')]:
   a=x[(x.method==m)&(x.reference==ref)].iloc[0];v.append(pm(a.paired_rejection_difference,3)+' ('+pct(a.paired_mcse,3)+')')
  rows.append(row([s.paper_id,'U' if score=='unadjusted' else 'A']+v))
g_tables.append(longtable(r'Primary paired rejection-rate differences with Monte Carlo standard errors (percentage points).','tab:supp-primary-paired','lllll',r'ID & Score & Base R--RT & Recon. R--RT & Recon. CRT--base CRT',rows,r'Entries are signed paired differences, with paired Monte Carlo standard errors in parentheses. The mapping from P01--P32 to settings is in Supplementary Table~\ref{tab:supp-primary-design}. R is the ordinary Gaussian SIGA-R. A displayed zero standard error means that the recorded paired decisions have zero empirical variance; it is not a proof of identical procedures.',font=r'\footnotesize',tabsep='6pt'))
# Variance summaries segregate null/power; report signed raw separately in supplied CSV.
rows=[]
for group in ['main','controls','sample_size','pbc95']:
 for F in sorted(VD[VD.group==group].F.unique()):
  for outcome in ['continuous','binary']:
   for score in ['unadjusted','adjusted']:
    for ev in ['type1','power']:
     x=VD[(VD.group==group)&(VD.F==F)&(VD.outcome==outcome)&(VD.score==score)&(VD.evaluation==ev)]
     if x.empty:continue
     vals=[]
     for c in ['original','refined']:
      a=x[x.construction==c].R_bias_ratio_of_means.abs(); vals += [pct(a.median(),2),pct(a.max(),2)]
     rows.append(row([{'main':'Primary','controls':'Sharp null','sample_size':'Larger $n$','pbc95':'$.95$'}[group],F,'Cont.' if outcome=='continuous' else 'Binary','U' if score=='unadjusted' else 'A','Null' if ev=='type1' else 'Power']+vals))
g_tables.append(longtable(r'Absolute relative discrepancies in mean conditional reference variance (\%).','tab:supp-variance-summary','lrlllrrrr',r'Family & $F$ & Outcome & Score & Target & \multicolumn{2}{c}{Base} & \multicolumn{2}{c}{Reconstructed}\\\cmidrule(lr){6-7}\cmidrule(lr){8-9} & & & & & Median & Maximum & Median & Maximum',rows,r'The discrepancy for each scenario and score is $100|\overline{\widehat V_R}/\overline{V^*}-1|$, where $V^*$ is the sample variance of 4,999 regenerated statistics and the bars average over outcome trials. Summaries give equal weight to each scenario--score combination within a row, not to each outcome trial. Null and power scenarios are kept separate. These diagnostics are ratios of means, not means of trial-level ratios. Values are displayed to two decimal places; the precision displayed is not a bound on allocation-calibration or outcome-simulation error.',font=r'\footnotesize',tabsep='4pt'))
# Lattice results, all relevant conditions in both runs.
rows=[]
for s in SC[(SC.outcome=='binary')&SC.role.str.startswith('superiority')].sort_values(['batch','index'],ascending=[False,True]).itertuples():
 x=RR[(RR.batch==s.batch)&(RR['index']==s.index)&(RR.score=='unadjusted')].set_index('method')
 vals=[pct(x.loc[m,'rate'],3 if s.evaluation=='type1' else 2) for m in ['RT','original_R_normal','refined_R_normal','original_R_lattice','refined_R_lattice','original_CRT_lattice','refined_CRT_lattice']]
 rows.append(row([s.paper_id,'Null' if s.evaluation=='type1' else 'Power']+vals))
g_tables.append(table(r'Unadjusted binary superiority: ordinary and lattice-based reference approximations (\%).','tab:supp-lattice-results','llrrrrrrr',r'ID & Target & RT & Normal $R_0$ & Normal $R_1$ & Lattice $R_0$ & Lattice $R_1$ & Lattice $C_0$ & Lattice $C_1$',rows,r'The lattice C columns approximate the rescaled, raw-tie-inclusive CRT tail as defined in Appendix~F; they are not the simulated CRT results in the main tables. The ordinary Gaussian S procedures remain separate. All rows use $b=0$ and unadjusted scores. Subscripts 0 and 1 denote base and reconstructed covariances.',font=r'\scriptsize',tabsep='3.1pt'))
# Sparsity counts: empty, one arm only, and their total. No duplication by score.
rows=[]
for group in ['main','controls','sample_size','pbc95']:
 for F,n in sorted(set(map(tuple,VD[VD.group==group][['F','n']].values))):
  x=VD[(VD.group==group)&(VD.F==F)&(VD.n==n)].drop_duplicates(['id'])
  fr=lambda a:f'{a.min():.2f}--{a.max():.2f}' if a.min()!=a.max() else f'{a.min():.2f}'
  rows.append(row([{'main':'Primary','controls':'Sharp null','sample_size':'Larger sample','pbc95':'Stronger coin'}[group],F,n,len(x),fr(2**F-x.mean_represented),fr(x.mean_represented-x.mean_both_arms),fr(x.mean_zero_filled_strata),fr(100*x.mean_zero_filled_profile_mass)]))
g_tables.append(table(r'Strata with at least one unrepresented treatment arm, including empty strata.','tab:supp-sparsity','lrrrllll',r'Family & $F$ & $n$ & Scenarios & Empty & One arm only & Total & Mass (\%)',rows,r'Entries summarise trial-average counts within scenarios; ranges are over scenario averages. Empty strata contain no participant. One-arm-only strata contain participants but none in the other arm. Total is their sum; all such strata use $\widehat d_s(b)=0$. Mass is the total design-profile probability of empty and one-arm-only strata combined. Displayed endpoints may not add exactly after rounding.',font=r'\footnotesize',tabsep='3pt'))
# Estimated variance ratios are diagnostics, not estimates of a null sampling
# variance under alternatives and not exact limiting variance ratios.
rows=[]
for F in [2,5]:
 for outcome in ['continuous','binary']:
  for score in ['unadjusted','adjusted']:
   x=VD[(VD.batch=='primary')&(VD.F==F)&(VD.outcome==outcome)&(VD.score==score)&(VD.evaluation=='type1')]
   vals=[]
   for c in ['original','refined']:
    z=x[x.construction==c];a=z.mean_estimated_R/z.mean_estimated_S
    vals.append(f'{a.min():.4f}--{a.max():.4f}')
   z=x[x.construction=='refined'];a=z.mean_reference_variance/z.empirical_sampling_variance
   vals.append(f'{a.min():.4f}--{a.max():.4f}')
   rows.append(row([F,'Cont.' if outcome=='continuous' else 'Binary','U' if score=='unadjusted' else 'A']+vals))
g_tables.append(table(r'Variance ratios across primary null scenarios.','tab:supp-primary-ratios','rlllll',r'$F$ & Outcome & Score & Estimated, base & Estimated, recon. & Empirical',rows,r'Estimated entries are ranges of $\overline{\widehat V_R}/\overline{\widehat V_S}$ over four null scenarios per row. Empirical entries use $\overline{V^*}/\widehat{\operatorname{Var}}(T_n)$, with sampling variance computed across outcome trials. These are finite-sample diagnostics, not evaluations of $v_R/v_S$ from known limiting covariances. They are not means of the trial-level ratios $\widehat\lambda_n^2$. Alternatives are excluded.',font=r'\footnotesize',tabsep='3pt'))
# Show all four stronger-coin null scenarios, not only the most favourable.
rows=[]
for sid in ['A13','A15','A17','A19']:
 for score in ['unadjusted','adjusted']:
  x=RR[(RR.paper_id==sid)&(RR.score==score)].set_index('method')
  vals=[]
  for m in ['RT','original_CRT_data','refined_CRT_data','original_CRT_noise','refined_CRT_noise']:
   vals.append(pct(x.loc[m,'rate'],3)+' ('+pct(x.loc[m,'mcse'],3)+')')
  rows.append(row([sid,'U' if score=='unadjusted' else 'A']+vals))
g_tables.append(table(r'Stronger-coin null rejection rates and exploratory noise subtraction (\%).','tab:supp-noise-null','lllllll',r'ID & Score & RT & CRT base & CRT recon. & Noise base & Noise recon.',rows,r'All scenarios use $F=5,n=400,p_{\rm bc}=.95$, two-sided level $.05$ and 100,000 trials. Entries are rates with marginal Monte Carlo standard errors in parentheses, not paired errors for the differences. Noise denotes the exploratory quadratic subtraction in Remark~\ref{rem:supp-noise}; it is not the primary CRT. The four null scenarios and both scores are all included. Original outcome simulations and test definitions are unchanged.',font=r'\footnotesize',tabsep='3pt'))


DEST=ROOT/'audit'/'regenerated_tables'
DEST.mkdir(exist_ok=True)
import hashlib, json
expected=json.loads((ROOT/'audit'/'table_expected_hashes.json').read_text())
checks=[]
for i,t in enumerate(main_tables,1):
    if i==3:
        continue  # Optional derived summary was not part of the 11 archived full tables.
    name=f'previous_main_{i:02}.tex'
    assert hashlib.sha256(t.encode()).hexdigest()==expected[name], name
    (DEST/name).write_text(t)
    checks.append(name)
for i,t in enumerate(g_tables,1):
    name=f'supplement_G_{i:02}.tex'
    assert hashlib.sha256(t.encode()).hexdigest()==expected[name], name
    (DEST/name).write_text(t)
    checks.append(name)
assert len(checks)==11 and set(checks)==set(expected)
print('Verified',len(checks),'table exports by SHA256 against the frozen aggregate results.')
print('Generated table fragments in audit/regenerated_tables; no journal manuscript source required.')
