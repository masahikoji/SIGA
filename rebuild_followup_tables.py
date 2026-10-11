from pathlib import Path
from decimal import Decimal, ROUND_HALF_UP
import pandas as pd
ROOT=Path(__file__).resolve().parent
R=ROOT/'results/report'
# Numeric helpers use the recorded integer rejections and deterministic decimal rounding.
rates=pd.read_csv(R/'rates.csv'); diag=pd.read_csv(R/'diagnostic_rates.csv'); pairs=pd.read_csv(R/'paired_comparisons.csv'); var=pd.read_csv(R/'variance_diagnostics.csv'); eff=pd.read_csv(R/'effect_decomposition.csv'); aud=pd.read_csv(R/'audit_paired.csv')
order=['G_C_n200_null','G_C_n200_power','G_C_n1000_null','G_C_n1000_power','A13_diagnostic']
labels={order[0]:'200, null',order[1]:'200, alternative',order[2]:'1000, null',order[3]:'1000, alternative',order[4]:'400, A13 null'}
def fmt(x,n=3):
    d=Decimal(str(float(x))).quantize(Decimal('1').scaleb(-n),rounding=ROUND_HALF_UP)
    if d==0:d=abs(d)
    return f'{d:.{n}f}'
def sel(df,id,score,construction,method=None):
    a=df[(df.id==id)&(df.score==score)&(df.construction==construction)]
    if method is not None:a=a[a.method==method]
    assert len(a)==1,(id,score,construction,method,len(a))
    return a.iloc[0]
def rate(id,score,cons,method):
    row=sel(rates if method in ['RT','SIGA_S','SIGA_R','CRT_plugin'] else diag,id,score,cons,method)
    p=Decimal(int(row.rejections))*100/Decimal(int(row.available));n=2 if 'power' in id else 3
    return f'{p.quantize(Decimal("1").scaleb(-n),rounding=ROUND_HALF_UP):.{n}f}'
def table(caption,label,spec,header,rows,note):
    return '\n'.join([r'\begin{table}[!htbp]',r'\centering',r'\caption{'+caption+'}',r'\label{'+label+'}',r'\begin{threeparttable}',r'\begin{singlespace}',r'\small',r'\setlength{\tabcolsep}{3.5pt}',r'\renewcommand{\arraystretch}{1.1}',r'\begin{tabular}{'+spec+'}',r'\toprule',header,r'\midrule',*rows,r'\bottomrule',r'\end{tabular}',r'\begin{tablenotes}[flushleft]\footnotesize',r'\item[] '+note,r'\end{tablenotes}',r'\end{singlespace}',r'\end{threeparttable}',r'\end{table}'])
rows=[]
for id in order:
 for score,st in [('unadjusted','U'),('adjusted','A')]:
    vals=[labels[id],st,rate(id,score,'none','RT')]
    vals +=[rate(id,score,c,m) for c in ['base','reconstructed'] for m in ['SIGA_S','SIGA_R','CRT_plugin']]
    rows.append(' & '.join(vals)+r' \\')
t1=table(r'Follow-up rejection probabilities (\%).','tab:supp-focused-rates','@{}llrrrrrrr@{}',r'& & & \multicolumn{3}{c}{Base} & \multicolumn{3}{c}{Reconstructed}\\'+'\n'+r'Setting ($n$, role) & Score & RT & S & R & CRT & S & R & CRT\\',rows,r'S and R denote SIGA-S and SIGA-R; U and A denote unadjusted and adjusted scores. Each null uses 100,000 trials and each alternative 10,000. The two-factor tests have one-sided level 2.5\%; A13 has two-sided level 5\%. All analyses use 4,999 reference assignments where required. The A13 rows are a fresh run, not replacements for Appendix~G.')
rows=[]
for id in order:
 for score,st in [('unadjusted','U'),('adjusted','A')]:
    vals=[labels[id],st]
    for ma,mb in [('SIGA_R','RT'),('CRT_plugin','RT'),('CRT_plugin','SIGA_S')]:
        q=pairs[(pairs.id==id)&(pairs.score==score)&(pairs.construction=='reconstructed')&(pairs.method_A==ma)&(pairs.method_B==mb)].iloc[0]
        vals.append(fmt(100*q.difference)+' ('+fmt(100*q.paired_mcse)+')')
    rows.append(' & '.join(vals)+r' \\')
t2=table(r'Paired differences in follow-up rejection rates (percentage points).','tab:supp-focused-paired','@{}llrrr@{}',r'Setting ($n$, role) & Score & SIGA-R $-$ RT & CRT $-$ RT & CRT $-$ SIGA-S\\',rows,r'Values are differences in rejection proportions followed by paired Monte Carlo standard errors. SIGA-S, SIGA-R and CRT use reconstructed covariances. Comparisons use the same trials and reference assignments; all other recorded comparisons are supplied with the code.')
rows=[]
metrics=[(r'$\overline{\widehat V_{\rm S}}/\widehat{\Var}(T_n)$','VS_over_sampling_variance_at_null',4),(r'$\overline{V^*}/\widehat{\Var}(T_n)$','MC_over_sampling_variance_at_null',4),(r'$\overline{\widehat V_{\rm R}^{\rm raw}}/\overline{\widehat V_{\rm S}}$','VR_plugin_raw_over_VS',4),(r'Median $\widehat\lambda$','VR_plugin_lambda_q50',4),(r'95th percentile of $\widehat\lambda$','VR_plugin_lambda_q95',4)]
for lab,field,dec in metrics:
    rows.append(' & '.join([lab]+[fmt(sel(var,'A13_diagnostic',score,'reconstructed')[field],dec) for score in ['unadjusted','adjusted']])+r' \\')
for lab,method in [(r'CRT using true $\bm d=0$ (\%)','CRT_oracle_d'),(r'Exploratory noise-subtracted CRT (\%)','CRT_noise')]:
    rows.append(' & '.join([lab]+[rate('A13_diagnostic',sc,'reconstructed',method) for sc in ['unadjusted','adjusted']])+r' \\')
for lab,component in [(r'Quadratic noise / mean $\widehat V_{\rm S}$ (\%)','quadratic_noise_component'),(r'Estimated subtraction / mean $\widehat V_{\rm S}$ (\%)','estimated_noise_subtraction'),(r'Singleton expected noise / mean $\widehat V_{\rm S}$ (\%)','singleton_expected_noise')]:
    vals=[]
    for sc in ['unadjusted','adjusted']:
        q=eff[(eff.id=='A13_diagnostic')&(eff.score==sc)&(eff.construction=='reconstructed')&(eff.component==component)].iloc[0]
        vals.append(fmt(100*q.ratio_to_mean_VS))
    rows.append(' & '.join([lab]+vals)+r' \\')
t3=table(r'A13 follow-up variance and effect-estimation diagnostics.','tab:supp-focused-A13','@{}lrr@{}',r'Quantity & Unadjusted & Adjusted\\',rows,r'All entries use 100,000 trials and reconstructed covariance. $V^*$ is the sample variance of the regenerated reference statistics. The true-effect and noise-subtracted procedures are diagnostics, not replacements for the principal CRT. The variance ratios are ratios of means; the quantiles describe the trial-specific scale factor.')
rows=[]
for score,st in [('unadjusted','U'),('adjusted','A')]:
 for method,lab in [('RT','RT'),('SIGA_S','SIGA-S'),('SIGA_R','SIGA-R'),('CRT_plugin','CRT')]:
    vals=[st,lab]
    for code in ['independent_calibration_01','independent_calibration_02','high_reference']:
        q=aud[(aud.score==score)&(aud.construction=='reconstructed')&(aud.method==method)&(aud.audit==code)].iloc[0]
        vals.append(fmt(100*q.difference)+' ('+fmt(100*q.paired_mcse)+')')
    rows.append(' & '.join(vals)+r' \\')
t4=table(r'A13 sensitivity: paired rejection-rate changes (percentage points).','tab:supp-focused-sensitivity','@{}llrrr@{}',r'Score & Method & Calibration 2 & Calibration 3 & 9,999 references\\',rows,r'Each change, with paired Monte Carlo standard error in parentheses, is relative to calibration 1 and 4,999 references on the same first 2,000 trials. SIGA uses reconstructed covariances. Calibration comparisons retain the reference paths; the 4,999 reference draws form the initial subset of the 9,999 draws. These are not additional independent outcome trials.')

out=ROOT/'derived_tables'
out.mkdir(exist_ok=True)
for i,t in enumerate([t1,t2,t3],15):
    (out/f'table_{i}.tex').write_text(t+'\n')
print('Rebuilt Tables 15--17 from the supplied CSVs:',out)
