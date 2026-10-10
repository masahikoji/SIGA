#!/usr/bin/env python3
"""Generate/check the compact primary power table using integer rejection counts."""
from pathlib import Path
from decimal import Decimal
import argparse, csv, json, re

ROOT = Path(__file__).resolve().parents[1]
METHODS = ['RT', 'refined_S_normal', 'refined_R_normal', 'refined_CRT_data']
START = '% BEGIN GENERATED POWER SUMMARY'
END = '% END GENERATED POWER SUMMARY'

def build():
    with (ROOT/'audit'/'primary_rejection_rates.csv').open(newline='') as f:
        source = list(csv.DictReader(f))
    rows, audit = [], []
    for role, title in [('superiority_power', r'Two-sided superiority ($\alpha=.05$)'),
                        ('noninferiority_power', r'Upper-tailed non-inferiority ($\alpha=.025$)')]:
        # Original code uses 'noninferiority_power'; fail rather than guess a match.
        if role not in {r['role'] for r in source}:
            raise ValueError(f'Missing exact role {role}; found {sorted({r["role"] for r in source})}')
        rows.append(r'\multicolumn{7}{@{}l}{\textit{'+title+r'}}\\')
        for F in [2,5]:
            for outcome in ['continuous','binary']:
                for score in ['unadjusted','adjusted']:
                    n = 200 if F==2 else 400
                    values=[]
                    for method in METHODS:
                        sub=[r for r in source if r['role']==role and int(r['F'])==F
                             and r['outcome']==outcome and r['score']==score and r['method']==method]
                        assert len(sub)==2, (role,F,outcome,score,method,len(sub))
                        assert {r['direction'] for r in sub}=={'first_factor','interaction'}
                        assert all(int(r['n_trials'])==10000 and int(r['n'])==n for r in sub)
                        pcts=[Decimal(r['rejections'])*100/Decimal(r['n_trials']) for r in sub]
                        cell=f'{min(pcts):.2f}--{max(pcts):.2f}'
                        values.append(cell)
                        audit.append({'role':role,'F':F,'n':n,'outcome':outcome,'score':score,'method':method,
                                      'source_indices':[int(r['index']) for r in sub],
                                      'rejections':[int(r['rejections']) for r in sub],
                                      'cell':cell})
                    vals=[f'$({F},{n})$', 'Continuous' if outcome=='continuous' else 'Binary',
                          'U' if score=='unadjusted' else 'A']+values
                    rows.append(' & '.join(vals)+r' \\')
        rows.append(r'\addlinespace[5pt]')
    assert len(audit)==64
    notes=r'''U and A denote unadjusted scores and scores adjusted for factor main effects. Each range spans both prespecified effect directions at the displayed design, outcome and testing objective. SIGA-S, SIGA-R and CRT use reconstructed covariances. All designs use $p_{\rm bc}=.80$ and 10,000 trials per alternative. The corresponding null calibration is shown in Figure~\ref{fig:primary-null}; power differences are not adjusted to equal test size. Ranges are descriptive, not uncertainty intervals. Exact effects and all base/reconstructed rates are in Supplementary Appendix~G.'''
    table=(r'''\begin{table}[!htbp]
\centering
\caption{Primary power by outcome, adjustment and testing objective (\%).}
\label{tab:main-power-summary}
\begin{threeparttable}
\begin{singlespace}\small
\setlength{\tabcolsep}{3.5pt}
\renewcommand{\arraystretch}{1.15}
\begin{tabular}{@{}llcrrrr@{}}
\toprule
$(F,n)$ & Outcome & Score & RT & SIGA-S & SIGA-R & CRT\\
\midrule
'''+ '\n'.join(rows[:-1]) + '\n'+r'''\bottomrule
\end{tabular}
\begin{tablenotes}[flushleft]\small
\item[] '''+notes+r'''
\end{tablenotes}
\end{singlespace}
\end{threeparttable}
\end{table}''')
    return table,audit

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--write',action='store_true')
    args=p.parse_args()
    table, data=build()
    target=ROOT/'audit'/'regenerated_tables'/'primary_power_summary.tex'
    target.parent.mkdir(exist_ok=True)
    if args.write or not target.exists():
        target.write_text(table+'\n')
    else:
        assert target.read_text()==table+'\n', 'Regenerated standalone summary has changed'
    (ROOT/'audit'/'power_summary_cells.json').write_text(json.dumps(data,indent=2)+'\n')
    print('PASS: 64 power-range cells, each from two exact rejection counts; all 32 scenario-score power comparisons retained.')

if __name__=='__main__': main()
