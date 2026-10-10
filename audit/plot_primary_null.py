#!/usr/bin/env python3
"""Plot all primary null rows, using supplied rates and Wilson intervals."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
df=pd.read_csv(ROOT/'data/primary/rejection_rates.csv')
df=df[df.evaluation.eq('type1')].copy()
methods=['RT','refined_S_normal','refined_R_normal','refined_CRT_data']
labels=['RT','SIGA-S, reconstructed','SIGA-R, reconstructed','CRT, reconstructed']
indices=sorted(df['index'].unique())
keys=[(i,s) for test in ['superiority_null','noninferiority_null'] for i in indices if df[df['index'].eq(i)].role.iloc[0]==test for s in ['unadjusted','adjusted']]
fig,ax=plt.subplots(figsize=(8.0,9.1))
y=np.arange(len(keys))
for j,(m,label) in enumerate(zip(methods,labels)):
    rows=pd.DataFrame([df[(df['index']==i)&(df.score==s)&(df.method==m)].iloc[0] for i,s in keys])
    rate=100*(rows.rate.to_numpy()-rows.alpha.to_numpy())
    err=np.array([100*(rows.rate-rows.ci_lower),100*(rows.ci_upper-rows.rate)])
    ax.errorbar(rate,y+(j-1.5)*.19,xerr=err,fmt=['o','s','^','D'][j],markersize=3,linewidth=.8,capsize=1.5,label=label)
ax.axvline(0,linestyle='--',linewidth=.8)
ax.set_yticks(y)
ax.set_yticklabels([f'P{i+1:02d}  '+('U' if s=='unadjusted' else 'A') for i,s in keys],fontsize=8)
ax.invert_yaxis()
ax.set_ylim(len(keys)-.3,-2.0)
ax.set_xlim(-1.12,.67)
ax.set_xlabel('Rejection probability minus nominal level (percentage points)',fontsize=10)
ax.set_ylabel('Scenario and score',fontsize=10)
ax.tick_params(axis='x',labelsize=9)
ax.grid(axis='x',linestyle=':',linewidth=.5,alpha=.5)
ax.legend(loc='upper left',bbox_to_anchor=(0,1.075),ncol=2,frameon=False,fontsize=9)
fig.subplots_adjust(left=.14,right=.98,bottom=.085,top=.925)
(ROOT/'figures').mkdir(exist_ok=True)
fig.savefig(ROOT/'figures/primary_null_rates.pdf')
fig.savefig(ROOT/'figures/primary_null_rates.png',dpi=160)
plotdata=df[df.method.isin(methods)].copy()
plotdata['paper_id']=plotdata['index'].map(lambda x:f'P{x+1:02d}')
plotdata.to_csv(ROOT/'audit/primary_null_figure_data.csv',index=False)
print('Figure saved. Rows:',len(plotdata),'settings:',len(keys))
