#!/usr/bin/env python3
"""Rebuild all manuscript computational-timing values from the two frozen M3 Ultra benchmark batches.

Run from any directory: python audit/verify_complete_benchmark.py
Only source files inside this manuscript package are used.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import json

root=Path(__file__).resolve().parents[1]
lookup={'RT':'RT','original_R_normal':'SIGA-R','refined_R_normal':'SIGA-R',
        'original_S_normal':'SIGA-S','refined_S_normal':'SIGA-S',
        'original_CRT_data':'CRT','refined_CRT_data':'CRT'}
indexes={'main':[0,1,8,9,16,17,24,25], 'additional':[4,5,8,9,12,13,14,15]}
get_case=lambda x: ('(2,200) p=.80' if x.F==2 and x.n==200 else
                    '(5,400) p=.80' if x.F==5 and x.n==400 and x.pbc==.8 else
                    '(2,1000) p=.80' if x.F==2 and x.n==1000 else
                    '(5,2000) p=.80' if x.F==5 and x.n==2000 else
                    '(5,400) p=.95' if x.F==5 and x.n==400 and x.pbc==.95 else 'UNKNOWN')
case_order=['(2,200) p=.80','(5,400) p=.80','(2,1000) p=.80','(5,2000) p=.80','(5,400) p=.95']
methods=['RT','SIGA-R','SIGA-S','CRT']
records=[];crossovers=[];spreads=[];cals=set();num_validations=0
for batch in indexes:
 d=root/'timing'/f'standalone_{batch}'
 assert d.is_dir(), d
 summ=pd.read_csv(d/'timing_summary.csv'); reps=pd.read_csv(d/'timing_repetitions.csv')
 val=pd.read_csv(d/'implementation_validation.csv'); cal=pd.read_csv(d/'calibration_times.csv')
 scope=json.loads((d/'scope.json').read_text()); plan=json.loads((d/'source_plan.json').read_text()); status=json.loads((d/'STATUS.json').read_text())
 assert scope['source_plan_hash']==plan['plan_hash']==status['source_plan_hash']
 assert scope['source_signature']==plan['source_signature']
 assert scope['scenario_indices']==indexes[batch] and scope['one_process'] and scope['numerical_threads']==1
 assert (scope['repetitions'],scope['calibration_repetitions'],scope['trials_per_timing_block'],scope['reference_draws'],scope['calibration_triples'])==(5,3,100,4999,100000)
 assert len(summ)==len(val)==112 and len(reps)==640
 assert status['complete'] and not status['production_simulations_changed'] and status['max_implementation_p_difference']==0.
 assert val['max_abs_p_difference'].max()==0
 assert (summ.groupby(['index','score']).size()==7).all()
 assert (reps.groupby(['index','score','method']).size()==5).all()
 assert (reps['calls']>=100).all() and (reps['elapsed_seconds']>=1).all()
 assert np.allclose(reps['seconds_per_analysis'],reps['elapsed_seconds']/reps['calls'],rtol=5e-13)
 med_cal=cal.groupby('calibration_key')['seconds'].median()
 assert (cal.groupby('calibration_key')['replicate'].apply(lambda x:sorted(x.tolist())==[1,2,3])).all()
 assert (cal['triples']==100000).all()
 cals.update(med_cal.keys())
 times=reps.pivot_table(index=['index','score','method'],values='seconds_per_analysis',aggfunc='median')['seconds_per_analysis']
 costs=reps[reps.method=='score_preparation'].groupby(['index','score'])['seconds_per_analysis'].median()
 for _,r in summ.iterrows():
  i=int(r['index']);sc=plan['scenarios'][i]
  assert r['id']==sc['id'] and r['F']==sc['F'] and r['n']==sc['n'] and r['outcome']==sc['outcome'] and r['role']==sc['role']
  assert r['pbc']==sc['pbc'] and r['method'] in lookup
  assert abs(r['median_analysis_seconds']-times.loc[(i,r.score,r.method)])<1e-10
  assert abs(r['common_score_seconds']-float(costs.loc[(i,r.score)]))<1e-10
  assert abs(r['common_bundle_setup_seconds']-(float(med_cal.loc[sc['calibration_key']]) if r['method']!='RT' else 0))<1e-10
  if r['method']!='RT':
   rt_time=times.loc[(i,r.score,'RT')]
   if 'normal' in r['method']:
    assert rt_time>r['median_analysis_seconds']
    crossovers.append(r['common_bundle_setup_seconds']/(rt_time-r['median_analysis_seconds']))
  projected=(100000*(r['median_analysis_seconds']+r['common_score_seconds'])+r['common_bundle_setup_seconds'])/60
  assert abs(r['projected_100000_analysis_minutes']-100000*r['median_analysis_seconds']/60)<1e-8
  assert abs(r['projected_target_with_score_and_setup_minutes']-(r['target_trials']*(r['median_analysis_seconds']+r['common_score_seconds'])+r['common_bundle_setup_seconds'])/60)<1e-8
  records.append({'batch':batch,'index':i,'F':r['F'],'n':r['n'],'pbc':r['pbc'],'case':get_case(r),'outcome':r['outcome'],'score':r.score,'role':r.role,'method':r.method,'procedure':lookup[r.method],'100k_projected_minutes':projected})
 v=reps[reps.method!='score_preparation'].groupby(['index','score','method'])['seconds_per_analysis'].agg(['min','median','max'])
 spreads.extend(((v['max']-v['min'])/v['median']).tolist())
 num_validations+=len(val)
assert len(records)==224 and num_validations==224 and len(cals)==5 and set(x['case'] for x in records)==set(case_order)
assert [int(np.floor(min(crossovers))+1),int(np.floor(max(crossovers))+1)]==[78,102]
tab=pd.DataFrame(records).groupby(['case','procedure'])['100k_projected_minutes'].agg(['min','max','count'])
# Compare the recomputed projections with the separately recorded merged benchmark export.
archived=pd.read_csv(root/'timing'/'combined_benchmark_computed.csv')
assert len(archived)==224 and len(records)==224
recal=pd.DataFrame(records)
for rec in recal.to_dict('records'):
    matches=archived[(archived['run']==('primary' if rec['batch']=='main' else rec['batch'])) &
                     (archived['index']==rec['index']) &
                     (archived['score']==rec['score']) &
                     (archived['method']==rec['method'])]
    assert len(matches)==1, rec
    value=float(matches.iloc[0]['projected_100k_with_score_and_setup_minutes'])
    assert abs(value-rec['100k_projected_minutes'])<1e-8, rec
print('Verified',len(records),'method-score-scenario combinations,',num_validations,'unchanged-code p-value checks')
print('Distinct calibration settings:',len(cals),'five-repeat dispersion median:',np.median(spreads),'max:',max(spreads))
print('Gaussian crossover integers:',int(np.floor(min(crossovers))+1),'through',int(np.floor(max(crossovers))+1))
for case in case_order:
 print('\n'+case)
 for m in methods:
  r=tab.loc[(case,m)]
  print(f' {m:7}: {r["min"]:8.4f}--{r["max"]:8.4f} projected min / 100000 analyses ({int(r["count"])} configurations)')
print('\nBenchmark projections, repetition records and timing scope verified without manuscript files.')
