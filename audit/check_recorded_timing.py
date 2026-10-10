#!/usr/bin/env python3
"""Validate and reproduce journal-article timing numbers; no trials are rerun."""
from collections import defaultdict
from pathlib import Path
import csv
import json
import math

ROOT=Path(__file__).resolve().parent.parent
SOURCE=ROOT/'timing'

def read(name):
    with (SOURCE/name).open(newline='',encoding='utf-8') as stream:
        return list(csv.DictReader(stream))

scenarios=read('timing_by_scenario.csv')
runs=read('timing_by_run.csv')
calibrations=read('calibration_times.csv')
assert len(scenarios)==52 and len(runs)==2 and len(calibrations)==7
assert all(s['complete']=='True' and int(s['trials_recorded'])==int(s['trials_expected']) for s in scenarios)
assert {r['run']:int(r['scenarios']) for r in runs}=={'main_d012_20261009':32,'additional_d012_20261009':20}
assert len({c['calibration_key'] for c in calibrations})==5
assert all(int(c['calibration_replicates'])==100000 for c in calibrations)
assert sum(int(s['trials_recorded']) for s in scenarios)==3040000
assert sum(int(r['trials_recorded']) for r in runs)==3040000

by_design=defaultdict(list)
for s in scenarios:
    key=(int(s['F']),int(s['sample_size']),float(s['p_bc']))
    by_design[key].append(s)
    n=int(s['trials_recorded'])
    t0=float(s['construction_seconds_sum'])
    t1=float(s['reference_seconds_sum'])
    t2=float(s['tail_seconds_sum'])
    assert n>0 and min(t0,t1,t2)>=0
    assert math.isclose(float(s['measured_seconds_sum']),t0+t1+t2,abs_tol=2e-5)
    assert math.isclose(float(s['reference_fraction_of_measured']),t1/(t0+t1+t2),abs_tol=1e-8)

result={'source_files':[p.name for p in SOURCE.glob('*times.csv')]+['timing_by_scenario.csv'],
        'run_checks':{}, 'design_reference_milliseconds_per_trial':{},
        'calibration_seconds_range':[round(min(float(c['measured_calibration_seconds']) for c in calibrations),6),round(max(float(c['measured_calibration_seconds']) for c in calibrations),6)],
        'distinct_calibration_settings':5, 'recorded_calibrations':7,
        'all_trials':3040000,'all_reference_sequences':3040000*4999}

for r in runs:
    mine=[s for s in scenarios if s['run']==r['run']]
    assert sum(int(s['trials_recorded']) for s in mine)==int(r['trials_recorded'])
    fields={'construction':'construction_seconds_sum','reference':'reference_seconds_sum','tail':'tail_seconds_sum'}
    for stage,field in fields.items():
        assert math.isclose(sum(float(s[field]) for s in mine),float(r[stage+'_measured_worker_seconds']),abs_tol=0.0002)
    assert math.isclose(float(r['reference_measured_worker_seconds'])/float(r['measured_worker_seconds_sum']),float(r['reference_fraction_of_measured']),abs_tol=1e-8)
    result['run_checks'][r['run']]={'trials':int(r['trials_recorded']),
        'reference_share_percent':round(100*float(r['reference_fraction_of_measured']),6),
        'construction_worker_seconds':float(r['construction_measured_worker_seconds']),
        'reference_worker_seconds':float(r['reference_measured_worker_seconds']),
        'tail_worker_seconds':float(r['tail_measured_worker_seconds']),
        'total_worker_seconds':float(r['measured_worker_seconds_sum'])}

for key, group in sorted(by_design.items()):
    n=sum(int(s['trials_recorded']) for s in group)
    ref=sum(float(s['reference_seconds_sum']) for s in group)
    label=f'F={key[0]},n={key[1]},pbc={key[2]:.2f}'
    result['design_reference_milliseconds_per_trial'][label]={
        'trials':n,'reference_milliseconds_mean':round(1000*ref/n,6),
        'range_by_scenario_ms':[round(min(float(s['reference_milliseconds_per_trial']) for s in group),6),round(max(float(s['reference_milliseconds_per_trial']) for s in group),6)]}

result['reference_fraction_all_percent']=round(100*sum(float(r['reference_measured_worker_seconds']) for r in runs)/sum(float(r['measured_worker_seconds_sum']) for r in runs),6)
result['worker_hours_all']=round(sum(float(r['measured_worker_seconds_sum']) for r in runs)/3600,6)
result['worker_hours_reference']=round(sum(float(r['reference_measured_worker_seconds']) for r in runs)/3600,6)

expected={'F=2,n=200,pbc=0.80':15.7,'F=5,n=400,pbc=0.80':34.8,'F=2,n=1000,pbc=0.80':69.0,'F=5,n=2000,pbc=0.80':151.5,'F=5,n=400,pbc=0.95':30.9}
for label,value in expected.items():
    assert abs(result['design_reference_milliseconds_per_trial'][label]['reference_milliseconds_mean']-value)<0.055,(label,value)
assert abs(result['run_checks']['main_d012_20261009']['reference_share_percent']-97.25)<0.01
assert abs(result['run_checks']['additional_d012_20261009']['reference_share_percent']-98.54)<0.01
path=SOURCE/'timing_audit.json'
path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(result,indent=2))
print('ALL CHECKS PASSED')
