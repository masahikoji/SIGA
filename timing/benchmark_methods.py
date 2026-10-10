#!/usr/bin/env python3
"""Read-only, standalone-method benchmark for the existing SIGA d012 program.

No production plan, calibration, checkpoint or analysis module is changed.
The current plan is loaded and checked by its own frozen-code/environment checks.
All timings use a single process and one numerical-library thread.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

# These must precede numpy/scipy/numba imports.
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
os.environ.setdefault('NUMBA_CACHE_DIR', str(Path.home()/'.cache'/'siga_timing_numba'))
os.environ.setdefault('PYTHONPYCACHEPREFIX', str(Path.home()/'.cache'/'siga_timing_python'))
sys.pycache_prefix = os.environ['PYTHONPYCACHEPREFIX']

METHODS = ('RT', 'original_S_normal', 'refined_S_normal', 'original_R_normal',
           'refined_R_normal', 'original_CRT_data', 'refined_CRT_data')


def save_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+'\n')


def save_csv(path, rows):
    if not rows:
        raise ValueError('No rows to save: '+str(path))
    with Path(path).open('w', newline='') as handle:
        writer=csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def choose_indices(plan, text):
    if text == 'all':
        return list(range(len(plan['scenarios'])))
    indices=list(dict.fromkeys(int(s) for s in text.split(',')))
    if not indices or min(indices)<0 or max(indices)>=len(plan['scenarios']):
        raise ValueError('Invalid zero-based scenario indices.')
    return indices


@dataclass
class Prepared:
    observation: object
    score: object
    scores: object
    statistic: float
    counts: object
    boundary: float
    tail: str


def prepare(obs, scenario, score_name):
    obs.validate()
    if len(scenario['boundaries']) != 1:
        raise ValueError('This benchmark supports one-boundary superiority/NI only.')
    b=float(scenario['boundaries'][0])
    w=obs.y-b*obs.assignment
    if score_name == 'unadjusted':
        score=w-w.mean()
    elif score_name == 'adjusted':
        c=np.column_stack((np.ones(len(w)),obs.patterns[obs.stratum]))
        score=w-c@np.linalg.lstsq(c,w,rcond=None)[0]
        score-=score.mean()
    else:
        raise ValueError(score_name)
    counts=np.bincount(obs.stratum,minlength=len(obs.patterns)).astype(float)
    return Prepared(obs,score,np.ascontiguousarray(score[:,None]),
                    float((2*obs.assignment-1)@score/2),counts,b,scenario['tails'][0])


def variances(prep, calibration, method):
    if method.startswith('original_'):
        sq=np.sqrt(prep.counts)
        omega=calibration['gamma']*sq[:,None]*sq[None,:]
    else:
        L=core.contrast_matrix(prep.observation.patterns)
        _,omega,_=core.covariance_constructions(prep.counts,calibration['gamma'],calibration['xi'],L)
    sv=core.sampling_variance(prep.score,prep.observation.stratum,prep.counts,omega)['variance']
    if method.endswith('_S_normal'):
        return sv, None
    dhat,_,_=core.effect_deviations(prep.observation,prep.boundary)
    rv=core.reference_variance(sv,prep.counts,dhat,calibration['psi'])['variance']
    return sv,rv


def references(prep, scenario, seed, trial_id, reference_draws):
    # Fresh reference generation for each call; do NOT reuse precomputed RT draws.
    rng=design.rng_for(seed,'reference',scenario['id'],int(trial_id))
    result=np.empty(reference_draws)
    for begin in range(0,reference_draws,256):
        m=min(256,reference_draws-begin)
        result[begin:begin+m]=allocation.reference_batch(
            prep.observation.stratum,prep.observation.patterns,scenario['pbc'],
            prep.scores,rng.random((m,scenario['n'])))[:,0]
    return result


def evaluate(prep, scenario, seed, trial_id, calibration, method, reference_draws):
    if method == 'RT':
        ref=references(prep,scenario,seed,trial_id,reference_draws)
        return core.mc_pvalue(ref,prep.statistic,prep.tail)
    sv,rv=variances(prep,calibration,method)
    if method.endswith('_S_normal'):
        return core.normal_pvalue(prep.statistic,sv,prep.tail)
    if method.endswith('_R_normal'):
        return core.normal_pvalue(prep.statistic,rv,prep.tail)
    if method.endswith('_CRT_data'):
        ref=references(prep,scenario,seed,trial_id,reference_draws)
        scale=float(np.sqrt(rv/sv)) if sv>0 and rv>0 else 1.0
        return core.mc_pvalue(ref,prep.statistic,prep.tail,scale,True)
    raise ValueError(method)


def measure(function, size, min_seconds):
    calls=0;checksum=0.0
    wall_start=time.perf_counter();cpu_start=time.process_time()
    while True:
        for j in range(size):
            value=function(j,calls)
            checksum+=float(value)
            calls+=1
        wall=time.perf_counter()-wall_start
        if wall>=min_seconds:
            break
    cpu=time.process_time()-cpu_start
    if not np.isfinite(checksum):
        raise FloatingPointError('Nonfinite benchmark output.')
    return dict(calls=calls,elapsed_seconds=wall,process_cpu_seconds=cpu,
                seconds_per_analysis=wall/calls,cpu_seconds_per_analysis=cpu/calls)


def benchmark(args):
    global np, core, design, allocation, experiment
    code=Path(args.code).expanduser().resolve();run=Path(args.run).expanduser().resolve()
    if not (code/'siga_validation'/'experiment.py').is_file():
        raise FileNotFoundError('Expected the existing SIGA code at '+str(code))
    sys.path.insert(0,str(code))
    import numpy as np
    from siga_validation import core, design, allocation, experiment
    plan=experiment.load_plan(run)  # includes frozen source and environment checks
    if any(m not in experiment.METHODS for m in METHODS):
        raise RuntimeError('Unexpected implementation method names.')
    indices=choose_indices(plan,args.scenarios)
    out=Path(args.out).expanduser().resolve()
    if out==run or run in out.parents or out==code or code in out.parents:
        raise ValueError('Use a separate output directory, outside the production run and code.')
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('Use a new, empty benchmark output directory: '+str(out))
    out.mkdir(parents=True,exist_ok=True)
    shutil.copy2(run/'plan.json',out/'source_plan.json')
    B=int(plan['configuration']['reference'])
    scenarios=[plan['scenarios'][i] for i in indices]
    scope=dict(created_utc=datetime.now(timezone.utc).isoformat(),kind='standalone_method_timing',
               code_path=str(code),source_run=str(run),source_plan_hash=plan['plan_hash'],
               source_signature=experiment.source_signature(),numerical_environment=experiment.numerical_environment(),
               platform=platform.platform(),processor=platform.processor(),logical_cpu_count=os.cpu_count(),
               one_process=True,numerical_threads=1,
               reference_draws=B,calibration_triples=int(plan['configuration']['calibration']),
               scenario_indices=indices,trials_per_timing_block=args.trials,repetitions=args.repeats,
               minimum_seconds_per_block=args.min_seconds,calibration_repetitions=args.calibration_repeats,
               generating_trial_start=2000000000,comparison_methods=list(METHODS),
               excluded=['data generation','disk IO','compilation warm-up'],
               score_preparation='timed separately and added for end-to-end projections',
               calibration_cost='common full three-sequence bundle; conservatively also charged to SIGA-S',
               note='Timing projections only. No operating-characteristic result is generated or replaced.')
    save_json(out/'scope.json',scope)
    # Load, but never overwrite, original calibration files.
    calibrations={}
    keys=list(dict.fromkeys(s['calibration_key'] for s in scenarios))
    for key in keys:
        cal=experiment.calibration_for(run,key)
        if str(cal['plan_hash'])!=plan['plan_hash']:
            raise RuntimeError('Calibration/plan mismatch.')
        if int(cal['B0'])!=int(plan['configuration']['calibration']):
            raise RuntimeError('Unexpected calibration budget.')
        calibrations[key]=cal
    # Warm up allocation, normal calculations and geometry; verify standalone
    # results against the unchanged original evaluation route on common trials.
    validations=[];prepared={};observations={}
    for s in scenarios:
        obs=[design.generate_observation(s,plan['seed'],2000000000+j) for j in range(args.trials)]
        observations[s['index']]=obs
        original=experiment.trial_evaluation(s,plan,2000000000,calibrations[s['calibration_key']])
        for k,score_name in enumerate(('unadjusted','adjusted')):
            items=[prepare(o,s,score_name) for o in obs];prepared[(s['index'],score_name)]=items
            for method in METHODS:
                got=evaluate(items[0],s,plan['seed'],2000000000,calibrations[s['calibration_key']],method,B)
                expected=float(original['p'][k,experiment.METHODS.index(method)])
                error=abs(got-expected)
                if error>5e-12:
                    raise AssertionError(f'Standalone mismatch {s["index"]}, {score_name}, {method}: {got} vs {expected}')
                validations.append(dict(index=s['index'],score=score_name,method=method,max_abs_p_difference=error))
    save_csv(out/'implementation_validation.csv',validations)
    print(f'Validated {len(validations)} standalone p-values against the unchanged original implementation.',flush=True)
    # Re-measure the FULL calibration bundle after warming compilation. Copies
    # of the frozen plan and fresh calibration files live only under the new output.
    cal_rows=[]
    for rep in range(args.calibration_repeats):
        temp=out/'calibration_remeasure'/f'repetition_{rep+1}'
        temp.mkdir(parents=True);shutil.copy2(run/'plan.json',temp/'plan.json')
        for key in keys:
            s=next(s for s in scenarios if s['calibration_key']==key)
            pat,pi=design.patterns_and_prob(s['F'],s['profile_law'])
            rng=np.random.default_rng(8910+rep)
            sid=np.searchsorted(np.cumsum(pi),rng.random((2,s['n']))).astype(np.int64)
            allocation.calibration_batch(sid,pat,s['pbc'],rng.random((2,3,s['n'])))
            message=experiment.calibrate_one((str(temp),key))
            with np.load(temp/'calibration'/f'{key}.npz',allow_pickle=False) as z:
                seconds=float(z['elapsed_seconds'])
            cal_rows.append(dict(calibration_key=key,F=s['F'],n=s['n'],pbc=s['pbc'],
                                 replicate=rep+1,triples=plan['configuration']['calibration'],seconds=seconds))
            print('Calibration benchmark:',message,flush=True)
    save_csv(out/'calibration_times.csv',cal_rows)
    timing_rows=[]
    for s in scenarios:
        cal=calibrations[s['calibration_key']]
        for score_name in ('unadjusted','adjusted'):
            items=prepared[(s['index'],score_name)]
            for rep in range(args.repeats):
                names=('score_preparation',)+METHODS
                shift=rep%len(names);names=names[shift:]+names[:shift]
                for method in names:
                    if method=='score_preparation':
                        fn=lambda j,c:prepare(observations[s['index']][j],s,score_name).statistic
                    else:
                        fn=lambda j,c,m=method:evaluate(items[j],s,plan['seed'],2100000000+rep*10000000+c,cal,m,B)
                    times=measure(fn,len(items),args.min_seconds)
                    timing_rows.append(dict(index=s['index'],id=s['id'],F=s['F'],n=s['n'],outcome=s['outcome'],
                        pbc=s['pbc'],role=s['role'],direction=s['direction'],score=score_name,method=method,
                        repetition=rep+1,reference_draws=B,**times))
            save_csv(out/'timing_repetitions.csv',timing_rows)
            print('Timed scenario',s['index'],score_name,flush=True)
    summaries=[]
    for s in scenarios:
        ck=s['calibration_key']
        cal_cost=float(np.median([r['seconds'] for r in cal_rows if r['calibration_key']==ck]))
        for score_name in ('unadjusted','adjusted'):
            candidates=[r for r in timing_rows if r['index']==s['index'] and r['score']==score_name]
            score_cost=float(np.median([r['seconds_per_analysis'] for r in candidates if r['method']=='score_preparation']))
            for method in METHODS:
                rows=[r for r in candidates if r['method']==method]
                values=np.array([r['seconds_per_analysis'] for r in rows])
                cost=float(np.median(values));setup=0. if method=='RT' else cal_cost
                target=int(s['outer'])
                summaries.append(dict(index=s['index'],id=s['id'],F=s['F'],n=s['n'],outcome=s['outcome'],pbc=s['pbc'],
                    role=s['role'],direction=s['direction'],score=score_name,method=method,reference_draws=B,
                    repetitions=args.repeats,median_analysis_seconds=cost,min_analysis_seconds=float(values.min()),
                    max_analysis_seconds=float(values.max()),common_score_seconds=score_cost,
                    common_bundle_setup_seconds=setup,target_trials=target,
                    projected_100000_analysis_minutes=cost*100000/60,
                    projected_100000_with_setup_minutes=(setup+100000*cost)/60,
                    projected_target_with_setup_minutes=(setup+target*cost)/60,
                    projected_target_with_score_and_setup_minutes=(setup+target*(cost+score_cost))/60))
    save_csv(out/'timing_summary.csv',summaries)
    save_json(out/'STATUS.json',dict(complete=True,kind='benchmark_only',source_plan_hash=plan['plan_hash'],
        scenarios=len(scenarios),summary_rows=len(summaries),production_simulations_changed=False,
        max_implementation_p_difference=max(r['max_abs_p_difference'] for r in validations)))
    print('Complete. Upload timing_summary.csv, calibration_times.csv and scope.json from:',out,flush=True)
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--code',required=True,help='Existing, unchanged SIGA_M3Ultra_d012_20261009 code directory')
    parser.add_argument('--run',default=str(Path.home()/'SIGA_runs'/'main_d012_20261009'))
    parser.add_argument('--out',default=str(Path.home()/'SIGA_runs'/('timing_d012_'+datetime.now().strftime('%Y%m%d_%H%M%S'))))
    parser.add_argument('--scenarios',default='0,1,8,9,16,17,24,25',help='Zero-based indices or all')
    parser.add_argument('--trials',type=int,default=30)
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--calibration-repeats',type=int,default=3)
    parser.add_argument('--min-seconds',type=float,default=.25)
    args=parser.parse_args()
    if args.trials<2 or args.repeats<1 or args.calibration_repeats<1 or args.min_seconds<=0:
        parser.error('Trials >=2, repetitions >=1 and minimum seconds >0 are required.')
    try:
        benchmark(args)
    except (ValueError,FileNotFoundError,RuntimeError,AssertionError,FileExistsError) as exc:
        parser.exit(1,'ERROR: '+str(exc)+'\n')

if __name__=='__main__':
    main()
