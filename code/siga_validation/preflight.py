"""Deterministic checks of the frozen design, not simulation results."""
from pathlib import Path
import numpy as np
from .experiment import load_plan, atomic_json
from .resimulation import BINARY_AMPLITUDE, work_estimate
from .report import write_csv


def audit_plan(root):
    root=Path(root); plan=load_plan(root)
    if plan.get('binary_amplitude') != BINARY_AMPLITUDE or plan.get('binary_feasibility') != 'error':
        raise ValueError('This is not a fixed d012, no-shrinking plan.')
    cfg=plan['configuration']
    if plan['profile']=='production':
        if (cfg['outer_null'],cfg['outer_power'],cfg['reference'])!=(100000,10000,4999):
            raise ValueError('Production replication counts are inconsistent.')
    rows=[]; summary=[]
    for sc in plan['scenarios']:
        pi=np.asarray(sc['profile_prob']); d=np.asarray(sc['dgp_deviations'])
        if not np.all(np.isfinite(d)) or abs(pi@d)>1e-12:
            raise ValueError('Invalid/uncentred DGP effect deviations.')
        expected=cfg['outer_null'] if sc['null_evaluation'] else cfg['outer_power']
        if sc['outer']!=expected: raise ValueError('Wrong scenario-specific replication count.')
        if sc['outcome']!='binary': continue
        p0=np.asarray(sc['control_risks']); p1=np.asarray(sc['treatment_risks'])
        amplitude=0.0 if sc['sharp_binary'] else BINARY_AMPLITUDE
        if sc['feasibility_scale']!=1.0 or abs(sc['actual_amplitude']-amplitude)>1e-12:
            raise ValueError('Binary heterogeneity has been changed.')
        np.testing.assert_allclose(p1-p0,sc['delta']+d,atol=1e-12,rtol=0)
        np.testing.assert_allclose(pi@(p1-p0),sc['delta'],atol=1e-12,rtol=0)
        if not np.all((p0>.02)&(p0<.98)&(p1>.02)&(p1<.98)):
            raise ValueError('A binary probability is outside (0.02,0.98).')
        if abs(np.max(np.abs((p1-p0)-sc['delta']))-amplitude)>1e-12:
            raise ValueError('The requested deviation from the generating mean is not preserved.')
        base=dict(index=sc['index'],F=sc['F'],n=sc['n'],direction=sc['direction'],
                  role=sc['role'],delta=sc['delta'],amplitude=amplitude,
                  feasibility_scale=sc['feasibility_scale'])
        summary.append(dict(**base,min_control=float(p0.min()),max_control=float(p0.max()),
                            min_treatment=float(p1.min()),max_treatment=float(p1.max()),
                            weighted_mean_deviation=float(pi@d)))
        for s in range(len(pi)):
            rows.append(dict(**base,stratum=s,profile_probability=float(pi[s]),
                             control_probability=float(p0[s]),treatment_probability=float(p1[s]),
                             deviation_from_generating_mean=float(d[s])))
    if rows:
        write_csv(root/'binary_probability_check.csv',rows)
        write_csv(root/'binary_probability_summary.csv',summary)
    report=dict(passed=True,plan_hash=plan['plan_hash'],binary_amplitude=BINARY_AMPLITUDE,
                sharp_null_amplitude=0.0,binary_scenarios=len(summary),
                automatic_shrinking=False,**work_estimate(plan))
    if summary:
        report['minimum_binary_treatment_probability']=min(s['min_treatment'] for s in summary)
        report['maximum_binary_treatment_probability']=max(s['max_treatment'] for s in summary)
    atomic_json(root/'preflight.json',report)
    return report
