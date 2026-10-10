"""Tests of the requested fixed-amplitude revision and execution safeguards."""
import json
import tempfile
from pathlib import Path
import numpy as np
from .design import make_scenario,stable_id
from .resimulation import make_resimulation_plan,freeze
from .preflight import audit_plan
from .experiment import load_plan


def test_fixed_binary_deviations():
    for suites in ('main','controls','sample_size','all'):
        plan=make_resimulation_plan(suites)
        for sc in plan['scenarios']:
            if sc['outcome']!='binary': continue
            pi=np.array(sc['profile_prob']);p0=np.array(sc['control_risks']);p1=np.array(sc['treatment_risks'])
            target=0. if sc['sharp_binary'] else .12
            assert sc['feasibility_scale']==1.
            np.testing.assert_allclose(np.max(np.abs(p1-p0-sc['delta'])),target,atol=1e-12)
            np.testing.assert_allclose(pi@(p1-p0),sc['delta'],atol=1e-12)
            assert np.all((p0>.02)&(p0<.98)&(p1>.02)&(p1<.98))


def test_no_silent_shrink():
    args=(2,200,'binary',.8,1.,.25,.15,'first_factor','superiority_power','test')
    try: make_scenario(*args,binary_policy='error')
    except ValueError as exc:
        assert 'Automatic shrinking and clipping are disabled' in str(exc)
    else: raise AssertionError('The infeasible 0.15 design was silently accepted.')
    valid=list(args);valid[6]=.12
    sc=make_scenario(*valid,binary_policy='error')
    assert sc['actual_amplitude']==.12 and sc['feasibility_scale']==1.


def test_preflight_outputs():
    with tempfile.TemporaryDirectory() as t:
        root=Path(t)/'run';freeze(root,make_resimulation_plan())
        result=audit_plan(root)
        assert result['passed'] and result['binary_scenarios']==16
        assert result['outer_trials']==1760000 and result['reference_paths']==8798240000
        assert (root/'binary_probability_check.csv').exists()
        assert result['maximum_binary_treatment_probability']<.98


def test_environment_change_rejected():
    with tempfile.TemporaryDirectory() as t:
        root=Path(t)/'run';freeze(root,make_resimulation_plan(profile='smoke'))
        p=json.loads((root/'plan.json').read_text());p.pop('plan_hash')
        p['numerical_environment']['numpy']='incompatible-test-version'
        p['plan_hash']=stable_id(p);(root/'plan.json').write_text(json.dumps(p))
        try: load_plan(root)
        except RuntimeError as exc: assert 'Numerical environment changed' in str(exc)
        else: raise AssertionError('Different numerical environment was accepted.')


def run_d012_tests():
    for test in (test_fixed_binary_deviations,test_no_silent_shrink,test_preflight_outputs,
                 test_environment_change_rejected):
        test();print('PASS',test.__name__,flush=True)
