from __future__ import annotations
import contextlib
import copy
import importlib.util
import inspect
import io
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from revision import engine
from revision.utils import rng_for,stable_hash,verify_vendor,load_plan
from revision.designs import make_scenarios,resolve_scenario,patterns_prob,maximum_ratio_direction
from revision.calibration import simulate_calibration,CovAccumulator,resolve_and_seal,calibration_jobs,calibrate_job,load_resolved,load_calibration
from revision.runner import initialise,_run_job
from revision.reporting import make_report,paired_summary,ratio_summary,load_scenario_results
from siga_validation.allocation import path_from_uniforms
from siga_validation.core import scores_and_variances,mc_pvalue,normal_pvalue,covariance_constructions,contrast_matrix


class FocusedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c2=simulate_calibration(2,40,.95,500,56,['test-cal'])
        v,cls.diag=maximum_ratio_direction(cls.c2['psi'],cls.c2['profile_prob'])
        s=make_scenarios('gap')[0];s['n']=40;s['calibration_key']='test'
        cls.sc=resolve_scenario(s,{'test':dict(direction=v.tolist())})
        cls.obs=__import__('revision.designs',fromlist=['generate_observation']).generate_observation(cls.sc,457,3)
        cls.entries,_=scores_and_variances(cls.obs,[cls.sc['boundary']],cls.c2,noise_adjustment=True)
        cls.refs=engine.reference_statistics(cls.obs,.95,cls.entries,789,cls.sc['id'],3,99)

    def test_01_vendor_bytes_unchanged(self):
        self.assertEqual(verify_vendor()['commit'],'76daf34906b365dc787e1cbeeb4114c89e558b0d')
    def test_02_scope(self):
        self.assertEqual(len(make_scenarios()),5)
        self.assertEqual(len(make_scenarios('gap')),4)
        self.assertEqual(len(make_scenarios('a13')),1)
        with self.assertRaises(ValueError):make_scenarios('core')
    def test_03_production_budget(self):
        with tempfile.TemporaryDirectory() as d:
            plan,work=initialise(Path(d)/'out','full','focused',123)
            self.assertEqual(work['outcome_trials'],320000)
            self.assertEqual(work['extra_reference_paths'],10000000)
            self.assertEqual(work['calibration_triples'],900000)
            self.assertEqual(len(calibration_jobs(plan,Path(d)/'out')),7)
            self.assertEqual(plan['configuration']['high_reference'],9999)
    def test_04_mean_zero_direction(self):
        v,di=maximum_ratio_direction(self.c2['psi'],self.c2['profile_prob'])
        self.assertAlmostEqual(self.c2['profile_prob']@v,0.,places=14)
        self.assertAlmostEqual(np.max(np.abs(v)),1.)
    def test_05_marginal_effects(self):
        for s in make_scenarios():
            v,_=maximum_ratio_direction(np.eye(2**s['F']),patterns_prob(s['F'])[1])
            sc=resolve_scenario(s,{s['calibration_key']:dict(direction=v.tolist())})
            self.assertAlmostEqual(np.array(sc['profile_prob'])@(np.array(sc['mu1'])-sc['mu0']),s['delta'],places=12)
    def test_06_fixed_power_shift(self):
        for s in make_scenarios('gap'):
            expected=0 if s['null_evaluation'] else 1.5/np.sqrt(s['n'])
            self.assertAlmostEqual(s['delta']-s['boundary'],expected)
    def test_07_a13_zero_deviation(self):
        sc=resolve_scenario(make_scenarios('a13')[0],{})
        np.testing.assert_array_equal(np.asarray(sc['mu1'])-sc['mu0'],np.zeros(32))
        self.assertEqual(sc['sigma'],.25);self.assertEqual(sc['eta'],.25)
    def test_08_independent_rng_streams(self):
        a=rng_for(2,'outer','a',1).random(20)
        np.testing.assert_array_equal(a,rng_for(2,'outer','a',1).random(20))
        self.assertFalse(np.array_equal(a,rng_for(2,'reference','a',1).random(20)))
        self.assertFalse(np.array_equal(rng_for(2,'calibration','analysis').random(20),rng_for(2,'calibration','direction').random(20)))
    def test_09_covariance_accumulator(self):
        x=np.random.default_rng(1).normal(size=(301,4));acc=CovAccumulator(4)
        for a in np.array_split(x,7): acc.add(a)
        np.testing.assert_allclose(acc.covariance(),np.cov(x,rowvar=False),atol=1e-14)
    def test_10_nested_reference_draws(self):
        refs=engine.reference_statistics(self.obs,.95,self.entries,789,self.sc['id'],3,49)
        np.testing.assert_array_equal(refs,self.refs[:49])
    def test_11_principal_pvalues_match_original(self):
        p,_=engine.feasible_from_entries(self.entries,self.refs,'upper')
        for i,e in enumerate(self.entries):
            for j,tag in enumerate(engine.UPSTREAM_TAGS):
                vs=e[tag]['S']['variance'];vr=e[tag]['R']['variance'];t=e['statistic']
                expected=[mc_pvalue(self.refs[:,i],t,'upper'),normal_pvalue(t,vs,'upper'),
                    normal_pvalue(t,vr,'upper'),mc_pvalue(self.refs[:,i],t,'upper',np.sqrt(vr/vs),True)]
                np.testing.assert_array_equal(p[i,j,:4],expected)
    def test_12_truth_is_diagnostic_only(self):
        p,m,e=engine.evaluate_with_truth(self.obs,self.entries,self.refs,self.c2,self.sc)
        fake=copy.deepcopy(self.sc);fake['mu1']=(np.array(fake['mu1'])+.03).tolist()
        pp,mm,ee=engine.evaluate_with_truth(self.obs,self.entries,self.refs,self.c2,fake)
        np.testing.assert_array_equal(p[:,:,:4],pp[:,:,:4])
        self.assertNotIn('sc',inspect.signature(engine.feasible_from_entries).parameters)
    def test_13_exact_correction_decomposition(self):
        _,metrics,eff=engine.evaluate_with_truth(self.obs,self.entries,self.refs,self.c2,self.sc)
        d=dict(zip(engine.EFFECT_METRICS,eff));e=self.entries[0]
        qhat=e['original']['R']['raw']-e['original']['S']['variance']
        self.assertAlmostEqual(qhat-d['oracle_correction'],d['missing_component']+d['linear_noise_component']+d['quadratic_noise_component'],places=10)
    def test_14_oracle_zero_CRT_equals_RT(self):
        sc=resolve_scenario(make_scenarios('a13')[0],{});sc['n']=40
        cal=simulate_calibration(5,40,.95,500,56,['test-cal5'])
        obs=__import__('revision.designs',fromlist=['generate_observation']).generate_observation(sc,88,1)
        en,_=scores_and_variances(obs,[0.],cal,noise_adjustment=True)
        refs=engine.reference_statistics(obs,.95,en,99,sc['id'],1,49)
        p,_,_=engine.evaluate_with_truth(obs,en,refs,cal,sc)
        np.testing.assert_array_equal(p[:,:,5],p[:,:,0])
    def test_15_no_comparator_or_MC_CRT(self):
        self.assertEqual(engine.PRINCIPAL,('RT','SIGA_S','SIGA_R','CRT_plugin'))
        self.assertNotIn('CRT_MC',engine.METHODS)
        self.assertFalse(hasattr(engine,'stratified_yys'))
    def test_16_paired_error(self):
        a=np.array([.01,.01,.9,.9]);b=np.array([.01,.9,.01,.9])
        out=paired_summary(a,b,.05)
        self.assertEqual(out['difference'],0.)
        self.assertAlmostEqual(out['paired_mcse'],np.std([0,1,-1,0],ddof=1)/2)
        self.assertEqual(out['A_only'],1);self.assertEqual(out['B_only'],1)
    def test_17_no_overwrite_and_plan_guard(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'out';initialise(root,'smoke','gap',8)
            with self.assertRaises(FileExistsError): initialise(root,'smoke','gap',8)
            plan=json.loads((root/'plan.json').read_text());plan['seed']=99
            (root/'plan.json').write_text(json.dumps(plan))
            with self.assertRaises(RuntimeError): load_plan(root)
    def test_18_reconstruction_preserves_balance(self):
        counts=np.array([10.,8.,12.,10.]);L=contrast_matrix(self.obs.patterns)
        _,ref,geom=covariance_constructions(counts,self.c2['gamma'],self.c2['xi'],L)
        self.assertFalse(geom['fallback'])
        np.testing.assert_allclose(L@ref@L.T,self.c2['xi'],atol=1e-10)
    def test_19_reference_tail_ties(self):
        x=np.array([-2.,-1.,0.,1.,2.])
        self.assertAlmostEqual(mc_pvalue(x,1.,'two',1.1,True),5/6)
        self.assertAlmostEqual(mc_pvalue(x,1.,'two',1.,True),mc_pvalue(x,1.,'two'))
    def test_20_allocation_matches_absolute_range(self):
        pat,_=patterns_prob(2);rng=np.random.default_rng(91);sid=rng.integers(0,4,100);u=rng.random(100)
        z=path_from_uniforms(sid,pat,.8,u);expected=[];overall=0;D=np.zeros((2,2),int)
        for s,x in zip(sid,u):
            plus=abs(overall+1)+sum(abs(D[f,pat[s,f]]+1) for f in range(2))
            minus=abs(overall-1)+sum(abs(D[f,pat[s,f]]-1) for f in range(2))
            p=.8 if plus<minus else .2 if plus>minus else .5
            a=1 if x<p else -1;expected.append(a);overall+=a
            for f in range(2):D[f,pat[s,f]]+=a
        np.testing.assert_array_equal(z,expected)

if __name__=='__main__':unittest.main()
