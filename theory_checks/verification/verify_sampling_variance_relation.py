#!/usr/bin/env python3
"""Exact algebraic comparison with Bugni, Canay and Shaikh (2018), eqs. 15-18.
This checks the variance identity on the overlapping unadjusted 1:1 setting.
It is not a test of the probabilistic assumptions or a Monte Carlo simulation.
"""
from pathlib import Path
import json
import sympy as sp

p, t, m, d, s0, s1, c = sp.symbols('p t m d sigma0_sq sigma1_sq cov01', real=True)
# Per-stratum contribution: our sign imbalance is twice their centred assignment.
sigma_d = 4*p*t
nu_q = p*(s0+s1+2*c)/4
delta_sq = p*(s0+s1-2*c+d*d)
sixteen_vs = 4*m*m*sigma_d + 4*nu_q + delta_sq
# BCS residual-outcome, stratum-effect and imbalance terms at target proportion 1/2.
bcs = 2*p*(s0+s1) + p*d*d + 16*p*t*m*m
difference = sp.expand(sixteen_vs-bcs)
assert difference == 0
assert sp.diff(sp.expand(sixteen_vs),c) == 0
result = {'target_allocation_fraction':'1/2','imbalance_scaling':'D_SIG A = 2 * D_BCS'.replace('SIG A','SIGA'),
          'variance_scale':'16 * v_S', 'symbolic_difference':str(difference),
          'dependence_on_potential_outcome_covariance':'0', 'all_exact_checks_passed':True}
Path('sampling_variance_relation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
