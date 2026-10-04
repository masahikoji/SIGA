#!/usr/bin/env python3
"""Recompute the two non-product witnesses and check with a second implementation."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import gzip, json
from verify_range_certificate import certificate, product_probabilities
from audit_certificate_independent import check

cases = {
 'main_correlated': (['.50','.40','.30','.20','.10'], '165433/773750'),
 'trial_rounded_correlated': (['.466','.592','.287','.154','.300'], '1911843769570152559/7556024433125000000'),
}
all_results={}; audit=[]
for name,(marginals,expected) in cases.items():
    ps=[Q(x) for x in marginals]; us=list(product((-1,1), repeat=len(ps)))
    mu=[2*p-1 for p in ps]
    probs=[p*(1+Q(1,10)*(u[0]-mu[0])*(u[1]-mu[1]))
           for p,u in zip(product_probabilities(marginals),us)]
    assert sum(probs)==1 and min(probs)>0
    for j in range(len(ps)):
        assert sum(p for p,u in zip(probs,us) if u[j]==1)==ps[j]
    result=certificate(probs,include_full=True)
    assert result['certificate_pass']
    assert result['minimum_zero_overall_slope']['value']==expected
    all_results[name]=result; audit.append(check(name,result))
with gzip.open('correlated_certificate_results.json.gz','wt') as f:
    json.dump(all_results,f,indent=2)
Path('correlated_independent_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
