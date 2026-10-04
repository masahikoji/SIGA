#!/usr/bin/env python3
"""Independent check of all product-box corners via conditional mean differences.
Does not import the primary checker, use its inverse, or reuse its chamber minima.
The coefficients are computed from independent-factor regression identities.
"""
from fractions import Fraction as Q
from itertools import product
from math import lcm, prod
from pathlib import Path
import json

F=5;us=list(product((-1,1),repeat=F))
sg=lambda x:(x>0)-(x<0)
chambers=[]
for c in product((-2,-1,0,1,2),repeat=F):
    f=tuple(sg(1+sum(sg(1+u[j]*c[j]) for j in range(F))) for u in us)
    vertices=list(product(*[(-1,1) if x==0 else (sg(x),) for x in c]))
    chambers.append((c,f,vertices))
faces=[(e,tuple(sg(sum(x*y for x,y in zip(e,u))) for u in us)) for e in product((-1,0,1),repeat=F) if any(e)]

def independently_check(corner):
    p=[Q(x) for x in corner];mu=[2*x-1 for x in p]
    probabilities=[prod(pj if uj==1 else 1-pj for pj,uj in zip(p,u)) for u in us]
    slopes=[]
    for j in range(F):
        slopes.append([Q(u[j],2)*prod(p[k] if u[k]==1 else 1-p[k] for k in range(F) if k!=j) for u in us])
    intercept=[probabilities[i]-sum(mu[j]*slopes[j][i] for j in range(F)) for i in range(len(us))]
    B=[intercept]+slopes;den=lcm(*(q.denominator for row in B for q in row))
    B=[[int(q*den) for q in row] for row in B];cache={}
    def coeff(f):
        if f not in cache:cache[f]=tuple(sum(a*b for a,b in zip(row,f)) for row in B)
        return cache[f]
    c0=None;ray=None;c1=None;nray=None
    for c,f,vertices in chambers:
        b=coeff(f)
        v=min(b[0]+sum(bj*x for bj,x in zip(b[1:],z)) for z in vertices)
        c0=v if c0 is None else min(c0,v)
        for j,cj in enumerate(c):
            if abs(cj)==2:
                v=sg(cj)*b[j+1];ray=v if ray is None else min(ray,v)
                if len(set(f))>1:
                    nray=v if nray is None else min(nray,v)
                if v==0: assert all(x==1 for x in f)
    for e,f in faces:
        b=coeff(f)
        for j,ej in enumerate(e):
            if ej:
                v=ej*b[j+1];c1=v if c1 is None else min(c1,v)
    return tuple(Q(x,den) for x in (c0,ray,c1,nray))

inp=json.loads(Path(__file__).with_name('product_box_certificates.json').read_text())
res={}
for name,r in inp.items():
    vals=[]
    for case in r['corners']:
        ans=independently_check(case['marginals'])
        expected=tuple(Q(case[k]['value']) for k in ('minimum_base','minimum_ray','minimum_zero_overall_slope','minimum_nonconstant_ray'))
        assert ans==expected,(name,case['marginals'],ans,expected)
        assert ans[0]>0 and ans[1]>=0 and ans[2]>0 and ans[3]>0
        vals.append(ans)
    res[name]={'corners_independently_verified':len(vals),'all_minima_match':True,
               'minima':[str(min(x[j] for x in vals)) for j in range(4)],
               'coefficient_method':'independent-factor conditional mean differences; exact rational arithmetic',
               'base_method':'explicit enumeration of every chamber vertex'}
Path(__file__).with_name('product_box_independent_audit.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
