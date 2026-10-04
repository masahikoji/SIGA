#!/usr/bin/env python3
"""Second implementation: SymPy normal equations and explicit chamber vertices.
This verifies the finite algebraic certificates, not stochastic convergence.
Run after verify_range_certificate.py --full --output certificate_results.json.
"""
from pathlib import Path
from itertools import product
import argparse, gzip, json, random
import sympy as sp

def sg(x):
    return int(x > 0) - int(x < 0)

def check(name, data):
    F=data['F']; us=list(product((-1,1),repeat=F))
    A=sp.Matrix([[1]*len(us)]+[[u[j] for u in us] for j in range(F)])
    p=sp.diag(*[sp.Rational(x) for x in data['profile_probabilities']])
    M=A*p*A.T
    B=M.inv(method='DM')*A*p
    assert B*A.T==sp.eye(F+1)
    cache={}
    def coeff(fs):
        if fs not in cache: cache[fs]=tuple(B*sp.Matrix(fs))
        return cache[fs]
    minbase=None; minray=None; minzero=None; nvertex=0; nonconstant=[]; structural=True
    for rec in data['positive_overall_chambers']:
        c=rec['chamber']
        fs=tuple(sg(1+sum(sg(1+u[j]*c[j]) for j in range(F))) for u in us)
        b=coeff(fs)
        assert b==tuple(sp.Rational(x) for x in rec['b'])
        endpoints=[(-1,1) if v==0 else (sg(v),) for v in c]
        vv=[b[0]+sum(b[j+1]*v[j] for j in range(F)) for v in product(*endpoints)]
        nvertex+=len(vv); val=min(vv)
        assert val==sp.Rational(rec['base'])
        minbase=val if minbase is None else min(minbase,val)
        for j,v in rec['rays']:
            val=sg(c[j])*b[j+1]
            assert val==sp.Rational(v)
            minray=val if minray is None else min(minray,val)
            if len(set(fs))>1: nonconstant.append(val)
            if val==0 and fs!=tuple([1]*len(us)): structural=False
    for rec in data['zero_overall_faces']:
        e=rec['sign_face']; fs=tuple(sg(sum(ej*uj for ej,uj in zip(e,u))) for u in us)
        b=coeff(fs)
        assert b==tuple(sp.Rational(x) for x in rec['d'])
        for j, val in rec['slopes']:
            val2=e[j]*b[j+1]; assert val2==sp.Rational(val)
            minzero=val2 if minzero is None else min(minzero,val2)
    assert str(minbase)==data['minimum_base']['value']
    assert str(minray)==data['minimum_ray']['value']
    assert str(minzero)==data['minimum_zero_overall_slope']['value']
    assert structural==data['all_zero_rays_structural_constant_sign']
    if nonconstant: assert str(min(nonconstant))==data['minimum_nonconstant_ray']['value']
    return {'case':name,'all_coefficients_and_inequalities_agree':True,
            'explicit_chamber_vertices_checked':nvertex,
            'normal_equations_verified':True,'structural_zero_rays_verified':structural,'c0':str(minbase),'c1':str(minzero),
            'minimum_ray':str(minray),'certificate_pass':bool(minbase>0 and minray>=0 and minzero>0)}

def integer_rule_check():
    rng=random.Random(20260920); count=0
    for F in (2,5):
        x=[0]*(F+1)
        for t in range(1500):
            u=[rng.choice((-1,1)) for _ in range(F)]
            active=[x[0]]+[(x[0]+u[j]*x[j+1])//2 for j in range(F)]
            assert all((x[0]+u[j]*x[j+1])%2==0 for j in range(F))
            score1=sum(abs(v+1) for v in active); score0=sum(abs(v-1) for v in active)
            f=sg(sg(x[0])+sum(sg(x[0]+u[j]*x[j+1]) for j in range(F)))
            assert sg(score1-score0)==f
            z=rng.choice((-1,1)); x=[v+z*a for v,a in zip(x,[1]+u)];count+=1
    return count

def two_factor_identity_check():
    rng=random.Random(102); cases=0
    us=list(product((-1,1),repeat=2)); A=sp.Matrix([[1,u,v] for u,v in us]).T
    for i in range(100):
        nums=[rng.randint(1,1000) for _ in us]; probs=[sp.Rational(v,sum(nums)) for v in nums]
        M=A*sp.diag(*probs)*A.T
        gamma=M.inv()*A*sp.diag(*probs)*sp.Matrix([u*v for u,v in us])
        iw=[1/v for v in probs]; w=[v/sum(iw) for v in iw]
        for j,(u,v) in enumerate(us):
            assert u*v-gamma[0]-gamma[1]*u-gamma[2]*v==4*w[j]*u*v
            assert u*v*gamma[0]+v*gamma[1]+u*gamma[2]==1-4*w[j]
        cases+=1
    return cases

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('certificate',type=Path)
    ap.add_argument('--output',type=Path,default=Path('independent_audit.json')); args=ap.parse_args()
    with (gzip.open(args.certificate, 'rt') if args.certificate.suffix=='.gz' else args.certificate.open()) as f:
        data=json.load(f)
    if 'F' in data: data={'provided_law':data}
    res={'implementation':'SymPy exact inverse and explicit finite chamber vertices; no import of original checker',
         'results':[check(k,v) for k,v in data.items()],
         'original_integer_rule_checks':integer_rule_check(),
         'arbitrary_joint_two_factor_projection_identity_checks':two_factor_identity_check()}
    args.output.write_text(json.dumps(res,indent=2)+'\n')
    print(json.dumps(res,indent=2))
if __name__=='__main__': main()
