#!/usr/bin/env python3
"""Exact rational certificates for the unchanged binary-factor range rule.

No Monte Carlo, floating-point inversion, or state-space truncation is used.
The checker tests all sign chambers and their recession rays.
Outputs JSON summaries; optionally all individual inequalities.
Python 3.10+; standard library only.
"""
from __future__ import annotations
from fractions import Fraction as Q
from itertools import product
import argparse
import json
from pathlib import Path
from typing import Iterable


def sign(x: int | Q) -> int:
    return (x > 0) - (x < 0)


def inverse(matrix: list[list[Q]]) -> list[list[Q]]:
    n = len(matrix)
    a = [list(row) + [Q(int(i == j)) for j in range(n)]
         for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            raise ValueError('The profile second-moment matrix is singular.')
        a[j], a[pivot] = a[pivot], a[j]
        scale = a[j][j]
        a[j] = [v / scale for v in a[j]]
        for i in range(n):
            if i != j:
                scale = a[i][j]
                a[i] = [v - scale * w for v, w in zip(a[i], a[j])]
    return [row[n:] for row in a]


def product_probabilities(marginals: Iterable[str]) -> list[Q]:
    p = [Q(v) for v in marginals]
    if not p or any(v <= 0 or v >= 1 for v in p):
        raise ValueError('All binary-factor probabilities must be in (0,1).')
    ans = []
    for u in product((-1, 1), repeat=len(p)):
        val = Q(1)
        for uj, pj in zip(u, p):
            val *= pj if uj == 1 else 1-pj
        ans.append(val)
    return ans


def certificate(probs: list[Q], include_full: bool = False) -> dict:
    n = len(probs)
    F = n.bit_length()-1
    if n != 2**F or F < 1 or sum(probs) != 1 or any(p <= 0 for p in probs):
        raise ValueError('Provide 2**F strictly positive probabilities summing to one.')
    us = list(product((-1, 1), repeat=F))
    aa = [(1,)+u for u in us]
    M = [[sum((p*a[i]*a[j] for p,a in zip(probs,aa)), Q(0))
          for j in range(F+1)] for i in range(F+1)]
    Mi = inverse(M)
    # Every coefficient vector is a rational matrix applied to an integer sign vector.
    K = [[p*sum((Mi[j][l]*a[l] for l in range(F+1)), Q(0))
          for p,a in zip(probs,aa)] for j in range(F+1)]
    # Clear denominators once. The chamber enumeration thereafter uses integers.
    from math import lcm
    den = lcm(*(x.denominator for row in K for x in row))
    Knum = [[int(x*den) for x in row] for row in K]
    def coefficient_numerators(f: tuple[int,...]) -> tuple[int,...]:
        return tuple(sum(x*y for x,y in zip(row,f)) for row in Knum)
    cached: dict[tuple[int,...], tuple[int,...]] = {}
    def coeff(f: tuple[int,...]) -> tuple[int,...]:
        if f not in cached:
            cached[f] = coefficient_numerators(f)
        return cached[f]
    bases=[]; rays=[]; boundary=[]; full=[]; full_zero=[]
    zero_ray_nonconstant=[]
    nonconstant_rays=[]
    for cats in product((-2,-1,0,1,2), repeat=F):
        f = tuple(sign(1+sum(sign(1+uj*cj) for uj,cj in zip(u,cats))) for u in us)
        b = coeff(f)
        base = b[0] + sum(-abs(b[j+1]) if c==0 else sign(c)*b[j+1]
                          for j,c in enumerate(cats))
        bases.append((base,cats))
        these=[]
        for j,c in enumerate(cats):
            if abs(c)==2:
                val=sign(c)*b[j+1]
                rays.append((val,cats,j))
                these.append((j,val))
                if len(set(f))>1:
                    nonconstant_rays.append((val,cats,j))
                    if val==0:
                        zero_ray_nonconstant.append((cats,j))
        if include_full:
            full.append({'chamber': cats, 'b': [str(Q(v,den)) for v in b],
                         'base': str(Q(base,den)),
                         'rays': [[j,str(Q(v,den))] for j,v in these]})
    for e in product((-1,0,1), repeat=F):
        if not any(e):
            continue
        f=tuple(sign(sum(ej*uj for ej,uj in zip(e,u))) for u in us)
        b=coeff(f)
        slopes=[]
        for j,ej in enumerate(e):
            if ej:
                boundary.append((ej*b[j+1],e,j))
                slopes.append((j,ej*b[j+1]))
        if include_full:
            full_zero.append({'sign_face':e, 'd':[str(Q(v,den)) for v in b],
                              'slopes':[[j,str(Q(v,den))] for j,v in slopes]})
    minbase=min(bases); minray=min(rays); minbd=min(boundary)
    def record(item: tuple) -> dict:
        return {'value':str(Q(item[0],den)), 'pattern':item[1],
                **({'coordinate_0_based':item[2]} if len(item)==3 else {})}
    result={
        'F':F, 'profile_order':'lexicographic (-1,+1)^F',
        'profile_probabilities':[str(p) for p in probs],
        'chambers_checked':len(bases), 'ray_inequalities_checked':len(rays),
        'zero_overall_inequalities_checked':len(boundary),
        'minimum_base':record(minbase), 'minimum_ray':record(minray),
        'minimum_zero_overall_slope':record(minbd),
        'minimum_nonconstant_ray':record(min(nonconstant_rays)) if nonconstant_rays else None,
        'zero_rays_in_nonconstant_chambers':zero_ray_nonconstant,
        'all_zero_rays_structural_constant_sign':len(zero_ray_nonconstant)==0,
        'certificate_pass': minbase[0]>0 and minray[0]>=0 and minbd[0]>0,
        'exact_arithmetic':True,
        'method':'rational Gauss-Jordan; denominator-cleared exhaustive integer chamber checks',
    }
    if include_full:
        result['positive_overall_chambers']=full
        result['zero_overall_faces']=full_zero
        result['design_second_moment']=[[str(x) for x in row] for row in M]
        result['inverse_design_second_moment']=[[str(x) for x in row] for row in Mi]
    return result


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('range_certificate_results.json'))
    parser.add_argument('--full', action='store_true', help='Include every positive-overall chamber.')
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--marginals', help='Comma-separated independent binary probabilities, e.g. 0.5,0.4')
    group.add_argument('--joint-file', type=Path, help='JSON list of rational probability strings in lexicographic (-1,+1)^F order')
    args=parser.parse_args()
    if args.marginals or args.joint_file:
        try:
            probs=(product_probabilities(args.marginals.split(',')) if args.marginals else
                   [Q(str(x)) for x in json.loads(args.joint_file.read_text())])
            r=certificate(probs,args.full)
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(r,indent=2)+'\n')
            print('PASS' if r['certificate_pass'] else 'FAIL')
            print('c0:',r['minimum_base'],'c1:',r['minimum_zero_overall_slope'])
            print('Saved',args.output)
            return
        except (ValueError,OSError,TypeError) as exc:
            parser.error(str(exc))
    cases={
      'main_two_factors':['0.50','0.40'],
      'main_five_factors':['0.50','0.40','0.30','0.20','0.10'],
      'trial_five_factors':['0.466','0.592','0.287','0.154','0.300'],
      'certificate_failure_not_instability':['0.01']*5,
    }
    results={}
    for name,marg in cases.items():
        r=certificate(product_probabilities(marg), args.full)
        r['marginal_probabilities']=marg
        results[name]=r
        print(name, 'PASS' if r['certificate_pass'] else 'FAIL',
              'c0='+r['minimum_base']['value'],
              'ray='+r['minimum_ray']['value'],
              'c1='+r['minimum_zero_overall_slope']['value'],
              'nonconstant_ray='+str(r['minimum_nonconstant_ray']),flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results,indent=2)+'\n')
    print('Saved',args.output)

if __name__=='__main__':
    main()
