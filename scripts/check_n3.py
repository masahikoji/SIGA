#!/usr/bin/env python3
"""Exact rational check of Appendix E's three-arrival counterexample.
Uses the original absolute-imbalance rule, not the formula being checked.
"""
from itertools import product
from fractions import Fraction as Q

def conditional_correlations(profiles,p):
    F=len(profiles[0]);vals=[Q(0),Q(0),Q(0)]
    for signs in product([-1,1],repeat=3):
        overall=0; margins=[[0,0] for _ in range(F)];prob=Q(1)
        for u,z in zip(profiles,signs):
            score={a:abs(overall+a)+sum(abs(margins[j][u[j]]+a) for j in range(F)) for a in [-1,1]}
            chance=Q(1,2) if score[-1]==score[1] else (p if score[z]<score[-z] else 1-p)
            prob*=chance;overall+=z
            for j in range(F): margins[j][u[j]]+=z
        for k,(i,j) in enumerate([(0,1),(0,2),(1,2)]):vals[k]+=prob*signs[i]*signs[j]
    return vals

def check(F,p):
    profiles=list(product([0,1],repeat=F));J=len(profiles);off=Q(0)
    for us in product(profiles,repeat=3):
        c=conditional_correlations(us,p)
        v=[(-1)**(F-sum(u)) for u in us]
        off+=sum(v[i]*v[j]*a*a for (i,j),a in zip([(0,1),(0,2),(1,2)],c))
    computed=1+Q(2,3*J**3)*off
    g=2*p-1
    gap={1:Q(4,3)*p*(1-p)*g*g,2:Q(1,6)*p*p*g*g,3:-Q(1,6)*p*(1-p)*g*g,4:-Q(1,32)*p*p*g*g}[F]
    assert computed==1+gap,(F,p,computed,1+gap)
    return str(computed)
if __name__=='__main__':
    for F in [1,2,3,4]:
        val=check(F,Q(4,5));print('Exact n=3, F='+str(F)+', p=4/5:',val)
