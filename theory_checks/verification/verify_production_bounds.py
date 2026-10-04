#!/usr/bin/env python3
"""Rigorous rational bounds for the production SWIFT-inspired probabilities.
Uses Machin's identity, alternating series and rational outward rounding.
The normal quantile is bracketed by verified CDF inequalities.
No floating point or approximate special-function call is used.
The formulas are those in production/swift_direct/simulation_swift_direct_siga_s_rt.R
at GitHub release v1.0.0 (lines containing AGE_*, NIHSS_*, and FACTOR_PREVALENCE_MASTER).
"""
from fractions import Fraction as Q
from math import factorial, isqrt
from pathlib import Path
import json

def arctan_bounds(x, N=90):
    s=sum(((-1)**k*x**(2*k+1)/Q(2*k+1) for k in range(N+1)),Q(0))
    t=s+(-1)**(N+1)*x**(2*N+3)/Q(2*N+3)
    return min(s,t),max(s,t)

a,b=arctan_bounds(Q(1,5));c,d=arctan_bounds(Q(1,239))
PI=(16*a-4*d,16*b-4*c)

def sqrt_bounds(x, digits=70):
    scale=10**digits
    n=isqrt(x.numerator*scale*scale//x.denominator)
    return Q(n,scale),Q(n+1,scale)
lo=sqrt_bounds(2*PI[0])[0];hi=sqrt_bounds(2*PI[1])[1]
C=(1/hi,1/lo)

def phi_bounds(z):
    if z<0:
        a,b=phi_bounds(-z); return 1-b,1-a
    if z==0:return Q(1,2),Q(1,2)
    assert z*z<36
    N=140
    s=Q(0); term=z
    for k in range(N+1):
        if k:term=-term*z*z*Q(2*k-1,2*k*(2*k+1))
        s+=term
    nxt=-term*z*z*Q(2*N+1,2*(N+1)*(2*N+3))
    A,B=min(s,s+nxt),max(s,s+nxt)
    assert A>0 and z*z/(2*(N+2))<1
    return Q(1,2)+C[0]*A,Q(1,2)+C[1]*B

def rnd(v,dec=12,upper=False):
    scale=10**dec;n=(v.numerator*scale)//v.denominator
    if upper and Q(n,scale)<v:n+=1
    return Q(n,scale)

def strdec(v,dec=12):return format(float(v),f'.{dec}f')
qlo=Q('0.67448975019');qhi=Q('0.67448975020')
assert phi_bounds(qlo)[1]<Q(3,4)<phi_bounds(qhi)[0]

def cdf_standardised(value,mean,numerator):
    # sd = numerator / (2 * normal_quantile(.75)); hence z = 2*q*(value-mean)/numerator.
    zlo=2*qlo*(value-mean)/numerator;zhi=2*qhi*(value-mean)/numerator
    a,b=min(zlo,zhi),max(zlo,zhi)
    return phi_bounds(a)[0],phi_bounds(b)[1]

def tail(cut,mean,sdnum,low,high):
    l=cdf_standardised(low,mean,sdnum);h=cdf_standardised(high,mean,sdnum)
    cc=cdf_standardised(cut,mean,sdnum)
    numerator=(h[0]-cc[1],h[1]-cc[0]);denominator=(h[0]-l[1],h[1]-l[0])
    assert numerator[0]>0 and denominator[0]>0
    return rnd(numerator[0]/denominator[1]),rnd(numerator[1]/denominator[0],upper=True)

nihss=tail(Q('17.5'),Q(17),Q(20)-Q(201*13+207*12,408),Q(5),Q(30))
age=tail(Q('69.5'),Q('72.5'),Q('16.5'),Q(18),Q(100))
vals=[nihss,age,(Q(117,408),)*2,(Q(63,408),)*2,(Q(3,10),)*2]
box=[('0.465','0.467'),('0.591','0.593'),('0.286','0.288'),('0.153','0.155'),('0.299','0.301')]
res={
 'method':'Exact rational alternating-series enclosures; 70-digit square-root bounds; verified normal-quantile bracket',
 'source':'masahikoji/SIGA v1.0.0 production/swift_direct/simulation_swift_direct_siga_s_rt.R',
 'normal_075_quantile_bracket':[str(qlo),str(qhi)],
 'probabilities':[]}
for name,ab,(l,u) in zip(['NIHSS_gt17','Age_ge70','ICA_related','Tandem','ASPECTS_4_7'],vals,box):
    assert Q(l)<ab[0]<=ab[1]<Q(u)
    res['probabilities'].append({'name':name,'exact_lower':str(ab[0]),'exact_upper':str(ab[1]),'display_lower':strdec(ab[0]),'display_upper':strdec(ab[1]),'certified_box':[l,u]})
if __name__=='__main__':
    dest=Path(__file__).with_name('production_probability_bounds.json')
    dest.write_text(json.dumps(res,indent=2)+'\n')
    print(json.dumps(res,indent=2))
