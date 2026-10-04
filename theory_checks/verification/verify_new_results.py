#!/usr/bin/env python3
"""Exact finite checks for the added identities; no operating-characteristic simulation.

These checks complement the written proofs and are not a formal proof assistant.
The CSV is a copy of an existing GitHub release result, not newly generated data.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import csv, hashlib, json
import numpy as np
BASE=Path(__file__).resolve().parent

def covariance(rows: list[tuple[F,list[int]]]) -> list[list[F]]:
    d=len(rows[0][1]); mu=[sum(w*x[j] for w,x in rows) for j in range(d)]
    return [[sum(w*x[j]*x[k] for w,x in rows)-mu[j]*mu[k] for k in range(d)] for j in range(d)]

def pair_options(s: tuple[int,int], r: F):
    for z in product((-1,1),repeat=2):
        q=F(1,2) if s[0]==s[1] else (r if z[1]==-z[0] else 1-r)
        yield z,F(1,2)*q

def check_counterexample(r: F):
    subsets=[c for m in range(1,4) for c in combinations(range(3),m)]
    rows=[]
    for s in product((0,1),repeat=2):
        options=list(pair_options(s,r))
        for sample in product(options,repeat=3):
            signs=[t[0] for t in sample]
            prob=F(1,4)
            for _,w in sample: prob*=w
            reward=[s.count(j)-1 for j in range(2)]
            for c in subsets:
                for j in range(2):
                    reward.append(sum(np.prod([signs[a][i] for a in c]).item() for i in range(2) if s[i]==j))
            rows.append((prob,reward))
    assert sum(w for w,x in rows)==1
    cov=covariance(rows); a=2*r-1
    for bi,B in enumerate(subsets):
        start=2+2*bi; t=(-a)**len(B)/4
        assert cov[start][start]/2==F(1,2)
        assert cov[start+1][start+1]/2==F(1,2)
        assert cov[start][start+1]/2==t
        for j in range(len(cov)):
            if j not in (start,start+1): assert cov[start][j]==0 and cov[start+1][j]==0
    # Unit effect direction (1,-1), divided by the score's factor 16.
    assert (2*(a*a/4)*(-1))/16 == -a*a/32
    return {'r':str(r),'configurations':len(rows),'gap_for_d_1_minus1':str(-a*a/32)}

def check_covariance_identity():
    signs=list(product((-1,1),repeat=4)); s=(0,1,0,1)
    weights=[F(1+sum(z)**2+(z[0]*z[1]+z[2]*z[3])**2) for z in signs]
    weights=[w/sum(weights) for w in weights]
    C=covariance(list(zip(weights,[list(z) for z in signs])))
    rows=[]
    for i,x in enumerate(signs):
        for j,y in enumerate(signs):
            rows.append((weights[i]*weights[j],[sum(x[t]*y[t] for t in range(4) if s[t]==u) for u in range(2)]))
    got=covariance(rows)
    expected=[[sum(C[i][j]**2 for i in range(4) for j in range(4) if s[i]==u and s[j]==v) for v in range(2)] for u in range(2)]
    assert got==expected
    return {'configurations':len(rows),'passed':True}

def check_efron_drift():
    checks=0
    for p,u in [(F(3,5),F(6,5)),(F(4,5),F(2)),(F(19,20),F(3))]:
        q=1-p; r=q*u+p/u
        assert 1<u<p/q and r<1
        for d in range(-12,13):
            pp=F(1,2) if d==0 else (q if d>0 else p)
            lhs=pp*u**abs(d+1)+(1-pp)*u**abs(d-1)
            rhs=u if d==0 else r*u**abs(d)
            assert lhs==rhs; checks+=1
    return {'exact_integer_state_checks':checks,'passed':True}

def check_spectrum():
    pi=np.array([.10,.20,.30,.40]); D=np.diag(pi)
    raw=np.array([[1,.2,.3,.4],[.2,1,.2,.1],[.3,.2,1,.2],[.4,.1,.2,1.]])
    Psi=np.diag(np.sqrt(pi))@raw@np.diag(np.sqrt(pi))
    B=np.linalg.qr(np.sqrt(pi)[:,None],mode='complete')[0][:,1:]
    Q=np.linalg.qr(pi[:,None],mode='complete')[0][:,1:]
    R=B.T@raw@B; M=Q.T@D@Q; S=Q.T@Psi@Q
    vals,vec=np.linalg.eigh(M); W=vec@np.diag(vals**-.5)@vec.T
    a=np.linalg.eigvalsh(R); b=np.linalg.eigvalsh(W@S@W)
    assert np.max(np.abs(a-b))<1e-12
    return {'max_absolute_eigenvalue_difference':float(np.max(np.abs(a-b))),'passed':True}

def check_archive():
    data=(BASE/'design_pair_path_diagnostics.csv').read_bytes()
    blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    assert blob=='27a39178656bb285340891193c64112e5e5d21b2',blob
    rows=list(csv.DictReader(data.decode().splitlines()))
    assert len(rows)==6
    ms=(BASE.parent/'supplement.tex').read_text()
    for row in rows:
        for key in ['analysis_min_ratio','analysis_max_ratio']:
            token=f'{float(row[key]):.6f}'
            assert token in ms or token[1:] in ms,(key,token)
    return {'source':'https://github.com/masahikoji/SIGA/blob/v1.0.0/data/production_results/raw/pair_path_supplemental/design_pair_path_diagnostics.csv','git_blob_sha':blob,'sha256':hashlib.sha256(data).hexdigest(),'rows':6,'values_match_manuscript':True,'covariance_matrices_recomputed':False}

if __name__=='__main__':
    result={'counterexample':[check_counterexample(r) for r in (F(3,5),F(4,5),F(19,20))], 'finite_covariance_identity':check_covariance_identity(),'efron_drift':check_efron_drift(),'spectral_coordinates':check_spectrum(),'archived_diagnostics':check_archive(),'new_trial_simulations':0}
    (BASE/'new_results_checks.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
