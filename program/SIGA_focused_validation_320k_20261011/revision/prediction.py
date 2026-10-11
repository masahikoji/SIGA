"""Null-only formula diagnostics with finite-n, simulated allocation covariances.

This is not an estimate of a known limiting covariance without finite-n error.
DGP truth is used only in these explicitly labelled diagnostic calculations.
"""
from __future__ import annotations
import numpy as np
from scipy.special import ndtr, ndtri
from .designs import patterns_prob
from .calibration import load_calibration


def null_predictions(root,plan,resolved):
    rows=[]
    for sc in resolved['scenarios']:
        if not sc['null_evaluation']:
            continue
        pat,pi=patterns_prob(sc['F'])
        cal=load_calibration(root,sc['calibration_key'])
        d=np.asarray(sc['mu1'])-np.asarray(sc['mu0'])-sc['boundary']
        q=(np.asarray(sc['mu1'])+np.asarray(sc['mu0'])-sc['boundary'])/2
        noise=float(pi@((np.asarray(sc['var0'])+np.asarray(sc['var1']))/2))
        for adjusted in (False,True):
            C=np.column_stack((np.ones(len(pi)),pat)) if adjusted else np.ones((len(pi),1))
            beta=np.linalg.solve(C.T@(pi[:,None]*C),C.T@(pi*q))
            m=q-C@beta
            core=float(m@cal['sigma_D_finite']@m)+noise
            vs=(core+float(pi@(d*d))/4)/4
            vr=(core+float(d@cal['psi']@d)/4)/4
            if vs<=0 or vr<=0: raise ValueError('Invalid diagnostic variance.')
            z=ndtri(1-sc['alpha']/(2 if sc['tail']=='two' else 1))
            size=float((2 if sc['tail']=='two' else 1)*ndtr(-z*np.sqrt(vr/vs)))
            rows.append(dict(id=sc['id'],score='adjusted' if adjusted else 'unadjusted',
                n=sc['n'],F=sc['F'],pbc=sc['pbc'],alpha=sc['alpha'],
                sampling_variance_per_n=vs,reference_variance_per_n=vr,
                reference_to_sampling_ratio=vr/vs,predicted_RT_null_rate=size,
                scope='First-order null formula with independent finite-n allocation covariance estimates; not exact finite-sample size or exact limiting covariance.'))
    return rows
