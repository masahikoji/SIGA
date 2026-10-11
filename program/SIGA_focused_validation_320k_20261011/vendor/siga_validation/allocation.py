"""The original equal-weight, overall-plus-marginal range rule.

All random numbers are generated outside compiled functions, from independent
NumPy PCG64DXSM streams. Compiled and noncompiled paths implement the same rule.
No outcome information enters allocation or covariance simulation.
"""
from __future__ import annotations
import numpy as np
from numba import njit

@njit(cache=True, nogil=True)
def path_from_uniforms(sid, patterns, pbc, uniforms):
    n = len(sid)
    F = patterns.shape[1]
    marginal = np.zeros((F, 2), dtype=np.int64)
    overall = 0
    out = np.empty(n, dtype=np.int8)
    for i in range(n):
        # For integer imbalances, |d+1|-|d-1| = 2 sign(d).
        vote = int(overall > 0) - int(overall < 0)
        s = sid[i]
        for f in range(F):
            current = marginal[f, patterns[s, f]]
            vote += int(current > 0) - int(current < 0)
        prob_plus = 1.0-pbc if vote > 0 else (pbc if vote < 0 else 0.5)
        z = 1 if uniforms[i] < prob_plus else -1
        out[i] = z
        overall += z
        for f in range(F):
            marginal[f, patterns[s, f]] += z
    return out

@njit(cache=True, nogil=True)
def calibration_batch(sids, patterns, pbc, uniforms):
    m, n = sids.shape
    J = len(patterns)
    U = np.zeros((m, J))
    D = np.zeros((m, J))
    G01 = np.zeros((m, J))
    G02 = np.zeros((m, J))
    treated = np.zeros(m, dtype=np.int64)
    for b in range(m):
        sid = sids[b]
        z0 = path_from_uniforms(sid, patterns, pbc, uniforms[b, 0])
        z1 = path_from_uniforms(sid, patterns, pbc, uniforms[b, 1])
        z2 = path_from_uniforms(sid, patterns, pbc, uniforms[b, 2])
        counts = np.zeros(J, dtype=np.int64)
        for i in range(n):
            s = sid[i]
            counts[s] += 1
            D[b, s] += z0[i]
            G01[b, s] += z0[i]*z1[i]/np.sqrt(n)
            G02[b, s] += z0[i]*z2[i]/np.sqrt(n)
            treated[b] += int(z0[i] == 1)
        for s in range(J):
            if counts[s] > 0:
                U[b, s] = D[b, s]/np.sqrt(counts[s])
    return U, D, G01, G02, treated

@njit(cache=True, nogil=True)
def reference_batch(sid, patterns, pbc, scores, uniforms):
    """All score columns use the same regenerated assignment sequences."""
    b, n = uniforms.shape
    L = scores.shape[1]
    out = np.zeros((b, L))
    for j in range(b):
        z = path_from_uniforms(sid, patterns, pbc, uniforms[j])
        for i in range(n):
            for l in range(L):
                out[j, l] += 0.5*z[i]*scores[i, l]
    return out
