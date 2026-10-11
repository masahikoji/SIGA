"""Data-only inference and covariance construction.

Functions here do not accept a generating treatment effect, outcome model,
potential outcomes, true stratum effects, or true sampling/reference variances.
The calibration object contains allocation-only Monte Carlo estimates.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import ndtr

TIE_TOL = 1e-12

@dataclass(frozen=True)
class Observation:
    y: np.ndarray
    assignment: np.ndarray
    stratum: np.ndarray
    patterns: np.ndarray

    def validate(self):
        n = len(self.y)
        if n < 2 or self.assignment.shape != (n,) or self.stratum.shape != (n,):
            raise ValueError("Inconsistent observation arrays.")
        if not np.isfinite(self.y).all() or not np.isin(self.assignment, [0, 1]).all():
            raise ValueError("Non-finite outcome or nonbinary assignment.")
        if not np.isin(self.patterns, [0, 1]).all():
            raise ValueError("Allocation factors must be coded 0/1.")
        if (self.stratum < 0).any() or (self.stratum >= len(self.patterns)).any():
            raise ValueError("Invalid stratum indices.")


def psd_roundoff(a):
    """Only remove numerical negative eigenvalues; never fit to error rates."""
    a = (a + a.T)*0.5
    w, v = np.linalg.eigh(a)
    tol = 1e-10*max(1.0, float(np.max(np.abs(w))))
    if w.min() < -tol:
        raise ValueError("Substantially indefinite input covariance.")
    w = np.maximum(w, 0.0)
    return (v*w)@v.T


def contrast_matrix(patterns):
    return np.column_stack((np.ones(len(patterns)), 2*patterns-1)).T.astype(float)


def covariance_constructions(counts, gamma, xi, L):
    """Return original/refined Omega and geometry diagnostics.

    If the balancing matrix is rank deficient, the whole original construction
    is retained. A pseudoinverse is NOT used to claim unattainable Xi matching.
    Empty joint strata are otherwise allowed.
    """
    counts = np.asarray(counts, dtype=float)
    if np.any(counts < 0) or counts.sum() < 1:
        raise ValueError("Invalid counts.")
    J = len(counts)
    sq = np.sqrt(counts)
    original = gamma*sq[:, None]*sq[None, :]
    B = sq[:, None]*L.T
    K = B.T@B
    eig = np.linalg.eigvalsh(K)
    # Vanishing theoretical threshold; floating roundoff is checked separately.
    scale = max(1.0, float(eig[-1]))
    threshold = min(1e-12, 1.0/counts.sum())*scale
    full_rank = bool(eig[0] > threshold)
    if not full_rank:
        return original, original.copy(), dict(fallback=True, balance_error=np.nan)
    Q = np.linalg.solve(K, B.T).T
    P = Q@B.T
    W = np.eye(J)-P
    refined_gamma = W@gamma@W.T + Q@xi@Q.T
    refined = refined_gamma*sq[:, None]*sq[None, :]
    refined = (refined+refined.T)*0.5
    error = float(np.max(np.abs(L@refined@L.T-xi)))
    return original, refined, dict(fallback=False, balance_error=error)


def sampling_variance(score, sid, counts, Omega):
    n = len(score)
    J = len(counts)
    active = counts > 0
    means = np.zeros(J)
    sums = np.bincount(sid, weights=score, minlength=J)
    means[active] = sums[active]/counts[active]
    e = score-means[sid]
    trace = float(np.sum(np.diag(Omega)[active]/counts[active]))
    denom = n-int(active.sum())
    kraw = (n-trace)/denom if denom else 0.0
    kappa = max(kraw, 0.0)
    between = 0.25*float(means@Omega@means)
    within = 0.25*kappa*float(e@e)
    value = between+within
    if value < -1e-8 or not np.isfinite(value):
        raise FloatingPointError("Invalid sampling variance.")
    return dict(variance=max(value, 0.0), kappa=kappa, kappa_active=kraw < 0,
                means=means, residual=e, between=between, within=within)


def effect_deviations(obs, boundary):
    """Within-stratum arm-mean difference minus tested b, not generating Delta."""
    J = len(obs.patterns)
    means, sizes, vm = [], [], []
    for arm in (0, 1):
        use = obs.assignment == arm
        count = np.bincount(obs.stratum[use], minlength=J)
        sy = np.bincount(obs.stratum[use], weights=obs.y[use], minlength=J)
        sy2 = np.bincount(obs.stratum[use], weights=obs.y[use]**2, minlength=J)
        m = sy/np.maximum(count, 1)
        vv = np.zeros(J)
        ok = count > 1
        vv[ok] = np.maximum((sy2[ok]-sy[ok]**2/count[ok])/(count[ok]-1), 0)/count[ok]
        means.append(m); sizes.append(count); vm.append(vv)
    both = (sizes[0] > 0) & (sizes[1] > 0)
    d = np.where(both, means[1]-means[0]-boundary, 0.0)
    return d, np.where(both, vm[0]+vm[1], 0.0), both


def reference_variance(vs, counts, d, psi, var_d=None):
    n = int(counts.sum())
    M = psi-np.diag(counts/n)
    quad = float(d@M@d)
    if var_d is not None:
        quad -= float(np.diag(M)@var_d)
    raw = vs+n*quad/16
    value = max(raw, vs/n)
    return dict(variance=max(value, 0.0), raw=raw, floor=raw < vs/n,
                correction=n*quad/16)


def scores_and_variances(obs, boundaries, calibration, noise_adjustment=False):
    """Data-only analysis, exactly the same dhat for both constructions."""
    obs.validate()
    n = len(obs.y)
    J = len(obs.patterns)
    counts = np.bincount(obs.stratum, minlength=J).astype(float)
    L = contrast_matrix(obs.patterns)
    O0, O1, geom = covariance_constructions(counts, calibration['gamma'], calibration['xi'], L)
    C = np.column_stack((np.ones(n), obs.patterns[obs.stratum]))
    out = []
    for b in boundaries:
        W = obs.y-b*obs.assignment
        unadj = W-W.mean()
        adj = W-C@np.linalg.lstsq(C, W, rcond=None)[0]
        adj -= adj.mean()
        dhat, vd, both = effect_deviations(obs, b)
        for name, score in [('unadjusted', unadj), ('adjusted', adj)]:
            entry = dict(boundary=b, score_name=name, score=score,
                         statistic=float((2*obs.assignment-1)@score/2),
                         counts=counts, both_arms=int(both.sum()), dhat=dhat)
            for tag, O in [('original', O0), ('refined', O1)]:
                sv = sampling_variance(score, obs.stratum, counts, O)
                rv = reference_variance(sv['variance'], counts, dhat, calibration['psi'])
                entry[tag] = dict(S=sv, R=rv)
                if noise_adjustment:
                    entry[tag]['R_noise'] = reference_variance(sv['variance'], counts, dhat, calibration['psi'], vd)
            out.append(entry)
    return out, geom


def normal_pvalue(t, variance, tail):
    if not variance > 0:
        return 1.0 if abs(t) <= TIE_TOL else 0.0
    z = t/np.sqrt(variance)
    if tail == 'two': return float(2*ndtr(-abs(z)))
    if tail == 'upper': return float(ndtr(-z))
    if tail == 'lower': return float(ndtr(z))
    raise ValueError(tail)


def mc_pvalue(samples, observed, tail, scale=1.0, retain_raw_ties=False):
    if not np.isfinite(scale) or scale <= 0:
        return 1.0
    t = float(observed)
    if tail == 'two':
        exceed = np.abs(samples)+TIE_TOL >= scale*abs(t)
        tie = np.abs(np.abs(samples)-abs(t)) <= TIE_TOL
    elif tail == 'upper':
        exceed = samples+TIE_TOL >= scale*t
        tie = np.abs(samples-t) <= TIE_TOL
    elif tail == 'lower':
        exceed = samples-TIE_TOL <= scale*t
        tie = np.abs(samples-t) <= TIE_TOL
    else: raise ValueError(tail)
    if retain_raw_ties: exceed = exceed | tie
    return float((1+np.count_nonzero(exceed))/(len(samples)+1))


def _normal_cell_mass(lo, hi):
    # Survival-tail subtraction avoids catastrophic cancellation at large z.
    return np.where(lo >= 0, ndtr(-lo)-ndtr(-hi), ndtr(hi)-ndtr(lo))


def lattice_pvalue(t, variance, n, successes, weights, scale=1.0, raw_ties=False):
    """b=0, unadjusted binary, inclusive absolute-value tail only.

    Use REFERENCE variance for the corrected RT target; scale=lambda and add
    raw-scale ties as a union, not twice. This never jitters data or uses mid-p.
    """
    if not variance > 0: return 1.0 if abs(t) <= TIE_TOL else 0.0
    sd = np.sqrt(variance)
    prob = 0.0
    for k in np.flatnonzero(weights > 0):
        lo, hi = max(0, k-(n-successes)), min(k, successes)
        x = np.arange(lo, hi+1, dtype=float)
        values = x-k*successes/n
        masses = _normal_cell_mass((values-.5)/sd, (values+.5)/sd)
        total = float(masses.sum())
        if total <= 0: raise FloatingPointError("Zero mass in feasible lattice.")
        use = np.abs(values)+TIE_TOL >= scale*abs(t)
        if raw_ties: use |= np.abs(np.abs(values)-abs(t)) <= TIE_TOL
        prob += float(weights[k])*float(masses[use].sum())/total
    return float(np.clip(prob, 0.0, 1.0))
