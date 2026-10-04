#!/usr/bin/env python3
"""Reproduces every entry of main-text Table 3 and Supplementary Tables 11-12 (sign of the variance
difference under the stochastic range method).  Deterministic forward recursion over the pair of
contrast states of two independent repetitions (sparsedp.py); no Monte Carlo.
Memory: the three-factor runs need about 3 GB; each block prints as it finishes.
Usage: python3 reproduce_sign_tables.py [f1|f2|f3|n3|p1|all]
"""
import sys, itertools, subprocess, numpy as np
from sparsedp import psi_sparse, restricted_eigs, make_profiles
what = sys.argv[1] if len(sys.argv) > 1 else "all"

def f1():
    print("== One factor (Table 3 rows 1-6, Supplementary Table 11); tolerance 1e-16 ==")
    for pplus, p in [(.5,.6),(.5,.7),(.5,.8),(.3,.8),(.1,.8),(.5,.95),(.3,.95),(.1,.95),(.5,.99),(.5,.999)]:
        out, probs = psi_sparse(1, [pplus], p, 80, L=82, tol=1e-16)
        v = np.array([probs[1], -probs[0]])
        r = {n: (v @ out[n][0] @ v) / (v @ np.diag(probs) @ v) for n in out}
        print(f"pi+={pplus} p={p}: " + " ".join(f"n={n}:{r[n]:.4f}" for n in (10,20,40,60,80)) +
              f" extrap={2*r[80]-r[40]:.4f} discarded={1-out[80][1]:.1e} min_n={min(r.values()):.4f}", flush=True)
    for pplus in (.5,.3,.1):
        out, probs = psi_sparse(1, [pplus], 1.0, 60, L=6, tol=0.0)
        v = np.array([probs[1], -probs[0]])
        r = [(v @ out[n][0] @ v) / (v @ np.diag(probs) @ v) for n in out]
        print(f"pi+={pplus} p=1: max |ratio-1| over n<=60 = {max(abs(x-1) for x in r):.1e}", flush=True)

def f2():
    print("== Two factors (Table 3 rows 7-8, Supplementary Table 12) ==")
    for p, nmax, tol in [(.8, 40, 1e-13), (.95, 30, 1e-14)]:
        out, probs = psi_sparse(2, [.5,.4], p, nmax, L=nmax+2, tol=tol)
        for n in (10,20,30,40):
            if n in out: print(f"(.5,.4) p={p} n={n}: eig={np.round(restricted_eigs(out[n][0],probs),4)} discarded={1-out[n][1]:.1e}", flush=True)
    out, probs = psi_sparse(2, [.5,.5], .8, 30, L=32, tol=1e-13)
    profs, _ = make_profiles(2, [.5,.5])
    W = np.column_stack([np.array([np.prod([u[f] for f in A]) if A else 1.0 for u in profs])/2 for A in [(),(0,),(1,),(0,1)]])
    D = W.T @ out[30][0] @ W * 4
    print(f"uniform p=.8 n=30 Walsh eigenvalues (|A|=0,1,1,2): {np.round(np.diag(D),4)} max offdiag={np.abs(D-np.diag(np.diag(D))).max():.1e}")

def f3():
    print("== Three factors (Table 3 row 9, Supplementary Table 12) ==")
    profs, _ = make_profiles(3, [.5,.5,.5])
    sets = [A for k in range(4) for A in itertools.combinations(range(3), k)]
    W = np.column_stack([np.array([np.prod([u[f] for f in A]) if A else 1.0 for u in profs])/np.sqrt(8) for A in sets])
    out, probs = psi_sparse(3, [.5,.5,.5], .8, 9, L=11, tol=0.0)
    for n in range(3, 10):
        D = np.diag(W.T @ out[n][0] @ W * 8)
        by = {k: np.mean([D[i] for i, A in enumerate(sets) if len(A) == k]) for k in (1,2,3)}
        print(f"uniform p=.8 n={n}: |A|=1 {by[1]:.4f} |A|=2 {by[2]:.4f} |A|=3 {by[3]:.4f} (no state discarded)", flush=True)
    out, probs = psi_sparse(3, [.5,.5,.5], .8, 10, L=12, tol=1e-11)
    D = np.diag(W.T @ out[10][0] @ W * 8)
    by = {k: np.mean([D[i] for i, A in enumerate(sets) if len(A) == k]) for k in (1,2,3)}
    print(f"uniform p=.8 n=10: |A|=1 {by[1]:.4f} |A|=2 {by[2]:.4f} |A|=3 {by[3]:.4f} discarded={1-out[10][1]:.1e}", flush=True)
    out, probs = psi_sparse(3, [.5,.4,.3], .8, 8, L=10, tol=1e-12)
    for n in (3,4,6,8): print(f"(.5,.4,.3) p=.8 n={n}: smallest eigenvalue {restricted_eigs(out[n][0],probs)[0]:.4f} discarded={1-out[n][1]:.1e}")

def n3():
    print("== n = 3 closed forms (Proposition 4): exact rational enumeration ==")
    subprocess.run([sys.executable, "exact_bruteforce.py", "3"]); subprocess.run([sys.executable, "exact_bruteforce.py", "4"])

def p1():
    print("== Deterministic one-factor rule (Corollary 4): symbolic Psi_D ==")
    subprocess.run([sys.executable, "exact_p1_F1.py"])

for name, fn in [("f1", f1), ("f2", f2), ("f3", f3), ("n3", n3), ("p1", p1)]:
    if what in ("all", name): fn()
