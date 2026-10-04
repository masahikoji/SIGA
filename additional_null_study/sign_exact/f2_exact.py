import numpy as np, time, sys, json
from sparsedp import psi_sparse, restricted_eigs
prev=[float(s) for s in sys.argv[1].split(',')]; p=float(sys.argv[2]); nmax=int(sys.argv[3]); tol=float(sys.argv[4])
t=time.time()
out,probs=psi_sparse(2,prev,p,nmax,L=nmax+2,tol=tol)
print(f"prev={prev} p={p} tol={tol} time={time.time()-t:.0f}s states={out[nmax][2]} lost_mass={1-out[nmax][1]:.2e}")
for n in [10,20,30,40]:
    if n in out: print(f"  n={n}: eig={np.round(restricted_eigs(out[n][0],probs),5)} lost={1-out[n][1]:.1e}")
