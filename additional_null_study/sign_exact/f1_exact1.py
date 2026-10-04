import numpy as np, sys, json, time
from sparsedp import psi_sparse
pplus=float(sys.argv[1]); p=float(sys.argv[2]); tol=float(sys.argv[3]); nmax=int(sys.argv[4])
t=time.time()
res,probs=psi_sparse(1,[pplus],p,nmax,L=nmax+2,tol=tol)
v=np.array([probs[1],-probs[0]])
r={n: float((v@res[n][0]@v)/(v@np.diag(probs)@v)) for n in res}
ex=2*r[nmax]-r[nmax//2]
print(f"pi+={pplus} p={p} tol={tol}: r_10={r[10]:.4f} r_20={r[20]:.4f} r_40={r[40]:.4f} r_{nmax}={r[nmax]:.4f} extrap={ex:.4f} lost_mass={1-res[nmax][1]:.2e} states={res[nmax][2]} min_n={min(r.values()):.4f} t={time.time()-t:.0f}s", flush=True)
json.dump({"ratio":r,"lost":{n:1-res[n][1] for n in res}},open(f"f1x_{pplus}_{p}.json","w"))
