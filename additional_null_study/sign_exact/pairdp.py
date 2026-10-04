import numpy as np, itertools, sys
from collections import defaultdict

def sgn(t): return (t>0)-(t<0)

def make_profiles(F, prev):
    # independent binary factors; prev[f] = P(u_f = +1)
    profs=[]; probs=[]
    for u in itertools.product([-1,1], repeat=F):
        pr=1.0
        for f in range(F): pr*= prev[f] if u[f]==1 else 1-prev[f]
        profs.append(u); probs.append(pr)
    return profs, np.array(probs)

def fval(x, u):
    # x = (x0, x1..xF) contrast state ; u profile
    s = sgn(x[0]) + sum(sgn(x[0]+u[f]*x[f+1]) for f in range(len(u)))
    return sgn(s)

def psi_n(F, prev, p, nmax, tol=1e-15, verbose=False):
    profs, probs = make_profiles(F, prev)
    J=len(profs); eta=2*p-1
    zero=tuple([0]*(F+1))
    # state: (x1, x2) ; store prob and M (J-vector)
    states={ (zero,zero): (1.0, np.zeros(J)) }
    S=np.zeros((J,J))
    out={}
    for n in range(1,nmax+1):
        new=defaultdict(lambda:[0.0, np.zeros(J)])
        for (x1,x2),(P,M) in states.items():
            for j,u in enumerate(profs):
                pu=probs[j]
                f1=fval(x1,u); f2=fval(x2,u)
                pz1={1:(1-eta*f1)/2, -1:(1+eta*f1)/2}
                pz2={1:(1-eta*f2)/2, -1:(1+eta*f2)/2}
                for z1 in (1,-1):
                    for z2 in (1,-1):
                        pr=pu*pz1[z1]*pz2[z2]
                        if pr==0: continue
                        xi=z1*z2
                        nx1=tuple(x1[k]+z1*(1 if k==0 else u[k-1]) for k in range(F+1))
                        nx2=tuple(x2[k]+z2*(1 if k==0 else u[k-1]) for k in range(F+1))
                        key=(nx1,nx2)
                        ent=new[key]
                        ent[0]+=P*pr
                        ent[1]+= pr*(M + P*xi*np.eye(J)[j])
                        # second moment update
                        S[:,j]+=pr*xi*M; S[j,:]+=pr*xi*M; S[j,j]+=pr*P
        # prune
        states={k:(v[0],v[1]) for k,v in new.items() if v[0]>tol}
        mass=sum(v[0] for v in states.values())
        Psi=S/n
        out[n]=(Psi.copy(), mass, len(states))
        if verbose: print(n, len(states), mass)
    return out, probs

def restricted_eigs(Psi, probs):
    J=len(probs); Pi=np.diag(probs)
    sq=np.sqrt(probs)
    # orthonormal basis of complement of sqrt(pi)
    Q,_=np.linalg.qr(np.column_stack([sq, np.eye(J)[:,:J-1]]))
    B=Q[:,1:]
    R=B.T@np.diag(1/sq)@Psi@np.diag(1/sq)@B
    return np.sort(np.linalg.eigvalsh((R+R.T)/2))

if __name__=="__main__":
    F=int(sys.argv[1]); prev=[float(s) for s in sys.argv[2].split(',')]; p=float(sys.argv[3]); nmax=int(sys.argv[4])
    out,probs=psi_n(F,prev,p,nmax)
    for n in sorted(out):
        Psi,mass,ns=out[n]
        ev=restricted_eigs(Psi,probs)
        if n%max(1,nmax//20)==0 or n<=6:
            print(f"n={n:4d} states={ns:7d} mass={mass:.12f} restricted eig(R_pi,n)=", np.round(ev,6))
