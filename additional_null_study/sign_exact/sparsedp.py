import numpy as np, itertools

def make_profiles(F, prev):
    profs=[]; probs=[]
    for u in itertools.product([-1,1], repeat=F):
        pr=1.0
        for f in range(F): pr*= prev[f] if u[f]==1 else 1-prev[f]
        profs.append(u); probs.append(pr)
    return profs, np.array(probs)

class Codec:
    def __init__(self, ncoord, L):
        self.L=L; self.base=2*L+1; self.nc=ncoord
        self.mult=np.array([self.base**k for k in range(ncoord)], dtype=np.int64)
    def encode(self, coords):  # coords: (nc, N) int arrays in [-L,L]
        return ((coords+self.L).astype(np.int64)*self.mult[:,None]).sum(axis=0)
    def decode(self, codes):
        out=np.empty((self.nc, codes.size), dtype=np.int64)
        c=codes.copy()
        for k in range(self.nc):
            out[k]=c % self.base - self.L
            c//=self.base
        return out

def fvals(x, u):
    # x: (d, N) coords, d=F+1
    s=np.sign(x[0]).astype(np.int64)
    for f in range(len(u)): s+=np.sign(x[0]+u[f]*x[f+1])
    return np.sign(s)

def psi_sparse(F, prev, p, nmax, L=20, tol=1e-13, report=None):
    profs, probs = make_profiles(F, prev); J=len(profs); eta=2*p-1
    d=F+1; cod=Codec(2*d, L)
    codes=cod.encode(np.zeros((2*d,1),dtype=np.int64))
    P=np.array([1.0]); M=np.zeros((1,J))
    S=np.zeros((J,J)); out={}
    for n in range(1,nmax+1):
        X=cod.decode(codes); x1=X[:d]; x2=X[d:]
        aggC=None; aggP=None; aggM=None
        for j,u in enumerate(profs):
            allcodes=[]; allP=[]; allM=[]
            pu=probs[j]; a=np.array([1]+list(u))
            f1=fvals(x1,u); f2=fvals(x2,u)
            for z1 in (1,-1):
                for z2 in (1,-1):
                    w=pu*(1-eta*z1*f1)/2*(1-eta*z2*f2)/2
                    xi=z1*z2
                    PW=P*w
                    MW=M*w[:,None]
                    mj=MW.sum(axis=0)
                    S[:,j]+=xi*mj; S[j,:]+=xi*mj; S[j,j]+=PW.sum()
                    MW[:,j]+=xi*PW
                    shift=np.concatenate([z1*a, z2*a])
                    newX=X+shift[:,None]
                    ok=(np.abs(newX)<=L).all(axis=0) & (PW>0)
                    allcodes.append(cod.encode(newX[:,ok])); allP.append(PW[ok]); allM.append(MW[ok])
            if aggC is not None: allcodes.append(aggC); allP.append(aggP); allM.append(aggM)
            C=np.concatenate(allcodes); Pv=np.concatenate(allP); Mv=np.concatenate(allM)
            uniq, inv = np.unique(C, return_inverse=True)
            aggP=np.bincount(inv, weights=Pv, minlength=uniq.size)
            aggM=np.zeros((uniq.size,J))
            for jj in range(J): aggM[:,jj]=np.bincount(inv, weights=Mv[:,jj], minlength=uniq.size)
            aggC=uniq
            del C,Pv,Mv,inv
        keep=aggP>tol
        codes=aggC[keep]; P=aggP[keep]; M=aggM[keep]
        out[n]=(S/n, P.sum(), codes.size)
        if report and (n in report): print(f"   step {n}: states={codes.size} mass={P.sum():.8f}", flush=True)
    return out, probs

def restricted_eigs(Psi, probs):
    J=len(probs); sq=np.sqrt(probs)
    Q,_=np.linalg.qr(np.column_stack([sq, np.eye(J)[:,:J-1]]))
    B=Q[:,1:]
    R=B.T@np.diag(1/sq)@Psi@np.diag(1/sq)@B
    return np.sort(np.linalg.eigvalsh((R+R.T)/2))
