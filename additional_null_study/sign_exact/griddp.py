import numpy as np, itertools

def sgn(a): return np.sign(a).astype(np.int8)

def make_profiles(F, prev):
    profs=[]; probs=[]
    for u in itertools.product([-1,1], repeat=F):
        pr=1.0
        for f in range(F): pr*= prev[f] if u[f]==1 else 1-prev[f]
        profs.append(u); probs.append(pr)
    return profs, np.array(probs)

def psi_grid(F, prev, p, nmax, L):
    """Deterministic DP over pair contrast states (x^(1),x^(2)), each coordinate in [-L,L].
    Mass that would leave the grid is dropped; the retained mass is reported.
    Returns Psi_n = Var(G_n)/n with G_n = sum_i h(S_i) Z_i^(1) Z_i^(2)."""
    profs, probs = make_profiles(F, prev); J=len(profs); eta=2*p-1
    d=F+1; side=2*L+1
    shape=(side,)*(2*d)
    P=np.zeros(shape); M=np.zeros((J,)+shape)
    c=(L,)*(2*d); P[c]=1.0
    S=np.zeros((J,J))
    ax=np.arange(-L,L+1)
    grids=np.meshgrid(*([ax]*(2*d)), indexing='ij')
    x1=grids[:d]; x2=grids[d:]
    fvals={}
    for j,u in enumerate(profs):
        for copy,x in ((1,x1),(2,x2)):
            s=sgn(x[0]).astype(np.int16)
            for f in range(F): s=s+sgn(x[0]+u[f]*x[f+1])
            fvals[(j,copy)]=np.sign(s).astype(np.int8)
    out={}
    sumaxes=tuple(range(1,2*d+1))
    for n in range(1,nmax+1):
        newP=np.zeros(shape); newM=np.zeros((J,)+shape)
        for j,u in enumerate(profs):
            pu=probs[j]
            a=[1]+list(u)
            f1=fvals[(j,1)]; f2=fvals[(j,2)]
            for z1 in (1,-1):
                for z2 in (1,-1):
                    w=pu*(1-eta*z1*f1)/2*(1-eta*z2*f2)/2
                    xi=z1*z2
                    shift=[z1*a[k] for k in range(d)]+[z2*a[k] for k in range(d)]
                    PW=P*w
                    MW=M*w
                    mj=MW.sum(axis=sumaxes)
                    S[:,j]+=xi*mj; S[j,:]+=xi*mj; S[j,j]+=PW.sum()
                    src=[]; dst=[]
                    for sh in shift:
                        if sh>=0: src.append(slice(0,side-sh)); dst.append(slice(sh,side))
                        else: src.append(slice(-sh,side)); dst.append(slice(0,side+sh))
                    src=tuple(src); dst=tuple(dst)
                    newP[dst]+=PW[src]
                    addM=MW; addM[j]+=xi*PW
                    newM[(slice(None),)+dst]+=addM[(slice(None),)+src]
        P,M=newP,newM
        out[n]=(S/n, P.sum())
    return out, probs

def restricted_eigs(Psi, probs):
    J=len(probs); sq=np.sqrt(probs)
    Q,_=np.linalg.qr(np.column_stack([sq, np.eye(J)[:,:J-1]]))
    B=Q[:,1:]
    R=B.T@np.diag(1/sq)@Psi@np.diag(1/sq)@B
    return np.sort(np.linalg.eigvalsh((R+R.T)/2))
