# Exact rational computation of Psi_n = Var(G_n)/n for small n by full enumeration (two copies).
from fractions import Fraction as Fr
import itertools, sys
def sgn(t): return (t>0)-(t<0)
def fval(x,u):
    s=sgn(x[0])+sum(sgn(x[0]+u[f]*x[f+1]) for f in range(len(u)))
    return sgn(s)
def psi_exact(F, prev, p, n):
    profs=list(itertools.product([-1,1],repeat=F)); J=len(profs)
    pu={u: Fr(1) for u in profs}
    for u in profs:
        for f in range(F): pu[u]*= prev[f] if u[f]==1 else 1-prev[f]
    eta=2*p-1
    def pz(x,u,z):
        f=fval(x,u); return (1-eta*z*f)/2
    # enumerate: recursive over steps, carrying pair state and G vector; accumulate E[G G^T]
    S=[[Fr(0)]*J for _ in range(J)]
    def rec(i,x1,x2,G,pr):
        if i==n:
            for a in range(J):
                if G[a]==0: continue
                for b in range(J): S[a][b]+=pr*G[a]*G[b]
            return
        for j,u in enumerate(profs):
            for z1 in (1,-1):
                p1=pz(x1,u,z1)
                if p1==0: continue
                for z2 in (1,-1):
                    p2=pz(x2,u,z2)
                    if p2==0: continue
                    G2=list(G); G2[j]+=z1*z2
                    y1=tuple(x1[k]+z1*(1 if k==0 else u[k-1]) for k in range(F+1))
                    y2=tuple(x2[k]+z2*(1 if k==0 else u[k-1]) for k in range(F+1))
                    rec(i+1,y1,y2,G2,pr*pu[u]*p1*p2)
    zero=tuple([0]*(F+1)); rec(0,zero,zero,[0]*J,Fr(1))
    Psi=[[S[a][b]/n for b in range(J)] for a in range(J)]
    return Psi, profs, pu
if __name__=="__main__":
    F=3; prev=[Fr(1,2)]*3; p=Fr(4,5); n=int(sys.argv[1]) if len(sys.argv)>1 else 3
    Psi,profs,pu=psi_exact(F,prev,p,n)
    # Walsh ratio for the three-way interaction: v_s = u1 u2 u3 ; ratio = v^T Psi v / v^T Pi v
    for A in [(0,),(0,1),(0,1,2)]:
        v=[Fr(int(1)) for _ in profs]
        for s,u in enumerate(profs):
            for f in A: v[s]*=u[f]
        num=sum(v[a]*Psi[a][b]*v[b] for a in range(len(profs)) for b in range(len(profs)))
        den=sum(v[a]*v[a]*pu[profs[a]] for a in range(len(profs)))
        r=num/den
        print(f"n={n} F=3 uniform p=4/5, Walsh set {A}: exact ratio = {r} = {float(r):.6f}")
