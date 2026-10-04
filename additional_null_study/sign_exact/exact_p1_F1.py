import sympy as sp, itertools
def sgn(t): return (t>0)-(t<0)
def fval(x,u):
    s=sgn(x[0])+sum(sgn(x[0]+u[f]*x[f+1]) for f in range(len(u)))
    return sgn(s)
q=sp.symbols('q', positive=True)      # q = pi_+ = P(U=+1)
profs=[(-1,),(1,)]; pu={(-1,):1-q,(1,):q}
J=2; idx={(-1,):0,(1,):1}
half=sp.Rational(1,2)
def step_dist(x,u):
    f=fval(x,u)
    if f==0: return {1:half,-1:half}
    return {-f:sp.Integer(1)}
# pair states reachable from (0,0)
zero=(0,0)
o=(zero,zero)
states=[o]; seen={o}; trans={}
k=0
while k<len(states):
    v=states[k]; k+=1
    x1,x2=v; trans[v]=[]
    for u in profs:
        d1=step_dist(x1,u); d2=step_dist(x2,u)
        for z1,p1 in d1.items():
            for z2,p2 in d2.items():
                y1=(x1[0]+z1, x1[1]+z1*u[0]); y2=(x2[0]+z2, x2[1]+z2*u[0])
                w=(y1,y2); pr=pu[u]*p1*p2
                r=sp.zeros(J,1); r[idx[u]]=z1*z2
                trans[v].append((w,pr,r))
                if w not in seen: seen.add(w); states.append(w)
print("pair states:",len(states))
others=[v for v in states if v!=o]
ind={v:i for i,v in enumerate(others)}
N=len(others)
# m: hitting times
A=sp.eye(N); b=sp.ones(N,1)
for v in others:
    for (w,pr,r) in trans[v]:
        if w!=o: A[ind[v],ind[w]]-=pr
m=A.LUsolve(b); m=m.applyfunc(sp.simplify)
# g: expected reward to o (J-vectors)
G=sp.zeros(N,J); bg=sp.zeros(N,J)
for v in others:
    for (w,pr,r) in trans[v]:
        bg[ind[v],:]+= (pr*r).T
g=A.LUsolve(bg); g=g.applyfunc(sp.simplify)
def gvec(w): return sp.zeros(J,1) if w==o else g[ind[w],:].T
# Q: second moment matrices (vectorise JxJ -> 4 columns)
bq=sp.zeros(N,J*J)
for v in others:
    for (w,pr,r) in trans[v]:
        gw=gvec(w)
        mat=r*r.T + r*gw.T + gw*r.T
        bq[ind[v],:]+= pr*sp.Matrix([mat[i,j] for i in range(J) for j in range(J)]).T
Q=A.LUsolve(bq); Q=Q.applyfunc(sp.simplify)
def Qmat(w): return sp.zeros(J,J) if w==o else sp.Matrix(J,J,list(Q[ind[w],:]))
Etau=1+sum(pr*(0 if w==o else m[ind[w]]) for (w,pr,r) in trans[o])
EW=sp.zeros(J,1); EWW=sp.zeros(J,J)
for (w,pr,r) in trans[o]:
    gw=gvec(w); EW+=pr*(r+gw); EWW+=pr*(r*r.T+r*gw.T+gw*r.T+Qmat(w))
Etau=sp.simplify(Etau); EW=EW.applyfunc(sp.simplify); EWW=EWW.applyfunc(sp.simplify)
Psi=(EWW/Etau).applyfunc(sp.factor)
print("E tau =",Etau); print("E W =",EW.T); print("Psi_D =",Psi)
Pi=sp.diag(1-q,q)
v=sp.Matrix([q,-(1-q)])  # pi^T v = (1-q)q - q(1-q)=0 with ordering (-1,+1)
ratio=sp.factor(sp.simplify((v.T*Psi*v)[0]/(v.T*Pi*v)[0]))
print("restricted ratio r(q) =",ratio)
print("r(q)-1 =",sp.factor(sp.simplify(ratio-1)))
for val in [sp.Rational(1,2),sp.Rational(1,4),sp.Rational(1,10),sp.Rational(1,100)]:
    print(val, sp.nsimplify(ratio.subs(q,val)), float(ratio.subs(q,val)))
# also eigen/trace info
print("Psi_D - Pi =", (Psi-Pi).applyfunc(sp.factor))
