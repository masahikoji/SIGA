import itertools
def sgn(t): return (t>0)-(t<0)
def fval(x,u):
    s=sgn(x[0])+sum(sgn(x[0]+u[f]*x[f+1]) for f in range(len(u)))
    return sgn(s)
def reachable(F, cap=200000, maxnorm=60):
    profs=list(itertools.product([-1,1],repeat=F))
    zero=tuple([0]*(F+1)); seen={zero}; frontier=[zero]; overflow=False
    while frontier:
        nxt=[]
        for x in frontier:
            for u in profs:
                f=fval(x,u)
                zs=[-f] if f!=0 else [1,-1]   # p=1: deterministic unless tie
                for z in zs:
                    y=tuple(x[k]+z*(1 if k==0 else u[k-1]) for k in range(F+1))
                    if max(abs(c) for c in y)>maxnorm: overflow=True; continue
                    if y not in seen:
                        seen.add(y); nxt.append(y)
        frontier=nxt
        if len(seen)>cap: return seen, True
    return seen, overflow
for F in [1,2,3]:
    S,ov=reachable(F)
    print(F, len(S), "overflow" if ov else "finite", "max|x| =", max(max(abs(c) for c in x) for x in S))
    if F==1: print(sorted(S))
