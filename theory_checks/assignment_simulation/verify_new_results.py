"""Checks for the new allocation-only computations and exact covariance identity."""
from pathlib import Path
from fractions import Fraction as Q
import itertools,json,subprocess,hashlib,tempfile
import numpy as np
ROOT=Path(__file__).resolve().parent
checks={}
exe=ROOT/'simulate_sign_products'
with tempfile.TemporaryDirectory() as folder:
 folder=Path(folder)
 args=['5','400','64','0.8','72635']
 subprocess.run([str(exe),*args,'1',str(folder/'one.bin')],check=True,capture_output=True)
 subprocess.run([str(exe),*args,'4',str(folder/'four.bin')],check=True,capture_output=True)
 assert (folder/'one.bin').read_bytes()==(folder/'four.bin').read_bytes()
 checks['thread_count_invariance']=True
 # Recompile the generator using direct absolute candidate scores at every step.
 cpp=(ROOT/'simulate_sign_products.cpp').read_text()
 for state in ['st0','st1','st2']:cpp=cpp.replace(f'assign({state},',f'reference({state},')
 (folder/'reference.cpp').write_text(cpp)
 subprocess.run(['g++','-O3','-std=c++17','-fopenmp',str(folder/'reference.cpp'),'-o',str(folder/'reference')],check=True)
 subprocess.run([str(folder/'reference'),*args,'2',str(folder/'reference.bin')],check=True,capture_output=True)
 assert (folder/'one.bin').read_bytes()==(folder/'reference.bin').read_bytes()
 checks['whole_stream_absolute_rule_equivalence']=True

# Exact enumeration for two independent strata, sizes 2 and 3, with rational
# label-symmetric allocation laws. Products in distinct repetitions are sampled
# independently under the same fixed profile sequence.
states=list(itertools.product([-1,1],repeat=5));S=[0,0,1,1,1];n=5
law=[]
for z in states:
 p0=(Q(1)+Q(-1,3)*z[0]*z[1])/4
 p1=(Q(1)+Q(1,5)*z[2]*z[3]+Q(-1,7)*z[3]*z[4])/8
 law.append(p0*p1)
assert sum(law)==1 and min(law)>0
C=[[sum(p*z[i]*z[j] for p,z in zip(law,states)) for j in range(n)] for i in range(n)]
assert all(C[i][j]==0 for i in range(n) for j in range(n) if S[i]!=S[j])
product=[[Q(0) for _ in range(2)] for _ in range(2)]
for z,pz in zip(states,law):
 for w,pw in zip(states,law):
  sums=[sum(z[i]*w[i] for i in range(n) if S[i]==s) for s in range(2)]
  for s in range(2):
   for t in range(2):product[s][t]+=pz*pw*sums[s]*sums[t]/n
formula=[[sum(C[i][j]**2 for i in range(n) for j in range(n) if S[i]==s and S[j]==t)/n for t in range(2)] for s in range(2)]
assert product==formula
excess=[[product[s][t]-(Q(S.count(s),n) if s==t else 0) for t in range(2)] for s in range(2)]
assert excess[0][1]==0 and excess[1][0]==0 and excess[0][0]>=0 and excess[1][1]>=0
checks['stratified_covariance_exact_enumeration']={'configurations':len(states)**2,'excess_matrix':[[str(x) for x in row] for row in excess]}

# Recompute projected second/fourth moments without einsum or covariance calls.
# Each row is one independent triple; the within-triple products stay grouped.
summaries=json.loads((ROOT/'split_results.json').read_text())
covs=json.loads((ROOT/'directions_and_covariances.json').read_text())
for r in summaries:
 tag=r['setting'];n=r['n'];d=np.array(covs[tag]['direction'])
 a=np.fromfile(ROOT/f'eval_{tag}.bin',dtype='<i2').reshape(-1,2,32)
 first=a[:,0,:]@d;second=a[:,1,:]@d
 moments=(first*first+second*second)/(2*n)
 ratio=float(sum(moments)/len(moments))
 var=float(sum((moments-ratio)**2)/(len(moments)-1))
 se=(var/len(moments))**.5
 assert abs(ratio-r['independent_ratio'])<5e-13
 assert abs(se-r['mcse'])<5e-14
 assert r['direction_norm']>0.999999999999
 checks[f'independent_summary_{tag}']={'ratio':ratio,'mcse':se}
checks['passed']=True
(ROOT/'NEW_CHECKS.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
