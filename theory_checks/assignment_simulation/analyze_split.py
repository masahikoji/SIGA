"""Independent selection/evaluation for a finite-n covariance ratio.
Covariance selection uses averaged centred sample covariances, as in archived code.
Evaluation exploits the exactly zero product mean; the SE clusters both products
within each independently generated triple. No iid assumption is made for the
individual components of G or for the two products within a triple.
"""
from pathlib import Path
import numpy as np
import json,csv,hashlib
from scipy.linalg import null_space, eigh
ROOT=Path(__file__).resolve().parent
fp=np.array([.5,.4,.3,.2,.1]); F=len(fp); J=1<<F
patterns=((np.arange(J)[:,None]>>np.arange(F))&1)
pi=np.prod(np.where(patterns,fp,1-fp),axis=1)
Bpi=null_space(np.sqrt(pi)[None,:]); W=Bpi/np.sqrt(pi)[:,None]
assert np.allclose(W.T@np.diag(pi)@W,np.eye(J-1))
assert np.max(np.abs(pi@W))<1e-14
records=[];covs={}
for tag,n,p in [('n400',400,.8),('n2000',2000,.8),('control',400,.5)]:
    train=np.fromfile(ROOT/f'train_{tag}.bin',dtype='<i2').reshape(-1,2,J).astype(float)/np.sqrt(n)
    evaluate=np.fromfile(ROOT/f'eval_{tag}.bin',dtype='<i2').reshape(-1,2,J).astype(float)/np.sqrt(n)
    assert train.shape[0]==200000 and evaluate.shape[0]==100000
    Psi=.5*(np.cov(train[:,0,:],rowvar=False,ddof=1)+np.cov(train[:,1,:],rowvar=False,ddof=1))
    lam,vec=eigh(W.T@Psi@W)
    direction=W@vec[:,0]; direction/=np.sqrt(np.dot(pi,direction**2))
    assert abs(pi@direction)<1e-13
    # Direction sign has no effect; orient for a canonical saved file.
    if direction[np.argmax(np.abs(direction))]<0:direction=-direction
    vals=np.einsum('bks,s->bk',evaluate,direction)
    squared=.5*np.sum(vals**2,axis=1)
    estimate=squared.mean();se=squared.std(ddof=1)/np.sqrt(len(squared))
    evalcov=.5*(np.cov(evaluate[:,0,:],rowvar=False,ddof=1)+np.cov(evaluate[:,1,:],rowvar=False,ddof=1))
    # Check independent calculation from summed raw second moments.
    second=.5*(evaluate[:,0,:].T@evaluate[:,0,:]+evaluate[:,1,:].T@evaluate[:,1,:])/len(evaluate)
    assert abs(direction@second@direction-estimate)<1e-12
    covs[tag]={'training_covariance':Psi.tolist(),'evaluation_covariance':evalcov.tolist(),
              'training_spectrum':lam.tolist(),'direction':direction.tolist(),
              'pi':pi.tolist(),'profiles':patterns.tolist()}
    record={'setting':tag,'F':F,'n':n,'pbc':p,'selection_replicates':len(train),
          'evaluation_replicates':len(evaluate),'in_sample_min':float(lam[0]),
          'independent_ratio':float(estimate),'mcse':float(se),
          'ci95_lower':float(estimate-1.959963984540054*se),
          'ci95_upper':float(estimate+1.959963984540054*se),
          'direction_mean':float(pi@direction),'direction_norm':float(pi@(direction**2)),
          'evaluation_sample_variance_ratio':float(direction@evalcov@direction),
          'evaluation_selected_min':float(eigh(W.T@evalcov@W,eigvals_only=True)[0]),
          'projection_mean01':float(vals[:,0].mean()),'projection_mean02':float(vals[:,1].mean()),
          'squared_product_correlation':float(np.corrcoef(vals[:,0]**2,vals[:,1]**2)[0,1]),
          'squared_projection_99_9pct':float(np.quantile(squared,.999))}
    for stage in ['train','eval']:
        record[stage+'_sha256']=hashlib.sha256((ROOT/f'{stage}_{tag}.bin').read_bytes()).hexdigest()
    records.append(record)
    print(json.dumps(record,indent=2))
(ROOT/'split_results.json').write_text(json.dumps(records,indent=2)+'\n')
(ROOT/'directions_and_covariances.json').write_text(json.dumps(covs,indent=2)+'\n')
with (ROOT/'split_results.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=records[0].keys());w.writeheader();w.writerows(records)
