import gzip, numpy as np, pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from pathlib import Path
P=Path('/mnt/data'); rng=np.random.default_rng(171213)
raw=P/'files/GSM5220920_HC1.counts.tsv.gz'
qc=pd.read_csv(P/'GSE171213_Ham_QC_Hucre_Metrikleri_v2.csv',encoding='utf-8-sig')
qc=qc[(qc['ornek']=='HC1') & (qc['QC_MT30'].astype(str).str.lower()=='true')].copy()
with gzip.open(raw,'rt') as f: header=f.readline().rstrip('\n').split('\t')
cols=header[1:]; selected=set(qc['CellIndex']); inds=[i for i,c in enumerate(cols) if c in selected]
assert len(inds)==len(qc), (len(inds),len(qc))
# Only selected QC cells, preserve raw gene count rows
rows=[]; ci=[]; vals=[]; genes=[]
with gzip.open(raw,'rt') as f:
 next(f)
 for gi,line in enumerate(f):
  parts=line.rstrip('\n').split('\t'); genes.append(parts[0]); arr=np.asarray(parts[1:],dtype=np.float32)[inds]
  nz=np.flatnonzero(arr)
  if len(nz): rows.extend([gi]*len(nz)); ci.extend(nz.tolist()); vals.extend(arr[nz].tolist())
X=sparse.coo_matrix((np.asarray(vals,dtype=np.float32),(rows,ci)),shape=(len(genes),len(inds))).tocsr().T.tocsr()
assert X.shape[0]==len(qc)
# Standard library-size normalization and log1p; synthetic doublets from raw count sums
lib=np.asarray(X.sum(axis=1)).ravel(); assert np.all(lib>0)
S=X[rng.integers(0,len(inds),len(inds))]+X[rng.integers(0,len(inds),len(inds))]
A=sparse.vstack([X,S]).tocsr()
lib_all=np.asarray(A.sum(axis=1)).ravel(); A=sparse.diags(1e4/lib_all)@A
A.data=np.log1p(A.data)
# Variable gene selection using sparse moments; fit SVD on observed + simulated in common embedding
mean=np.asarray(A.mean(axis=0)).ravel(); sq=np.asarray(A.power(2).mean(axis=0)).ravel(); var=np.maximum(sq-mean**2,0)
idx=np.argsort(var)[-2000:]; A=A[:,idx]
Z=TruncatedSVD(n_components=30,random_state=171213).fit_transform(A)
Z=Z/(Z.std(axis=0)+1e-8)
knn=NearestNeighbors(n_neighbors=41,algorithm='auto',n_jobs=2).fit(Z)
dist,neighbors=knn.kneighbors(Z[:len(inds)])
score=(neighbors[:,1:]>=len(inds)).mean(axis=1)
result=qc[['ornek','CellIndex','CellName','n_genes_raw','n_UMI_raw','mito_pct_raw']].copy()
score_by_index=dict(zip(np.array(cols)[inds],score)); result['simulated_neighbor_fraction']=result.CellIndex.map(score_by_index)
result['rank_percentile']=result.simulated_neighbor_fraction.rank(pct=True)
result.to_csv(P/'GSE171213_HC1_Doublet_Skorlari_Pilot_v1.csv',index=False)
summary=pd.DataFrame([{'ornek':'HC1','QC30_cells':len(inds),'simulated_doublets':len(inds),'genes_raw':len(genes),'nonzero_entries':X.nnz,'embedding_components':30,'variable_genes':2000,'neighbors':40,'seed':171213,'median_simulated_neighbor_fraction':float(np.median(score)),'p95_simulated_neighbor_fraction':float(np.quantile(score,.95)),'doublets_removed':0,'classification':'not calibrated; pilot score only'}]);summary.to_csv(P/'GSE171213_HC1_Doublet_Yontem_Ozeti_v1.csv',index=False)
print(summary.to_string(index=False))
