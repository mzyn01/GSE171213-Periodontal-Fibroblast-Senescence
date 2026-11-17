import os
os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['OMP_NUM_THREADS']='2'
import glob,gzip,time
import numpy as np,pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score
B='/mnt/data'; a=pd.read_csv(B+'/GSE171213_43438_Ortak_Transkriptom_Kumeleri_v3.csv',low_memory=False)
hvg=pd.read_csv(B+'/GSE171213_Kumeleme_HVG_v3.csv').gene.astype(str).tolist()
assert len(a)==43438 and len(hvg)==2500
parts=[]; indices=[]; samples=[]; audit=[]
for sample,g in a.groupby('ornek',sort=False):
 f=glob.glob(f'{B}/files/GSM*_{sample}.counts.tsv.gz');assert len(f)==1
 with gzip.open(f[0],'rt') as h: header=h.readline().rstrip('\n').split('\t')
 colmap={v:i for i,v in enumerate(header)}; keep=g.CellIndex.astype(str).tolist();assert all(x in colmap for x in keep)
 pos=[colmap[x] for x in keep]
 genes=[]; blocks=[]
 with gzip.open(f[0],'rt') as h:
  for chunk in pd.read_csv(h,sep='\t',chunksize=200,low_memory=False):
   z=chunk.iloc[:,0].astype(str).isin(hvg)
   if z.any():
    genes.extend(chunk.loc[z].iloc[:,0].astype(str).tolist())
    blocks.append(sparse.csr_matrix(chunk.loc[z].iloc[:,pos].to_numpy(dtype=np.float32)))
 mat=sparse.vstack(blocks,format='csr'); loc={v:i for i,v in enumerate(genes)}
 assert all(x in loc for x in hvg),(sample,len(genes))
 # normalize using original full-library sizes, avoiding HVG-only library bias
 mat=mat[[loc[x] for x in hvg]].T.tocsr()
 lib=g.n_UMI_raw.to_numpy(dtype=np.float32);assert (lib>0).all()
 mat.data *= np.repeat(1e4/lib,np.diff(mat.indptr));np.log1p(mat.data,out=mat.data)
 parts.append(mat);indices.extend(g.index);samples.extend([sample]*len(g));audit.append((sample,len(g),mat.nnz))
 print('READ',sample,len(g),flush=True)
X=sparse.vstack(parts,format='csr');del parts
assert X.shape==(43438,2500)
a=a.loc[indices].reset_index(drop=True)
Z=TruncatedSVD(n_components=31,n_iter=4,random_state=42).fit_transform(X)[:,1:31].astype(np.float32)
Z=(Z-Z.mean(0))/(Z.std(0)+1e-6)
# sample-centered latent representation: sensitivity analysis, NOT validated batch correction
sample=np.asarray(samples);Zc=Z.copy()
for s in np.unique(sample):
 m=sample==s;Zc[m]-=Z[m].mean(0)
Zc=(Zc-Zc.mean(0))/(Zc.std(0)+1e-6)
km=MiniBatchKMeans(n_clusters=24,random_state=42,n_init=5,batch_size=2048,max_iter=150)
alt=km.fit_predict(Zc)
a['batch_centered_cluster_v6']=alt
# cluster stability and fibroblast marker enrichment using precomputed per-cell independent marker readouts
res=[]
for c,g in a.groupby('batch_centered_cluster_v6'):
 counts=g.ornek.value_counts(normalize=True)
 res.append({'cluster':c,'cells':len(g),'dominant_sample':counts.index[0],'dominant_sample_fraction':counts.iloc[0], 'fibro_score_mean':g.Fibroblast_log_score.mean(),'fibro_two_marker_fraction':(g.Fibroblast_genes_detected>=2).mean(),'old_fibro_candidate_fraction':g.fibro_high_confidence_candidate.mean(),'old_cluster10_fraction':(g.transkriptom_kume==10).mean()})
r=pd.DataFrame(res)
# cluster-level provisional fibro identity; not external reference
r['fibro_cluster_candidate']=(r.fibro_score_mean>.15)&(r.fibro_two_marker_fraction>=.25)
a['centered_fibro_candidate_v6']=a.batch_centered_cluster_v6.isin(r.loc[r.fibro_cluster_candidate,'cluster']) & (a.Fibroblast_genes_detected>=2) & (a.Fibroblast_log_score>a[[k+'_log_score' for k in ['Immune','Epithelial','Endothelial','Pericyte','Plasma','Mast']]].max(axis=1))
old=a.fibro_high_confidence_candidate.astype(bool)
new=a.centered_fibro_candidate_v6.astype(bool)
metrics={'cells':len(a),'old_fibro':int(old.sum()),'centered_fibro':int(new.sum()),'overlap':int((old&new).sum()),'old_retention_fraction':float((old&new).sum()/old.sum()),'cluster_ARI':float(adjusted_rand_score(a.transkriptom_kume,alt)),'method':'sample-centered SVD sensitivity only; no external reference atlas; no validated batch correction'}
a[['ornek','CellIndex','CellName','transkriptom_kume','batch_centered_cluster_v6','fibro_high_confidence_candidate','centered_fibro_candidate_v6']].to_csv(B+'/GSE171213_Batch_Merkezleme_Hucre_v6.csv',index=False)
r.to_csv(B+'/GSE171213_Batch_Merkezleme_Kume_v6.csv',index=False)
pd.DataFrame([metrics]).to_csv(B+'/GSE171213_Batch_Merkezleme_Ozet_v6.csv',index=False)
print('DONE',metrics,flush=True)
