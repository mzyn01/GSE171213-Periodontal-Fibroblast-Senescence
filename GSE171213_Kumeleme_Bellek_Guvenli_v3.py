"""GSE171213: memory-safe per-sample sparse checkpoints, then joint analysis.
Requires GSE171213_43438_Korunan_Hucre_Anotasyon_v1.csv and 12 raw counts files.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['OMP_NUM_THREADS']='2'
import glob, gzip, json, hashlib
import numpy as np,pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import MiniBatchKMeans
BASE='/mnt/data'; PREFIX=BASE+'/GSE171213_'; CACHE=BASE+'/gse171213_sparse_cache';os.makedirs(CACHE,exist_ok=True)
a=pd.read_csv(PREFIX+'43438_Korunan_Hucre_Anotasyon_v1.csv',low_memory=False)
assert len(a)==43438 and not a[['ornek','CellIndex']].duplicated().any()
parts=[]; genesets=[]; orders=[]; manifest=[]
for sample,group in a.groupby('ornek',sort=False):
    f=glob.glob(f'{BASE}/files/GSM*_{sample}.counts.tsv.gz');assert len(f)==1
    keep=group.CellIndex.astype(str).tolist()
    with gzip.open(f[0],'rt') as h: header=h.readline().rstrip('\n').split('\t')
    colmap={s:i for i,s in enumerate(header)}
    assert len(set(keep))==len(keep) and all(k in colmap for k in keep)
    # pandas processes genes in bounded row blocks, never accumulating Python lists of individual nonzeros
    matrices=[];genes=[]
    with gzip.open(f[0],'rt') as h:
        for chunk in pd.read_csv(h,sep='\t',chunksize=300,low_memory=False):
            names=chunk.iloc[:,0].astype(str).tolist();genes.extend(names)
            values=chunk[keep].to_numpy(dtype=np.float32,copy=True)
            assert np.isfinite(values).all() and (values>=0).all()
            matrices.append(sparse.csr_matrix(values))
    mat=sparse.vstack(matrices,format='csr');del matrices
    assert len(genes)==mat.shape[0] and len(set(genes))==len(genes)
    sparse.save_npz(f'{CACHE}/{sample}.npz',mat)
    pd.Series(genes).to_csv(f'{CACHE}/{sample}_genes.csv',index=False,header=False)
    parts.append((sample,group.index.to_numpy()));genesets.append(set(genes))
    manifest.append({'sample':sample,'cells':len(keep),'genes':len(genes),'nnz':int(mat.nnz),'source':os.path.basename(f[0])})
    print('CHECKPOINT',sample,mat.shape,mat.nnz,flush=True)
common=sorted(set.intersection(*genesets));assert len(common)>10000
mats=[]
for sample,ix in parts:
    mat=sparse.load_npz(f'{CACHE}/{sample}.npz')
    genes=pd.read_csv(f'{CACHE}/{sample}_genes.csv',header=None).iloc[:,0].astype(str).tolist()
    lookup={g:i for i,g in enumerate(genes)}
    mats.append(mat[[lookup[g] for g in common]].T.tocsr());orders.extend(ix.tolist())
X=sparse.vstack(mats,format='csr');del mats
assert X.shape==(43438,len(common))
a=a.loc[orders].reset_index(drop=True)
lib=np.asarray(X.sum(axis=1)).ravel();assert np.all(lib>0)
Y=X.copy();Y.data*=np.repeat((1e4/lib).astype(np.float32),np.diff(Y.indptr));np.log1p(Y.data,out=Y.data)
mean=np.asarray(Y.mean(axis=0)).ravel();var=np.asarray(Y.power(2).mean(axis=0)).ravel()-mean**2
valid=np.array([not g.upper().startswith(('MT-','RPS','RPL')) for g in common]);var[~valid]=-np.inf
hvg=np.argsort(var)[-2500:][::-1]; Z=TruncatedSVD(n_components=35,random_state=42,n_iter=5).fit_transform(Y[:,hvg])[:,1:31]
Z=((Z-Z.mean(0))/(Z.std(0)+1e-6)).astype(np.float32)
labels=MiniBatchKMeans(n_clusters=24,random_state=42,batch_size=2048,n_init=5,max_iter=150).fit_predict(Z)
a['transkriptom_kume']=labels
fib=['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA','FAP','THY1','DPT','COL6A1']
neg={'Immune':['PTPRC','LST1','CD3D','MS4A1'],'Epithelial':['EPCAM','KRT14','KRT5','KRT19'],'Endothelial':['PECAM1','VWF','KDR'],'Pericyte':['RGS5','MCAM','ACTA2','TAGLN'],'Plasma':['MZB1','JCHAIN'],'Mast':['TPSAB1','KIT']}
lookup={g:i for i,g in enumerate(common)}
for name,gs in [('Fibroblast',fib)]+list(neg.items()):
    ids=[lookup[g] for g in gs if g in lookup];sub=Y[:,ids]
    a[name+'_log_score']=np.asarray(sub.mean(axis=1)).ravel() if ids else 0
    a[name+'_genes_detected']=np.asarray((sub>0).sum(axis=1)).ravel() if ids else 0
stats=[]
for c,q in a.groupby('transkriptom_kume'):
    stats.append({'cluster':int(c),'cells':len(q),'fibro_score':q.Fibroblast_log_score.mean(),'fibro_two_marker_fraction':(q.Fibroblast_genes_detected>=2).mean(),'negative_max_score':max(q[n+'_log_score'].mean() for n in neg),**{n+'_score':q[n+'_log_score'].mean() for n in neg}})
s=pd.DataFrame(stats);s['fibro_cluster_candidate']=(s.fibro_score>1.5*s.negative_max_score)&(s.fibro_two_marker_fraction>=.25)&(s.fibro_score>=.15)
a['fibro_cluster_candidate']=a.transkriptom_kume.map(s.set_index('cluster').fibro_cluster_candidate)
a['fibro_high_confidence_candidate']=a.fibro_cluster_candidate&(a.Fibroblast_genes_detected>=2)&(a.Fibroblast_log_score>a[[n+'_log_score' for n in neg]].max(axis=1))
a.to_csv(PREFIX+'43438_Ortak_Transkriptom_Kumeleri_v3.csv',index=False)
s.to_csv(PREFIX+'Ortak_Kume_Belirtec_Ozeti_v3.csv',index=False)
pd.DataFrame(manifest).to_csv(PREFIX+'Kumeleme_Girdi_Denetimi_v3.csv',index=False)
pd.DataFrame({'gene':np.array(common)[hvg],'variance':var[hvg]}).to_csv(PREFIX+'Kumeleme_HVG_v3.csv',index=False)
np.savez_compressed(PREFIX+'Kumeleme_Latent_v3.npz',embedding=Z,cluster=labels)
summary={'cells':len(a),'common_genes':len(common),'clusters':len(s),'fibro_cluster_cells':int(a.fibro_cluster_candidate.sum()),'high_confidence_fibro_candidates':int(a.fibro_high_confidence_candidate.sum()),'method':'exploratory pooled SVD plus MiniBatchKMeans; not batch corrected or reference validated'}
with open(PREFIX+'Kumeleme_Durum_v3.json','w') as h:json.dump(summary,h,indent=2)
print('DONE',summary,flush=True)
