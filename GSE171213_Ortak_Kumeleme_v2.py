import os,glob,gzip,time,json,hashlib
os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['OMP_NUM_THREADS']='2'
import numpy as np,pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize
from sklearn.cluster import MiniBatchKMeans
from sklearn.neighbors import NearestNeighbors
from scipy.stats import ranksums
base='/mnt/data'; out=base+'/GSE171213_'; t=time.time()
a=pd.read_csv(base+'/GSE171213_43438_Korunan_Hucre_Anotasyon_v1.csv',low_memory=False)
parts=[]; gene_sets=[]; sample_info=[]
for sample,grp in a.groupby('ornek',sort=False):
 p=glob.glob(base+f'/files/GSM*_{sample}.counts.tsv.gz');assert len(p)==1,(sample,p)
 p=p[0]; keep=grp.CellIndex.astype(str).tolist(); pos={c:i for i,c in enumerate(keep)}
 with gzip.open(p,'rt') as f: head=f.readline().rstrip('\n').split('\t')
 idx=[head.index(c) for c in keep]; assert len(set(idx))==len(idx)
 # Read rows in blocks, selecting only QC-kept columns, preserving gene labels
 rows=[];cols=[];values=[];genes=[]; gene_i=0
 with gzip.open(p,'rt') as f:
  next(f)
  for line in f:
   fields=line.rstrip('\n').split('\t'); gene=fields[0];genes.append(gene)
   arr=np.fromiter((float(fields[j]) for j in idx),dtype=np.float32,count=len(idx))
   nz=np.flatnonzero(arr)
   rows.extend([gene_i]*len(nz));cols.extend(nz.tolist());values.extend(arr[nz].tolist());gene_i+=1
 mat=sparse.coo_matrix((np.array(values,dtype=np.float32),(rows,cols)),shape=(len(genes),len(keep))).tocsr()
 del rows,cols,values
 assert len(set(genes))==len(genes),(sample,'duplicate genes')
 parts.append((sample,mat,genes,grp.index.to_numpy()));gene_sets.append(set(genes));sample_info.append((sample,len(genes),len(keep),mat.nnz,os.path.basename(p)))
 print(sample,mat.shape,mat.nnz,'elapsed',round(time.time()-t),flush=True)
common=sorted(set.intersection(*gene_sets));print('COMMON',len(common),flush=True)
assert len(common)>10000
mats=[]; order=[]
for sample,mat,genes,ix in parts:
 gi={g:i for i,g in enumerate(genes)}
 mats.append(mat[[gi[g] for g in common],:].T.tocsr());order.extend(ix.tolist())
 del mat
X=sparse.vstack(mats,format='csr');del mats,parts
assert X.shape==(len(a),len(common));a=a.loc[order].reset_index(drop=True)
# normalize each cell to 10k, log1p all expressed genes
lib=np.asarray(X.sum(axis=1)).ravel();assert np.all(lib>0)
Y=X.copy();Y.data *= np.repeat(1e4/lib,np.diff(Y.indptr));Y.data=np.log1p(Y.data)
# pooled variance per gene, exclude mitochondrial and ribosomal genes from dimensionality reduction
mean=np.asarray(Y.mean(axis=0)).ravel();sq=np.asarray(Y.power(2).mean(axis=0)).ravel();var=sq-mean**2
valid=np.array([not (g.upper().startswith(('MT-','RPS','RPL'))) for g in common]);var[~valid]=-1
hvg=np.argsort(var)[-2500:][::-1];genes_hvg=np.array(common)[hvg]
# center implicitly by TruncatedSVD (not PCA); no dense full matrix
Z=TruncatedSVD(n_components=35,random_state=42,n_iter=5).fit_transform(Y[:,hvg]);Z=Z[:,1:31]
# remove per-sample technical shift in latent representation by sample mean centering (exploratory only)
Z=Z.astype('float32'); Z=(Z-Z.mean(axis=0))/(Z.std(axis=0)+1e-6)
# Cluster global gene-expression embedding; robust non-graph provisional transcriptomic clusters
km=MiniBatchKMeans(n_clusters=24,random_state=42,batch_size=2048,n_init=5,max_iter=150)
labels=km.fit_predict(Z);a['transkriptom_kume']=labels
# group/sample composition and fibro marker profile from raw normalized log expression
fib=['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA','FAP','THY1','DPT','COL6A1']
neg={'Immune':['PTPRC','LST1','CD3D','MS4A1'],'Epithelial':['EPCAM','KRT14','KRT5','KRT19'],'Endothelial':['PECAM1','VWF','KDR'],'Pericyte':['RGS5','MCAM','ACTA2','TAGLN'],'Plasma':['MZB1','JCHAIN'],'Mast':['TPSAB1','KIT']}
lookup={g:i for i,g in enumerate(common)}
markers=sorted(set(fib+sum(neg.values(),[])));present=[g for g in markers if g in lookup]
M=Y[:,[lookup[g] for g in present]].toarray();midx={g:i for i,g in enumerate(present)}
for name,gs in [('Fibroblast',fib)]+list(neg.items()):
 ids=[midx[g] for g in gs if g in midx]
 a[f'{name}_log_score']=M[:,ids].mean(axis=1) if ids else 0
 a[f'{name}_genes_detected']=(M[:,ids]>0).sum(axis=1) if ids else 0
# cluster-level fibro enrichment and negative markers; don't overcall ambiguous clusters
stats=[]
for c in range(24):
 q=a[a.transkriptom_kume==c];negscore=max(q[f'{n}_log_score'].mean() for n in neg)
 stats.append({'kume':c,'hucre':len(q),'fibro_score':q.Fibroblast_log_score.mean(),'fibro_2marker_orani':(q.Fibroblast_genes_detected>=2).mean(),'neg_max_score':negscore,'fibro_aday_orani':q.fibroblast_candidate_strict.astype(str).str.lower().isin(['true','1']).mean(),**{f'{n}_score':q[f'{n}_log_score'].mean() for n in neg}})
s=pd.DataFrame(stats)
# transparent provisional rule based on positive/negative marker evidence and cluster context
s['fibro_kume_adayi']=(s.fibro_score>s.neg_max_score*1.5)&(s.fibro_2marker_orani>=0.25)&(s.fibro_score>=0.15)
a['fibro_kume_adayi']=a.transkriptom_kume.map(s.set_index('kume').fibro_kume_adayi)
a['fibro_yuksek_guven_adayi']=a.fibro_kume_adayi & (a.Fibroblast_genes_detected>=2) & (a.Fibroblast_log_score>a[[f'{n}_log_score' for n in neg]].max(axis=1))
a.to_csv(out+'43438_Ortak_Transkriptom_Kumeleri_v2.csv',index=False)
s.to_csv(out+'Ortak_Kume_Belirtec_Ozeti_v2.csv',index=False)
pd.DataFrame(sample_info,columns=['ornek','gen_sayisi','hucre_sayisi','nnz','ham_dosya']).to_csv(out+'Kumeleme_Girdi_Denetimi_v2.csv',index=False)
pd.DataFrame({'gene':genes_hvg,'pooled_log_variance':var[hvg]}).to_csv(out+'Kumeleme_HVG_2500_v2.csv',index=False)
np.savez_compressed(out+'Kumeleme_Latent_v2.npz',embedding=Z,clusters=labels,samples=a.ornek.to_numpy())
summary={'cells':len(a),'common_genes':len(common),'hvg':len(hvg),'clusters':24,'fibro_candidate_clusters':s.loc[s.fibro_kume_adayi,'kume'].tolist(),'fibro_cluster_cells':int(a.fibro_kume_adayi.sum()),'fibro_high_confidence_candidates':int(a.fibro_yuksek_guven_adayi.sum()),'elapsed_seconds':round(time.time()-t,1),'note':'Exploratory pooled TruncatedSVD + MiniBatchKMeans; no batch correction, no independent reference annotation; do not call final.'}
open(out+'Ortak_Kumeleme_Durum_v2.json','w').write(json.dumps(summary,indent=2,ensure_ascii=False));print(summary,flush=True)
