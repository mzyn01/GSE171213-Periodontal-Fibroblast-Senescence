import os,gzip,glob,zipfile,time
import numpy as np,pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score,adjusted_rand_score
from sklearn.neighbors import NearestNeighbors
from scipy.stats import mannwhitneyu
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT='/mnt/data/GSE171213_Fibro_AltKume_v1';os.makedirs(OUT,exist_ok=True)
with zipfile.ZipFile('/mnt/data/GSE171213_Harmony_Sonuc_Paketi.zip') as z:
 meta=pd.read_csv(z.open('GSE171213_Harmony_Hucre_Etiketleri.csv'))
meta=meta.loc[meta.harmony_kume24.astype(str)=='4'].copy().reset_index(drop=True)
print('Selected',len(meta),flush=True)
matrices=[];gene_sets=[];indices=[];qc=[]
for sample in ['HC1','HC2','HC3','HC4','PD1','PD2','PD3','PD4','PD5','PDT1','PDT2','PDT3']:
 t=time.time();path=glob.glob('/mnt/data/files/*_'+sample+'.counts.tsv.gz')[0]
 sel=meta.index[meta.ornek==sample].to_numpy(); cells=meta.loc[sel,'CellIndex'].astype(str).tolist()
 # Header names from compressed input
 with gzip.open(path,'rt') as h:header=h.readline().rstrip('\n').split('\t')
 missing=set(cells)-set(header);assert not missing,(sample,missing)
 cols=['gene']+cells;df=pd.read_csv(path,sep='\t',usecols=cols,compression='gzip',dtype={c:'int32' for c in cells})
 assert df.gene.is_unique
 genes=df.gene.astype(str).to_numpy();arr=sparse.csr_matrix(df[cells].to_numpy(dtype=np.float32).T)
 matrices.append(arr);gene_sets.append(genes);indices.append(sel)
 qc.append({'ornek':sample,'fibro_hucre':len(sel),'ham_gen':len(genes),'okuma_sn':round(time.time()-t,1)})
 print(sample,len(sel),len(genes),round(time.time()-t,1),flush=True)
common=sorted(set.intersection(*(set(g) for g in gene_sets)));print('common genes',len(common),flush=True)
parts=[]
for a,genes in zip(matrices,gene_sets):
 lookup=pd.Index(genes).get_indexer(common);assert (lookup>=0).all();parts.append(a[:,lookup])
X=sparse.vstack(parts,format='csr');order=np.concatenate(indices);meta=meta.iloc[order].reset_index(drop=True)
assert X.shape==(1343,len(common));pd.DataFrame(qc).to_csv(OUT+'/Ham_Girdi_Denetimi.csv',index=False)
# remove genes seen in fewer than 5 cells, normalize by total counts across common genes
umi=np.asarray(X.sum(axis=1)).ravel();assert (umi>0).all()
Y=X.astype(np.float32).copy();Y.data*=np.repeat(1e4/umi,np.diff(Y.indptr));Y.data=np.log1p(Y.data)
express=np.asarray((X>0).sum(axis=0)).ravel();valid=express>=5
means=np.asarray(Y.mean(axis=0)).ravel();sq=np.asarray(Y.power(2).mean(axis=0)).ravel();variance=np.maximum(0,sq-means**2)
# high variance genes without technical covariate correction, exploratory
v=np.where(valid)[0];hvg=v[np.argsort(variance[v])[-min(2000,len(v)):]]
P=TruncatedSVD(n_components=25,random_state=42).fit_transform(Y[:,hvg]);P=(P-P.mean(0))/(P.std(0)+1e-6)
res=[];fits={}
for k in range(2,7):
 labels=KMeans(n_clusters=k,n_init=15,random_state=42).fit_predict(P)
 sil=silhouette_score(P,labels, sample_size=min(1100,len(P)),random_state=42)
 other=KMeans(n_clusters=k,n_init=15,random_state=77).fit_predict(P)
 ari=adjusted_rand_score(labels,other)
 sizes=np.bincount(labels,minlength=k)
 res.append({'k':k,'silhouette':sil,'seed_ARI':ari,'min_cluster':sizes.min(),'max_cluster':sizes.max()});fits[k]=labels
pd.DataFrame(res).to_csv(OUT+'/K_Secim_Denetimi.csv',index=False)
# prefer stable, not tiny clusters and reasonable silhouette
r=pd.DataFrame(res);c=r[(r.min_cluster>=35)&(r.seed_ARI>=0.75)];k=int(c.sort_values(['silhouette','seed_ARI'],ascending=False).iloc[0].k) if len(c) else 3
lab=fits[k];meta['fibro_alt_kume']=lab;meta['grup']=meta.ornek.str.extract(r'^(PDT|PD|HC)')[0]
meta.to_csv(OUT+'/Fibroblast_1343_AltKume_Etiketleri.csv',index=False)
ct=pd.crosstab([meta.grup,meta.ornek],meta.fibro_alt_kume);ct.to_csv(OUT+'/Ornek_AltKume_Dagilimi.csv')
# unbiased gene marker rankings: mean log normalized enrichment, expression fraction
markers=[]
for j in range(k):
 inside=(lab==j);outside=~inside
 m1=np.asarray(Y[inside].mean(axis=0)).ravel();m0=np.asarray(Y[outside].mean(axis=0)).ravel()
 d1=np.asarray((X[inside]>0).mean(axis=0)).ravel();d0=np.asarray((X[outside]>0).mean(axis=0)).ravel()
 delta=m1-m0
 score=delta*np.sqrt(np.maximum(0,d1-d0))
 top=np.argsort(score)[-30:][::-1]
 for rank,g in enumerate(top,1):markers.append({'alt_kume':j,'sira':rank,'gen':common[g],'log_norm_fark':round(float(delta[g]),4),'pozitif_oran_kume':round(float(d1[g]),4),'pozitif_oran_diger':round(float(d0[g]),4)})
pd.DataFrame(markers).to_csv(OUT+'/AltKume_Ayirici_Genler.csv',index=False)
# gene panels, exploratory not senescence diagnostic
panels={'fibroblast':['COL1A1','COL1A2','DCN','LUM','PDGFRA','COL3A1'], 'senesans_hucre_dongusu':['CDKN1A','CDKN2A','CDKN2B','GLB1','GDF15'], 'SASP':['IL6','IL8','CXCL1','CXCL2','CCL2','MMP1','MMP3','SERPINE1'], 'proliferasyon':['MKI67','TOP2A','PCNA'], 'inflamasyon':['NFKBIA','TNFAIP3','ICAM1']}
lookup={g:i for i,g in enumerate(common)};rows=[]
for name,genes in panels.items():
 avail=[g for g in genes if g in lookup];ix=[lookup[g] for g in avail]
 vals=np.asarray(Y[:,ix].mean(axis=1)).ravel() if ix else np.full(len(meta),np.nan)
 meta[name+'_skor']=vals
 for j in range(k):
  mask=lab==j
  rows.append({'panel':name,'alt_kume':j,'n':mask.sum(),'gen_sayisi':len(avail),'mevcut_genler':';'.join(avail),'ortalama_skor':float(np.mean(vals[mask])),'medyan_skor':float(np.median(vals[mask]))})
pd.DataFrame(rows).to_csv(OUT+'/Senesans_SASP_Panel_AltKume.csv',index=False)
meta.to_csv(OUT+'/Fibroblast_1343_AltKume_Skorlari.csv',index=False)
# sample level scores (avoid pseudoreplication)
meta.groupby(['grup','ornek']).agg(n=('CellIndex','size'),**{name+'_ortalama':(name+'_skor','mean') for name in panels}).reset_index().to_csv(OUT+'/Ornek_Bazli_Skorlar.csv',index=False)
# diagnostic UMAP, non-batch-corrected subset
try:
 import umap
 embedding=umap.UMAP(n_neighbors=20,min_dist=0.3,random_state=42,n_jobs=1).fit_transform(P)
except Exception:
 embedding=P[:,:2]
fig,axes=plt.subplots(1,3,figsize=(16,4.8))
for ax,col,title in zip(axes,[lab,meta.grup,meta.senesans_hucre_dongusu_skor],['Fibroblast alt kumeleri','HC / PD / PDT','Senesans aday panel skoru']):
 if title=='HC / PD / PDT':
  for group in ['HC','PD','PDT']:
   mask=(meta.grup==group).to_numpy();ax.scatter(embedding[mask,0],embedding[mask,1],s=3,alpha=.65,label=group)
  ax.legend(markerscale=3)
 else:
  sc=ax.scatter(embedding[:,0],embedding[:,1],c=col,s=3,cmap='viridis',alpha=.75)
  fig.colorbar(sc,ax=ax,shrink=.7)
 ax.set_title(title);ax.set_xticks([]);ax.set_yticks([])
fig.tight_layout();fig.savefig(OUT+'/Fibroblast_AltKume_UMAP_300dpi.png',dpi=300);plt.close(fig)
with open(OUT+'/Turkce_Analiz_Raporu.md','w') as f:
 f.write(f'# GSE171213 – Fibroblast alt küme analizi (keşifsel)\n\nHam sayımlardan 12 örnekteki 1.343 Harmony küme 4 hücresi çıkarıldı. {len(common)} ortak gen; {len(hvg)} yüksek varyanslı gen, 25 boyutlu SVD ve KMeans ile k={k} alt küme seçildi.\n\n')
 f.write('## Küme seçimi\n'+r.to_markdown(index=False)+'\n\n')
 f.write('## Alt küme sayıları\n'+meta.groupby(['grup','fibro_alt_kume']).size().unstack(fill_value=0).to_markdown()+'\n\n')
 f.write('## Yorum sınırları\n- Bu analiz fibroblast alt kümelerini **keşifsel** olarak tanımlar; bağımsız atlas eşleştirmesi ve biyolojik fonksiyon doğrulaması henüz yapılmadı.\n- Senesans/SASP skorları yalnızca aday gen panelleridir; hücresel senesans kanıtı değildir. Özellikle CDKN2A gibi düşük ekspresyonlu genlerin sıfırları biyolojik negatiflik olarak yorumlanamaz.\n- Alt küme belirteçleri log-normalize ortalama farkı ve pozitiflik farkına göre sıralandı; çoklu test düzeltmeli diferansiyel ekspresyon iddiası yoktur.\n- PDT hücre sayısı düşüktür. Örnekler biyolojik tekrar kabul edilmeli, hücreler bağımsız tekrar olarak test edilmemelidir.\n- Önceki Harmony kümesi 4 seçimine bağlıdır ve yeni alt kümeleme ham sayım normalizasyonundan yapılmıştır; örnek kaynaklı batch etkisi kalabilir.\n')
print('FINISHED',k,meta.groupby(['grup','fibro_alt_kume']).size().to_dict(),flush=True)
