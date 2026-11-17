import zipfile,io,numpy as np,pandas as pd
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from scipy.stats import mannwhitneyu
P='/mnt/data/'; z=zipfile.ZipFile(P+'GSE171213_Harmony_Sonuc_Paketi.zip'); q=zipfile.ZipFile(P+'GSE171213_Harmony_Colab_Girdi.zip')
a=pd.read_csv(z.open('GSE171213_Harmony_Hucre_Etiketleri.csv')); old=pd.read_csv(q.open('GSE171213_Harmony_Colab_Metadata.csv')); X=np.load(io.BytesIO(q.read('GSE171213_Harmony_Colab_Embedding.npy'))); Y=np.load(io.BytesIO(z.read('GSE171213_Harmony_Duzeltilmis_Latent.npy')))
assert len(a)==len(old)==len(X)==len(Y)==43438
assert (a['ornek'].values==old['ornek'].values).all()
grp=a.ornek.str.extract(r'^(HC|PDT|PD)')[0].to_numpy(); samp=a.ornek.to_numpy(); fib=a.fibro_high_confidence_candidate.astype(str).str.lower().isin(['true','1']).to_numpy(); cl=a.harmony_kume24.to_numpy()
rng=np.random.default_rng(42); ix=np.sort(rng.choice(len(a),min(8000,len(a)),replace=False)); out=[]
for name,E in [('once',X),('harmony_sonra',Y)]:
 E=E.astype('float32'); nn=NearestNeighbors(n_neighbors=31,algorithm='auto',n_jobs=2).fit(E); idx=nn.kneighbors(E[ix],return_distance=False)[:,1:]
 out.append(dict(durum=name,ornek_silhouette=silhouette_score(E[ix],samp[ix],sample_size=4000,random_state=42),grup_silhouette=silhouette_score(E[ix],grp[ix],sample_size=4000,random_state=42),komsu_ayni_ornek=float((samp[idx]==samp[ix,None]).mean()),komsu_ayni_grup=float((grp[idx]==grp[ix,None]).mean()),fibro_silhouette=silhouette_score(E[ix],fib[ix],sample_size=4000,random_state=42)))
pd.DataFrame(out).to_csv(P+'GSE171213_Biyolojik_Korunum_Metrikleri_v8.csv',index=False)
# preexisting expression marker scores, not Harmony-modified
b=pd.read_csv(P+'GSE171213_43438_Ortak_Transkriptom_Kumeleri_v3.csv'); assert (b.ornek.values==a.ornek.values).all() and (b.CellIndex.astype(str).values==a.CellIndex.astype(str).values).all()
markers=['Fibroblast_log_score','Fibroblast_genes_detected','Immune_log_score','Epithelial_log_score','Endothelial_log_score','Pericyte_log_score','Plasma_log_score','Mast_log_score']; tab=[]
for k in sorted(np.unique(cl)):
 mask=cl==k; d={'harmony_kume24':int(k),'hucre':int(mask.sum()),'fibro_aday':int((mask&fib).sum()),'fibro_aday_orani':float(fib[mask].mean()),'ornek_sayisi':int(a.ornek[mask].nunique()),'grup_sayisi':int(pd.Series(grp[mask]).nunique())}
 for m in markers:d[m+'_ort']=float(b.loc[mask,m].mean())
 tab.append(d)
pd.DataFrame(tab).to_csv(P+'GSE171213_Harmony_Kume_Marker_Korunum_v8.csv',index=False)
# sample-level fibro and cluster4 representation, group comparisons descriptive only
c=a.assign(grup=grp,fibro=fib,kume4=(cl==4)).groupby(['ornek','grup']).agg(toplam=('fibro','size'),fibro_aday=('fibro','sum'),kume4=('kume4','sum')).reset_index();c['fibro_yuzde']=100*c.fibro_aday/c.toplam;c['kume4_yuzde']=100*c.kume4/c.toplam
c.to_csv(P+'GSE171213_Harmony_Biyolojik_Ornek_Ozeti_v8.csv',index=False)
print(pd.DataFrame(out).to_string(index=False));print(c.to_string(index=False));print('cluster4',pd.DataFrame(tab).query('harmony_kume24==4').to_string(index=False))
