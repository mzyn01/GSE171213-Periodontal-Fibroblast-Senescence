import os
os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['OMP_NUM_THREADS']='2'
import numpy as np,pandas as pd
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
from scipy.stats import chi2_contingency
from pathlib import Path
P=Path('/mnt/data'); a=pd.read_csv(P/'GSE171213_43438_Ortak_Transkriptom_Kumeleri_v3.csv',low_memory=False)
z=np.load(P/'GSE171213_Kumeleme_Latent_v3.npz')['embedding'].astype('float32'); assert len(a)==len(z)==43438
assert (a.transkriptom_kume.to_numpy()==np.load(P/'GSE171213_Kumeleme_Latent_v3.npz')['cluster']).all()
# Existing SVD embedding is grouped by sample in source checkpoint; no row reorder.
# Sample centering within HC/PD/PDT preserves each group centroid, but may remove real sample biology.
groups=a.ornek.str.extract(r'^(HC|PD|PDT)')[0].to_numpy(); samples=a.ornek.to_numpy(); original=a.transkriptom_kume.to_numpy()
centered=z.copy()
for g in np.unique(groups):
    gi=np.where(groups==g)[0]; gm=z[gi].mean(0)
    for s in np.unique(samples[gi]):
        ix=np.where(samples==s)[0]; sm=z[ix].mean(0)
        centered[ix]-=(sm-gm)

def cramers_v(lab):
    t=pd.crosstab(samples,lab);chi2=chi2_contingency(t,correction=False)[0]
    return float(np.sqrt(chi2/(len(lab)*min(t.shape[0]-1,t.shape[1]-1))))
def concentration(lab):
    t=pd.crosstab(lab,samples); return float((t.max(axis=1)/t.sum(axis=1)).median()),int(((t.max(axis=1)/t.sum(axis=1))>.5).sum())
rows=[]; labels={};
for strength in [0.,.5,1.]:
    x=z+strength*(centered-z)
    k=MiniBatchKMeans(n_clusters=24,random_state=42,n_init=8,batch_size=4096,max_iter=220).fit_predict(x)
    labels[strength]=k
    median,nmajor=concentration(k)
    rows.append(dict(strength=strength,method='within-group sample-centering sensitivity; NOT Harmony/scVI',n_cells=len(k),n_clusters=len(set(k)),cramers_v=cramers_v(k),median_cluster_top_sample_share=median,clusters_majority_one_sample=nmajor,ARI_vs_original=adjusted_rand_score(original,k),NMI_vs_original=normalized_mutual_info_score(original,k)))
    a[f'cluster_centering_{int(strength*100)}']=k
pd.DataFrame(rows).to_csv(P/'GSE171213_Batch_Duyarlilik_Metrikleri_v7.csv',index=False)
# Gene markers scores computed independently of the corrected embedding, but from same source raw counts.
neg=['Immune','Epithelial','Endothelial','Pericyte','Plasma','Mast']
marker=a['Fibroblast_log_score'].to_numpy(); marker2=a['Fibroblast_genes_detected'].to_numpy()>=2
maxneg=a[[f'{x}_log_score' for x in neg]].max(axis=1).to_numpy()
strict=marker2&(marker>maxneg)&(marker>=.15)
stat=[]
for strength,k in labels.items():
    for c in np.unique(k):
        ix=k==c
        stat.append(dict(strength=strength,cluster=int(c),cells=int(ix.sum()),strict_marker_cells=int(strict[ix].sum()),strict_marker_fraction=float(strict[ix].mean()),mean_fibro_score=float(marker[ix].mean()),mean_negative_max=float(maxneg[ix].mean()),top_sample_share=float(pd.Series(samples[ix]).value_counts(normalize=True).iloc[0])))
s=pd.DataFrame(stat);s.to_csv(P/'GSE171213_Batch_Duyarlilik_Kume_Marker_v7.csv',index=False)
a[['ornek','CellIndex','CellName','transkriptom_kume','cluster_centering_0','cluster_centering_50','cluster_centering_100','fibro_high_confidence_candidate']].to_csv(P/'GSE171213_Batch_Duyarlilik_Hucre_v7.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
print('strict markers',int(strict.sum()),'original candidate',int(a.fibro_high_confidence_candidate.sum()))
