import gzip,glob,zipfile,os
import numpy as np,pandas as pd
from scipy.stats import mannwhitneyu
ROOT='/mnt/data'; OUT=ROOT+'/GSE171213_Fibroblast_Gen_Dogrulama_v1'
os.makedirs(OUT,exist_ok=True)
with zipfile.ZipFile(ROOT+'/GSE171213_Harmony_Sonuc_Paketi.zip') as z:
 df=pd.read_csv(z.open('GSE171213_Harmony_Hucre_Etiketleri.csv'))
qc=pd.read_csv(ROOT+'/GSE171213_43438_Korunan_Hucre_Anotasyon_v1.csv',usecols=['ornek','CellIndex','n_UMI_raw'])
df=df.merge(qc,on=['ornek','CellIndex'],how='left',validate='one_to_one')
assert len(df)==43438 and df.n_UMI_raw.notna().all()
panels={
 'Fibroblast':['COL1A1','COL1A2','DCN','LUM','PDGFRA','COL3A1','FBLN1','DPT','THY1','PDPN','MFAP5'],
 'Immune':['PTPRC','LST1','CD3D','CD79A','MS4A1','LYZ'],
 'Epithelial':['EPCAM','KRT14','KRT5','KRT19','KRT17'],
 'Endothelial':['PECAM1','VWF','KDR','CLDN5','EMCN'],
 'Pericyte_SMC':['RGS5','CSPG4','MCAM','ACTA2','MYH11'],
 'Osteoblast':['BGLAP','RUNX2','SP7'],
 'Proliferation':['MKI67','TOP2A']}
genes=list(dict.fromkeys(sum(panels.values(),[])))
values=np.zeros((len(df),len(genes)),dtype=np.float32)
found={g:0 for g in genes}; sample_coverage=[]
for sample in df.ornek.unique():
 path=glob.glob(ROOT+'/files/*_'+sample+'.counts.tsv.gz')
 assert len(path)==1,(sample,path)
 ids=np.flatnonzero(df.ornek.to_numpy()==sample)
 selected=df.iloc[ids].CellIndex.astype(str).tolist()
 with gzip.open(path[0],'rt') as f:
  header=next(f).rstrip('\n\r').split('\t')
  pos={c:i for i,c in enumerate(header[1:])}
  assert len(set(selected))==len(selected) and all(c in pos for c in selected)
  ix=np.array([pos[c] for c in selected],dtype=np.int32)
  present=set()
  for line in f:
   gene=line.split('\t',1)[0]
   if gene not in found:continue
   chunks=line.rstrip('\n\r').split('\t',1)[1]
   arr=np.fromstring(chunks,sep='\t',dtype=np.float32)
   assert len(arr)==len(header)-1,(sample,gene,len(arr),len(header))
   values[ids,genes.index(gene)]=arr[ix]
   found[gene]+=1;present.add(gene)
 sample_coverage.append({'ornek':sample,'hucre':len(ids),'marker_mevcut':len(present),'marker_toplam':len(genes)})
 print(sample,len(ids),len(present),flush=True)
# UMI-normalize log1p counts per 10k
norm=np.log1p(values/df.n_UMI_raw.to_numpy(dtype=np.float32)[:,None]*10000)
cluster=df.harmony_kume24.astype(str).eq('4').to_numpy()
assert cluster.sum()==1343
rows=[]
for j,g in enumerate(genes):
 a=norm[cluster,j];b=norm[~cluster,j]
 det_a=(values[cluster,j]>0).mean();det_b=(values[~cluster,j]>0).mean()
 rows.append({'gen':g,'panel':next(k for k,v in panels.items() if g in v),'ornek_kapsami':found[g],
  'fibro_kume4_pozitif_yuzde':round(det_a*100,2),'diger_pozitif_yuzde':round(det_b*100,2),
  'fibro_kume4_log1p_CPM_ortalama':round(float(a.mean()),4),'diger_log1p_CPM_ortalama':round(float(b.mean()),4),
  'pozitif_oran_farki_yuzde_puan':round((det_a-det_b)*100,2),
  'log_ortalama_farki':round(float(a.mean()-b.mean()),4)})
res=pd.DataFrame(rows)
res.to_csv(OUT+'/Gen_Belirtec_Karsilastirmasi.csv',index=False)
# per-sample panel score and contamination
srows=[]
for sample in df.ornek.unique():
 idx=(df.ornek.to_numpy()==sample)&cluster
 for panel,gg in panels.items():
  avail=[genes.index(g) for g in gg if found[g]>0]
  srows.append({'ornek':sample,'grup':sample.rstrip('12345'),'panel':panel,'fibro_hucre':int(idx.sum()),
   'panel_log1p_CPM_ortalama':round(float(norm[np.ix_(idx,avail)].mean()),4) if idx.sum() and avail else np.nan,
   'en_az_bir_marker_pozitif_yuzde':round(float((values[np.ix_(idx,avail)]>0).any(axis=1).mean()*100),2) if idx.sum() and avail else np.nan})
pd.DataFrame(srows).to_csv(OUT+'/Ornek_Bazli_Panel_Dogrulamasi.csv',index=False)
pd.DataFrame(sample_coverage).to_csv(OUT+'/Ham_Veri_Kapsam_Denetimi.csv',index=False)
# summary panel group comparisons
summary=[]
for panel,gg in panels.items():
 avail=[genes.index(g) for g in gg if found[g]>0]
 a=norm[np.ix_(cluster,avail)].mean(axis=1)
 b=norm[np.ix_(~cluster,avail)].mean(axis=1)
 summary.append({'panel':panel,'fibro_kume4_skor_ortalama':float(a.mean()),'diger_skor_ortalama':float(b.mean()),'fark':float(a.mean()-b.mean()),'gen_sayisi':len(avail)})
pd.DataFrame(summary).to_csv(OUT+'/Panel_Skor_Ozeti.csv',index=False)
# save only cluster 4 gene expression for downstream
subset=df.loc[cluster,['ornek','CellIndex','CellName','harmony_kume24','n_UMI_raw']].reset_index(drop=True)
for j,g in enumerate(genes):subset[g]=values[cluster,j]
subset.to_csv(OUT+'/Fibro_Kume4_Ham_Marker_Sayimlari.csv.gz',index=False,compression='gzip')
print('DONE',OUT,flush=True)
