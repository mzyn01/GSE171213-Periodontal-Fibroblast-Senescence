import os,gzip,glob
import numpy as np,pandas as pd
from scipy.stats import wilcoxon,kruskal
from statsmodels.stats.multitest import multipletests
out='/mnt/data/GSE171213_Senesans_Genisletilmis_v3'
meta=pd.read_csv('/mnt/data/GSE171213_Fibro_AltKume_v1/Fibroblast_1343_AltKume_Skorlari.csv')
panels={
 'cell_cycle_arrest':['CDKN1A','CDKN2A','CDKN2B','TP53','RB1','GADD45A','GADD45B','BTG2'],
 'senescence_support':['GLB1','GDF15','SERPINE1','IGFBP3','IGFBP5','IGFBP7','HMGB1','LMNB1'],
 'SASP_inflammatory':['IL6','IL1A','IL1B','CXCL1','CXCL2','CXCL3','CXCL8','CCL2','CCL7','CCL8','TNF','TGFB1','MMP1','MMP3','MMP9','MMP10','MMP12','CSF2'],
 'SASP_ECM':['SERPINE1','MMP2','MMP14','TIMP1','TIMP2','COL1A1','COL3A1','FN1','POSTN'],
 'proliferation':['MKI67','TOP2A','PCNA','MCM2','MCM5','MCM6','UBE2C','CENPF','CDK1'],
 'stress_inflammation':['NFKBIA','TNFAIP3','ICAM1','JUN','FOS','ATF3','DDIT3','HIF1A'],
 'fibroblast':['COL1A1','COL1A2','DCN','LUM','PDGFRA','COL3A1']}
allgenes=set(sum(panels.values(),[])); rows=[]; data={}; coverage=[]
for sample,g in meta.groupby('ornek',sort=False):
 f=glob.glob('/mnt/data/files/*_'+sample+'.counts.tsv.gz')[0]
 indices=g.index.to_numpy(); cellids=g.CellIndex.astype(str).tolist()
 with gzip.open(f,'rt') as fh:
  hdr=fh.readline().rstrip('\n').split('\t'); loc={v:i for i,v in enumerate(hdr[1:])}
  col=np.array([loc[c]+1 for c in cellids]); found=set()
  for line in fh:
   gene=line.split('\t',1)[0]
   if gene not in allgenes:continue
   vals=np.fromstring(line.split('\t',1)[1],sep='\t',dtype=np.float32)
   data[gene,sample]=vals[col-1];found.add(gene)
 coverage.append({'ornek':sample,'n_hucre':len(g),'mevcut_gen':len(found),'eksik_gen':';'.join(sorted(allgenes-found))})
 for gene in allgenes:
  if gene in found:
   for ix,v in zip(indices,data[gene,sample]):rows.append((ix,gene,v))
print('extracted',len(rows),flush=True)
long=pd.DataFrame(rows,columns=['ix','gene','count']); long=long.merge(meta[['ornek','grup','fibro_alt_kume']],left_on='ix',right_index=True)
# each gene per sample/subcluster detection, sample paired tests on per-gene detection
agg=long.groupby(['ornek','grup','fibro_alt_kume','gene'],as_index=False).agg(n=('count','size'),detect=('count',lambda x:np.mean(x>0)),mean_count=('count','mean'))
agg.to_csv(out+'/Ornek_AltKume_Gen_Deteksiyon.csv',index=False)
results=[]
for gene,g in agg.groupby('gene'):
 pv=g.pivot(index='ornek',columns='fibro_alt_kume',values='detect').dropna()
 if len(pv)<3:continue
 diff=(pv[1]-pv[0]).values
 try:p=wilcoxon(diff,zero_method='wilcox').pvalue if np.any(diff!=0) else 1.
 except ValueError:p=1.
 results.append({'gen':gene,'n_eslesmis_ornek':len(diff),'alt0_ortalama_deteksiyon':pv[0].mean(),'alt1_ortalama_deteksiyon':pv[1].mean(),'ortalama_fark':diff.mean(),'pozitif_ornek':int((diff>0).sum()),'p_wilcoxon':p})
res=pd.DataFrame(results).sort_values('p_wilcoxon');res['q_BH']=multipletests(res.p_wilcoxon,method='fdr_bh')[1]
res.to_csv(out+'/Gen_Bazli_Eslesmis_Dogrulama.csv',index=False)
# panel summary mean detection over genes available for each sample/subcluster, paired difference
panelrows=[]
for name,genes in panels.items():
 z=agg[agg.gene.isin(genes)].groupby(['ornek','grup','fibro_alt_kume'],as_index=False).agg(panel_detect=('detect','mean'),n_gene=('gene','nunique'))
 z['panel']=name;panelrows.append(z)
pan=pd.concat(panelrows,ignore_index=True);pan.to_csv(out+'/Panel_Ornek_AltKume.csv',index=False)
pr=[]
for name,g in pan.groupby('panel'):
 pv=g.pivot(index='ornek',columns='fibro_alt_kume',values='panel_detect').dropna()
 diff=pv[1]-pv[0]
 try:p=wilcoxon(diff).pvalue if np.any(diff!=0) else 1.
 except ValueError:p=1.
 pr.append({'panel':name,'n_ornek':len(diff),'alt0_ortalama':pv[0].mean(),'alt1_ortalama':pv[1].mean(),'fark':diff.mean(),'pozitif_ornek':int((diff>0).sum()),'p_wilcoxon':p})
pr=pd.DataFrame(pr);pr['q_BH']=multipletests(pr.p_wilcoxon,method='fdr_bh')[1];pr.to_csv(out+'/Panel_Eslesmis_Testler.csv',index=False)
pd.DataFrame(coverage).to_csv(out+'/Girdi_Gen_Kapsami.csv',index=False)
print(pr.to_string(index=False),flush=True)
print('selected genes',res[res.gen.isin(['CDKN1A','CDKN2A','LMNB1','MKI67','TOP2A','PCNA','IL6','CXCL8','SERPINE1','GDF15'])].to_string(index=False),flush=True)
