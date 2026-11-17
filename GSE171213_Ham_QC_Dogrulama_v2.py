"""GSE171213 raw counts QC, version 2. Python 3; pandas numpy matplotlib. No pilot CSV used in calculation."""
from pathlib import Path
import gzip,hashlib,json
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('/mnt/data'); IN=ROOT/'files'
SAMPLES=['HC1','HC2','HC3','HC4','PD1','PD2','PD3','PD4','PD5','PDT1','PDT2','PDT3']
records=[]; summaries=[]
for k,s in enumerate(SAMPLES):
 stem=f'GSM{5220920+k}_{s}'
 counts=IN/f'{stem}.counts.tsv.gz'; barcode=IN/f'{stem}.cellname.list.txt.gz'
 if not counts.exists() or not barcode.exists():raise FileNotFoundError(stem)
 b=pd.read_csv(barcode,sep='\t',dtype=str)
 if not {'CellIndex','CellName'}.issubset(b.columns):raise ValueError(f'{s}: barcode columns {b.columns.tolist()}')
 with gzip.open(counts,'rt') as f: head=f.readline().rstrip('\n\r').split('\t')
 idx=head[1:]
 if len(idx)!=len(b) or idx!=b.CellIndex.tolist():raise ValueError(f'{s}: column barcode order mismatch')
 if len(idx)!=len(set(idx)) or b.CellName.duplicated().any():raise ValueError(f'{s}: duplicate barcode')
 umi=np.zeros(len(idx),dtype=np.int64); ng=np.zeros(len(idx),dtype=np.int32); mt=np.zeros(len(idx),dtype=np.int64)
 genes=0;mtgenes=0
 for chunk in pd.read_csv(counts,sep='\t',chunksize=200,compression='gzip',dtype={ 'gene':str }):
  if chunk.columns[0]!='gene':raise ValueError(f'{s}: invalid gene col')
  names=chunk['gene'].astype(str)
  arr=chunk.iloc[:,1:].to_numpy(dtype=np.int64)
  if (arr<0).any():raise ValueError('negative count')
  umi+=arr.sum(axis=0);ng+=(arr>0).sum(axis=0).astype(np.int32)
  is_mt=names.str.upper().str.startswith('MT-').to_numpy()
  if is_mt.any():mt+=arr[is_mt].sum(axis=0);mtgenes+=int(is_mt.sum())
  genes+=len(chunk)
 pct=np.divide(mt*100,umi,out=np.zeros(len(idx),dtype=float),where=umi>0)
 group='HC' if s.startswith('HC') else 'PDT' if s.startswith('PDT') else 'PD'
 z=pd.DataFrame({'ornek':s,'grup':group,'CellIndex':idx,'CellName':b.CellName,'n_genes_raw':ng,'n_UMI_raw':umi,'mito_pct_raw':np.round(pct,8)})
 for cut in (15,20,25,30):z[f'QC_MT{cut}']=((ng>=300)&(ng<=6000)&(umi>0)&(pct<cut))
 records.append(z)
 rec={'ornek':s,'grup':group,'ham_hucre':len(idx),'gen_sayisi':genes,'mt_gen_sayisi':mtgenes,'medyan_mito_pct':float(np.median(pct)),'medyan_gen':float(np.median(ng)),'medyan_UMI':float(np.median(umi)),'sifir_UMI':int((umi==0).sum()),'counts_sha256':hashlib.sha256(counts.read_bytes()).hexdigest(),'barcode_sha256':hashlib.sha256(barcode.read_bytes()).hexdigest()}
 for cut in (15,20,25,30):rec[f'MT{cut}_hucre']=int(z[f'QC_MT{cut}'].sum())
 summaries.append(rec)
 print(s,len(idx),genes,round(np.median(pct),2),rec['MT30_hucre'],flush=True)
allcells=pd.concat(records,ignore_index=True); summary=pd.DataFrame(summaries)
allcells.to_csv(ROOT/'GSE171213_Ham_QC_Hucre_Metrikleri_v2.csv',index=False,encoding='utf-8-sig')
summary.to_csv(ROOT/'GSE171213_Ham_QC_Ornek_Ozeti_v2.csv',index=False,encoding='utf-8-sig')
summary.groupby('grup',sort=False)[['ham_hucre','MT15_hucre','MT20_hucre','MT25_hucre','MT30_hucre']].sum().to_csv(ROOT/'GSE171213_Ham_QC_Grup_Ozeti_v2.csv',encoding='utf-8-sig')
fig,ax=plt.subplots(figsize=(11,5));x=np.arange(len(summary));w=.19
for i,c in enumerate((15,20,25,30)):ax.bar(x+(i-1.5)*w,100*summary[f'MT{c}_hucre']/summary.ham_hucre,w,label=f'MT <%{c}')
ax.set_xticks(x,summary.ornek);ax.set_ylabel('Ham hücrelerin korunan yüzdesi');ax.set_title('Ham sayım matrislerinden yeniden hesaplanan QC eşik duyarlılığı');ax.legend(ncol=4);fig.tight_layout();fig.savefig(ROOT/'GSE171213_Ham_QC_Esik_Duyarliligi_v2.png',dpi=300);plt.close(fig)
with open(ROOT/'GSE171213_Ham_QC_Calisma_Ozeti_v2.json','w') as f:json.dump({'ham_hucre':int(summary.ham_hucre.sum()),'QC_15':int(summary.MT15_hucre.sum()),'QC_20':int(summary.MT20_hucre.sum()),'QC_25':int(summary.MT25_hucre.sum()),'QC_30':int(summary.MT30_hucre.sum()),'n_samples':len(summary),'mt_definition':'genes with symbol prefix MT- (case insensitive)','n_genes_definition':'number of genes with counts > 0 per cell','qc':'300 <= n_genes <= 6000 and UMI > 0 and mito_pct < threshold'},f,ensure_ascii=False,indent=2)
print('TOTAL',summary[['ham_hucre','MT15_hucre','MT20_hucre','MT25_hucre','MT30_hucre']].sum().to_dict())
