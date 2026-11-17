import os,glob,time,itertools,zipfile
import numpy as np,pandas as pd
from scipy.stats import wilcoxon,spearmanr
from statsmodels.stats.multitest import multipletests
OUT='/mnt/data/GSE171213_Pseudobulk_Referans_v4';os.makedirs(OUT,exist_ok=True)
meta=pd.read_csv('/mnt/data/GSE171213_Fibro_AltKume_v1/Fibroblast_1343_AltKume_Etiketleri.csv')
samples=['HC1','HC2','HC3','HC4','PD1','PD2','PD3','PD4','PD5','PDT1','PDT2','PDT3']
common=set.intersection(*[set(pd.read_csv('/mnt/data/gse171213_sparse_cache/'+s+'_genes.csv',header=None).iloc[:,0].astype(str)) for s in samples]);genes=sorted(common);idx={g:i for i,g in enumerate(genes)}
print('Common',len(genes),flush=True)
# columns: whole fibroblast sample and each subcluster; sums of original UMI counts (not averaged single-cell expression)
records=[];matrix=[]; audit=[]
for s in samples:
 t=time.time(); sub=meta[meta.ornek==s];cells=sub.CellIndex.astype(str).tolist();p=glob.glob('/mnt/data/files/*_'+s+'.counts.tsv.gz')[0]
 # select gene column and only target fibroblast cell columns
 df=pd.read_csv(p,sep='\t',compression='gzip',usecols=['gene']+cells)
 assert df.gene.is_unique and set(cells).issubset(df.columns)
 df=df.set_index('gene').reindex(genes).fillna(0)
 A=df[cells].to_numpy(dtype=np.int32)
 assert (A>=0).all()
 for cl in ['all',0,1]:
  sel=np.ones(len(cells),dtype=bool) if cl=='all' else sub.fibro_alt_kume.to_numpy()==cl
  x=A[:,sel].sum(axis=1,dtype=np.int64)
  matrix.append(x)
  records.append({'ornek':s,'grup':'PDT' if s.startswith('PDT') else 'PD' if s.startswith('PD') else 'HC','alt_kume':str(cl),'n_hucre':int(sel.sum()),'library_UMI':int(x.sum())})
 audit.append({'ornek':s,'hucre':len(cells),'gen':len(df),'sn':round(time.time()-t,1)})
 print(s,len(cells),round(time.time()-t,1),flush=True)
counts=np.stack(matrix).T
m=pd.DataFrame(records);pd.DataFrame(audit).to_csv(OUT+'/Ham_Girdi_Denetimi.csv',index=False)
m.to_csv(OUT+'/Pseudobulk_Kutuphane_Ozeti.csv',index=False)
pd.DataFrame(counts,index=genes,columns=[f'{r.ornek}_{r.alt_kume}' for r in m.itertuples()]).to_csv(OUT+'/Pseudobulk_Ham_Sayimlar.csv.gz',compression='gzip',index_label='gene')
# external literature-grounded targeted panels: NOT official full SenMayo or CellAge downloads
panels={
 'arrest_p21_p16':['CDKN1A','CDKN2A','CDKN2B','GADD45A','GADD45B','BTG2'],
 'SASP_secretory':['IL6','CXCL8','CXCL1','CXCL2','CCL2','CCL7','MMP1','MMP3','MMP10','SERPINE1','IGFBP7','GDF15'],
 'stress_immediate_early':['FOS','JUN','JUNB','ATF3','DUSP1','HSPA1A','HSPA1B','DDIT3','HMOX1'],
 'inflammation_NFkB':['NFKBIA','TNFAIP3','ICAM1','CXCL8','IL6','CCL2'],
 'proliferation':['MKI67','TOP2A','PCNA','MCM2','MCM5','TYMS','UBE2C'],
 'ECM_fibroblast':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA']
}
# log CPM of pseudobulk profiles; summarize per sample and subcluster, and per-sample gene ranks
lib=counts.sum(axis=0).astype(float)
logcpm=np.log2(1+counts/lib[None,:]*1e6)
# Within sample use normalized panel gene expression and matched differences; exploratory targeted score
p_rows=[];coverage=[]
for name,gs in panels.items():
 found=[g for g in gs if g in idx];coverage.append({'panel':name,'planlanan':len(gs),'bulunan':len(found),'eksik':';'.join(sorted(set(gs)-set(found)))})
 z=logcpm[[idx[g] for g in found],:].mean(axis=0)
 for j,r in m.iterrows():p_rows.append({**r.to_dict(),'panel':name,'mean_log2_1p_CPM':float(z[j])})
pd.DataFrame(coverage).to_csv(OUT+'/Referans_Panel_Kapsami.csv',index=False)
score=pd.DataFrame(p_rows);score.to_csv(OUT+'/Pseudobulk_Panel_Skorlari.csv',index=False)
# per-gene pairwise sample-level within-group and pooled comparisons; require >=10 cells per subcluster
pairs=[]
for s in samples:
 i0=m.index[(m.ornek==s)&(m.alt_kume=='0')][0];i1=m.index[(m.ornek==s)&(m.alt_kume=='1')][0]
 if m.loc[i0,'n_hucre']>=10 and m.loc[i1,'n_hucre']>=10:
  pairs.append((s,i0,i1))
print('Qualified paired samples',len(pairs),[x[0] for x in pairs],flush=True)
if pairs:
 delta=np.stack([logcpm[:,i1]-logcpm[:,i0] for s,i0,i1 in pairs],axis=0)
 # exact two-sided sign-flip test on mean paired delta, gene-wise
 perms=np.array(list(itertools.product([-1.,1.],repeat=len(pairs))),dtype=np.float32)
 stat=np.abs(delta.mean(axis=0));null=np.abs(perms@delta/len(pairs))
 p=(null>=stat[None,:]-1e-10).mean(axis=0)
 q=multipletests(p,method='fdr_bh')[1]
 dd=pd.DataFrame({'gene':genes,'mean_delta_log2_1p_CPM':delta.mean(axis=0),'paired_p':p,'paired_q':q,'n_pairs':len(pairs)})
 dd.to_csv(OUT+'/AltKume_Pseudobulk_Gen_Testleri.csv',index=False)
# exact permutation test of sample-level pseudobulk group contrasts on all fibroblast cells
allix=m.index[m.alt_kume=='all'].to_numpy();a=m.loc[allix].reset_index(drop=True)
L=logcpm[:,allix];raw=counts[:,allix]
valid=((raw/lib[allix][None,:]*1e6)>1).sum(axis=1)>=3
rows=[];summary=[]
for ga,gb in [('HC','PD'),('HC','PDT'),('PD','PDT')]:
 ia=np.where(a.grup.to_numpy()==ga)[0];ib=np.where(a.grup.to_numpy()==gb)[0]
 z=L[:,np.r_[ia,ib]];n1=len(ia);n2=len(ib)
 obs=z[:,n1:].mean(axis=1)-z[:,:n1].mean(axis=1)
 perm=np.array(list(itertools.combinations(range(n1+n2),n1)))
 mask=np.zeros((len(perm),n1+n2),dtype=np.float64)
 mask[np.arange(len(perm))[:,None],perm]=1
 weights=(1-mask)/n2-mask/n1
 null=np.abs(weights@z.T)
 p=np.maximum((null>=np.abs(obs)[None,:]-1e-7).mean(axis=0),1/len(perm))
 p[~valid]=np.nan
 q=np.full(len(genes),np.nan);q[valid]=multipletests(p[valid],method='fdr_bh')[1]
 rows.extend({'karsilastirma':ga+'_vs_'+gb,'gene':genes[j],'log2_1p_CPM_fark_'+gb+'_eksi_'+ga:float(obs[j]),'permutation_p':float(p[j]) if valid[j] else np.nan,'BH_q':float(q[j]) if valid[j] else np.nan,'filtre_gecer':bool(valid[j]),'n_'+ga:n1,'n_'+gb:n2} for j in range(len(genes)))
 summary.append({'karsilastirma':ga+'_vs_'+gb,'n_'+ga:n1,'n_'+gb:n2,'perm_sayisi':len(perm),'gen_test':int(valid.sum()),'q_0_05':int(np.nansum(q<.05)),'min_p':float(np.nanmin(p))})
pd.DataFrame(rows).to_csv(OUT+'/Grup_Pseudobulk_Gen_Testleri.csv.gz',index=False,compression='gzip')
pd.DataFrame(summary).to_csv(OUT+'/Grup_Test_Ozeti.csv',index=False)
# exact sample-level permutation of targeted scores (all fibro sample) and matched subtype score
pr=[]
for panel in panels:
 sub=score[(score.panel==panel)&(score.alt_kume=='all')]
 for ga,gb in [('HC','PD'),('HC','PDT'),('PD','PDT')]:
  v1=sub.loc[sub.grup==ga,'mean_log2_1p_CPM'].to_numpy();v2=sub.loc[sub.grup==gb,'mean_log2_1p_CPM'].to_numpy()
  v=np.r_[v1,v2];n1=len(v1);diff=v2.mean()-v1.mean()
  null=[]
  for ids in itertools.combinations(range(len(v)),n1):
   mask=np.zeros(len(v),bool);mask[list(ids)]=True;null.append(v[~mask].mean()-v[mask].mean())
  p=np.mean(np.abs(null)>=abs(diff)-1e-10)
  pr.append({'panel':panel,'karsilastirma':ga+'_vs_'+gb,'fark_log2_1p_CPM':diff,'p_exact':p,'n1':len(v1),'n2':len(v2)})
 # qualified paired
 sub=score[(score.panel==panel)&(score.alt_kume.isin(['0','1']))]
 diffs=[]
 for s,i0,i1 in pairs:
  vv=sub[sub.ornek==s].set_index('alt_kume')['mean_log2_1p_CPM'];diffs.append(vv['1']-vv['0'])
 if diffs:
  d=np.array(diffs);null=np.abs(np.array(list(itertools.product([-1.,1.],repeat=len(d))))@d/len(d))
  p=np.mean(null>=abs(d.mean())-1e-10)
  pr.append({'panel':panel,'karsilastirma':'alt1_vs_alt0_paired_min10','fark_log2_1p_CPM':d.mean(),'p_exact':p,'n1':len(d),'n2':len(d)})
paneltests=pd.DataFrame(pr);paneltests['BH_q']=multipletests(paneltests.p_exact,method='fdr_bh')[1];paneltests.to_csv(OUT+'/Panel_Ornek_Duzeyi_Testleri.csv',index=False)
with open(OUT+'/Turkce_Bilimsel_Rapor.md','w') as f:
 f.write('# GSE171213 — bağımsız imza kontrolü ve pseudobulk v4\n\n')
 f.write('## Analiz kapsamı\n12 örnek, 1.343 fibroblast, 16.790 ortak gen. Orijinal UMI sayımları her örnek için ve fibroblast alt küme için toplandı. İstatistiksel tekrar **örnek** düzeyindedir.\n\n')
 f.write('## Referans sınırlılığı\nSenMayo ve CellAge yayımlanmış bağımsız referanslardır; bu çalıştırmada bunların **tam, orijinal gen listesi indirilemedi**. Bu nedenle kullanılan paneller literatürle uyumlu, **önceden tanımlanmış aday gen panelleridir**, resmi SenMayo/CellAge puanları ya da bağımsız atlas doğrulaması olarak sunulamaz.\n\n')
 f.write('## Yöntem\n- 12 örnekte ortak bulunan genler üzerinden örnek ve alt küme bazında **ham UMI toplama** yapıldı.\n- log2(1+CPM) ile keşifsel örnek düzeyi karşılaştırmalar yapıldı. Bu işlem edgeR/DESeq2 negatif binom GLM yerine geçmez.\n- Gruplar arasında tam permütasyon testi; alt kümeler arasında her alt kümede en az 10 hücresi olan örneklerde eşleştirilmiş tam işaret değiştirme testi kullanıldı.\n- Gen testleri BH-FDR düzeltmesiyle raporlandı. Panel testlerinde BH düzeltmesi uygulandı.\n- HC 4, PD 5, PDT 3 biyolojik örnek; PDT toplam 28 fibroblast içerir ve örnek başına sayım düşük olduğundan PDT analizleri özellikle kırılgandır.\n\n')
 f.write('## Grup karşılaştırmaları\n'+pd.DataFrame(summary).to_markdown(index=False)+'\n\n')
 f.write('## Panel karşılaştırmaları\n'+paneltests.to_markdown(index=False)+'\n\n')
 f.write('## Yorum\n- Örnek temelli sonuçlar hücreleri bağımsız tekrar saymaz.\n- Küçük örneklem, alt kümelerde hücre dengesizliği ve imza özgüllüğü senesans tanısı koymayı engeller.\n- SASP, inflamasyon ve stres imzaları örtüşebilir; hücre döngüsü baskılanması bağımsız olarak doğrulanmadıkça senesan fibroblast olarak kesin etiketleme yapılmaz.\n')
print('DONE',flush=True)
