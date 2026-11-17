import pandas as pd,numpy as np,itertools,os
from scipy.stats import spearmanr
root='/mnt/data'; out=root+'/GSE171213_Referans_OnKontrol_v5';os.makedirs(out,exist_ok=True)
d=pd.read_csv(root+'/GSE171213_Pseudobulk_Referans_v4/Pseudobulk_Ham_Sayimlar.csv.gz').set_index('gene')
libs=d.sum(axis=0); x=np.log2(1+d.div(libs,axis=1)*1e6)
# Curated diagnostic controls, NOT original SenMayo or CellAge gene lists
panels={
 'senescence_arrest':['CDKN1A','CDKN2A','CDKN2B','GADD45A','GADD45B','BTG2','TP53','RB1'],
 'SASP_secretory':['IL6','CXCL8','CCL2','CXCL1','CXCL2','MMP1','MMP3','SERPINE1','IGFBP7','GDF15','IL1A','IL1B'],
 'stress_response':['JUN','FOS','ATF3','DDIT3','HSPA1A','HSPA1B','HSP90AA1','DNAJB1','NFKBIA','IER3'],
 'proliferation':['MKI67','TOP2A','PCNA','MCM2','MCM5','CDK1','CCNB1','UBE2C']}
coverage=[]; scores={}
for k,genes in panels.items():
 found=[g for g in genes if g in x.index]; missing=[g for g in genes if g not in x.index]
 coverage.append({'panel':k,'requested':len(genes),'available':len(found),'missing':';'.join(missing),'available_genes':';'.join(found)})
 scores[k]=x.loc[found].mean(axis=0)
sc=pd.DataFrame(scores);sc.index.name='sample_cluster';sc.to_csv(out+'/Pseudobulk_Kontrol_Panel_Skorlari.csv')
pd.DataFrame(coverage).to_csv(out+'/Panel_Gen_Kapsami.csv',index=False)
lib=pd.read_csv(root+'/GSE171213_Pseudobulk_Referans_v4/Pseudobulk_Kutuphane_Ozeti.csv');print('LIB',lib.columns.tolist());print(lib.head(3).to_string(index=False))
# Paired only when both subclusters have >=10 cells
counts={}
for _,r in lib.iterrows():
 keys=r.to_dict();
 for col in ['ornek','alt_kume','n_hucre']:
  if col not in keys:break
 else: counts[f"{r['ornek']}_{r['alt_kume']}"]=int(r['n_hucre'])
valid=[s for s in ['HC1','HC2','HC3','HC4','PD1','PD2','PD3','PD4','PD5','PDT1','PDT2','PDT3'] if counts.get(s+'_0',0)>=10 and counts.get(s+'_1',0)>=10]
print('valid',valid)
rows=[];deltas={}
for panel in panels:
 v=np.array([sc.loc[s+'_1',panel]-sc.loc[s+'_0',panel] for s in valid]);deltas[panel]=v
 obs=abs(v.mean()); perms=np.array([abs(np.mean(v*np.array(signs))) for signs in itertools.product([-1,1],repeat=len(v))]);p=(perms>=obs-1e-12).mean()
 rows.append({'panel':panel,'n_samples':len(valid),'mean_delta_log2_1p_CPM':v.mean(),'median_delta':np.median(v),'n_positive':sum(v>0),'p_exact_two_sided':p})
result=pd.DataFrame(rows);
from statsmodels.stats.multitest import multipletests
result['BH_q']=multipletests(result.p_exact_two_sided,method='fdr_bh')[1];result.to_csv(out+'/Eslesmis_Panel_Duyarlilik.csv',index=False)
pd.DataFrame({'sample':valid,**deltas}).to_csv(out+'/Ornek_Bazli_Panel_Farklari.csv',index=False)
cor=[]
for a,b in itertools.combinations(panels,2):
 rho,p=spearmanr(deltas[a],deltas[b]);cor.append({'panel_a':a,'panel_b':b,'rho':rho,'p':p})
pd.DataFrame(cor).to_csv(out+'/Panel_Fark_Korelasyonlari.csv',index=False)
print(result.to_string(index=False));print(pd.DataFrame(cor).to_string(index=False))
