import pandas as pd,numpy as np,zipfile,io,os,itertools
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
base='/mnt/data/';out=base+'GSE171213_Orijinal_Referans_Dogrulama_v6';os.makedirs(out,exist_ok=True)
with zipfile.ZipFile(base+'cellAge.zip') as z: ca=pd.read_csv(io.BytesIO(z.read('cellage3.tsv')),sep='\t')
with zipfile.ZipFile(base+'cellSignatures.zip') as z: cs=pd.read_csv(io.BytesIO(z.read('signatures1.csv')),sep=';')
from openpyxl import load_workbook
w=load_workbook(base+'41467_2022_32552_MOESM4_ESM.xlsx',read_only=True,data_only=True)
sm=[str(row[0]).strip().upper() for row in list(w['human'].values)[1:] if row[0]]
d=pd.read_csv(base+'GSE171213_Pseudobulk_Referans_v4/Pseudobulk_Ham_Sayimlar.csv.gz').set_index('gene');d=d.groupby(level=0).sum();lib=pd.read_csv(base+'GSE171213_Pseudobulk_Referans_v4/Pseudobulk_Kutuphane_Ozeti.csv');good=lib[(lib.alt_kume.astype(str).isin(['0','1'])) & (lib.n_hucre>=10)].groupby('ornek').alt_kume.nunique();samples=good[good==2].index.tolist();print('samples',samples)
# CPM normalization
cpm=np.log2(1+ d.div(d.sum(axis=0),axis=1)*1e6)
ca['Gene symbol']=ca['Gene symbol'].astype(str).str.upper();cs['gene_symbol']=cs['gene_symbol'].astype(str).str.upper();
sets={'SenMayo_125':(set(sm),set()),'CellAge_induces':(set(ca.loc[ca['Senescence Effect']=='Induces','Gene symbol']),set()),'CellAge_inhibits':(set(ca.loc[ca['Senescence Effect']=='Inhibits','Gene symbol']),set()),'CellAge_directional':(set(ca.loc[ca['Senescence Effect']=='Induces','Gene symbol']),set(ca.loc[ca['Senescence Effect']=='Inhibits','Gene symbol'])),'Signatures_overexpressed':(set(cs.loc[(cs.ovevrexp>cs.underexp)&(cs.p_value<.05),'gene_symbol']),set()),'Signatures_underexpressed':(set(cs.loc[(cs.underexp>cs.ovevrexp)&(cs.p_value<.05),'gene_symbol']),set()),'Signatures_directional':(set(cs.loc[(cs.ovevrexp>cs.underexp)&(cs.p_value<.05),'gene_symbol']),set(cs.loc[(cs.underexp>cs.ovevrexp)&(cs.p_value<.05),'gene_symbol']))}
results=[];sample_rows=[];coverage=[]
for name,(pos,neg) in sets.items():
 pos=pos&set(d.index);neg=neg&set(d.index); allgenes=sets[name][0]|sets[name][1]
 coverage.append(dict(panel=name,reference_genes=len(allgenes),measured_genes=len(allgenes&set(d.index)),coverage_pct=100*len(allgenes&set(d.index))/len(allgenes),positive_measured=len(pos),negative_measured=len(neg)))
 diffs=[]
 for s in samples:
  v=(cpm.loc[list(pos),s+'_1'].mean()-cpm.loc[list(pos),s+'_0'].mean()) if pos else 0
  if neg:v-=cpm.loc[list(neg),s+'_1'].mean()-cpm.loc[list(neg),s+'_0'].mean()
  diffs.append(float(v));sample_rows.append(dict(sample=s,panel=name,difference=float(v)))
 vals=np.array(diffs);obs=abs(vals.mean());perms=np.array([abs(np.mean(vals*np.array(sign))) for sign in itertools.product([-1,1],repeat=len(vals))]);p=(perms>=obs-1e-10).mean()
 results.append(dict(panel=name,n_samples=len(samples),mean_difference=vals.mean(),median_difference=np.median(vals),positive_samples=int((vals>0).sum()),permutation_p=p))
res=pd.DataFrame(results);res['BH_q']=multipletests(res.permutation_p,method='fdr_bh')[1];res.to_csv(out+'/Referans_Test_Sonuclari.csv',index=False);pd.DataFrame(coverage).to_csv(out+'/Referans_Kapsami.csv',index=False);pd.DataFrame(sample_rows).to_csv(out+'/Ornek_Bazli_Referans_Farklari.csv',index=False)
print(res.to_string(index=False));print(pd.DataFrame(coverage).to_string(index=False))
# Evaluate correlations with stress from previous control
prev=pd.read_csv(base+'GSE171213_Referans_OnKontrol_v5/Ornek_Bazli_Panel_Farklari.csv').set_index('sample');wide=pd.DataFrame(sample_rows).pivot(index='sample',columns='panel',values='difference');rows=[]
for panel in wide.columns:
 for other in ['stress_response','SASP_secretory']:
  common=wide.index.intersection(prev.index);rho,p=spearmanr(wide.loc[common,panel],prev.loc[common,other]);rows.append(dict(reference=panel,comparator=other,spearman_rho=rho,p=p,n=len(common)))
pd.DataFrame(rows).to_csv(out+'/Stres_SASP_Korelasyonlari.csv',index=False)
with open(out+'/Turkce_Rapor.md','w') as f:
 f.write('# Özgün senesans referansları ile doğrulama\n\nKaynaklar: SenMayo Nature Communications 2022 ek veri insan 125 gen; HAGR CellAge cellage3.tsv; HAGR cellSignatures signatures1.csv.\n\nYöntem: Ham pseudobulk UMI matrisi, log2(1+CPM), iki alt kümede en az 10 hücresi olan 7 örnek, eşleştirilmiş iki yönlü eksiksiz işaret değiştirme testi (128 permütasyon), 7 panel için BH-FDR. CellAge yönlü paneli senesansı indükleyen genler eksi baskılayan genler; gen başına eşit ağırlık. SenMayo SASP ağırlıklı olup bağımsız yaşlanma kanıtı sayılmaz.\n\n'+res.to_markdown(index=False)+'\n\n## Gen kapsamı\n\n'+pd.DataFrame(coverage).to_markdown(index=False)+'\n\nYorum: Referans imzaları arasında gen örtüşmesi vardır ve testler bağımsız değildir. Pseudobulk panel skoru hücre düzeyinde yaşlanmayı doğrudan göstermez. Stres ve SASP korelasyonları ayrışma hakkında keşifsel göstergedir, nedensellik veya bağımsız etki testi değildir.\n')
with zipfile.ZipFile(base+'GSE171213_Orijinal_Referans_Dogrulama_v6.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
 for fn in os.listdir(out):z.write(out+'/'+fn,fn)
