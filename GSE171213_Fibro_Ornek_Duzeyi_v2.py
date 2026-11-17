import pandas as pd,numpy as np
from scipy.stats import mannwhitneyu, kruskal, wilcoxon, spearmanr
from itertools import combinations
from pathlib import Path
p=Path('/mnt/data/GSE171213_Fibro_AltKume_v1'); out=Path('/mnt/data/GSE171213_Fibro_Ornek_Duzeyi_v2');out.mkdir(exist_ok=True)
d=pd.read_csv(p/'Fibroblast_1343_AltKume_Skorlari.csv'); assert len(d)==1343
panels=['senesans_hucre_dongusu_skor','SASP_skor','inflamasyon_skor','proliferasyon_skor']
a=d.groupby(['grup','ornek','fibro_alt_kume']).agg(n=('CellName','size'),**{x:(x,'mean') for x in panels}).reset_index()
a.to_csv(out/'Ornek_AltKume_Panel_Ortalamalari.csv',index=False)
rows=[]
for x in panels:
 for g in ['HC','PD','PDT']:
  z=d[d.grup==g].groupby('ornek')[x].mean().values
  rows.append(dict(panel=x,grup=g,n_ornek=len(z),ortalama=float(z.mean()),medyan=float(np.median(z)),min=float(z.min()),max=float(z.max())))
summary=pd.DataFrame(rows);summary.to_csv(out/'Grup_Bazli_Tanimlayicilar.csv',index=False)
res=[]
for x in panels:
 samples={g:d[d.grup==g].groupby('ornek')[x].mean().values for g in ['HC','PD','PDT']}
 H,pv=kruskal(*samples.values())
 res.append(dict(panel=x,karsilastirma='HC-PD-PDT',test='Kruskal-Wallis',n1=4,n2=5,istatistik=H,p_degeri=pv))
 for g1,g2 in [('HC','PD'),('HC','PDT'),('PD','PDT')]:
  v1,v2=samples[g1],samples[g2]
  u,pv=mannwhitneyu(v1,v2,alternative='two-sided',method='exact')
  res.append(dict(panel=x,karsilastirma=f'{g1}-{g2}',test='Mann-Whitney U exact',n1=len(v1),n2=len(v2),istatistik=u,p_degeri=pv))
pd.DataFrame(res).to_csv(out/'Ornek_Duzeyi_Kesifsel_Testler.csv',index=False)
paired=[]
for x in panels:
 pivot=a.pivot(index='ornek',columns='fibro_alt_kume',values=x).dropna()
 diff=pivot[1]-pivot[0]
 # two-sided exact sign-flip permutation of mean differences
 n=len(diff); vals=diff.to_numpy(); obs=abs(vals.mean()); allmeans=np.array([abs(np.mean(vals*np.array([1 if (bits>>i)&1 else -1 for i in range(n)]))) for bits in range(1<<n)])
 pexact=(allmeans>=obs-1e-12).mean()
 paired.append(dict(panel=x,n_eslesmis_ornek=n,ortalama_fark_alt1_eksi_alt0=diff.mean(),medyan_fark=diff.median(),pozitif_ornek_sayisi=int((diff>0).sum()),iki_yonlu_isaret_degistirme_p=pexact))
pd.DataFrame(paired).to_csv(out/'AltKume_Eslesmis_Ornek_Testleri.csv',index=False)
count=pd.crosstab(d.ornek,d.fibro_alt_kume); count['grup']=[x.rstrip('12345') for x in count.index];count['alt1_orani']=count[1]/(count[0]+count[1]);count.to_csv(out/'AltKume_Ornek_Oranlari.csv')
with open(out/'Turkce_Ornek_Duzeyi_Raporu.md','w') as f:
 f.write('# GSE171213 – Fibroblast senesans/SASP örnek düzeyi keşifsel denetim\n\n')
 f.write('Birim: biyolojik örnek (HC n=4, PD n=5, PDT n=3); 1.343 hücre bağımsız tekrar kabul edilmedi. Analiz önceki ham sayım temelli hücre skorlarını kullanır.\n\n')
 f.write('## Gruplar arası testler\n\n'+pd.DataFrame(res).round(5).to_markdown(index=False)+'\n\n')
 f.write('## Alt küme 1–0: örnek içinde eşleştirilmiş farklar\n\n'+pd.DataFrame(paired).round(5).to_markdown(index=False)+'\n\n')
 f.write('## Kısıtlar\n\nPDT hücreleri yalnızca 28 ve bazı örneklerde 6–11 hücredir. PD4 yalnızca 8 hücredir. Panel skorları keşifseldir; özellikle SASP panelinde gen kapsamı eksiktir. Hücre tipi alt kümeleri silhouette açısından zayıf ayrılmıştır. Testler çoklu karşılaştırma düzeltmesi olmadan keşifseldir; nedensel veya kesin senesans sınıflandırması yapılmaz.\n')
print(pd.DataFrame(res).to_string(index=False));print('\nPAIRED\n',pd.DataFrame(paired).to_string(index=False))
