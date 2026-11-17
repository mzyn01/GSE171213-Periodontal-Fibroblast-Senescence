from pathlib import Path
import pandas as pd, numpy as np
from scipy.stats import spearmanr, mannwhitneyu
from docx import Document
from docx.shared import Pt
P=Path('/mnt/data')
a=pd.read_csv(P/'GSE171213_HC1_Doublet_Skorlari_Pilot_v1.csv'); b=pd.read_csv(P/'GSE171213_QC_Karar_Hucre_Denetimi_v3.csv'); b=b[(b.ornek=='HC1') & (b.QC_MT30==True)].copy()
assert a.CellName.is_unique and b.barcode.is_unique
x=a.merge(b,left_on='CellName',right_on='barcode',how='left',validate='one_to_one',indicator=True)
assert len(x)==2794
matched=int((x._merge=='both').sum()); missing=len(x)-matched
for v in ['fibro_candidate','doublet_flag']:
 x[v]=x[v].fillna(False).astype(bool)
s=x.simulated_neighbor_fraction
for cut in [.5,.6,.7,.8]:x[f'skor_ge_{str(cut).replace(".","")}']=s>=cut
markers=['fibro_score','immune_score','epithelial_score','endothelial_score','pericyte_score']
x['positive_marker_groups']=(x[markers].fillna(0)>0).sum(axis=1)
x['mixed_marker_groups_2plus']=x.positive_marker_groups>=2
for col in ['n_UMI_raw_x','n_genes_raw_x']:
 assert col in x.columns
rows=[]
for cut in [.5,.6,.7,.8]:
 high=x[s>=cut]; low=x[s<cut]
 rows.append(dict(threshold=cut,high_n=len(high),high_pct=100*len(high)/len(x),high_median_UMI=high.n_UMI_raw_x.median(),low_median_UMI=low.n_UMI_raw_x.median(),high_median_genes=high.n_genes_raw_x.median(),low_median_genes=low.n_genes_raw_x.median(),high_mixed_2plus_pct=100*high.loc[high._merge=='both','mixed_marker_groups_2plus'].mean(),low_mixed_2plus_pct=100*low.loc[low._merge=='both','mixed_marker_groups_2plus'].mean(),high_pilot_doublet_flag_n=int(high.doublet_flag.sum()),high_fibro_candidate_n=int(high.fibro_candidate.sum())))
t=pd.DataFrame(rows)
cor=[]
for v in ['n_UMI_raw_x','n_genes_raw_x','mito_pct_raw_x','positive_marker_groups']:
 sub=x[x._merge=='both'] if v=='positive_marker_groups' else x
 r,p=spearmanr(sub.simulated_neighbor_fraction,sub[v]);cor.append(dict(variable=v,spearman_rho=r,p_value=p))
c=pd.DataFrame(cor)
cols=['CellIndex','CellName','n_genes_raw_x','n_UMI_raw_x','mito_pct_raw_x','simulated_neighbor_fraction','rank_percentile','fibro_candidate','doublet_flag']+markers+['positive_marker_groups','mixed_marker_groups_2plus']
x[cols+['_merge']].to_csv(P/'GSE171213_HC1_Doublet_Degerlendirme_Hucre_v2.csv',index=False)
t.to_csv(P/'GSE171213_HC1_Doublet_Esik_Duyarlilik_v2.csv',index=False)
c.to_csv(P/'GSE171213_HC1_Doublet_Korelasyon_v2.csv',index=False)
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axs=plt.subplots(1,2,figsize=(11,4))
axs[0].hist(s,bins=np.arange(-.0125,1.025,.025));axs[0].set(xlabel='Yapay doublet komşuluk oranı',ylabel='Hücre sayısı',title='HC1 skor dağılımı')
axs[1].scatter(np.log10(x.n_UMI_raw_x),s,s=7,alpha=.35);axs[1].set(xlabel='log10(UMI)',ylabel='Yapay doublet komşuluk oranı',title='UMI ile ilişki')
fig.tight_layout();fig.savefig(P/'GSE171213_HC1_Doublet_Degerlendirme_v2.png',dpi=300);plt.close(fig)
d=Document();sec=d.sections[0];sec.header.paragraphs[0].text='GSE171213  |  HC1 DOUBLEt DENETİMİ'
d.add_heading('HC1 doublet skoru: değerlendirme ve karar',0)
d.add_paragraph('15 Kasım 2026 • Örnek: HC1 • QC30 hücreleri: 2.794')
d.add_heading('Amaç ve kapsam',level=1)
d.add_paragraph('Önceki yapay doublet komşuluk skorlarını, QC metriklerini ve pilot hücre tipi panel skorlarını barkod bazında birleştirerek teknik ve biyolojik tutarlılık açısından denetlemek. Bu çalışma yeni bir bağımsız doublet algoritması veya deneysel doğrulama değildir.')
d.add_heading('Veri bütünlüğü',level=1)
d.add_paragraph(f'HC1 skor tablosundaki {len(a):,} hücrenin {matched:,} tanesi pilot anotasyon denetim tablosuyla barkod üzerinden eşleşti; {missing:,} QC30 hücresi pilot panelde yok. Çift eşleşme: 0.'.replace(',','.'))
d.add_heading('Skor eşiklerinin duyarlılığı',level=1)
tab=d.add_table(rows=1,cols=5);tab.style='Light Shading Accent 1'
for cell,h in zip(tab.rows[0].cells,['Eşik','Hücre','%','Yüksek skor medyan UMI','Yüksek skor medyan gen']):cell.text=h
for _,r in t.iterrows():
 for cell,val in zip(tab.add_row().cells,[f'{r.threshold:.2f}',str(int(r.high_n)),f'{r.high_pct:.1f}',f'{r.high_median_UMI:.0f}',f'{r.high_median_genes:.0f}']):cell.text=val
d.add_heading('İlişki denetimi',level=1)
for _,r in c.iterrows():d.add_paragraph(f"{r.variable}: Spearman ρ={r.spearman_rho:.3f}; p={r.p_value:.3g}",style='List Bullet')
d.add_heading('Yorum ve sınırlamalar',level=1)
d.add_paragraph('Skor, gerçek hücrelerin yakın komşuları arasında yapay çift hücre oranıdır; kalibre edilmiş doublet olasılığı değildir. Skor için seçilen 0,50–0,80 eşikleri yalnızca duyarlılık amaçlıdır, önceden doğrulanmış karar sınırları değildir. Pilot anotasyon panelinde bulunmayan 797 hücrede marker karışımı değerlendirilemez; marker temelli oranlar yalnızca eşleşen hücrelerde hesaplanmıştır. Pilot hücre tipi skorlarının sıfır veya pozitif olması güvenilir hücre tipi anotasyonu kanıtlamaz; iki pozitif panel gerçek karışık transkriptomu da yansıtabilir.')
d.add_paragraph('Önceki pilot yöntem 2.794 yapay doublet oluşturmuş, log-normalizasyon, 2.000 yüksek varyanslı gen, 30 bileşenli TruncatedSVD ve 40 komşu kullanmıştır. Bu uygulamada yapay/gerçek komşu sayısı için teorik kalibrasyon, gerçek doublet prevalansı, hücre tipi çiftleri, bağımsız tekrar/tohum stabilitesi veya deneysel etiket bulunmamaktadır.')
d.add_heading('Karar',level=1)
d.add_paragraph('HC1 için skorların barkod ve QC anotasyonuyla birleştirilmesi, dağılım ve eşik duyarlılık denetimi tamamlandı. Ancak doublet sınıflandırması doğrulanmadı. Şimdilik 0 hücre çıkarılacak; HC1 doublet filtresi nihai olarak kapatılmayacak. Bir sonraki zorunlu iş, bağımsız/kalibre edilmiş yöntemle tekrar ve biyolojik belirteç çapraz kontrolüdür.')
d.add_heading('Tekrar üretilebilirlik',level=1)
d.add_paragraph('Girdiler: GSE171213_HC1_Doublet_Skorlari_Pilot_v1.csv; GSE171213_QC_Karar_Hucre_Denetimi_v3.csv. Çıktılar: hücre denetimi, eşik duyarlılığı, korelasyon CSV, 300 dpi şekil ve bu rapor. Kod: GSE171213_HC1_Doublet_Degerlendirme_v2.py.')
style=d.styles['Normal'];style.font.name='Aptos';style.font.size=Pt(9)
d.save(P/'GSE171213_HC1_Doublet_Degerlendirme_Raporu_v2.docx')
print(t.to_string(index=False));print(c.to_string(index=False));print('matched',matched,'missing',missing)
