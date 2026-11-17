from pathlib import Path
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from docx import Document
from docx.shared import Cm,Pt
ROOT=Path('/mnt/data'); samples=['HC1','HC2','HC3','HC4','PD1','PD2','PD3','PD4','PD5','PDT1','PDT2','PDT3']
cols={'HC1':('n_genes','n_counts','pct_mito'),'HC2':('n_genes','n_counts','pct_mito'),'HC3':('n_genes','n_UMI','pct_mito'),'HC4':('n_gen','n_UMI','mito_yuzde'),'PD1':('n_genes','n_counts','pct_mito'),'PD2':('n_gen','n_UMI','mitokondri_yuzde'),'PD3':('nGene','nUMI','mito_pct'),'PD4':('n_genes','n_counts','pct_mito'),'PD5':('n_genes','n_UMI','mito_pct'),'PDT1':('n_genes','n_UMI','mito_pct'),'PDT2':('nGenes','nUMI','mito_pct'),'PDT3':('saptanan_gen','UMI','mitokondri_yuzde')}
rows=[]; cellframes=[]
for s in samples:
 d=pd.read_csv(ROOT/f'{s}_Barkod_QC_v1.csv'); a,b,c=cols[s]; z=pd.DataFrame({'ornek':s,'grup':('HC' if s.startswith('HC') else 'PDT' if s.startswith('PDT') else 'PD'),'gen':pd.to_numeric(d[a],errors='coerce'),'UMI':pd.to_numeric(d[b],errors='coerce'),'mito_yuzde':pd.to_numeric(d[c],errors='coerce')});z['hucre_sirasi']=np.arange(1,len(z)+1)
 for mt in [15,20,25,30]:z[f'QC_MT{mt}']=(z.gen.between(300,6000)&(z.UMI>0)&(z.mito_yuzde<mt))
 cellframes.append(z)
 rec={'ornek':s,'grup':z.grup.iloc[0],'ham_hucre':len(z),'eksik_metrik':int(z[['gen','UMI','mito_yuzde']].isna().any(axis=1).sum()),'medyan_gen':float(z.gen.median()),'medyan_UMI':float(z.UMI.median()),'medyan_mito_yuzde':float(z.mito_yuzde.median()),'mito_100_ustu':int((z.mito_yuzde>100).sum())}
 for mt in [15,20,25,30]:rec[f'QC_MT{mt}_hucre']=int(z[f'QC_MT{mt}'].sum());rec[f'QC_MT{mt}_oran']=round(float(z[f'QC_MT{mt}'].mean()*100),2)
 rows.append(rec)
summary=pd.DataFrame(rows); cells=pd.concat(cellframes,ignore_index=True)
summary.to_csv(ROOT/'GSE171213_Ortak_QC_Ornek_Ozeti_v1.csv',index=False,encoding='utf-8-sig');cells.to_csv(ROOT/'GSE171213_Ortak_QC_Hucre_Metrikleri_v1.csv',index=False,encoding='utf-8-sig')
group=summary.groupby('grup',sort=False).agg(ornek_sayisi=('ornek','count'),ham_hucre=('ham_hucre','sum'),MT15=('QC_MT15_hucre','sum'),MT20=('QC_MT20_hucre','sum'),MT25=('QC_MT25_hucre','sum'),MT30=('QC_MT30_hucre','sum')).reset_index();group.to_csv(ROOT/'GSE171213_Ortak_QC_Grup_Ozeti_v1.csv',index=False,encoding='utf-8-sig')
fig,ax=plt.subplots(figsize=(12,5));x=np.arange(len(summary));width=.19
for i,mt in enumerate([15,20,25,30]):ax.bar(x+(i-1.5)*width,summary[f'QC_MT{mt}_oran'],width,label=f'MT < %{mt}')
ax.set_xticks(x,summary.ornek);ax.set_ylabel('QC sonrası korunan hücre (%)');ax.set_ylim(0,105);ax.legend(ncol=4,loc='upper right');ax.set_title('GSE171213: Ortak gen aralığında (300–6000) mitokondri eşik duyarlılığı');fig.tight_layout();fig.savefig(ROOT/'GSE171213_Ortak_QC_Esik_Duyarliligi_v1.png',dpi=300);plt.close(fig)
fig,ax=plt.subplots(figsize=(11,4));ax.scatter(range(len(summary)),summary.medyan_mito_yuzde);ax.set_xticks(range(len(summary)),summary.ornek);ax.set_ylabel('Ham hücrelerde medyan MT (%)');ax.set_title('Örnekler arası mitokondri yükü');ax.axhline(15,linestyle='--',alpha=.5);ax.axhline(30,linestyle='--',alpha=.5);fig.tight_layout();fig.savefig(ROOT/'GSE171213_Ortak_QC_Mitokondri_v1.png',dpi=300);plt.close(fig)
doc=Document();sec=doc.sections[0];sec.header.paragraphs[0].text='GSE171213 | Tekrarlanabilir ortak QC – Aşama 1'
doc.add_heading('GSE171213 – Ortak kalite kontrol ve eşik duyarlılığı',0)
doc.add_paragraph('Türkçe yöntem ve ara sonuç raporu | Sürüm 1 | 15 Kasım 2026')
doc.add_heading('Kapsam ve amaç',1);doc.add_paragraph('HC1–HC4 (sağlıklı), PD1–PD5 (periodontitis), PDT1–PDT3 (tedavi sonrası) olmak üzere 12 örneğin daha önce dışa aktarılmış hücre düzeyi QC tabloları tek şemada birleştirildi. Bu adım yeni bir hücre anotasyonu, doublet doğrulaması veya diferansiyel gen ekspresyonu analizi değildir.')
doc.add_heading('Tekrarlanabilir yöntem',1);doc.add_paragraph('Girdi: her örneğin *_Barkod_QC_v1.csv dosyası. Her dosyadaki gen, UMI ve mitokondri yüzdesi sütunları örnek bazlı açık eşleme sözlüğüyle birleştirildi. Aynı gen sayısı aralığı (300–6000), UMI>0 ve dört ayrı MT eşiği (<%15, <%20, <%25, <%30) uygulandı. Yüzde, ham hücre sayısı paydasına göre hesaplandı. Script: GSE171213_ortak_qc_v1.py. Girdiler değiştirilmeden bırakıldı.')
doc.add_heading('Temel bulgular',1);doc.add_paragraph(f'Toplam ham hücre: {len(cells):,}. Eksik QC metriği bulunan hücre: {summary.eksik_metrik.sum():,}.'.replace(',','.'))
for mt in [15,20,25,30]:doc.add_paragraph(f'MT <%{mt}: {int(summary[f"QC_MT{mt}_hucre"].sum()):,} hücre ({100*summary[f"QC_MT{mt}_hucre"].sum()/len(cells):.1f}%).'.replace(',','.'))
t=doc.add_table(rows=1, cols=5);t.style='Light Shading Accent 1'
for cell,label in zip(t.rows[0].cells,['Örnek','Ham','MT<15','MT<25','MT<30']):cell.text=label
for _,r in summary.iterrows():
 cc=t.add_row().cells
 for c,v in zip(cc,[r.ornek,int(r.ham_hucre),int(r.QC_MT15_hucre),int(r.QC_MT25_hucre),int(r.QC_MT30_hucre)]):c.text=str(v)
doc.add_picture(str(ROOT/'GSE171213_Ortak_QC_Esik_Duyarliligi_v1.png'),width=Cm(16))
doc.add_heading('Yorum ve önceden tanımlanacak kararlar',1);doc.add_paragraph('Mitokondri eşiklerine duyarlılık örnekler arasında değişir. Eşiği yalnızca en çok hücreyi koruduğu için seçmek uygun değildir. Nihai QC için dağılımlar, hücre tipleri, platform özellikleri, olası hasarlı hücreler ve çift hücre olasılığı birlikte incelenmelidir. Bu raporda hiçbir eşik nihai olarak ilan edilmemiştir.')
doc.add_heading('Sınırlamalar ve sonraki adım',1);doc.add_paragraph('Bu analiz, önceki pilot analizlerin QC dışa aktarımlarına dayanır; ham sayım matrislerinden tüm metrikleri bağımsız olarak yeniden hesaplayıp doğrulamaz. Hücre/donör kimlik eşlemesi, biyolojik eşleşmiş tedavi örnekleri, örnekler arası entegrasyon, batch etkileri, güvenilir fibroblast anotasyonu, doublet algoritması ve donör düzeyinde pseudobulk henüz tamamlanmadı. Bir sonraki adım ham matrislerden ortak QC doğrulaması, çift hücre taraması ve hücre tipi anotasyonudur.')
doc.add_paragraph('Üretilen dosyalar: örnek özeti CSV, grup özeti CSV, hücre metrikleri CSV, iki 300 dpi grafik, bu Word raporu ve yeniden çalıştırılabilir Python kodu.')
doc.save(ROOT/'GSE171213_Ortak_QC_Asama1_Rapor_v1.docx')
print(summary[['ornek','ham_hucre','medyan_mito_yuzde','QC_MT15_hucre','QC_MT25_hucre','QC_MT30_hucre']].to_string(index=False));print('TOTAL',len(cells),[int(summary[f'QC_MT{m}_hucre'].sum()) for m in [15,20,25,30]])
