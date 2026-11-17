import pandas as pd, numpy as np, matplotlib.pyplot as plt, json
from docx import Document
from docx.shared import Inches
from pathlib import Path
p=Path('/mnt/data'); q=pd.read_csv(p/'GSE171213_Ham_QC_Hucre_Metrikleri_v2.csv'); a=pd.read_csv(p/'GSE171213_Barkod_Uzlastirma_Hucre_Esleme_v1.csv'); q['barcode']=q.CellName.astype(str); a['barcode']=a.barcode.astype(str); a=a.drop(columns=['QC_MT15','QC_MT20','QC_MT25','QC_MT30','_merge'],errors='ignore')
assert len(q)==68542 and len(a)==43194
assert not q.duplicated(['ornek','barcode']).any() and not a.duplicated(['ornek','barcode']).any()
m=a.merge(q[['ornek','barcode','n_genes_raw','n_UMI_raw','mito_pct_raw','QC_MT15','QC_MT20','QC_MT25','QC_MT30']],on=['ornek','barcode'],validate='one_to_one',indicator='merge_check')
assert (m['merge_check']=='both').all()
for c in ['fibro_candidate','doublet_flag','QC_MT15','QC_MT20','QC_MT25','QC_MT30']:
 m[c]=m[c].astype(str).str.lower().eq('true')
for c in ['fibro_score','immune_score','epithelial_score','endothelial_score','pericyte_score']:
 m[c]=pd.to_numeric(m[c],errors='coerce')
rows=[]
for s,g in q.groupby('ornek',sort=False):
 x=m[m.ornek==s]; z={'ornek':s,'grup':g.grup.iloc[0],'ham_hucre':len(g),'pilot_panel':len(x),'pilot_fibro_adayi':int(x.fibro_candidate.sum()),'pilot_doublet_isaretli':int(x.doublet_flag.sum())}
 for t in [15,20,25,30]:
  key=f'QC_MT{t}'; y=x[x[key]]; f=y[y.fibro_candidate]
  z.update({f'QC{t}_tum':int(g[key].astype(str).str.lower().eq('true').sum()),f'QC{t}_panel':len(y),f'QC{t}_fibro':len(f),f'QC{t}_fibro_doublet_flag':int(f.doublet_flag.sum())})
 f=x[x.fibro_candidate]; z['fibro_skor_sifir']=int((f.fibro_score<=0).sum()); z['fibro_skor_pozitif']=int((f.fibro_score>0).sum()); z['fibro_immunden_yuksek']=int((f.fibro_score>f.immune_score).sum());rows.append(z)
out=pd.DataFrame(rows);out.to_csv(p/'GSE171213_QC_Karar_Ornek_Duyarlilik_v3.csv',index=False)
summary=[]
for t in [15,20,25,30]:
 key=f'QC_MT{t}'; y=m[m[key]]; f=y[y.fibro_candidate]; summary.append({'MT_esik':t,'ham_QC_hucre':int(q[key].astype(str).str.lower().eq('true').sum()),'panelde_QC_hucre':len(y),'fibro_adayi':len(f),'fibro_doublet_flag':int(f.doublet_flag.sum()),'fibro_skor_sifir':int((f.fibro_score<=0).sum()),'fibro_skor_pozitif':int((f.fibro_score>0).sum())})
sumdf=pd.DataFrame(summary);sumdf.to_csv(p/'GSE171213_QC_Karar_Toplam_Duyarlilik_v3.csv',index=False)
m[['ornek','grup','barcode','QC_MT15','QC_MT20','QC_MT25','QC_MT30','fibro_candidate','doublet_flag','fibro_score','immune_score','epithelial_score','endothelial_score','pericyte_score','n_genes_raw','n_UMI_raw','mito_pct_raw']].to_csv(p/'GSE171213_QC_Karar_Hucre_Denetimi_v3.csv',index=False)
fig,ax=plt.subplots(figsize=(8,4.8)); ax.plot(sumdf.MT_esik,sumdf.ham_QC_hucre,'o-',label='Tüm QC hücreleri');ax.plot(sumdf.MT_esik,sumdf.panelde_QC_hucre,'s-',label='Pilot paneldeki QC hücreleri');ax.plot(sumdf.MT_esik,sumdf.fibro_adayi,'^-',label='Pilot fibroblast adayları');ax.set(xlabel='Mitokondri eşiği (%)',ylabel='Hücre sayısı',title='QC eşik duyarlılığı ve pilot fibroblast adayları');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(p/'GSE171213_QC_Karar_Duyarlilik_v3.png',dpi=300);plt.close(fig)
D=Document();D.add_heading('GSE171213 | QC karar ve pilot doublet denetimi',0);D.add_paragraph('Sürüm v3 · 15 Kasım 2026 · Keşifsel karar denetimi; nihai QC veya doublet çağrısı değildir.')
D.add_heading('Kaynak ve doğrulama',1);D.add_paragraph('Girdi: ham matrislerden daha önce doğrulanmış 68.542 hücrelik QC CSV ve 43.194 hücrelik barkodları uzlaştırılmış pilot panel CSV. Eşleştirme örnek + tam barkod ile bire-bir yapıldı; yinelenen anahtar ve eşleşmeyen hücre kontrolü sıfır hata verdi. Ham sayım matrisleri bu sürümde tekrar işlenmedi. Önceki standart panel CSV dosyasında PD2 dahil bazı örneklerde barkod alanları eksikti; bu nedenle uzlaştırılmış barkod tablosu kullanıldı.')
D.add_heading('Eşik duyarlılığı',1);tab=D.add_table(rows=1,cols=5);tab.style='Light Shading Accent 1';
for c,v in zip(tab.rows[0].cells,['MT (%)','QC toplam','Panel QC','Fibro adayı','Fibro doublet işaretli']):c.text=v
for _,r in sumdf.iterrows():
 cells=tab.add_row().cells
 for c,v in zip(cells,[int(r.MT_esik),int(r.ham_QC_hucre),int(r.panelde_QC_hucre),int(r.fibro_adayi),int(r.fibro_doublet_flag)]):c.text=str(v)
D.add_picture(str(p/'GSE171213_QC_Karar_Duyarlilik_v3.png'),width=Inches(5.8))
D.add_heading('Anotasyon ve doublet metodolojik uyarıları',1);D.add_paragraph(f"Pilot fibroblast adayı toplamı {int(m.fibro_candidate.sum()):,}; bunların {int(((m.fibro_candidate)&(m.fibro_score<=0)).sum()):,} tanesinde fibro_score sıfır veya daha küçük. Bu bulgu, küme temelli aday etiketlerinin tek hücre belirteç düzeyinde doğrulanmadığını gösterir; fibroblast hücre sayıları biyolojik olarak kesin kabul edilemez.")
D.add_paragraph('Pilot doublet_flag, doğrulanmış Scrublet/scDblFinder/DoubletFinder sonucu değildir. Bu bayraklarla hücre silinmedi. Gerçek doublet analizi için ham sayım matrisi, örnek bazında beklenen doublet oranı, simüle çiftler ve bağımsız belirteç denetimi gereklidir.')
D.add_paragraph('Mitokondri eşiği %15–30 arasında duyarlılık analiziyle gösterildi. Nihai eşik, hücre tipi kompozisyonu ve örnek kalitesi incelenmeden seçilmedi. Hücreleri bağımsız biyolojik tekrar olarak kabul etmeyin.')
D.add_heading('Açık işler ve karar',1)
for s in ['Ham matris tabanlı gerçek doublet skorlaması: TAMAMLANMADI.','QC eşiği ve örnek bazlı aykırı örnek kararları: SABİTLENMEDİ.','Ortak hücre tipi anotasyonu ve fibroblast doğrulaması: TAMAMLANMADI.','Pseudobulk ve tedavi etkisi istatistikleri: YAPILMADI.'] : D.add_paragraph(s,style='List Bullet')
D.add_paragraph('Karar: QC duyarlılık ve pilot etiket tutarlılık denetimi tamamlandı; biyolojik QC/doublet aşaması kapatılamaz.')
D.save(p/'GSE171213_QC_Karar_Denetim_Raporu_v3.docx');print(sumdf.to_string(index=False));print('fibro zero',int(((m.fibro_candidate)&(m.fibro_score<=0)).sum()))
