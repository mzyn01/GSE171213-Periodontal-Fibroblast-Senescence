from pathlib import Path
import pandas as pd,numpy as np,json
from docx import Document
from docx.shared import Cm
R=Path('/mnt/data')
raw=pd.read_csv(R/'GSE171213_Ham_QC_Hucre_Metrikleri_v2.csv');old=pd.read_csv(R/'GSE171213_Ortak_QC_Hucre_Metrikleri_v1.csv')
raw['hucre_sirasi']=raw.groupby('ornek').cumcount()+1
m=raw.merge(old,on=['ornek','hucre_sirasi'],validate='one_to_one',suffixes=('_ham','_pilot'))
assert len(m)==len(raw)==len(old)==68542
assert (m.grup_ham==m.grup_pilot).all()
for a,b in [('n_genes_raw','gen'),('n_UMI_raw','UMI')]: assert (m[a]==m[b]).all(),a
max_delta=float((m.mito_pct_raw-m.mito_yuzde).abs().max());assert max_delta<0.001,max_delta
for cut in [15,20,25,30]:assert (m[f'QC_MT{cut}_ham']==m[f'QC_MT{cut}_pilot']).all(),cut
rec=pd.read_csv(R/'GSE171213_Barkod_Uzlastirma_Hucre_Esleme_v1.csv')
print('reconcile fields',rec[['ornek','cell','barcode']].head().to_dict('records'))
rawkey=raw[['ornek','CellIndex','CellName']].rename(columns={'CellIndex':'cell','CellName':'barcode_raw'})
matched=rec.merge(rawkey,on=['ornek','cell'],how='left',validate='many_to_one',indicator='raw_merge')
assert len(matched)==len(rec)
barcode_match=(matched.barcode.astype(str)==matched.barcode_raw.astype(str))
res=pd.DataFrame({'ornek':matched.ornek,'cell':matched.cell,'barcode_match':barcode_match,'raw_found':matched.raw_merge.eq('both')})
res.to_csv(R/'GSE171213_Ham_QC_Barkod_Dogrulama_v2.csv',index=False)
summary=raw.groupby(['grup','ornek'],sort=False).agg(ham=('CellIndex','size'),mt30=('QC_MT30','sum')).reset_index()
summary['barcode_match']=summary.ornek.map(res.groupby('ornek').barcode_match.sum()).astype(int)
summary['pilot_count']=summary.ornek.map(res.groupby('ornek').size()).astype(int)
summary['pilot_barcode_unmatched']=summary.pilot_count-summary.barcode_match
summary.to_csv(R/'GSE171213_Ham_QC_Karsilastirma_Ozeti_v2.csv',index=False)
result={'raw_cells':len(raw),'pilot_cells':len(old),'barcode_panel_cells':len(rec),'barcode_matched':int(barcode_match.sum()),'barcode_unmatched':int((~barcode_match).sum()),'gene_mismatch':int((m.n_genes_raw!=m.gen).sum()),'UMI_mismatch':int((m.n_UMI_raw!=m.UMI).sum()),'mito_max_absolute_difference':max_delta,'qc_mismatches':{str(c):int((m[f'QC_MT{c}_ham']!=m[f'QC_MT{c}_pilot']).sum()) for c in [15,20,25,30]},'MT30':int(raw.QC_MT30.sum()),'pilot_QC30_overlap':int(rec.QC_MT30.sum()),'QC30_outside_panel':int(raw.QC_MT30.sum()-rec.QC_MT30.sum()),'panel_outside_QC30':int((~rec.QC_MT30).sum())}
(R/'GSE171213_Ham_QC_Dogrulama_Kontrol_v2.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
D=Document();D.add_heading('GSE171213 | Ham sayım matrislerinden QC doğrulaması',0);D.add_paragraph('12 örnek | 15 Kasım 2026 | Sürüm 2')
D.add_heading('Amaç ve girdi',1);D.add_paragraph('12 örneğin orijinal counts.tsv.gz ve cellname.list.txt.gz dosyaları kullanılarak gen sayısı, toplam UMI ve MT- ile başlayan genlerin toplam UMI içindeki oranı sıfırdan hesaplandı. Her ham girdinin SHA-256 özeti örnek tablosunda kayıtlıdır.')
D.add_heading('Yöntem',1);D.add_paragraph('Gen sayısı: pozitif sayımı bulunan gen sayısı; UMI: tüm genlerin sayım toplamı; MT%: MT- önekli genlerin toplamı / tüm UMI ×100. Gen 300–6000 (sınırlar dahil), UMI>0, MT <%15/20/25/30 (katı küçük). Ham matris sütunları barkod CellIndex sırasıyla birebir denetlendi. Karşılaştırma örnek ve hücre sırası üzerinden yapıldı; ayrıca pilot panel barkodları ham CellIndex ve CellName ile doğrulandı.')
D.add_heading('Bağımsız yeniden hesaplama sonuçları',1);D.add_paragraph(f"Toplam {result['raw_cells']:,} hücre; MT<%30 koşulunda {result['MT30']:,} hücre. Gen uyuşmazlığı {result['gene_mismatch']}; UMI uyuşmazlığı {result['UMI_mismatch']}; maksimum MT yüzdesi mutlak farkı {max_delta:.2g}.".replace(',','.'))
D.add_paragraph(f"Pilot paneldeki {len(rec):,} hücre kimliğinin ham barkodla tam eşleşeni {result['barcode_matched']:,}; eşleşmeyen {result['barcode_unmatched']:,}.".replace(',','.'))
D.add_paragraph('MT eşiklerine göre hamdan hesaplanan ve önceki QC dışa aktarımı arasındaki hücre düzeyinde uyumsuzluklar: '+str(result['qc_mismatches']))
D.add_heading('Örnek düzeyi denetim',1);t=D.add_table(rows=1,cols=5);t.style='Light Shading Accent 1'
for c,x in zip(t.rows[0].cells,['Örnek','Ham','MT<30','Pilot panel','Barkod eşleşen']):c.text=x
for _,r in summary.iterrows():
 for c,x in zip(t.add_row().cells,[r.ornek,r.ham,r.mt30,r.pilot_count,r.barcode_match]):c.text=str(x)
D.add_picture(str(R/'GSE171213_Ham_QC_Esik_Duyarliligi_v2.png'),width=Cm(15.5))
D.add_heading('Yorum ve kapsam sınırı',1);D.add_paragraph('Bu aşama ham QC metriklerinin ve barkod kimliklerinin doğrulanmasıdır. Eşik duyarlılığı gösterilmiştir, ancak biyolojik hücre tipine göre eşik uygunluğu, doublet tespiti, ortak anotasyon, batch etkisi ve donor düzeyinde karşılaştırmalar henüz tamamlanmadı. MT<%30 nihai tercih olarak sabitlenmemiştir.')
D.add_heading('Tekrarlanabilirlik',1);D.add_paragraph('Kodlar: GSE171213_Ham_QC_Dogrulama_v2.py ve GSE171213_Ham_QC_Karsilastirma_v2.py. Çıktılar: hücre ve örnek CSV, grup özeti CSV, barkod denetimi CSV, kontrol JSON, 300 dpi PNG ve bu rapor. Kaynak dosyaların SHA-256 değerleri GSE171213_Ham_QC_Ornek_Ozeti_v2.csv içindedir.')
D.save(R/'GSE171213_Ham_QC_Dogrulama_Raporu_v2.docx');print(json.dumps(result,indent=2));print(summary.to_string(index=False))
