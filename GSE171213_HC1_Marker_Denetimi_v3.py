import gzip,csv,numpy as np,pandas as pd
from pathlib import Path
P=Path('/mnt/data'); root=P/'files'; base='GSM5220920_HC1'
sets={'fibroblast':['COL1A1','COL1A2','DCN','LUM','PDGFRA','COL3A1'],'immune':['PTPRC','LST1','TYROBP','CD3D','MS4A1'],'epithelial':['EPCAM','KRT14','KRT5','KRT19'],'endothelial':['PECAM1','VWF','KDR','CLDN5'],'pericyte':['RGS5','ACTA2','MCAM','PDGFRB']}
allgenes=set(sum(sets.values(),[])); selected={}
with gzip.open(root/(base+'.counts.tsv.gz'),'rt') as f:
 r=csv.reader(f,delimiter='\t'); header=next(r); idx={x:i for i,x in enumerate(header[1:])}
 for row in r:
  if row[0] in allgenes:selected[row[0]]=np.array(row[1:],dtype=np.float32)
qc=pd.read_csv(P/'GSE171213_HC1_Doublet_Skorlari_Pilot_v1.csv'); n=len(qc)
res=qc.copy(); res['pilot_annotation_present']=res.CellIndex.isin(pd.read_csv(P/'HC1_Hucre_Panel_Skorlari_v1.csv').cell_index)
for typ,genes in sets.items():
 found=[g for g in genes if g in selected]; mat=np.vstack([selected[g][[idx[x] for x in res.CellIndex]] for g in found]) if found else np.zeros((0,n))
 res[typ+'_detected_markers']=(mat>0).sum(axis=0) if len(mat) else 0
 res[typ+'_total_counts']=mat.sum(axis=0) if len(mat) else 0
 res[typ+'_evidence']=res[typ+'_detected_markers']>=2
res['lineages_with_2markers']=sum(res[k+'_evidence'].astype(int) for k in sets)
res['mixed_lineage_candidate']=res.lineages_with_2markers>=2
res['high_doublet_score_060']=res.simulated_neighbor_fraction>=0.6
res['high_score_and_mixed']=res.high_doublet_score_060 & res.mixed_lineage_candidate
res['high_score_and_fibro_evidence']=res.high_doublet_score_060 & res.fibroblast_evidence
res.to_csv(P/'GSE171213_HC1_Marker_Denetimi_Hucre_v3.csv',index=False)
summary={'QC30_hucre':n,'pilot_anotasyonu_olmayan':int((~res.pilot_annotation_present).sum()),'ham_panelde_bulunan_gen':len(selected),'panel_gen_toplam':len(allgenes),'iki_ve_uzeri_fibro_belirteci':int(res.fibroblast_evidence.sum()),'iki_ve_uzeri_soy_belirteci_karisimi':int(res.mixed_lineage_candidate.sum()),'skor_ge_060':int(res.high_doublet_score_060.sum()),'skor_ge_060_ve_karisik':int(res.high_score_and_mixed.sum()),'skor_ge_060_ve_fibro':int(res.high_score_and_fibro_evidence.sum()),'cikarilan_hucre':0}
pd.DataFrame([summary]).to_csv(P/'GSE171213_HC1_Marker_Denetimi_Ozet_v3.csv',index=False)
from docx import Document
d=Document();d.add_heading('GSE171213 / HC1 — Ham belirteç ve doublet denetimi',0);d.add_paragraph('HC1 QC30 kapsamındaki tüm hücreler için ham counts.tsv.gz dosyasından hücre soyuna özgü belirteçler yeniden okundu. Pilot anotasyon tablosuna bağımlı kalınmadı.')
d.add_heading('Doğrulanan sayılar',1)
for k,v in summary.items():d.add_paragraph(f'{k}: {v}',style='List Bullet')
d.add_heading('Yöntem ve karar',1);d.add_paragraph('Fibroblast: COL1A1, COL1A2, DCN, LUM, PDGFRA, COL3A1; immün: PTPRC, LST1, TYROBP, CD3D, MS4A1; epitel: EPCAM, KRT14, KRT5, KRT19; endotel: PECAM1, VWF, KDR, CLDN5; perisit: RGS5, ACTA2, MCAM, PDGFRB. Her soy için en az iki belirtecin pozitif olması tarama kanıtı sayıldı; bu bir nihai hücre tipi etiketi değildir.')
d.add_paragraph('Yapay-doublet komşuluk skoru >=0,60 ile iki farklı soydan en az ikişer belirteç saptanması birlikte raporlandı. Skor kalibrasyonu ve beklenen doublet oranı belirlenmediğinden hiçbir hücre silinmedi. Homotipik doubletler yalnızca soy karışımıyla yakalanamaz. HC1 kesin doublet sınıflandırması açık kalmaktadır.')
d.add_paragraph('Pilot anotasyonu olmayan hücrelerin tamamı aynı ham belirteç kontrolünden geçti. Marker yokluğu, düşük yakalama derinliği nedeniyle hücre soyunu dışlamaz.')
d.save(P/'GSE171213_HC1_Marker_Denetimi_Rapor_v3.docx')
print(summary)
