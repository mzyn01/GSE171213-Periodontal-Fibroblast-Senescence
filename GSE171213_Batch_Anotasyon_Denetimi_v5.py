import pandas as pd,numpy as np
from scipy.stats import chi2_contingency, fisher_exact
from pathlib import Path
P=Path('/mnt/data');d=pd.read_csv(P/'GSE171213_43438_Ortak_Transkriptom_Kumeleri_v3.csv'); assert len(d)==43438 and d[['ornek','CellName']].duplicated().sum()==0
ct=pd.crosstab(d.transkriptom_kume,d.ornek);chi,p,df,ex=chi2_contingency(ct);v=np.sqrt(chi/(len(d)*min(ct.shape[0]-1,ct.shape[1]-1)))
base=d.ornek.value_counts(normalize=True);rows=[]
for k,g in d.groupby('transkriptom_kume'):
 f=g.ornek.value_counts(normalize=True);top=f.idxmax(); expected=base[top]; H=-(f*np.log(f)).sum()/np.log(12)
 rows.append(dict(kume=k,n=len(g),en_baskin_ornek=top,en_baskin_ornek_orani=f.max(),beklenen_ornek_orani=expected,zenginlesme_kati=f.max()/expected,normalize_entropi=H,ornek_sayisi=g.ornek.nunique(),fibro_iki_belirtec_orani=(g.Fibroblast_genes_detected>=2).mean(),fibro_skor_ort=g.Fibroblast_log_score.mean(),karsi_soy_en_yuksek_ort=max(g[x].mean() for x in ['Immune_log_score','Epithelial_log_score','Endothelial_log_score','Pericyte_log_score','Plasma_log_score','Mast_log_score'])))
a=pd.DataFrame(rows).sort_values('kume');a.to_csv(P/'GSE171213_Batch_Kume_Denetimi_v5.csv',index=False)
# Cross-check cluster identity against separately computed original marker detection columns, NOT independent reference
q=d[d.transkriptom_kume==10].copy(); cols=['ornek','Fibroblast_genes_detected','fibroblast_detected_markers','fibroblast_candidate_strict','fibro_high_confidence_candidate','Fibroblast_log_score','Immune_log_score','Epithelial_log_score','Endothelial_log_score','Pericyte_log_score'];
summary=q.groupby('ornek').agg(n=('CellName','size'),fibro_ge2=('Fibroblast_genes_detected',lambda x:(x>=2).sum()),pilot_ge2=('fibroblast_detected_markers',lambda x:(x>=2).sum()),strict=('fibroblast_candidate_strict','sum'),high_conf=('fibro_high_confidence_candidate','sum'),fibro_mean=('Fibroblast_log_score','mean'),immune_mean=('Immune_log_score','mean'),epithelial_mean=('Epithelial_log_score','mean'),endothelial_mean=('Endothelial_log_score','mean'),pericyte_mean=('Pericyte_log_score','mean')).reset_index();summary.to_csv(P/'GSE171213_Fibroblast_Kume10_Capraz_Denetim_v5.csv',index=False)
# quantitative batch metrics, not corrected
with open(P/'GSE171213_Batch_Anotasyon_Degerlendirme_v5.md','w') as f:
 f.write('# GSE171213 – Batch etkisi ve fibroblast anotasyonu denetimi (v5)\n\n')
 f.write('## Analiz kapsamı\n43.438 hücre, 12 örnek, 24 keşifsel KMeans kümesi. Mevcut hücre bazlı küme etiketleri ve önceden hesaplanan belirteç alanları yeniden analiz edildi. Bu işlem batch düzeltmesi veya bağımsız referans anotasyonu değildir.\n\n')
 f.write(f'## Örnek-küme ilişkisi\nKi-kare={chi:.1f}; serbestlik derecesi={df}; p={p:.3g}; Cramér V={v:.3f}. Hücreler bağımsız biyolojik tekrar olmadığı için p değeri yalnızca tanımlayıcıdır; örnek bazlı çıkarım için kullanılamaz.\n\n')
 f.write('Örnek ağırlıkları eşit olmadığından yalnızca baskın örnek yüzdesi yerine, beklenen örnek payına göre zenginleşme ve normalize Shannon entropisi de hesaplandı.\n\n')
 f.write('### En fazla örnek yoğunlaşması gösteren kümeler\n\n'+a.sort_values('zenginlesme_kati',ascending=False).head(8)[['kume','n','en_baskin_ornek','en_baskin_ornek_orani','zenginlesme_kati','normalize_entropi']].to_markdown(index=False,floatfmt='.3f')+'\n\n')
 f.write('## Fibroblast kümesi 10\n'+summary.to_markdown(index=False,floatfmt='.3f')+'\n\n')
 f.write('Fibroblast kimliği için iki önceden hesaplanmış belirteç alanı çapraz kontrol edildi. Bunlar aynı ham veriden türetilmiştir; bağımsız yöntem veya dış referans doğrulaması sayılamaz.\n\n')
 f.write('## Karar ve kalan doğrulamalar\nBatch ile küme ilişkisi ölçüldü; fakat bu ilişki teknik batch etkisi ile gerçek örnek/hastalık bileşimi farkını tek başına ayıramaz. Batch-düzeltilmiş komşuluk/kümeleme, bağımsız referans temelli hücre tipi ataması, gen düzeyinde küme marker zenginleşmesi ve kararlılık karşılaştırması tamamlanmadan fibroblast sayısı nihai kabul edilmemelidir.\n')
print('done',len(q),'Cramers V',v,'cluster10 sample',q.ornek.nunique());print(a.sort_values('zenginlesme_kati',ascending=False).head(3).to_string(index=False))
