# GSE171213 — 705 hücre kimlik denetimi

Kaynak dosyalar: `GSE171213_Ham_QC_Hucre_Metrikleri_v2.csv`, `GSE171213_Doublet_44143_Hucre_Karar_v5.csv`, `GSE171213_43438_Korunan_Hucre_Anotasyon_v1.csv`, `GSE171213_Doublet_12_Ornek_Rapor_v5.md`.

- Ham 68.542; QC30 44.143; yüksek doublet şüphesiyle çıkarılan 705; kalan 43.438.
- QC30 hücrelerinin tamamı doublet karar tablosunda mevcut; 43.438 korunan hücrenin tamamı nihai anotasyon tablosunda mevcut.
- Tüm (ornek, CellIndex) birleşik kimlikleri üç dosyada benzersiz.
- Çıkarma koşulu: `simulated_neighbor_fraction >= 0.60` VE `lineages_with_2markers >= 2`. 705/705 çıkarılan hücre koşulu sağlıyor; korunanların hiçbiri sağlamıyor.
- 705 hücrenin hiçbiri nihai anotasyon dosyasında yok.
- Bu, operasyonel yüksek-şüphe filtrelemesidir; kalibre edilmiş doublet olasılığı ya da kesin doublet tanısı değildir.
- Sentetik doublet taraması: örnek bazında 1:1 sentetik doublet, 10.000'e normalize edilmiş log1p ifade, 2.000 değişken gen, 30 TruncatedSVD bileşeni, 40 komşuda sentetik komşu oranı; 5 hücre soyuna ait belirteçler değerlendirilmiş.
- Bu denetim hücre kimliği zincirini doğrular; doublet filtresinin biyolojik duyarlılık/özgüllüğünü bağımsız olarak doğrulamaz.
