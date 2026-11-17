# GSE171213 — Kod ve yazılım sürümü denetimi

Arşivlenen kod/not defteri: 23 dosya. Toplam arşivlenen kaynak: 26 dosya.

## Doğrulanabilenler
- Harmony not defterinde `harmonypy` ve `scanpy` importları, `random_state=42`, `max_iter_harmony=20` ve KMeans parametreleri mevcut.
- GSEA kodu `GSE171213_GO_Reactome_v7/GO_Reactome_Colab.py` arşivlendi.
- Pseudobulk ve fibroblast alt küme analiz betikleri arşivlendi.
- Kaynak dosyaların SHA-256 değerleri `SHA256_Dosya_Manifestosu.csv` dosyasında.

## Henüz doğrulanamayan
- Özgün Google Colab oturumunun Python ve paket sürümlerine ait tam `pip freeze` / `environment.yml` kaydı mevcut dosyalarda saptanamadı.
- Bugünkü çalışma ortamındaki paket sürümleri özgün analiz ortamının sürümleri değildir; bu nedenle rapora aktarılmadı.
- `Orijinal_Colab_Surum_Kayit_Kodu.py` özgün oturum halen erişilebiliyorsa orada çalıştırılmalıdır. Oturum yeniden oluşturulmuşsa çıktı yalnızca yeni ortamı belgeler.
- Bu arşiv kaynak kodu ve seçili yöntem kayıtlarını içerir; ham veri ve bütün ara çıktılar pakete dahil edilmemiştir. Dolayısıyla tek başına eksiksiz yeniden çalıştırılabilir bir veri paketi değildir.

## Metne önerilen ifade
“Analiz betikleri, Harmony not defteri ve seçili denetim kayıtları sürüm kontrollü bir ek arşivde saklandı. Arşivlenen dosyaların SHA-256 sağlama toplamları kaydedildi. Özgün analiz ortamına ait tüm yazılım sürümleri geriye dönük olarak doğrulanamadığından, bunlar kesin sürüm bilgisi olarak raporlanmadı.”
