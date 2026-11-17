# Bu kodu ÖZGÜN Colab çalışma ortamında çalıştırın.
import sys, platform, importlib.metadata as md, json
packages = ['numpy','pandas','scipy','scikit-learn','scanpy','anndata','harmonypy','gseapy','umap-learn','matplotlib','seaborn']
versions = {}
for package in packages:
    try: versions[package] = md.version(package)
    except md.PackageNotFoundError: versions[package] = 'kurulu değil / bulunamadı'
record = {'python':sys.version, 'platform':platform.platform(), 'packages':versions}
print(json.dumps(record,ensure_ascii=False,indent=2))
with open('GSE171213_Orijinal_Ortam_Surumleri.json','w',encoding='utf-8') as f:
    json.dump(record,f,ensure_ascii=False,indent=2)
