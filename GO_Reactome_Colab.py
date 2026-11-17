# GSE171213 fibroblast alt-kume 1 vs 0: exploratory ranked pathway analysis
# Run in Google Colab with internet, after uploading this folder's input CSV.
import sys, subprocess
subprocess.check_call([sys.executable,'-m','pip','-q','install','gseapy>=1.1.0'])
import pandas as pd, numpy as np, gseapy as gp
from pathlib import Path
out=Path('GO_Reactome_sonuclar');out.mkdir(exist_ok=True)
df=pd.read_csv('AltKume_Pseudobulk_Gen_Testleri.csv')
# Gene-level p-values are discrete with 7 paired samples; avoid transforming to pseudo-z.
# Rank by signed mean log2(1+CPM) differences; descriptive exploratory GSEA only.
r=df[['gene','mean_delta_log2_1p_CPM']].dropna().copy()
r=r[~r.gene.duplicated()]
r=r.sort_values('mean_delta_log2_1p_CPM',ascending=False)
r.columns=['gene','score']
r.to_csv(out/'Kullanilan_Gen_Siralamasi.csv',index=False)
for library in ['GO_Biological_Process_2023','Reactome_2022']:
    try:
        # Enrichr's curated library is fetched from online service by GSEApy.
        gmt=gp.get_library(name=library,organism='Human')
        matched={k:[g for g in v if g in set(r.gene)] for k,v in gmt.items()}
        matched={k:v for k,v in matched.items() if 10<=len(v)<=500}
        print(library, 'gene sets',len(matched))
        res=gp.prerank(rnk=r,gene_sets=matched,min_size=10,max_size=500,
                       permutation_num=1000,seed=20261116,threads=2,
                       outdir=None,verbose=False)
        tab=res.res2d.copy();tab.to_csv(out/(library+'_GSEA.csv'),index=False)
        print(tab[['Term','NES','NOM p-val','FDR q-val']].head(12).to_string(index=False))
    except Exception as e:
        print('ERROR',library,repr(e))
print('Done. Note: exploratory cell-state comparison; pathway q-values do not imply clinical group differences.')
