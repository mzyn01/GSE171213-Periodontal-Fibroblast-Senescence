# GSE171213 — Periodontal Fibroblast Transcriptional Programs

## Overview
This repository contains archived Python scripts and selected audit records for an exploratory reanalysis of the public single-cell RNA-seq dataset **GSE171213**. The analysis examines periodontal fibroblast heterogeneity, extracellular-matrix-associated pathways, and senescence/SASP-associated transcriptional signatures.

**Interpretation:** The two computational fibroblast subclusters are weakly separated (silhouette approximately 0.0744). The results do **not** establish two discrete biological fibroblast types or prove cellular senescence. Post-treatment fibroblasts are sparsely represented (28 cells), so treatment-related conclusions are not supported.

## Data source
- GEO accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171213
- Study groups: HC (4 samples), PD (5 samples), PDT (3 samples).
- Raw count matrices should be obtained from GEO; raw sequencing/count data are **not redistributed** here.

## Archived analysis components
The uploaded scripts cover the following stages (filenames may differ slightly across analysis revisions):

1. Raw-count QC and QC decision audits: `GSE171213_ortak_qc_v1.py`, `GSE171213_Ham_QC_Dogrulama_v2.py`, `GSE171213_QC_Karar_Denetimi_v3.py`.
2. Candidate doublet assessment and audit: `GSE171213_HC1_Doublet_Pilot_v1.py`, `GSE171213_HC1_Doublet_Degerlendirme_v2.py` and associated audit records.
3. Global clustering and integration: `GSE171213_Ortak_Kumeleme_v2.py`, `GSE171213_Harmony_Colab_Analiz.ipynb`.
4. Fibroblast identification and subclustering: `GSE171213_Fibroblast_Gen_Dogrulama_v1.py`, `GSE171213_Fibro_AltKume_Analiz_v1.py`.
5. Sample-level and senescence-signature analyses: `GSE171213_Fibro_Ornek_Duzeyi_v2.py`, `GSE171213_Pseudobulk_Referans_v4.py` and reference-check scripts.
6. Preranked GSEA: `GO_Reactome_Colab.py` (GO Biological Process 2023 and Reactome 2022).

Some scripts are intermediate analysis/audit versions, not a single end-to-end pipeline. **Do not assume the entire study can be reproduced by running every script sequentially without reconstructing inputs and paths.** Consult the manuscript Methods and script headers before execution.

## Main reported data flow
68,542 initial cells → 44,143 cells after QC → 43,438 cells after excluding 705 high-doublet-suspicion cells. The retained candidate fibroblast population comprises 1,343 cells, divided computationally into 540 and 803 cells.

## Code integrity
`SHA256_Dosya_Manifestosu.csv` contains SHA-256 hashes for the original archived files. **Important:** the manifest may record paths relative to the original archive (`Kaynaklar/`, `Denetim/`). If the files were uploaded into the GitHub repository root, those paths must be mapped to the new layout before running automated integrity checks. Hashes verify file contents, not the scientific validity of the analysis.

## Software environment and reproducibility limitations
The original Google Colab environment's complete package-version snapshot was not preserved. Do not represent versions from a later environment as the versions used for the original analysis. `Orijinal_Colab_Surum_Kayit_Kodu.py` can record versions only for the environment in which it is executed.

This archive does not include all intermediate data and outputs. Some scripts may require manually downloaded GEO matrices, locally prepared tables, and adapted file paths. An exact one-command reproduction workflow is **not yet available**.

## Statistical scope
The subcluster comparison and pathway enrichment are exploratory. Gene-level pseudobulk comparisons did not yield FDR-significant genes. Reactome collagen chain trimerization was enriched toward computational subcluster 0; senescence signatures did not provide FDR-significant, directionally consistent evidence for definitive senescence.

## Citation and archival DOI
A version-specific Zenodo DOI will be added **only after** a public release is deposited and the DOI is confirmed. Until then, cite the GEO accession and the repository URL:
https://github.com/mzyn01/GSE171213-Periodontal-Fibroblast-Senescence

## License
No license has been selected yet. A license should be chosen before public release, taking into account authorship and institutional requirements.


