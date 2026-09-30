# Ascending aortic aneurysm vs. acute type A dissection: integrative transcriptomic meta-analysis

Reproducible analysis accompanying the manuscript *"Ascending thoracic aortic aneurysm and acute type A aortic dissection are governed by largely distinct molecular programmes: an integrative meta-analysis of human aortic transcriptomes with single-cell validation"* (target journal: Journal of Translational Medicine).

## Layout

| Path | Content |
|---|---|
| `scripts/00_search_geo.py` | GEO DataSets search |
| `scripts/01_inspect_gse.py` | Inspect series metadata / files |
| `scripts/02_download.sh` | Download processed matrices, raw files and platform annotations |
| `scripts/03_preprocess.py` | Harmonised per-cohort pre-processing (symbol harmonisation, re-processing from raw where needed) |
| `scripts/04_de.py` | Per-cohort differential expression (limma-trend / PyDESeq2) + QC figures |
| `scripts/05_meta.py` | Random-effects meta-analysis, LOCO robustness, ATAA-vs-TAAD gene classes, direct-cohort validation, cohort concordance |
| `scripts/06_pathways.py` | GSEA (per cohort + meta-rankings), pathway classification, ORA |
| `scripts/07_singlecell.py`, `07b_sc_refine.py` | scRNA-seq integration (GSE155468 + GSE213740), annotation, marker sets, pseudo-bulk |
| `scripts/08_cellcomp.py` | ssGSEA cell-composition inference in bulk cohorts |
| `scripts/09_signature_ml.py` | Control-referenced signature scores, LOCO elastic-net classifier |
| `scripts/10_sc_validation.py` | Cell-type-resolved validation of gene classes |
| `scripts/11_ppi.py` | STRING PPI networks and hubs |
| `scripts/12_fig1_table1.py` | Figure 1 schematic and Table 1 |
| `scripts/13_export_docx.py` | Manuscript export to Word |
| `results/` | Result tables (DE, meta, pathways, cellcomp, ml, ppi, sc, sc_validation) |
| `figures/` | Main figures (PNG + PDF), QC and single-cell figures |
| `manuscript/` | `manuscript.md`, `manuscript_JTM.docx`, `Table1_cohorts.csv`, `README_zh.md` |

## Reproduce

```bash
pip install pandas numpy scipy statsmodels scikit-learn matplotlib seaborn gseapy pydeseq2 inmoose scanpy leidenalg igraph harmonypy python-docx networkx requests
bash scripts/02_download.sh                 # ~4 GB of GEO data (raw tars for GSE318877, GSE98770, GSE155468, GSE213740 downloaded separately, see script comments)
python scripts/03_preprocess.py && python scripts/04_de.py && python scripts/05_meta.py && python scripts/06_pathways.py
python scripts/07_singlecell.py && python scripts/07b_sc_refine.py && python scripts/08_cellcomp.py
python scripts/09_signature_ml.py && python scripts/10_sc_validation.py && python scripts/11_ppi.py && python scripts/12_fig1_table1.py && python scripts/13_export_docx.py
```

Raw data (`data/geo/`) and large intermediates are not committed (see `.gitignore`).
