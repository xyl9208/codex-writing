"""Single-cell validation: GSE155468 (ATAA, 8 vs 3) and GSE213740 (TAAD, 6 vs 3), ascending aorta.

Steps: per-sample QC -> downsample -> concatenate -> normalise/log -> HVG -> PCA -> Harmony (sample) -> Leiden -> marker-based annotation
Outputs: results/sc/celltype_proportions.csv, results/sc/markers_<celltype>.csv (aorta-specific marker sets for bulk deconvolution),
         results/sc/pseudobulk_celltype_sample.csv (mean log-normalised expression per cell type x sample), UMAP figures.
"""
import os, glob, gzip, warnings, numpy as np, pandas as pd, scipy.sparse as sp
import scanpy as sc, anndata as ad
warnings.filterwarnings("ignore")
sc.settings.n_jobs = 4; sc.settings.verbosity = 1
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC = os.path.join(ROOT, "data", "geo", "sc"); OUT = os.path.join(ROOT, "results", "sc"); os.makedirs(OUT, exist_ok=True)
FIG = os.path.join(ROOT, "figures", "sc"); os.makedirs(FIG, exist_ok=True)
MAX_CELLS = 4000; rng = np.random.default_rng(0)

def qc(a, sample, dataset, group):
    a.var_names_make_unique(); a.obs_names = [f"{sample}_{c}" for c in a.obs_names]
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
    a = a[(a.obs.n_genes_by_counts >= 200) & (a.obs.n_genes_by_counts <= 6000) & (a.obs.pct_counts_mt < 20)].copy()
    if a.n_obs > MAX_CELLS: a = a[rng.choice(a.obs_names, MAX_CELLS, replace=False)].copy()
    a.obs["sample"] = sample; a.obs["dataset"] = dataset; a.obs["group"] = group
    return a

adatas = []
# GSE155468: dense gene x cell tables per sample
for f in sorted(glob.glob(os.path.join(SC, "GSE155468", "*.txt.gz"))):
    sample = os.path.basename(f).split("_")[1].replace(".txt.gz", "")
    df = pd.read_csv(f, sep="\t", index_col=0)
    a = ad.AnnData(X=sp.csr_matrix(df.values.T.astype(np.float32)), obs=pd.DataFrame(index=df.columns.astype(str)), var=pd.DataFrame(index=df.index.astype(str)))
    del df
    adatas.append(qc(a, sample, "GSE155468", "ATAA" if sample.startswith("TAA") else "Control"))
    print("GSE155468", sample, adatas[-1].shape, flush=True)
# GSE213740: 10x mtx per sample
pref = sorted(set(os.path.basename(f).rsplit("_", 1)[0] for f in glob.glob(os.path.join(SC, "GSE213740", "*_matrix.mtx.gz"))))
for p in pref:
    sample = p.split("_", 1)[1].replace("_replicate_", "")
    a = sc.read_mtx(os.path.join(SC, "GSE213740", p + "_matrix.mtx.gz")).T
    feats = pd.read_csv(os.path.join(SC, "GSE213740", p + "_features.tsv.gz"), sep="\t", header=None)
    bcs = pd.read_csv(os.path.join(SC, "GSE213740", p + "_barcodes.tsv.gz"), sep="\t", header=None)
    a.var_names = feats[1].values; a.obs_names = bcs[0].values
    a = a[:, feats[2].values == "Gene Expression"].copy() if feats.shape[1] > 2 else a
    adatas.append(qc(a, sample, "GSE213740", "TAAD" if sample.startswith("AD") else "Control"))
    print("GSE213740", sample, adatas[-1].shape, flush=True)

adata = ad.concat(adatas, join="inner", merge="same"); adata.obs_names_make_unique()
sc.pp.filter_genes(adata, min_cells=10)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.pp.pca(adata, n_comps=40, mask_var="highly_variable")
import harmonypy as hm
ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "sample", max_iter_harmony=20, verbose=False)
Z = ho.Z_corr; Z = Z.T if Z.shape[0] == 40 else Z
adata.obsm["X_harmony"] = np.asarray(Z)
sc.pp.neighbors(adata, use_rep="X_harmony", n_neighbors=15)
sc.tl.leiden(adata, resolution=0.8, flavor="igraph", n_iterations=2, key_added="leiden")
sc.tl.umap(adata)

MARKERS = {
 "SMC": ["MYH11", "ACTA2", "CNN1", "TAGLN", "MYL9", "LMOD1"],
 "Fibroblast": ["DCN", "LUM", "PDGFRA", "FBLN1", "COL1A2", "C7"],
 "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"],
 "Macrophage": ["CD68", "CD163", "C1QA", "C1QB", "AIF1", "LYZ", "MS4A7"],
 "Monocyte": ["FCN1", "S100A12", "VCAN", "LYZ", "CD14"],
 "Neutrophil": ["S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2"],
 "T_cell": ["CD3D", "CD3E", "CD2", "TRAC", "IL7R"],
 "NK": ["NKG7", "GNLY", "KLRD1", "PRF1"],
 "B_cell": ["MS4A1", "CD79A", "CD79B", "BANK1"],
 "Plasma": ["MZB1", "JCHAIN", "IGKC", "XBP1"],
 "Mast": ["TPSAB1", "CPA3", "KIT", "MS4A2"],
 "Pericyte": ["RGS5", "KCNJ8", "ABCC9", "HIGD1B"],
 "Mesothelial": ["MSLN", "UPK3B", "KRT19"],
 "Proliferating": ["MKI67", "TOP2A", "CENPF"],
}
for ct, g in MARKERS.items():
    g = [x for x in g if x in adata.var_names]; sc.tl.score_genes(adata, g, score_name=f"score_{ct}")
score_cols = [f"score_{ct}" for ct in MARKERS]
cl_scores = adata.obs.groupby("leiden")[score_cols].mean()
cl_scores.columns = [c.replace("score_", "") for c in cl_scores.columns]
zs = (cl_scores - cl_scores.mean()) / cl_scores.std()
annot = zs.idxmax(axis=1)
# proliferating cells -> assign to second-best lineage label, keep note
for cl in annot.index:
    if annot[cl] == "Proliferating":
        annot[cl] = zs.loc[cl].drop("Proliferating").idxmax() + "_cycling"
adata.obs["celltype"] = adata.obs["leiden"].map(annot).astype(str)
adata.obs["celltype_major"] = adata.obs["celltype"].str.replace("_cycling", "")
cl_scores.to_csv(os.path.join(OUT, "cluster_marker_scores.csv")); annot.to_csv(os.path.join(OUT, "cluster_annotation.csv"))
print(adata.obs.celltype_major.value_counts(), flush=True)

# proportions per sample
prop = pd.crosstab(adata.obs["sample"], adata.obs["celltype_major"], normalize="index")
meta = adata.obs.drop_duplicates("sample").set_index("sample")[["dataset", "group"]]
prop = meta.join(prop); prop.to_csv(os.path.join(OUT, "celltype_proportions.csv"))
print(prop.groupby(["dataset", "group"]).mean(numeric_only=True).round(3).T.to_string(), flush=True)

# aorta-specific marker sets (from control cells only, both datasets) for bulk ssGSEA deconvolution
ctrl = adata[adata.obs.group == "Control"].copy()
keep = ctrl.obs.celltype_major.value_counts(); keep = keep[keep >= 30].index
ctrl = ctrl[ctrl.obs.celltype_major.isin(keep)].copy()
sc.tl.rank_genes_groups(ctrl, "celltype_major", method="wilcoxon", pts=True)
marker_sets = {}
for ct in keep:
    df = sc.get.rank_genes_groups_df(ctrl, ct)
    df = df[(df.logfoldchanges > 1.5) & (df.pvals_adj < 0.01) & (df.pct_nz_group > 0.3)]
    df = df[~df.names.str.startswith(("MT-", "RPL", "RPS"))].head(50)
    marker_sets[ct] = df.names.tolist(); df.to_csv(os.path.join(OUT, f"markers_{ct}.csv"), index=False)
with open(os.path.join(OUT, "aorta_celltype_markers.gmt"), "w") as fh:
    for ct, g in marker_sets.items(): fh.write("\t".join([ct, "scRNA_control_aorta"] + g) + "\n")

# pseudobulk mean log-normalised expression per cell type x sample (for signature-gene validation)
X = adata.X.tocsr(); rows = []
for (s, ct), idx in adata.obs.groupby(["sample", "celltype_major"]).indices.items():
    if len(idx) < 10: continue
    mu = np.asarray(X[idx].mean(axis=0)).ravel()
    rows.append(pd.Series(mu, index=adata.var_names, name=f"{s}|{ct}"))
pb = pd.DataFrame(rows).T; pb.to_csv(os.path.join(OUT, "pseudobulk_celltype_sample.csv"))
# pseudobulk counts (sum) per sample x cell type for DE-style comparisons
C = adata.layers["counts"].tocsr(); rows = []
for (s, ct), idx in adata.obs.groupby(["sample", "celltype_major"]).indices.items():
    if len(idx) < 10: continue
    rows.append(pd.Series(np.asarray(C[idx].sum(axis=0)).ravel(), index=adata.var_names, name=f"{s}|{ct}"))
pd.DataFrame(rows).T.to_csv(os.path.join(OUT, "pseudobulk_counts_celltype_sample.csv"))
adata.obs.to_csv(os.path.join(OUT, "cell_metadata.csv"))
np.save(os.path.join(OUT, "umap.npy"), adata.obsm["X_umap"])

# figures
sc.settings.figdir = FIG
sc.pl.umap(adata, color=["celltype_major"], legend_loc="on data", legend_fontsize=6, save="_celltype.png", show=False, frameon=False, title="Cell types (GSE155468 + GSE213740, Harmony)")
sc.pl.umap(adata, color=["dataset", "group"], save="_dataset_group.png", show=False, frameon=False, ncols=2)
adata.write_h5ad(os.path.join(OUT, "aorta_sc_integrated.h5ad"), compression="gzip")
print("done", adata.shape, flush=True)
