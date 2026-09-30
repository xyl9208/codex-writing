"""Disease-stage single-cell analysis: GSE222318 (TAAD acute n=4 / subacute n=3 / chronic n=2 / normal n=5)
and GSE189795 (acute TAAD n=5 vs control n=4; external acute replicate).

Both datasets are processed with the same pipeline as the atlas (QC -> <=4000 cells/sample -> normalise -> HVG -> PCA ->
Harmony(sample) -> Leiden -> marker-based annotation into the atlas categories).  Module scores are then computed at the
PATIENT x CELL-TYPE level (mean of per-gene z-scores relative to control samples of the same dataset, within cell type),
so that the statistical unit is the patient.
"""
import os, gzip, warnings, numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad
from scipy import stats
warnings.filterwarnings("ignore"); sc.settings.n_jobs = 4; sc.settings.verbosity = 0
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC = os.path.join(ROOT, "data", "geo", "sc"); OUT = os.path.join(ROOT, "results", "sc_stage"); os.makedirs(OUT, exist_ok=True); MOD = os.path.join(ROOT, "results", "modules")
MAX_CELLS = 4000; rng = np.random.default_rng(0)
modules = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(os.path.join(MOD, "modules.gmt"))}
mtab = pd.read_csv(os.path.join(MOD, "modules_table.csv")).set_index("module")
ORDER = ["SMC_contractile", "SMC_modulated", "Fibromyocyte", "Fibroblast", "Endothelial", "Macrophage", "Inflammatory_myeloid", "T_cell", "NK", "B_cell", "Plasma", "Mast"]
MARKERS = {"SMC": ["MYH11", "ACTA2", "CNN1", "TAGLN", "MYL9", "LMOD1"], "Fibromyocyte": ["LTBP2", "EFEMP1", "DKK3", "COL1A1", "BGN", "AEBP1"], "Fibroblast": ["DCN", "LUM", "PDGFRA", "FBLN1", "C7", "SERPINF1"],
           "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"], "Macrophage": ["CD68", "CD163", "C1QA", "C1QB", "AIF1", "MS4A7"], "Inflammatory_myeloid": ["S100A8", "S100A9", "S100A12", "FCN1", "CSF3R", "FCGR3B"],
           "T_cell": ["CD3D", "CD3E", "CD2", "TRAC", "IL7R"], "NK": ["NKG7", "GNLY", "KLRD1", "PRF1"], "B_cell": ["MS4A1", "CD79A", "CD79B", "BANK1"], "Plasma": ["MZB1", "JCHAIN", "IGKC", "XBP1"], "Mast": ["TPSAB1", "CPA3", "KIT", "MS4A2"],
           "Doublet": ["MKI67", "TOP2A"]}

def read_dense_csv_subsampled(path, sample_of_col, sep=",", chunk=1500):
    """Read a genes x cells dense CSV in row chunks, keeping <= MAX_CELLS randomly chosen cells per sample."""
    with gzip.open(path, "rt") as fh: header = fh.readline().rstrip("\n").split(sep)
    cells = [h.strip('"') for h in header[1:]]; samples = pd.Series([sample_of_col(c) for c in cells], index=cells)
    keep = []
    for s, idx in samples.groupby(samples).groups.items():
        idx = list(idx); keep += list(rng.choice(idx, min(MAX_CELLS, len(idx)), replace=False))
    keep = sorted(set(keep), key=cells.index); col_idx = [cells.index(c) + 1 for c in keep]
    blocks = []; genes = []
    for ch in pd.read_csv(path, sep=sep, usecols=[0] + col_idx, chunksize=chunk, index_col=0):
        genes += [str(g) for g in ch.index]; blocks.append(sp.csr_matrix(ch.values.astype(np.float32)))
    X = sp.vstack(blocks).T.tocsr()
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=keep), var=pd.DataFrame(index=genes)); a.var_names_make_unique(); a.obs["sample"] = samples.loc[keep].values
    return a

def process(a, dataset):
    a.var["mt"] = a.var_names.str.startswith("MT-"); sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
    a = a[(a.obs.n_genes_by_counts >= 200) & (a.obs.n_genes_by_counts <= 6000) & (a.obs.pct_counts_mt < 20)].copy()
    sc.pp.filter_genes(a, min_cells=10); a.layers["counts"] = a.X.copy(); sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, batch_key="sample"); sc.pp.pca(a, n_comps=40, mask_var="highly_variable")
    import harmonypy as hm
    Z = hm.run_harmony(a.obsm["X_pca"], a.obs, "sample", max_iter_harmony=20, verbose=False).Z_corr; Z = Z.T if Z.shape[0] == 40 else Z
    a.obsm["X_harmony"] = np.asarray(Z); sc.pp.neighbors(a, use_rep="X_harmony", n_neighbors=15); sc.tl.leiden(a, resolution=0.8, flavor="igraph", n_iterations=2)
    for ct, g in MARKERS.items(): sc.tl.score_genes(a, [x for x in g if x in a.var_names], score_name=f"s_{ct}")
    cl = a.obs.groupby("leiden")[[f"s_{ct}" for ct in MARKERS]].mean(); cl.columns = [c[2:] for c in cl.columns]
    lab = {}
    for c in cl.index:
        r = cl.loc[c].drop("Doublet"); best = r.idxmax()
        if best == "SMC": best = "SMC_contractile" if r["SMC"] > 1.2 * max(r["Fibromyocyte"], 0.01) else "SMC_modulated"
        if r.drop("SMC").max() > 0.4 and r["SMC"] > 0.8 and best in ("Macrophage", "Inflammatory_myeloid", "T_cell", "Endothelial"): best = "Doublet_like"
        lab[c] = best
    a.obs["celltype"] = a.obs.leiden.map(lab).astype(str); a = a[a.obs.celltype != "Doublet_like"].copy()
    a.obs["dataset"] = dataset; return a

def patient_celltype_module_scores(a, ctrl_label):
    """z per gene relative to control samples within each cell type (pseudo-bulk mean log-normalised expression), then module mean."""
    X = a.X.tocsr(); rows = []
    for (s, ct), idx in a.obs.groupby(["sample", "celltype"]).indices.items():
        if len(idx) < 20: continue
        rows.append(pd.Series(np.asarray(X[idx].mean(axis=0)).ravel(), index=a.var_names, name=f"{s}|{ct}"))
    pb = pd.DataFrame(rows).T; cols = pd.DataFrame([c.split("|") for c in pb.columns], columns=["sample", "celltype"], index=pb.columns)
    grp = a.obs.drop_duplicates("sample").set_index("sample")["group"]; cols["group"] = grp.loc[cols["sample"]].values
    out = []
    for ct in cols.celltype.unique():
        cc = cols.index[cols.celltype == ct]; ref = [c for c in cc if cols.loc[c, "group"] == ctrl_label]
        if len(ref) < 2: continue
        x = pb[cc]; mu = x[ref].mean(axis=1); sd = x[ref].std(axis=1).clip(lower=0.05); z = x.sub(mu, axis=0).div(sd, axis=0).clip(-10, 10)
        for k, v in modules.items():
            g = [y for y in v if y in z.index and pb.loc[y, cc].mean() > 0.05]  # genes detectably expressed in this cell type
            if len(g) < 5: continue
            sgn = 1 if mtab.loc[k, "direction"] == "up" else -1
            for c in cc: out.append(dict(sample=cols.loc[c, "sample"], celltype=ct, group=cols.loc[c, "group"], module=k, score=sgn * z.loc[g, c].mean(), n_genes=len(g)))
    return pd.DataFrame(out), pb

# ---------------------------------------------------------------- GSE222318
meta = pd.read_csv(os.path.join(SC, "GSE222318_Cell_metadata.csv.gz"), index_col=0)
samp_of = lambda c: meta.loc[c, "orig.ident"] if c in meta.index else c.split("_")[-1]
a1 = read_dense_csv_subsampled(os.path.join(SC, "GSE222318_Raw_gene_counts_matrix.csv.gz"), samp_of)
a1.obs = a1.obs.join(meta[["Group1", "Group2", "Position", "Sex"]]); a1.obs["group"] = a1.obs["Group2"].replace({"normal": "Control", "Normal": "Control"})
print("GSE222318 loaded:", a1.shape, a1.obs.groupby(["sample", "group"]).size().to_dict(), flush=True)
a1 = process(a1, "GSE222318"); print("GSE222318 processed:", a1.shape, a1.obs.celltype.value_counts().to_dict(), flush=True)
s1, pb1 = patient_celltype_module_scores(a1, "Control"); s1.to_csv(os.path.join(OUT, "GSE222318_patient_celltype_module_scores.csv"), index=False); pb1.to_csv(os.path.join(OUT, "GSE222318_pseudobulk_mean.csv"))
a1.obs.to_csv(os.path.join(OUT, "GSE222318_cell_metadata.csv")); pd.crosstab(a1.obs["sample"], a1.obs["celltype"], normalize="index").join(a1.obs.drop_duplicates("sample").set_index("sample")[["group"]]).to_csv(os.path.join(OUT, "GSE222318_celltype_proportions.csv"))
del a1

# ---------------------------------------------------------------- GSE189795 (dense normalised counts; treated as normalised expression -> log1p)
a2 = read_dense_csv_subsampled(os.path.join(SC, "GSE189795_AllSample.NormalizedCounts.txt.gz"), lambda c: c.split("_")[0], sep="\t")
a2.obs["group"] = np.where(a2.obs["sample"].str.startswith("AD"), "Acute", "Control")
print("GSE189795 loaded:", a2.shape, a2.obs.groupby(["sample", "group"]).size().to_dict(), flush=True)
# values are already normalised counts: skip normalize_total but keep the same downstream steps
a2.var["mt"] = a2.var_names.str.startswith("MT-"); sc.pp.calculate_qc_metrics(a2, qc_vars=["mt"], inplace=True, percent_top=None)
a2 = a2[(a2.obs.n_genes_by_counts >= 200)].copy(); sc.pp.filter_genes(a2, min_cells=10); a2.layers["counts"] = a2.X.copy(); sc.pp.log1p(a2)
sc.pp.highly_variable_genes(a2, n_top_genes=2000, batch_key="sample"); sc.pp.pca(a2, n_comps=40, mask_var="highly_variable")
import harmonypy as hm
Z = hm.run_harmony(a2.obsm["X_pca"], a2.obs, "sample", max_iter_harmony=20, verbose=False).Z_corr; Z = Z.T if Z.shape[0] == 40 else Z
a2.obsm["X_harmony"] = np.asarray(Z); sc.pp.neighbors(a2, use_rep="X_harmony", n_neighbors=15); sc.tl.leiden(a2, resolution=0.8, flavor="igraph", n_iterations=2)
for ct, g in MARKERS.items(): sc.tl.score_genes(a2, [x for x in g if x in a2.var_names], score_name=f"s_{ct}")
cl = a2.obs.groupby("leiden")[[f"s_{ct}" for ct in MARKERS]].mean(); cl.columns = [c[2:] for c in cl.columns]; lab = {}
for c in cl.index:
    r = cl.loc[c].drop("Doublet"); best = r.idxmax()
    if best == "SMC": best = "SMC_contractile" if r["SMC"] > 1.2 * max(r["Fibromyocyte"], 0.01) else "SMC_modulated"
    lab[c] = best
a2.obs["celltype"] = a2.obs.leiden.map(lab).astype(str); a2.obs["dataset"] = "GSE189795"
print("GSE189795 processed:", a2.shape, a2.obs.celltype.value_counts().to_dict(), flush=True)
s2, pb2 = patient_celltype_module_scores(a2, "Control"); s2.to_csv(os.path.join(OUT, "GSE189795_patient_celltype_module_scores.csv"), index=False); pb2.to_csv(os.path.join(OUT, "GSE189795_pseudobulk_mean.csv"))
a2.obs.to_csv(os.path.join(OUT, "GSE189795_cell_metadata.csv"))

# ---------------------------------------------------------------- stage comparison (patient level)
def hedges(x, y):
    nx_, ny = len(x), len(y)
    if nx_ < 2 or ny < 2: return np.nan, np.nan
    spv = np.sqrt(((nx_ - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / (nx_ + ny - 2)); g = (np.mean(x) - np.mean(y)) / spv * (1 - 3 / (4 * (nx_ + ny) - 9))
    return g, np.sqrt((nx_ + ny) / (nx_ * ny) + g ** 2 / (2 * (nx_ + ny)))
rows = []
for ds, s in [("GSE222318", s1), ("GSE189795", s2)]:
    for (ct, k), d in s.groupby(["celltype", "module"]):
        ctrl = d.loc[d.group == "Control", "score"].values
        for stage in [g for g in d.group.unique() if g != "Control"]:
            x = d.loc[d.group == stage, "score"].values
            if len(x) < 2 or len(ctrl) < 2: continue
            g, se = hedges(x, ctrl); rows.append(dict(dataset=ds, celltype=ct, module=k, stage=stage, n_stage=len(x), n_ctrl=len(ctrl), mean_diff=x.mean() - ctrl.mean(), g=g, se=se, p_mwu=stats.mannwhitneyu(x, ctrl).pvalue))
st = pd.DataFrame(rows); st.to_csv(os.path.join(OUT, "stage_module_effects.csv"), index=False)
pd.set_option("display.width", 250)
for ct in ["SMC_contractile", "SMC_modulated", "Fibromyocyte", "Fibroblast", "Endothelial", "Macrophage"]:
    sub = st[(st.dataset == "GSE222318") & (st.celltype == ct)].pivot(index="module", columns="stage", values="g")
    if len(sub): print(f"\nGSE222318 {ct}: Hedges' g vs control by stage\n", sub.round(2).to_string())
print("\nGSE189795 (acute vs control) SMC_contractile / Fibroblast:\n", st[(st.dataset == "GSE189795") & (st.celltype.isin(["SMC_contractile", "Fibroblast"]))].pivot(index="module", columns="celltype", values="g").round(2).to_string())
print("done", flush=True)
