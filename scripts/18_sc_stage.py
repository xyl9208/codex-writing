"""v2.1 disease-stage single-cell analysis.

GSE222318 (TAAD acute n=4 [intima-media] / subacute n=3 [intima-media] / chronic n=2 [whole wall] / normal n=5 [whole wall]):
  ESTIMATION ONLY (Position confounded with stage; n too small for any inferential test) -> patient x cell-type Tier 1
  programme scores C / I / R with label-free standardisation; effects (Hedges' g, bootstrap CI) for
  (i) each stage vs control, (ii) acute vs non-acute (subacute+chronic), (iii) acute vs subacute (same Position).
GSE189795 (acute TAAD n=5 vs control n=4): PRIMARY inference = R in contractile SMC, exact enumeration (126 assignments).
Both datasets are processed with the same pipeline as the atlas; cells are down-sampled to <= 4000 per sample.
Outputs include the patient / Position / cell-type / n_cells index so that 'patient-level' is verifiable."""
import os, sys, gzip, warnings, numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad
from scipy import stats
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from v21_lib import *
warnings.filterwarnings("ignore"); sc.settings.n_jobs = 4; sc.settings.verbosity = 0
SC = os.path.join(ROOT, "data", "geo", "sc"); OUT = os.path.join(ROOT, "results", "v21_stage"); os.makedirs(OUT, exist_ok=True)
MOD = load_tier1(); MAX_CELLS = 4000; rng = np.random.default_rng(0)
MARKERS = {"SMC": ["MYH11", "ACTA2", "CNN1", "TAGLN", "MYL9", "LMOD1"], "Fibromyocyte": ["LTBP2", "EFEMP1", "DKK3", "COL1A1", "BGN", "AEBP1"], "Fibroblast": ["DCN", "LUM", "PDGFRA", "FBLN1", "C7", "SERPINF1"],
           "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"], "Macrophage": ["CD68", "CD163", "C1QA", "C1QB", "AIF1", "MS4A7"], "Inflammatory_myeloid": ["S100A8", "S100A9", "S100A12", "FCN1", "CSF3R", "FCGR3B"],
           "T_cell": ["CD3D", "CD3E", "CD2", "TRAC", "IL7R"], "NK": ["NKG7", "GNLY", "KLRD1", "PRF1"], "B_cell": ["MS4A1", "CD79A", "CD79B", "BANK1"], "Plasma": ["MZB1", "JCHAIN", "IGKC", "XBP1"], "Mast": ["TPSAB1", "CPA3", "KIT", "MS4A2"]}

def read_dense_csv_subsampled(path, sample_of_col, sep=",", chunk=1500):
    with gzip.open(path, "rt") as fh: header = fh.readline().rstrip("\n").split(sep)
    cells = [h.strip('"') for h in header[1:]]; samples = pd.Series([sample_of_col(c) for c in cells], index=cells); keep = []
    for s, idx in samples.groupby(samples).groups.items(): idx = list(idx); keep += list(rng.choice(idx, min(MAX_CELLS, len(idx)), replace=False))
    keep = sorted(set(keep), key=cells.index); col_idx = [cells.index(c) + 1 for c in keep]; blocks = []; genes = []
    for ch in pd.read_csv(path, sep=sep, usecols=[0] + col_idx, chunksize=chunk, index_col=0): genes += [str(g) for g in ch.index]; blocks.append(sp.csr_matrix(ch.values.astype(np.float32)))
    a = ad.AnnData(X=sp.vstack(blocks).T.tocsr(), obs=pd.DataFrame(index=keep), var=pd.DataFrame(index=genes)); a.var_names_make_unique(); a.obs["sample"] = samples.loc[keep].values; return a

def annotate(a, normalise=True):
    a.var["mt"] = a.var_names.str.startswith("MT-"); sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
    a = a[(a.obs.n_genes_by_counts >= 200) & (a.obs.n_genes_by_counts <= 6000) & (a.obs.pct_counts_mt < 20)].copy(); sc.pp.filter_genes(a, min_cells=10); a.layers["counts"] = a.X.copy()
    if normalise: sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a); sc.pp.highly_variable_genes(a, n_top_genes=2000, batch_key="sample"); sc.pp.pca(a, n_comps=40, mask_var="highly_variable")
    import harmonypy as hm
    Z = hm.run_harmony(a.obsm["X_pca"], a.obs, "sample", max_iter_harmony=20, verbose=False).Z_corr; Z = Z.T if Z.shape[0] == 40 else Z; a.obsm["X_harmony"] = np.asarray(Z)
    sc.pp.neighbors(a, use_rep="X_harmony", n_neighbors=15); sc.tl.leiden(a, resolution=0.8, flavor="igraph", n_iterations=2)
    for ct, g in MARKERS.items(): sc.tl.score_genes(a, [x for x in g if x in a.var_names], score_name=f"s_{ct}")
    cl = a.obs.groupby("leiden")[[f"s_{ct}" for ct in MARKERS]].mean(); cl.columns = [c[2:] for c in cl.columns]; lab = {}
    for c in cl.index:
        r = cl.loc[c]; best = r.idxmax()
        if best == "SMC": best = "SMC_contractile" if r["SMC"] > 1.2 * max(r["Fibromyocyte"], 0.01) else "SMC_modulated"
        if r["SMC"] > 0.8 and best in ("Macrophage", "Inflammatory_myeloid", "T_cell", "Endothelial") and r[best] > 0.4: best = "Doublet_like"
        lab[c] = best
    a.obs["celltype"] = a.obs.leiden.map(lab).astype(str); return a[a.obs.celltype != "Doublet_like"].copy()

def stage_scores(a, ds):
    pb, ncells = celltype_pseudobulk(a, min_cells=20); cols = pd.DataFrame([c.split("|") for c in pb.columns], columns=["sample", "celltype"], index=pb.columns)
    sm = a.obs.drop_duplicates("sample").set_index("sample"); cols = cols.join(sm[[c for c in ["group", "Position", "Sex"] if c in sm.columns]], on="sample"); cols["n_cells"] = ncells; cols["dataset"] = ds
    out = []
    for ct in cols.celltype.unique():
        cc = cols.index[cols.celltype == ct]
        if len(cc) < 4: continue
        x = pb[cc]; x = x[(x > 0).mean(axis=1) >= 0.2]; ms, n = module_scores(x, MOD); pr = program_scores(ms); s = pd.concat([pr, ms], axis=1); s = s.join(cols.loc[cc, [c for c in ["sample", "group", "Position", "n_cells"] if c in cols.columns]]); s["celltype"] = ct; s["dataset"] = ds; out.append(s)
    return pd.concat(out), cols

# ---------------------------------------------------------------- GSE222318 (estimation only)
meta = pd.read_csv(os.path.join(SC, "GSE222318_Cell_metadata.csv.gz"), index_col=0)
a1 = read_dense_csv_subsampled(os.path.join(SC, "GSE222318_Raw_gene_counts_matrix.csv.gz"), lambda c: meta.loc[c, "orig.ident"] if c in meta.index else c.split("_")[-1])
a1.obs = a1.obs.join(meta[["Group2", "Position", "Sex"]]); a1.obs["group"] = a1.obs["Group2"].replace({"normal": "Control", "Normal": "Control"})
print("GSE222318 loaded:", a1.shape, flush=True); a1 = annotate(a1, normalise=True); print("GSE222318 annotated:", a1.shape, a1.obs.celltype.value_counts().to_dict(), flush=True)
s1, idx1 = stage_scores(a1, "GSE222318"); s1.to_csv(os.path.join(OUT, "GSE222318_patient_celltype_scores.csv"), index=False); idx1.to_csv(os.path.join(OUT, "GSE222318_pseudobulk_index.csv"))
a1.obs[["sample", "group", "Position", "Sex", "celltype", "n_genes_by_counts", "pct_counts_mt"]].to_csv(os.path.join(OUT, "GSE222318_cell_metadata.csv")); del a1
rows = []
for ct in s1.celltype.unique():
    d = s1[s1.celltype == ct]; ctrl = d[d.group == "Control"]
    for k in ["C", "I", "R"] + [m for m in MOD if m in d.columns]:
        for name, A, B in [("acute_vs_control", d[d.group == "Acute"], ctrl), ("subacute_vs_control", d[d.group == "Subacute"], ctrl), ("chronic_vs_control", d[d.group == "Chronic"], ctrl),
                           ("acute_vs_nonacute", d[d.group == "Acute"], d[d.group.isin(["Subacute", "Chronic"])]), ("acute_vs_subacute_samePosition", d[d.group == "Acute"], d[d.group == "Subacute"])]:
            x, y = A[k].dropna().values, B[k].dropna().values
            if len(x) < 2 or len(y) < 2: continue
            g, lo, hi = hedges_g(x, y); obs, p, na = exact_perm_test(x, y) if len(x) + len(y) <= 12 else (x.mean() - y.mean(), np.nan, np.nan)
            rows.append(dict(dataset="GSE222318", celltype=ct, measure=k, contrast=name, n_A=len(x), n_B=len(y), positions_A="/".join(sorted(A.Position.unique())), positions_B="/".join(sorted(B.Position.unique())), mean_diff=obs, g=g, ci_low=lo, ci_high=hi, p_exact_descriptive=p, n_assignments=na))
E1 = pd.DataFrame(rows); E1.to_csv(os.path.join(OUT, "GSE222318_stage_estimates.csv"), index=False)
pd.set_option("display.width", 250)
sub = E1[(E1.celltype == "SMC_contractile") & (E1.measure.isin(["C", "I", "R"]))].pivot_table(index="measure", columns="contrast", values="g").round(2); print("\nGSE222318 contractile SMC, Hedges' g (estimation only):\n", sub.to_string())
print(s1[s1.celltype == "SMC_contractile"][["sample", "group", "Position", "n_cells", "C", "I", "R"]].round(2).to_string(index=False))

# ---------------------------------------------------------------- GSE189795 (acute vs control; primary R in contractile SMC)
a2 = read_dense_csv_subsampled(os.path.join(SC, "GSE189795_AllSample.NormalizedCounts.txt.gz"), lambda c: c.split("_")[0], sep="\t")
a2.obs["group"] = np.where(a2.obs["sample"].str.startswith("AD"), "Acute", "Control"); print("GSE189795 loaded:", a2.shape, flush=True)
a2 = annotate(a2, normalise=False); print("GSE189795 annotated:", a2.shape, a2.obs.celltype.value_counts().to_dict(), flush=True)
s2, idx2 = stage_scores(a2, "GSE189795"); s2.to_csv(os.path.join(OUT, "GSE189795_patient_celltype_scores.csv"), index=False); idx2.to_csv(os.path.join(OUT, "GSE189795_pseudobulk_index.csv"))
a2.obs[["sample", "group", "celltype", "n_genes_by_counts"]].to_csv(os.path.join(OUT, "GSE189795_cell_metadata.csv"))
rows = []
for ct in s2.celltype.unique():
    d = s2[s2.celltype == ct]
    for k in ["C", "I", "R"] + [m for m in MOD if m in d.columns]:
        x, y = d.loc[d.group == "Acute", k].dropna().values, d.loc[d.group == "Control", k].dropna().values
        if len(x) < 2 or len(y) < 2: continue
        g, lo, hi = hedges_g(x, y); obs, p, na = exact_perm_test(x, y)
        rows.append(dict(dataset="GSE189795", celltype=ct, measure=k, n_acute=len(x), n_ctrl=len(y), mean_diff=obs, g=g, ci_low=lo, ci_high=hi, p_exact=p, n_assignments=na, primary=(ct == "SMC_contractile" and k == "R")))
E2 = pd.DataFrame(rows); m = ~E2.primary; E2["p_bh_exploratory"] = np.nan; E2.loc[m, "p_bh_exploratory"] = multipletests(E2.loc[m, "p_exact"], method="fdr_bh")[1]; E2.to_csv(os.path.join(OUT, "GSE189795_tests.csv"), index=False)
print("\nGSE189795 PRIMARY (R in contractile SMC, exact):"); print(E2[E2.primary][["n_acute", "n_ctrl", "mean_diff", "g", "ci_low", "ci_high", "p_exact", "n_assignments"]].round(3).to_string(index=False))
print(E2[(E2.celltype == "SMC_contractile") & (E2.measure.isin(["C", "I", "R"] + INJURY_SIX))][["measure", "mean_diff", "g", "ci_low", "ci_high", "p_exact", "p_bh_exploratory"]].round(3).to_string(index=False))
print("done", flush=True)
