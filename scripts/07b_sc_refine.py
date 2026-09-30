"""Refine single-cell annotation (marker-informed relabelling of Leiden clusters), recompute proportions,
aorta-specific marker sets, pseudo-bulk tables and figures."""
import os, warnings, numpy as np, pandas as pd, scanpy as sc
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "sc"); FIG = os.path.join(ROOT, "figures"); os.makedirs(os.path.join(FIG, "sc"), exist_ok=True)
a = sc.read_h5ad(os.path.join(OUT, "aorta_sc_integrated.h5ad"))
LABELS = {"0": "SMC_contractile", "1": "SMC_modulated", "2": "SMC_modulated", "3": "Fibromyocyte", "8": "Fibroblast", "16": "Fibroblast",
          "13": "Endothelial", "21": "Endothelial", "22": "Endothelial", "20": "Endothelial", "4": "Macrophage", "5": "Macrophage", "7": "Macrophage", "12": "Macrophage",
          "6": "Inflammatory_myeloid", "9": "T_cell", "10": "T_cell", "11": "T_cell", "15": "NK", "19": "B_cell", "14": "Plasma", "18": "Mast", "17": "Doublet_like"}
a.obs["celltype"] = a.obs["leiden"].astype(str).map(LABELS).astype("category")
a = a[a.obs.celltype != "Doublet_like"].copy()
a.obs["lineage"] = a.obs.celltype.astype(str).replace({"SMC_contractile": "SMC", "SMC_modulated": "SMC"})
print(a.obs.celltype.value_counts().to_string())
order = ["SMC_contractile", "SMC_modulated", "Fibromyocyte", "Fibroblast", "Endothelial", "Macrophage", "Inflammatory_myeloid", "T_cell", "NK", "B_cell", "Plasma", "Mast"]

# ---------------------------------------------------------------- proportions per sample and group tests
prop = pd.crosstab(a.obs["sample"], a.obs["celltype"], normalize="index")[order]
meta = a.obs.drop_duplicates("sample").set_index("sample")[["dataset", "group"]]
prop = meta.join(prop); prop.to_csv(os.path.join(OUT, "celltype_proportions.csv"))
rows = []
for ds, dis in [("GSE155468", "ATAA"), ("GSE213740", "TAAD")]:
    p = prop[prop.dataset == ds]
    for ct in order:
        d, c = p.loc[p.group == dis, ct], p.loc[p.group == "Control", ct]
        rows.append(dict(dataset=ds, disease=dis, celltype=ct, mean_disease=d.mean(), mean_control=c.mean(), log2_ratio=np.log2((d.mean() + 1e-3) / (c.mean() + 1e-3)), p_mannwhitney=stats.mannwhitneyu(d, c).pvalue))
pt = pd.DataFrame(rows); pt.to_csv(os.path.join(OUT, "celltype_proportion_tests.csv"), index=False)
print(pt.round(3).to_string(index=False))

# ---------------------------------------------------------------- marker sets (one-vs-rest, all cells; specificity filters)
sc.tl.rank_genes_groups(a, "celltype", method="wilcoxon", pts=True)
marker_sets = {}
for ct in order:
    df = sc.get.rank_genes_groups_df(a, ct)
    df = df[(df.logfoldchanges > 1.5) & (df.pvals_adj < 0.01) & (df.pct_nz_group > 0.3) & (df.pct_nz_reference < 0.25)]
    df = df[~df.names.str.startswith(("MT-", "RPL", "RPS", "MALAT1"))].head(50)
    marker_sets[ct] = df.names.tolist(); df.to_csv(os.path.join(OUT, f"markers_{ct}.csv"), index=False)
    print(ct, len(df), ", ".join(df.names[:8]))
with open(os.path.join(OUT, "aorta_celltype_markers.gmt"), "w") as fh:
    for ct, g in marker_sets.items():
        if len(g) >= 8: fh.write("\t".join([ct, "scRNA_aorta_GSE155468_GSE213740"] + g) + "\n")

# ---------------------------------------------------------------- pseudo-bulk tables
X = a.X.tocsr(); C = a.layers["counts"].tocsr(); rows_m = []; rows_c = []
for (s, ct), idx in a.obs.groupby(["sample", "celltype"]).indices.items():
    if len(idx) < 10: continue
    rows_m.append(pd.Series(np.asarray(X[idx].mean(axis=0)).ravel(), index=a.var_names, name=f"{s}|{ct}"))
    rows_c.append(pd.Series(np.asarray(C[idx].sum(axis=0)).ravel(), index=a.var_names, name=f"{s}|{ct}"))
pd.DataFrame(rows_m).T.to_csv(os.path.join(OUT, "pseudobulk_celltype_sample.csv")); pd.DataFrame(rows_c).T.to_csv(os.path.join(OUT, "pseudobulk_counts_celltype_sample.csv"))
a.obs.to_csv(os.path.join(OUT, "cell_metadata.csv"))
a.write_h5ad(os.path.join(OUT, "aorta_sc_integrated.h5ad"), compression="gzip")

# ---------------------------------------------------------------- figures: UMAP, dotplot, proportions
pal = dict(zip(order, sns.color_palette("tab20", len(order))))
fig, ax = plt.subplots(1, 3, figsize=(18, 5.6))
um = a.obsm["X_umap"]
for ct in order:
    idx = (a.obs.celltype == ct).values; ax[0].scatter(um[idx, 0], um[idx, 1], s=0.5, c=[pal[ct]], label=ct, rasterized=True)
    ax[0].annotate(ct, um[idx].mean(axis=0), fontsize=6.5, ha="center", weight="bold")
ax[0].set_title(f"Integrated ascending-aorta scRNA-seq (n={a.n_obs:,} cells; GSE155468 + GSE213740)", fontsize=9); ax[0].axis("off")
for i, (ds, dis) in enumerate([("GSE155468", "ATAA"), ("GSE213740", "TAAD")]):
    p = prop[prop.dataset == ds]; mp = p.groupby("group")[order].mean().loc[["Control", dis]]
    mp.plot(kind="bar", stacked=True, color=[pal[c] for c in order], ax=ax[i + 1], width=0.6, legend=(i == 1))
    ax[i + 1].set_ylabel("fraction of cells"); ax[i + 1].set_title(f"{ds}: {dis} (n={int((p.group==dis).sum())}) vs control (n={int((p.group=='Control').sum())})", fontsize=9); ax[i + 1].tick_params(axis="x", rotation=0)
    if i == 1: ax[i + 1].legend(fontsize=6, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig7_scRNA_overview.png"), dpi=200); fig.savefig(os.path.join(FIG, "Fig7_scRNA_overview.pdf")); plt.close(fig)
sc.settings.figdir = os.path.join(FIG, "sc")
canon = ["MYH11", "ACTA2", "CNN1", "LMOD1", "MGP", "SPARC", "TNFRSF11B", "COL1A1", "BGN", "LTBP2", "DCN", "LUM", "PDGFRA", "PECAM1", "VWF", "CD68", "C1QA", "CD163", "S100A8", "S100A9", "FCGR3B", "CD3E", "IL7R", "NKG7", "GNLY", "MS4A1", "CD79A", "MZB1", "JCHAIN", "TPSAB1", "KIT"]
canon = [g for g in canon if g in a.var_names]
sc.pl.dotplot(a, canon, groupby="celltype", categories_order=order, save="_canonical_markers.png", show=False, standard_scale="var")
print("done", a.shape)
