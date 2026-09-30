"""Cell-type-resolved single-cell validation of the bulk-derived gene classes.

1. Cellular source of each gene class: mean log-normalised expression of class genes per cell type (control cells), specificity.
2. Cell-type-level pseudo-bulk DE (disease vs control within GSE155468 [ATAA] and GSE213740 [TAAD]) per cell type:
   limma-trend on log-CPM of summed counts per sample x cell type. Concordance of TAAD-specific / ATAA-specific / shared
   gene directions with bulk meta-analysis, within each cell type.
3. Heatmap of top class genes: log2FC per cell type in ATAA (GSE155468) and TAAD (GSE213740).
"""
import os, json, warnings, numpy as np, pandas as pd, scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
META = os.path.join(ROOT, "results", "meta"); SCD = os.path.join(ROOT, "results", "sc"); OUT = os.path.join(ROOT, "results", "sc_validation"); os.makedirs(OUT, exist_ok=True); FIG = os.path.join(ROOT, "figures")
cmp = pd.read_csv(os.path.join(META, "ATAA_vs_TAAD_gene_comparison.csv"), index_col=0)
order = ["SMC_contractile", "SMC_modulated", "Fibromyocyte", "Fibroblast", "Endothelial", "Macrophage", "Inflammatory_myeloid", "T_cell", "NK", "B_cell", "Plasma", "Mast"]
pbm = pd.read_csv(os.path.join(SCD, "pseudobulk_celltype_sample.csv"), index_col=0)
pbc = pd.read_csv(os.path.join(SCD, "pseudobulk_counts_celltype_sample.csv"), index_col=0)
cm = pd.read_csv(os.path.join(SCD, "cell_metadata.csv"), index_col=0).drop_duplicates("sample").set_index("sample")[["dataset", "group"]]
cols = pd.DataFrame([c.split("|") for c in pbm.columns], columns=["sample", "celltype"], index=pbm.columns).join(cm, on="sample")

# ---------------------------------------------------------------- 1. cellular source of gene classes (control samples, mean of per-sample cell-type means)
ctrl_cols = cols.index[cols.group == "Control"]
src = pd.DataFrame({ct: pbm[ctrl_cols[cols.loc[ctrl_cols, "celltype"] == ct]].mean(axis=1) for ct in order})
src_rel = src.div(src.sum(axis=1) + 1e-9, axis=0)  # fraction of expression per cell type
rows = []
for cl in ["TAAD_specific", "ATAA_specific", "shared_concordant", "discordant"]:
    s = cmp[cmp["class"] == cl]; ref = s.lfc_TAAD if "TAAD" in cl else s.lfc_ATAA
    for direction, g in [("up", s.index[ref > 0]), ("down", s.index[ref < 0])]:
        g = g.intersection(src.index)
        if len(g) < 5: continue
        rel = src_rel.loc[g].mean(); rel.name = f"{cl}_{direction}"; rows.append(rel)
        top_ct = src.loc[g].idxmax(axis=1).value_counts(normalize=True)
        print(f"{cl} {direction} (n={len(g)}): dominant cell type of expression -> " + ", ".join(f"{k} {v*100:.0f}%" for k, v in top_ct.head(4).items()))
source = pd.DataFrame(rows); source.to_csv(os.path.join(OUT, "gene_class_cellular_source.csv"))

# ---------------------------------------------------------------- 2. cell-type pseudo-bulk DE
from inmoose.limma import lmFit, eBayes
def limma(expr, case_mask):
    design = pd.DataFrame({"i": 1.0, "case": case_mask.astype(float)}, index=expr.columns)
    fit = eBayes(lmFit(expr, design), trend=True)
    c = fit.coefficients.columns[1]
    return pd.DataFrame({"lfc": fit.coefficients[c], "t": fit.t[c], "p": fit.p_value[c]})
de = {}
for ds, dis in [("GSE155468", "ATAA"), ("GSE213740", "TAAD")]:
    for ct in order:
        cc = cols.index[(cols.dataset == ds) & (cols.celltype == ct)]
        grp = cols.loc[cc, "group"]
        if (grp == dis).sum() < 3 or (grp == "Control").sum() < 2: continue
        counts = pbc[cc]; counts = counts[(counts >= 5).sum(axis=1) >= 2]
        lcpm = np.log2(counts / counts.sum(axis=0) * 1e6 + 1)
        d = limma(lcpm, (grp == dis).values); d["dataset"] = ds; d["disease"] = dis; d["celltype"] = ct
        de[(dis, ct)] = d
        print(ds, ct, f"n={len(cc)} genes={len(d)} nominal p<0.05: {(d.p<0.05).sum()}", flush=True)
big = pd.concat([d.assign(gene=d.index) for d in de.values()]); big.to_csv(os.path.join(OUT, "celltype_pseudobulk_DE.csv"), index=False)

# concordance of bulk gene classes within cell types
rows = []
for cl in ["TAAD_specific", "ATAA_specific", "shared_concordant"]:
    s = cmp[cmp["class"] == cl]; ref = (s.lfc_TAAD if "TAAD" in cl else s.lfc_ATAA)
    for (dis, ct), d in de.items():
        g = s.index.intersection(d.index)
        # restrict to genes expressed in this cell type (top 60% by mean pseudo-bulk expression)
        if len(g) < 15: continue
        sign_ref = np.sign(ref[g]); conc = (np.sign(d.loc[g, "lfc"]) == sign_ref).mean() * 100
        rho = stats.spearmanr(ref[g], d.loc[g, "lfc"])[0]
        signed = (d.loc[g, "lfc"] * sign_ref); w = stats.wilcoxon(signed).pvalue if len(g) > 10 else np.nan
        rows.append(dict(gene_class=cl, disease_dataset=dis, celltype=ct, n_genes=len(g), direction_concordance=conc, rho=rho, median_signed_lfc=signed.median(), wilcoxon_p=w))
conc = pd.DataFrame(rows); conc.to_csv(os.path.join(OUT, "class_concordance_by_celltype.csv"), index=False)
pd.set_option("display.width", 220); print(conc.round(3).to_string(index=False))

# ---------------------------------------------------------------- 3. figure: heatmap of top genes (cell type x dataset)
fig, axes = plt.subplots(1, 3, figsize=(19, 7), gridspec_kw=dict(width_ratios=[1.2, 1.2, 1]))
for ax, cl in zip(axes[:2], ["TAAD_specific", "ATAA_specific"]):
    s = cmp[cmp["class"] == cl].copy(); ref = s.lfc_TAAD if "TAAD" in cl else s.lfc_ATAA
    s["absz"] = s.z_TAAD.abs() if "TAAD" in cl else s.z_ATAA.abs()
    top = pd.concat([s[ref > 0].sort_values("absz", ascending=False).head(15), s[ref < 0].sort_values("absz", ascending=False).head(10)]).index
    M = pd.DataFrame({f"{dis}|{ct}": de[(dis, ct)]["lfc"].reindex(top) for (dis, ct) in de if ct in order[:8]})
    M = M[[c for c in [f"{dis}|{ct}" for dis in ["ATAA", "TAAD"] for ct in order[:8]] if c in M.columns]]
    sns.heatmap(M, cmap="RdBu_r", center=0, vmin=-3, vmax=3, ax=ax, cbar_kws=dict(label="log2FC disease vs control (pseudo-bulk)", shrink=0.5), linewidths=0.2, linecolor="white")
    ax.set_xticklabels([c.replace("|", "\n") for c in M.columns], rotation=90, fontsize=7); ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)
    ax.axvline(sum(c.startswith("ATAA") for c in M.columns), color="k", lw=1.5); ax.set_title(f"{cl.replace('_','-')} genes (bulk meta): cell-type-resolved change\nleft: ATAA (GSE155468) | right: TAAD (GSE213740)", fontsize=9)
S = source.loc[[i for i in source.index if not i.startswith("discordant")], order]
sns.heatmap(S, cmap="viridis", ax=axes[2], cbar_kws=dict(label="fraction of class expression per cell type", shrink=0.5), annot=True, fmt=".2f", annot_kws=dict(size=6))
axes[2].set_xticklabels(order, rotation=90, fontsize=7); axes[2].set_yticklabels(S.index, fontsize=7, rotation=0); axes[2].set_title("Cellular source of gene classes (control aorta)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig8_sc_validation.png"), dpi=200); fig.savefig(os.path.join(FIG, "Fig8_sc_validation.pdf")); plt.close(fig)
