"""Per-cohort differential expression + QC figures.

RNA-seq count cohorts  -> pyDESeq2 (Wald test, design ~ group)
Array / FPKM cohorts   -> limma (inmoose port): lmFit + eBayes(trend=True)
Outputs results/de/<cohort>_de.csv with unified columns:
  gene, log2FC, se, stat, pvalue, padj, AveExpr, method
"""
import os, json, sys, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "results", "processed")
DE = os.path.join(ROOT, "results", "de"); os.makedirs(DE, exist_ok=True)
FIG = os.path.join(ROOT, "figures", "qc"); os.makedirs(FIG, exist_ok=True)
REG = json.load(open(os.path.join(PROC, "cohorts.json")))

def limma_de(expr, meta, case, ctrl, trend=True, min_expr=None):
    from inmoose.limma import lmFit, eBayes, topTable
    keep = meta["group"].isin([case, ctrl])
    expr = expr.loc[:, keep[keep].index]; meta = meta.loc[expr.columns]
    if min_expr is not None:  # drop genes never expressed (FPKM-derived matrices)
        expr = expr[(expr > min_expr).mean(axis=1) >= 0.25]
    expr = expr.dropna()
    design = pd.DataFrame({"Intercept": 1.0, "case": (meta["group"] == case).astype(float).values}, index=expr.columns)
    fit = lmFit(expr, design)
    fit = eBayes(fit, trend=trend, robust=False)
    from statsmodels.stats.multitest import multipletests
    coef_col = fit.coefficients.columns[1]  # inmoose names design columns column0/column1
    pv = fit.p_value[coef_col].values
    out = pd.DataFrame({"gene": fit.coefficients.index, "log2FC": fit.coefficients[coef_col].values,
                        "stat": fit.t[coef_col].values, "pvalue": pv,
                        "padj": multipletests(pv, method="fdr_bh")[1], "AveExpr": fit.Amean.values})
    out["se"] = out["log2FC"] / out["stat"]
    out["method"] = "limma-trend" if trend else "limma"
    return out.set_index("gene")

def deseq_de(counts, meta, case, ctrl):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    keep = meta["group"].isin([case, ctrl])
    counts = counts.loc[:, keep[keep].index]; meta = meta.loc[counts.columns].copy()
    n_min = meta["group"].value_counts().min()
    counts = counts[(counts >= 10).sum(axis=1) >= n_min]
    counts = counts.round().astype(int)
    meta["group"] = meta["group"].astype(str)
    dds = DeseqDataSet(counts=counts.T, metadata=meta[["group"]], design="~group", quiet=True, n_cpus=4)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["group", case, ctrl], quiet=True, n_cpus=4)
    st.summary()
    r = st.results_df
    out = pd.DataFrame({"log2FC": r["log2FoldChange"], "se": r["lfcSE"], "stat": r["stat"], "pvalue": r["pvalue"],
                        "padj": r["padj"], "AveExpr": np.log2(r["baseMean"] + 1)})
    out.index.name = "gene"; out["method"] = "DESeq2"
    return out

def qc_plots(name, expr, meta, de, case, ctrl):
    keep = meta["group"].isin([case, ctrl]); x = expr.loc[:, keep[keep].index].dropna()
    x = x.loc[x.var(axis=1).sort_values(ascending=False).index[:2000]]
    z = ((x.T - x.T.mean()) / (x.T.std() + 1e-9))
    pcs = PCA(n_components=2).fit(z); p = pcs.transform(z)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    cols = {ctrl: "#4C72B0", case: "#C44E52"}
    for g in [ctrl, case]:
        idx = (meta.loc[x.columns, "group"] == g).values
        ax[0].scatter(p[idx, 0], p[idx, 1], c=cols[g], label=f"{g} (n={idx.sum()})", s=40, edgecolor="k", linewidth=0.4)
    ax[0].set_xlabel(f"PC1 ({pcs.explained_variance_ratio_[0]*100:.1f}%)"); ax[0].set_ylabel(f"PC2 ({pcs.explained_variance_ratio_[1]*100:.1f}%)")
    ax[0].set_title(f"{name}: PCA (top 2000 variable genes)", fontsize=10); ax[0].legend(frameon=False, fontsize=8)
    d = de.dropna(subset=["padj"])
    sig = (d["padj"] < 0.05) & (d["log2FC"].abs() > 1)
    ax[1].scatter(d["log2FC"], -np.log10(d["pvalue"].clip(1e-300)), s=6, c="lightgrey", rasterized=True)
    ax[1].scatter(d.loc[sig & (d.log2FC > 0), "log2FC"], -np.log10(d.loc[sig & (d.log2FC > 0), "pvalue"].clip(1e-300)), s=6, c="#C44E52", rasterized=True)
    ax[1].scatter(d.loc[sig & (d.log2FC < 0), "log2FC"], -np.log10(d.loc[sig & (d.log2FC < 0), "pvalue"].clip(1e-300)), s=6, c="#4C72B0", rasterized=True)
    top = d[sig].assign(score=lambda t: -np.log10(t.pvalue.clip(1e-300)) * t.log2FC.abs()).sort_values("score", ascending=False).head(12)
    for g, r in top.iterrows(): ax[1].annotate(g, (r.log2FC, -np.log10(max(r.pvalue, 1e-300))), fontsize=6)
    ax[1].set_xlabel(f"log2FC ({case} vs {ctrl})"); ax[1].set_ylabel("-log10 P"); ax[1].set_title(f"{name}: {sig.sum()} DEGs (FDR<0.05, |log2FC|>1)", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, f"{name}_qc.png"), dpi=150); plt.close(fig)

def run(name):
    info = REG[name]
    expr = pd.read_csv(os.path.join(PROC, f"{name}_expr.csv"), index_col=0)
    meta = pd.read_csv(os.path.join(PROC, f"{name}_meta.csv"), index_col=0)
    meta.index = meta.index.astype(str); expr.columns = expr.columns.astype(str)
    if info["disease"] == "ATAA_vs_TAAD": case, ctrl = "TAAD", "ATAA"
    else: case, ctrl = info["disease"], "Control"
    if info.get("counts"):
        counts = pd.read_csv(os.path.join(PROC, f"{name}_counts.csv"), index_col=0); counts.columns = counts.columns.astype(str)
        de = deseq_de(counts, meta, case, ctrl)
    else:
        fpkm_like = "FPKM" in info["platform"] or "RNA-seq" in info["platform"]
        de = limma_de(expr, meta, case, ctrl, trend=True, min_expr=1.0 if fpkm_like else None)
    de["cohort"] = name; de["contrast"] = f"{case}_vs_{ctrl}"
    de.to_csv(os.path.join(DE, f"{name}_de.csv"))
    qc_plots(name, expr, meta, de, case, ctrl)
    n_sig = int(((de.padj < 0.05) & (de.log2FC.abs() > 1)).sum()); n_sig05 = int((de.padj < 0.05).sum())
    print(f"[{name}] {case} vs {ctrl}: genes tested={de.shape[0]}, FDR<0.05: {n_sig05}, FDR<0.05 & |LFC|>1: {n_sig}, method={de.method.iloc[0]}", flush=True)
    return dict(cohort=name, contrast=f"{case}_vs_{ctrl}", n_tested=de.shape[0], n_fdr05=n_sig05, n_fdr05_lfc1=n_sig, method=de.method.iloc[0])

if __name__ == "__main__":
    names = sys.argv[1:] or list(REG.keys())
    summ = [run(n) for n in names]
    pd.DataFrame(summ).to_csv(os.path.join(DE, "de_summary.csv"), index=False)
