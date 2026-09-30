"""Pathway-level comparison of ATAA and TAAD.

1. GSEA (gseapy prerank, 1000 permutations) on each cohort ranked by the DE statistic; gene sets: MSigDB Hallmark, KEGG (legacy), Reactome.
2. Disease-level GSEA on the Stouffer meta-z ranking (ATAA, TAAD), on the difference statistic z_diff (TAAD - ATAA) and on the direct contrast (GSE318877).
3. NES heatmap across cohorts; scatter of disease-level NES (ATAA vs TAAD).
4. ORA (hypergeometric, tested-gene background) of TAAD-specific / ATAA-specific / shared DEGs on GO-BP, KEGG, Reactome, TRRUST TFs.
"""
import os, json, warnings, numpy as np, pandas as pd, gseapy as gp
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DE = os.path.join(ROOT, "results", "de"); META = os.path.join(ROOT, "results", "meta"); GS = os.path.join(ROOT, "data", "genesets")
OUT = os.path.join(ROOT, "results", "pathways"); os.makedirs(OUT, exist_ok=True); FIG = os.path.join(ROOT, "figures")
REG = json.load(open(os.path.join(ROOT, "results", "processed", "cohorts.json"))); DESIGN = json.load(open(os.path.join(META, "design.json")))
GMT = {"Hallmark": os.path.join(GS, "h.all.v2024.1.Hs.symbols.gmt"), "KEGG": os.path.join(GS, "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"),
       "Reactome": os.path.join(GS, "c2.cp.reactome.v2024.1.Hs.symbols.gmt")}

def prerank(rnk, tag):
    res = []
    for lib, path in GMT.items():
        pr = gp.prerank(rnk=rnk, gene_sets=path, permutation_num=1000, min_size=15, max_size=500, seed=7, threads=4, outdir=None, verbose=False)
        r = pr.res2d.copy(); r["library"] = lib; r["ranking"] = tag; res.append(r)
    r = pd.concat(res); r = r.rename(columns={"Term": "term", "NES": "nes", "FDR q-val": "fdr", "NOM p-val": "pval"})
    r["nes"] = r.nes.astype(float); r["fdr"] = r.fdr.astype(float)
    r.to_csv(os.path.join(OUT, f"gsea_{tag}.csv"), index=False)
    return r[["term", "library", "nes", "fdr", "pval", "Lead_genes", "ranking"]]

def rank_from_de(name):
    d = pd.read_csv(os.path.join(DE, f"{name}_de.csv"), index_col=0).dropna(subset=["stat"])
    d = d[~d.index.duplicated()]
    return d["stat"].sort_values(ascending=False)

allres = {}
for name in REG:
    if os.path.exists(os.path.join(OUT, f"gsea_{name}.csv")):
        r = pd.read_csv(os.path.join(OUT, f"gsea_{name}.csv")); r = r.rename(columns={"Term": "term", "NES": "nes", "FDR q-val": "fdr", "NOM p-val": "pval"})
    else:
        r = prerank(rank_from_de(name), name)
    allres[name] = r; print(name, "GSEA done:", (r.fdr < 0.25).sum(), "sets FDR<0.25", flush=True)

# disease-level rankings
for dis in ["ATAA", "TAAD"]:
    m = pd.read_csv(os.path.join(META, f"{dis}_meta.csv"), index_col=0)
    allres[f"meta_{dis}"] = prerank(m["z_meta"].sort_values(ascending=False), f"meta_{dis}")
cmp = pd.read_csv(os.path.join(META, "ATAA_vs_TAAD_gene_comparison.csv"), index_col=0)
allres["meta_diff_TAADminusATAA"] = prerank(cmp["z_diff"].dropna().sort_values(ascending=False), "meta_diff_TAADminusATAA")
print("meta-level GSEA done", flush=True)

# ---------------------------------------------------------------- NES matrices
big = pd.concat(allres.values())
nes = big.pivot_table(index=["library", "term"], columns="ranking", values="nes")
fdr = big.pivot_table(index=["library", "term"], columns="ranking", values="fdr")
nes.to_csv(os.path.join(OUT, "NES_matrix.csv")); fdr.to_csv(os.path.join(OUT, "FDR_matrix.csv"))
ataa_cols = DESIGN["ATAA"]["cohorts"]; taad_cols = DESIGN["TAAD"]["cohorts"]
summary = pd.DataFrame({"NES_meta_ATAA": nes["meta_ATAA"], "FDR_meta_ATAA": fdr["meta_ATAA"], "NES_meta_TAAD": nes["meta_TAAD"], "FDR_meta_TAAD": fdr["meta_TAAD"],
                        "NES_diff": nes["meta_diff_TAADminusATAA"], "FDR_diff": fdr["meta_diff_TAADminusATAA"], "NES_direct_GSE318877": nes["GSE318877"], "FDR_direct_GSE318877": fdr["GSE318877"],
                        "mean_NES_ATAA_cohorts": nes[ataa_cols].mean(axis=1), "mean_NES_TAAD_cohorts": nes[taad_cols].mean(axis=1),
                        "n_ATAA_cohorts_FDR25_sameDir": ((fdr[ataa_cols] < 0.25) & (np.sign(nes[ataa_cols]).eq(np.sign(nes["meta_ATAA"]), axis=0))).sum(axis=1),
                        "n_TAAD_cohorts_FDR25_sameDir": ((fdr[taad_cols] < 0.25) & (np.sign(nes[taad_cols]).eq(np.sign(nes["meta_TAAD"]), axis=0))).sum(axis=1)})
def pclass(r):
    a = r.FDR_meta_ATAA < 0.05; t = r.FDR_meta_TAAD < 0.05
    if a and t: return "shared" if np.sign(r.NES_meta_ATAA) == np.sign(r.NES_meta_TAAD) else "opposite"
    if t and not a: return "TAAD_specific" if r.FDR_diff < 0.05 else "TAAD_only"
    if a and not t: return "ATAA_specific" if r.FDR_diff < 0.05 else "ATAA_only"
    return "ns"
summary["class"] = summary.apply(pclass, axis=1)
summary = summary.sort_values("NES_diff"); summary.to_csv(os.path.join(OUT, "pathway_summary.csv"))
print(summary["class"].value_counts().to_string())
from scipy import stats
ok = summary.dropna(subset=["NES_diff", "NES_direct_GSE318877"])
print("NES_diff (meta) vs direct GSE318877 NES: rho=%.3f p=%.1e (n=%d)" % (*stats.spearmanr(ok.NES_diff, ok.NES_direct_GSE318877)[:2], len(ok)))
hm = summary.loc["Hallmark"].dropna(subset=["NES_diff", "NES_direct_GSE318877"])
print("Hallmark only: rho=%.3f p=%.1e (n=%d)" % (*stats.spearmanr(hm.NES_diff, hm.NES_direct_GSE318877)[:2], len(hm)))

# ---------------------------------------------------------------- Figure: Hallmark NES heatmap across cohorts + disease-level scatter
H = nes.loc["Hallmark"]; Hf = fdr.loc["Hallmark"]
order_cols = ataa_cols + taad_cols + ["meta_ATAA", "meta_TAAD", "meta_diff_TAADminusATAA", "GSE318877"]
H = H[order_cols]; Hf = Hf[order_cols]
H.index = [i.replace("HALLMARK_", "").replace("_", " ").title() for i in H.index]; Hf.index = H.index
H = H.loc[H["meta_diff_TAADminusATAA"].sort_values().index]
fig, ax = plt.subplots(figsize=(11, 13))
sns.heatmap(H, cmap="RdBu_r", center=0, vmin=-3, vmax=3, ax=ax, cbar_kws=dict(label="NES", shrink=0.4), linewidths=0.3, linecolor="white")
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        if Hf.iloc[i, j] < 0.05: ax.text(j + 0.5, i + 0.5, "*", ha="center", va="center", fontsize=8, color="k")
labels = [f"{c}\n({REG[c]['disease']}, n={sum(REG[c]['groups'].values())})" if c in REG else c.replace("meta_", "meta ").replace("TAADminusATAA", "TAAD-ATAA") for c in H.columns]
ax.set_xticklabels(labels, rotation=90, fontsize=7); ax.set_yticklabels(ax.get_yticklabels(), fontsize=7.5)
ax.axvline(len(ataa_cols), color="k", lw=1.5); ax.axvline(len(ataa_cols) + len(taad_cols), color="k", lw=1.5); ax.axvline(len(ataa_cols) + len(taad_cols) + 2, color="k", lw=1)
ax.set_title("Hallmark GSEA: ATAA cohorts | TAAD cohorts | meta-rankings | direct TAAD-vs-ATAA (GSE318877)   (* FDR<0.05)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig3_hallmark_NES_heatmap.png"), dpi=200); fig.savefig(os.path.join(FIG, "Fig3_hallmark_NES_heatmap.pdf")); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(12, 5.2))
s = summary.loc["Hallmark"]
ax[0].scatter(s.NES_meta_ATAA, s.NES_meta_TAAD, c=np.where(s["class"].isin(["shared"]), "#7B2CBF", np.where(s["class"].str.startswith("TAAD"), "#C44E52", np.where(s["class"].str.startswith("ATAA"), "#2A9D8F", "lightgrey"))), s=30, edgecolor="k", linewidth=0.3)
for t, r in s.iterrows():
    if r["class"] != "ns" and (abs(r.NES_meta_ATAA) > 1.8 or abs(r.NES_meta_TAAD) > 1.8 or abs(r.NES_meta_ATAA - r.NES_meta_TAAD) > 1.5):
        ax[0].annotate(t.replace("HALLMARK_", "").replace("_", " ").title(), (r.NES_meta_ATAA, r.NES_meta_TAAD), fontsize=5.5)
ax[0].axhline(0, c="k", lw=0.5); ax[0].axvline(0, c="k", lw=0.5); ax[0].plot([-3, 3], [-3, 3], "k--", lw=0.5)
ax[0].set_xlabel("NES, ATAA vs control (meta ranking)"); ax[0].set_ylabel("NES, TAAD vs control (meta ranking)"); ax[0].set_title("Hallmark pathways: disease-level GSEA", fontsize=10)
ax[1].scatter(ok.NES_diff, ok.NES_direct_GSE318877, s=8, c="lightgrey", rasterized=True, label="all gene sets")
ax[1].scatter(hm.NES_diff, hm.NES_direct_GSE318877, s=22, c="#C44E52", edgecolor="k", linewidth=0.3, label="Hallmark")
ax[1].axhline(0, c="k", lw=0.5); ax[1].axvline(0, c="k", lw=0.5)
ax[1].set_xlabel("NES on meta difference ranking (TAAD - ATAA)"); ax[1].set_ylabel("NES, direct TAAD vs ATAA (GSE318877)")
ax[1].set_title("Pathway-level validation in the direct-contrast cohort\nrho=%.2f (all), %.2f (Hallmark)" % (stats.spearmanr(ok.NES_diff, ok.NES_direct_GSE318877)[0], stats.spearmanr(hm.NES_diff, hm.NES_direct_GSE318877)[0]), fontsize=10)
ax[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig3b_pathway_scatter.png"), dpi=200); fig.savefig(os.path.join(FIG, "Fig3b_pathway_scatter.pdf")); plt.close(fig)

# ---------------------------------------------------------------- ORA of gene classes
libs = {"GO_Biological_Process_2023": None, "KEGG_2021_Human": None, "Reactome_2022": None, "TRRUST_Transcription_Factors_2019": None, "MSigDB_Hallmark_2020": None}
os.makedirs(os.path.join(OUT, "enrichr_libs"), exist_ok=True)
for lib in libs:
    p = os.path.join(OUT, "enrichr_libs", lib + ".gmt")
    if not os.path.exists(p):
        d = gp.get_library(name=lib, organism="Human")
        with open(p, "w") as fh:
            for k, v in d.items(): fh.write("\t".join([k, "na"] + list(v)) + "\n")
    libs[lib] = p
background = cmp.index.tolist()
ora_all = []
for cl in ["TAAD_specific", "ATAA_specific", "shared_concordant", "TAAD_stronger", "ATAA_stronger"]:
    genes = cmp.index[cmp["class"] == cl].tolist()
    if len(genes) < 10: continue
    for direction in ["up", "down"]:
        ref = cmp.lfc_TAAD if "TAAD" in cl else cmp.lfc_ATAA
        g = [x for x in genes if (ref[x] > 0) == (direction == "up")]
        if len(g) < 10: continue
        for lib, p in libs.items():
            try:
                e = gp.enrich(gene_list=g, gene_sets=p, background=background, outdir=None, verbose=False).results
            except Exception as ex:
                print("ORA failed", cl, direction, lib, ex); continue
            e["class"] = cl; e["direction"] = direction; e["library"] = lib; e["n_genes"] = len(g); ora_all.append(e)
ora = pd.concat(ora_all); ora.to_csv(os.path.join(OUT, "ORA_gene_classes.csv"), index=False)
for (cl, d), s in ora[ora["Adjusted P-value"] < 0.05].groupby(["class", "direction"]):
    print(f"\n== ORA {cl} {d} (n={s.n_genes.iloc[0]}) top terms:")
    print(s.sort_values("Adjusted P-value").head(12)[["library", "Term", "Overlap", "Adjusted P-value"]].to_string(index=False))
