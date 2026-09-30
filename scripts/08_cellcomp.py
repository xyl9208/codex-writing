"""Cell-composition inference in bulk cohorts by ssGSEA with aorta-specific scRNA-derived marker sets.

For every bulk cohort: ssGSEA enrichment score per sample and cell type (gseapy.ssgsea) -> z-score within cohort ->
Hedges' g (disease vs control) per cell type -> random-effects meta per disease -> ATAA vs TAAD comparison.
Also applied to GSE318877 (TAAD vs ATAA direct) and to the scRNA pseudo-bulk (sanity check).
"""
import os, json, warnings, numpy as np, pandas as pd, gseapy as gp
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "results", "processed"); META = os.path.join(ROOT, "results", "meta"); SCD = os.path.join(ROOT, "results", "sc")
OUT = os.path.join(ROOT, "results", "cellcomp"); os.makedirs(OUT, exist_ok=True); FIG = os.path.join(ROOT, "figures")
REG = json.load(open(os.path.join(PROC, "cohorts.json"))); DESIGN = json.load(open(os.path.join(META, "design.json")))
GMT = os.path.join(SCD, "aorta_celltype_markers.gmt")
CTS = [l.split("\t")[0] for l in open(GMT)]

def ssgsea(expr):
    expr = expr[~expr.index.duplicated()].dropna()
    r = gp.ssgsea(data=expr, gene_sets=GMT, sample_norm_method="rank", outdir=None, min_size=5, threads=4, verbose=False)
    m = r.res2d.pivot(index="Term", columns="Name", values="NES").astype(float)
    return m.reindex(CTS)

def hedges_g(x, y):
    nx_, ny = len(x), len(y); sp = np.sqrt(((nx_ - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / (nx_ + ny - 2))
    g = (np.mean(x) - np.mean(y)) / sp * (1 - 3 / (4 * (nx_ + ny) - 9))
    se = np.sqrt((nx_ + ny) / (nx_ * ny) + g ** 2 / (2 * (nx_ + ny)))
    return g, se

rows = []; allscores = []
cohorts = {c: "ATAA" for c in DESIGN["ATAA"]["cohorts"]}; cohorts.update({c: "TAAD" for c in DESIGN["TAAD"]["cohorts"]}); cohorts["GSE318877"] = "direct"
for c, dis in cohorts.items():
    e = pd.read_csv(os.path.join(PROC, f"{c}_expr.csv"), index_col=0); m = pd.read_csv(os.path.join(PROC, f"{c}_meta.csv"), index_col=0)
    e.columns = e.columns.astype(str); m.index = m.index.astype(str)
    s = ssgsea(e); s = (s.sub(s.mean(axis=1), axis=0)).div(s.std(axis=1), axis=0)  # z within cohort
    st = s.T.join(m[["group"]]); st["cohort"] = c; allscores.append(st)
    case = "TAAD" if dis in ("TAAD", "direct") else "ATAA"; ctrl = "ATAA" if dis == "direct" else "Control"
    for ct in CTS:
        x, y = st.loc[st.group == case, ct].values, st.loc[st.group == ctrl, ct].values
        g, se = hedges_g(x, y); p = stats.mannwhitneyu(x, y).pvalue
        rows.append(dict(cohort=c, disease=dis, celltype=ct, g=g, se=se, p=p, n_case=len(x), n_ctrl=len(y)))
    print(c, "done", flush=True)
scores = pd.concat(allscores); scores.to_csv(os.path.join(OUT, "ssgsea_celltype_scores_z.csv"))
eff = pd.DataFrame(rows); eff.to_csv(os.path.join(OUT, "celltype_effects_per_cohort.csv"), index=False)

# random-effects meta per disease and cell type
def rem(g, se):
    w = 1 / se ** 2; yfix = np.sum(w * g) / np.sum(w); Q = np.sum(w * (g - yfix) ** 2); df = len(g) - 1
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w); tau2 = max(0, (Q - df) / C) if C > 0 else 0
    wr = 1 / (se ** 2 + tau2); y = np.sum(wr * g) / np.sum(wr); s = np.sqrt(1 / np.sum(wr))
    return y, s, 2 * stats.norm.sf(abs(y / s)), (Q - df) / Q * 100 if Q > 0 else 0
mrows = []
for dis in ["ATAA", "TAAD"]:
    for ct in CTS:
        sub = eff[(eff.disease == dis) & (eff.celltype == ct)]
        y, s, p, i2 = rem(sub.g.values, sub.se.values)
        mrows.append(dict(disease=dis, celltype=ct, g_meta=y, se_meta=s, p_meta=p, I2=max(i2, 0), n_cohorts=len(sub), frac_same_dir=float((np.sign(sub.g) == np.sign(y)).mean())))
mt = pd.DataFrame(mrows)
piv = mt.pivot(index="celltype", columns="disease", values=["g_meta", "se_meta", "p_meta"])
piv.columns = [f"{a}_{b}" for a, b in piv.columns]
piv["diff_TAAD_minus_ATAA"] = piv.g_meta_TAAD - piv.g_meta_ATAA
piv["z_diff"] = piv.diff_TAAD_minus_ATAA / np.sqrt(piv.se_meta_TAAD ** 2 + piv.se_meta_ATAA ** 2); piv["p_diff"] = 2 * stats.norm.sf(piv.z_diff.abs())
direct = eff[eff.disease == "direct"].set_index("celltype")[["g", "p"]].rename(columns={"g": "g_direct_GSE318877", "p": "p_direct_GSE318877"})
piv = piv.join(direct).loc[CTS]; piv.to_csv(os.path.join(OUT, "celltype_meta_summary.csv"))
print(piv.round(3).to_string())
print("Spearman rho between meta difference (TAAD-ATAA) and direct-cohort g:", stats.spearmanr(piv.diff_TAAD_minus_ATAA, piv.g_direct_GSE318877))

# scRNA sanity check: ssGSEA on pseudo-bulk vs observed proportions
pbc = pd.read_csv(os.path.join(SCD, "pseudobulk_counts_celltype_sample.csv"), index_col=0)
samples = sorted(set(c.split("|")[0] for c in pbc.columns))
pb = pd.DataFrame({s: pbc[[c for c in pbc.columns if c.startswith(s + "|")]].sum(axis=1) for s in samples})
lcpm = np.log2(pb / pb.sum(axis=0) * 1e6 + 1); lcpm = lcpm[(lcpm > 1).mean(axis=1) > 0.2]
ss = ssgsea(lcpm); prop = pd.read_csv(os.path.join(SCD, "celltype_proportions.csv"), index_col=0)
corr = {ct: stats.spearmanr(ss.loc[ct, prop.index], prop[ct])[0] for ct in CTS if ct in prop.columns}
print("ssGSEA score vs true proportion (scRNA pseudo-bulk, n=%d samples): " % len(prop), {k: round(v, 2) for k, v in corr.items()})
pd.Series(corr).to_csv(os.path.join(OUT, "ssgsea_validation_pseudobulk_spearman.csv"))

# figure: forest-style dot plot of meta effect sizes per cell type (ATAA vs TAAD) + per-cohort heatmap
fig, ax = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw=dict(width_ratios=[1, 1.5]))
yy = np.arange(len(CTS))
for dis, col, off in [("ATAA", "#2A9D8F", -0.15), ("TAAD", "#C44E52", 0.15)]:
    ax[0].errorbar(piv[f"g_meta_{dis}"], yy + off, xerr=1.96 * piv[f"se_meta_{dis}"], fmt="o", color=col, ecolor=col, capsize=2, label=f"{dis} vs control (meta)")
ax[0].scatter(piv.g_direct_GSE318877, yy, marker="D", facecolor="none", edgecolor="k", s=30, label="direct TAAD vs ATAA (GSE318877)")
ax[0].axvline(0, c="k", lw=0.6); ax[0].set_yticks(yy); ax[0].set_yticklabels(CTS, fontsize=8); ax[0].invert_yaxis()
ax[0].set_xlabel("Hedges' g of ssGSEA cell-type score"); ax[0].legend(fontsize=7, frameon=False, loc="lower right"); ax[0].set_title("Inferred cell-composition shifts", fontsize=10)
for i, ct in enumerate(CTS):
    if piv.loc[ct, "p_diff"] < 0.05: ax[0].text(ax[0].get_xlim()[1], i, "#", fontsize=8, va="center")
H = eff[eff.disease != "direct"].pivot(index="celltype", columns="cohort", values="g").loc[CTS, DESIGN["ATAA"]["cohorts"] + DESIGN["TAAD"]["cohorts"]]
P = eff[eff.disease != "direct"].pivot(index="celltype", columns="cohort", values="p").loc[CTS, H.columns]
sns.heatmap(H, cmap="RdBu_r", center=0, vmin=-3, vmax=3, ax=ax[1], cbar_kws=dict(label="Hedges' g (disease vs control)", shrink=0.6), linewidths=0.3, linecolor="white")
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        if P.iloc[i, j] < 0.05: ax[1].text(j + 0.5, i + 0.5, "*", ha="center", va="center", fontsize=9)
ax[1].set_xticklabels([f"{c}\n({REG[c]['disease']})" for c in H.columns], rotation=90, fontsize=7); ax[1].set_yticklabels(CTS, fontsize=8); ax[1].axvline(4, color="k", lw=1.5); ax[1].set_ylabel(""); ax[1].set_title("Per-cohort effects (* P<0.05)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig4_cell_composition.png"), dpi=200); fig.savefig(os.path.join(FIG, "Fig4_cell_composition.pdf")); plt.close(fig)
