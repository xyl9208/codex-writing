"""Cross-cohort meta-analysis and ATAA-vs-TAAD comparison at the gene level.

Design (bulk tissue, disease vs non-diseased aorta)
  ATAA : GSE26155, GSE140947, GSE235161, GSE202267           (4 cohorts, n = 115)
  TAAD : GSE52093, GSE153434, GSE98770, GSE294606, GSE267434, GSE147026   (6 cohorts, n = 76)
  excluded: GSE190635 (DE profile anti-correlated with all other TAAD cohorts, rho -0.37 to -0.68 -> sample-label inconsistency)
  Direct : GSE318877 (TAAD vs ATAA aortic media, micro-regional RNA-seq)
  Independent validation: pseudo-bulk of scRNA-seq (GSE155468 ATAA, GSE213740 TAAD) handled in 08_validation.py

Per disease: Stouffer weighted-Z (weights sqrt(n)) + DerSimonian-Laird random-effects log2FC.
Consensus DEG: meta-FDR < 0.05, same direction in >= 75 % of cohorts, nominal P < 0.05 (same direction) in >= 50 % of cohorts.
Leave-one-cohort-out (LOCO) robustness: fraction of consensus DEGs retained when each cohort is dropped.
Gene classes: shared / TAAD-specific / ATAA-specific / discordant; disease-difference z_diff = (LFC_T - LFC_A)/sqrt(se_T^2+se_A^2).
"""
import os, json, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DE = os.path.join(ROOT, "results", "de"); META = os.path.join(ROOT, "results", "meta"); os.makedirs(META, exist_ok=True)
FIG = os.path.join(ROOT, "figures"); os.makedirs(FIG, exist_ok=True)
REG = json.load(open(os.path.join(ROOT, "results", "processed", "cohorts.json")))
DESIGN = {"ATAA": dict(cohorts=["GSE26155", "GSE140947", "GSE235161", "GSE202267"]),
          "TAAD": dict(cohorts=["GSE52093", "GSE153434", "GSE98770", "GSE294606", "GSE267434", "GSE147026"]),
          "excluded": ["GSE190635"], "direct": ["GSE318877"]}
json.dump(DESIGN, open(os.path.join(META, "design.json"), "w"), indent=1)

def load(name):
    d = pd.read_csv(os.path.join(DE, f"{name}_de.csv"), index_col=0)
    d = d[~d.index.duplicated()].dropna(subset=["pvalue", "log2FC"])
    d["z"] = np.sign(d["log2FC"]) * stats.norm.isf(d["pvalue"].clip(1e-300) / 2)
    d["se"] = d["se"].abs().replace(0, np.nan)
    d["n"] = sum(REG[name]["groups"].values())
    return d

def meta(dfs, min_cohorts=None):
    """Stouffer + DL random effects over genes present in >= min_cohorts cohorts."""
    k = len(dfs); min_cohorts = min_cohorts or max(2, int(np.ceil(0.75 * k)))
    genes = pd.Index(sorted(set().union(*[set(d.index) for d in dfs])))
    Z = pd.DataFrame({d["cohort"].iloc[0]: d["z"].reindex(genes) for d in dfs})
    Y = pd.DataFrame({d["cohort"].iloc[0]: d["log2FC"].reindex(genes) for d in dfs})
    V = pd.DataFrame({d["cohort"].iloc[0]: (d["se"] ** 2).reindex(genes) for d in dfs})
    P = pd.DataFrame({d["cohort"].iloc[0]: d["pvalue"].reindex(genes) for d in dfs})
    W = pd.Series({d["cohort"].iloc[0]: np.sqrt(d["n"].iloc[0]) for d in dfs})
    present = Z.notna(); nk = present.sum(axis=1)
    keep = nk >= min_cohorts
    Z, Y, V, P, present, nk = Z[keep], Y[keep], V[keep], P[keep], present[keep], nk[keep]
    Wm = present * W
    zmeta = (Z.fillna(0) * Wm).sum(axis=1) / np.sqrt((Wm ** 2).sum(axis=1))
    p = 2 * stats.norm.sf(np.abs(zmeta))
    V = V.where(V.notna() & np.isfinite(V), np.nan)
    V = V.apply(lambda col: col.fillna(col.max()))
    w = (1 / V).where(present)
    yfix = (w * Y).sum(axis=1) / w.sum(axis=1)
    Q = (w * (Y.sub(yfix, axis=0)) ** 2).sum(axis=1); dfq = nk - 1
    C = w.sum(axis=1) - (w ** 2).sum(axis=1) / w.sum(axis=1)
    tau2 = ((Q - dfq) / C).clip(lower=0)
    wr = (1 / V.add(tau2, axis=0)).where(present)
    yrem = (wr * Y).sum(axis=1) / wr.sum(axis=1); se_rem = np.sqrt(1 / wr.sum(axis=1))
    sign_meta = np.sign(yrem)
    same = (np.sign(Y).eq(sign_meta, axis=0) & present).sum(axis=1) / nk
    nominal = ((P < 0.05) & np.sign(Y).eq(sign_meta, axis=0) & present).sum(axis=1) / nk
    out = pd.DataFrame({"z_meta": zmeta, "p_meta": p, "fdr_meta": multipletests(p, method="fdr_bh")[1], "lfc_rem": yrem, "se_rem": se_rem,
                        "I2": ((Q - dfq) / Q).clip(0, 1) * 100, "n_cohorts": nk, "frac_same_direction": same, "frac_nominal_same_dir": nominal})
    out["consensus_DEG"] = (out.fdr_meta < 0.05) & (out.frac_same_direction >= 0.75) & (out.frac_nominal_same_dir >= 0.5)
    out["direction"] = np.where(out.lfc_rem > 0, "up", "down")
    for c in Y.columns: out[f"lfc_{c}"] = Y[c]; out[f"p_{c}"] = P[c]
    return out

res = {}; loco = {}
for dis in ["ATAA", "TAAD"]:
    dfs = [load(c) for c in DESIGN[dis]["cohorts"]]
    m = meta(dfs).sort_values("p_meta"); m.to_csv(os.path.join(META, f"{dis}_meta.csv")); res[dis] = m
    deg = m[m.consensus_DEG]
    print(f"[{dis}] {len(dfs)} cohorts, n={sum(d['n'].iloc[0] for d in dfs)}; genes={len(m)}; consensus DEGs={len(deg)} (up {int((deg.direction=='up').sum())}, down {int((deg.direction=='down').sum())}); "
          f"median I2 (DEGs)={deg.I2.median():.1f}%; median I2 (all)={m.I2.median():.1f}%")
    # leave-one-cohort-out robustness
    rows = []
    for drop in DESIGN[dis]["cohorts"]:
        sub = [d for d in dfs if d["cohort"].iloc[0] != drop]
        ml = meta(sub, min_cohorts=max(2, int(np.ceil(0.75 * len(sub)))))
        kept = ml.consensus_DEG.reindex(deg.index).fillna(False)
        same = (np.sign(ml.lfc_rem.reindex(deg.index)) == np.sign(deg.lfc_rem))
        rows.append(dict(disease=dis, dropped=drop, n_DEG_without=int(ml.consensus_DEG.sum()), frac_consensus_retained=float(kept.mean()), frac_direction_retained=float(same.mean()),
                         rho_lfc=float(stats.spearmanr(ml.lfc_rem.reindex(deg.index).fillna(0), deg.lfc_rem)[0])))
    loco[dis] = pd.DataFrame(rows); print(loco[dis].round(3).to_string(index=False))
pd.concat(loco.values()).to_csv(os.path.join(META, "LOCO_robustness.csv"), index=False)

# ---------------------------------------------------------------- comparison ATAA vs TAAD
A, T = res["ATAA"], res["TAAD"]
common = A.index.intersection(T.index)
cmp = pd.DataFrame({"lfc_ATAA": A.loc[common, "lfc_rem"], "se_ATAA": A.loc[common, "se_rem"], "z_ATAA": A.loc[common, "z_meta"], "fdr_ATAA": A.loc[common, "fdr_meta"], "DEG_ATAA": A.loc[common, "consensus_DEG"], "I2_ATAA": A.loc[common, "I2"],
                    "lfc_TAAD": T.loc[common, "lfc_rem"], "se_TAAD": T.loc[common, "se_rem"], "z_TAAD": T.loc[common, "z_meta"], "fdr_TAAD": T.loc[common, "fdr_meta"], "DEG_TAAD": T.loc[common, "consensus_DEG"], "I2_TAAD": T.loc[common, "I2"]})
cmp["diff_lfc"] = cmp.lfc_TAAD - cmp.lfc_ATAA
cmp["z_diff"] = cmp.diff_lfc / np.sqrt(cmp.se_TAAD ** 2 + cmp.se_ATAA ** 2)
cmp["p_diff"] = 2 * stats.norm.sf(cmp.z_diff.abs()); cmp["fdr_diff"] = multipletests(cmp.p_diff, method="fdr_bh")[1]
def classify(r):
    a, t = r.DEG_ATAA, r.DEG_TAAD
    if a and t: return "shared_concordant" if np.sign(r.lfc_ATAA) == np.sign(r.lfc_TAAD) else "discordant"
    if t and not a: return "TAAD_specific" if (r.fdr_ATAA > 0.10 or np.sign(r.lfc_ATAA) != np.sign(r.lfc_TAAD)) else "TAAD_stronger"
    if a and not t: return "ATAA_specific" if (r.fdr_TAAD > 0.10 or np.sign(r.lfc_ATAA) != np.sign(r.lfc_TAAD)) else "ATAA_stronger"
    return "none"
cmp["class"] = cmp.apply(classify, axis=1)
cmp["class_strict"] = np.where(cmp["class"].isin(["TAAD_specific", "ATAA_specific"]) & (cmp.fdr_diff >= 0.05), cmp["class"] + "_ns", cmp["class"])
cmp = cmp.sort_values("p_diff"); cmp.to_csv(os.path.join(META, "ATAA_vs_TAAD_gene_comparison.csv"))
print("\nGene classes (consensus DEGs):"); print(cmp["class"].value_counts().to_string())
print("\nStrict classes (specific + significant disease difference FDR<0.05):"); print(cmp["class_strict"].value_counts().to_string())
r_all = stats.spearmanr(cmp.lfc_ATAA, cmp.lfc_TAAD); degs = cmp[cmp.DEG_ATAA | cmp.DEG_TAAD]; r_deg = stats.spearmanr(degs.lfc_ATAA, degs.lfc_TAAD)
r_z = stats.spearmanr(cmp.z_ATAA, cmp.z_TAAD)
print(f"\nmeta log2FC correlation ATAA vs TAAD: all genes rho={r_all[0]:.3f} (n={len(cmp)}); z_meta rho={r_z[0]:.3f}; union of DEGs rho={r_deg[0]:.3f} (n={len(degs)})")
# overlap statistics (hypergeometric) of consensus DEG sets
from scipy.stats import hypergeom
nA, nT, nAT, N = int(cmp.DEG_ATAA.sum()), int(cmp.DEG_TAAD.sum()), int((cmp.DEG_ATAA & cmp.DEG_TAAD).sum()), len(cmp)
p_over = hypergeom.sf(nAT - 1, N, nA, nT); exp = nA * nT / N
print(f"DEG overlap: ATAA {nA}, TAAD {nT}, overlap {nAT} (expected {exp:.1f}; hypergeometric P={p_over:.2e}); Jaccard={nAT/(nA+nT-nAT):.3f}")

# ---------------------------------------------------------------- validation with the direct contrast (GSE318877, TAAD vs ATAA)
direct = load("GSE318877")
cc = cmp.join(direct[["log2FC", "pvalue", "stat"]].rename(columns={"log2FC": "lfc_direct", "pvalue": "p_direct", "stat": "stat_direct"}), how="inner")
rho_all = stats.spearmanr(cc.diff_lfc, cc.lfc_direct); rho_z = stats.spearmanr(cc.z_diff, cc.stat_direct)
sig = cc[cc.fdr_diff < 0.05]
rho_sig = stats.spearmanr(sig.diff_lfc, sig.lfc_direct) if len(sig) > 10 else (np.nan, np.nan)
conc = (np.sign(sig.diff_lfc) == np.sign(sig.lfc_direct)).mean() * 100 if len(sig) else np.nan
print(f"\nDirect contrast GSE318877 (TAAD vs ATAA): rho(diff_lfc vs direct lfc, all genes)={rho_all[0]:.3f} p={rho_all[1]:.1e} n={len(cc)}; rho(z_diff vs direct stat)={rho_z[0]:.3f}; "
      f"genes with FDR_diff<0.05: n={len(sig)}, rho={rho_sig[0]:.3f}, direction concordance={conc:.1f}%")
val_rows = []
for cl in ["TAAD_specific", "ATAA_specific", "shared_concordant", "discordant"]:
    s = cc[cc["class"] == cl]
    if len(s) < 5: continue
    exp = np.sign(s.lfc_TAAD) if cl == "TAAD_specific" else (-np.sign(s.lfc_ATAA) if cl == "ATAA_specific" else np.sign(s.diff_lfc))
    agree = (np.sign(s.lfc_direct) == exp).mean() * 100; med = np.median(s.lfc_direct * exp)
    w = stats.wilcoxon(s.lfc_direct * exp) if len(s) > 10 else None
    val_rows.append(dict(cls=cl, n=len(s), sign_agreement=agree, median_signed_lfc=med, wilcoxon_p=w.pvalue if w else np.nan))
    print(f"    {cl}: n={len(s)}, sign agreement with direct contrast={agree:.1f}%, median signed direct lfc={med:.3f}" + (f", Wilcoxon p={w.pvalue:.1e}" if w else ""))
pd.DataFrame(val_rows).to_csv(os.path.join(META, "direct_contrast_validation_by_class.csv"), index=False)
cc.to_csv(os.path.join(META, "ATAA_vs_TAAD_with_direct_GSE318877.csv"))

# ---------------------------------------------------------------- cohort concordance heatmap (Figure 1 panel)
names = DESIGN["ATAA"]["cohorts"] + DESIGN["TAAD"]["cohorts"] + DESIGN["excluded"]
de = {n: load(n) for n in names}
M = pd.DataFrame(index=names, columns=names, dtype=float)
for a in names:
    for b in names:
        common_ = de[a].index.intersection(de[b].index); sa, sb = de[a].loc[common_], de[b].loc[common_]
        sel = (sa.padj < 0.05) | (sb.padj < 0.05)
        if sel.sum() < 200:
            top = pd.concat([sa.stat.abs(), sb.stat.abs()], axis=1).max(axis=1).sort_values(ascending=False).index[:1000]; sel = common_.isin(top)
        M.loc[a, b] = stats.spearmanr(sa.stat[sel], sb.stat[sel])[0]
M.to_csv(os.path.join(META, "cohort_concordance_rho.csv"))
lab = [f"{n} ({REG[n]['disease']}, n={sum(REG[n]['groups'].values())})" + (" *excluded" if n in DESIGN["excluded"] else "") for n in names]
fig, ax = plt.subplots(figsize=(8.2, 7))
sns.heatmap(M.astype(float), cmap="RdBu_r", center=0, vmin=-0.8, vmax=0.8, annot=True, fmt=".2f", annot_kws={"size": 7}, ax=ax, xticklabels=lab, yticklabels=lab, cbar_kws=dict(label="Spearman rho of DE statistics", shrink=0.6))
ax.tick_params(axis="x", labelsize=7, rotation=90); ax.tick_params(axis="y", labelsize=7)
ax.axhline(4, c="k", lw=1.5); ax.axvline(4, c="k", lw=1.5); ax.axhline(10, c="k", lw=1); ax.axvline(10, c="k", lw=1)
ax.set_title("Concordance of disease-vs-control DE profiles across cohorts", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig1b_cohort_concordance.png"), dpi=200); fig.savefig(os.path.join(FIG, "Fig1b_cohort_concordance.pdf")); plt.close(fig)
within = {"ATAA": M.loc[DESIGN["ATAA"]["cohorts"], DESIGN["ATAA"]["cohorts"]].values, "TAAD": M.loc[DESIGN["TAAD"]["cohorts"], DESIGN["TAAD"]["cohorts"]].values}
for k, v in within.items():
    iu = np.triu_indices(v.shape[0], 1); print(f"within-{k} cohort concordance: median rho={np.median(v[iu]):.2f} (range {v[iu].min():.2f} to {v[iu].max():.2f})")
between = M.loc[DESIGN["ATAA"]["cohorts"], DESIGN["TAAD"]["cohorts"]].values; print(f"between ATAA-TAAD cohort concordance: median rho={np.median(between):.2f} (range {between.min():.2f} to {between.max():.2f})")
wa = within["ATAA"][np.triu_indices(4, 1)]; wt = within["TAAD"][np.triu_indices(6, 1)]
print("Mann-Whitney within-TAAD vs within-ATAA rho: p=%.3g; within-TAAD vs between: p=%.3g" % (stats.mannwhitneyu(wt, wa).pvalue, stats.mannwhitneyu(wt, between.ravel()).pvalue))

# ---------------------------------------------------------------- summary figure (gene level)
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
cols = {"shared_concordant": "#7B2CBF", "TAAD_specific": "#C44E52", "ATAA_specific": "#2A9D8F", "discordant": "#F4A261", "TAAD_stronger": "#E5989B", "ATAA_stronger": "#95D5B2", "none": "#DDDDDD"}
for cl in ["none", "TAAD_stronger", "ATAA_stronger", "shared_concordant", "TAAD_specific", "ATAA_specific", "discordant"]:
    s = cmp[cmp["class"] == cl]
    ax[0].scatter(s.lfc_ATAA, s.lfc_TAAD, s=4 if cl == "none" else 9, c=cols[cl], label=f"{cl} ({len(s)})", rasterized=True, alpha=0.8)
ax[0].axhline(0, c="k", lw=0.5); ax[0].axvline(0, c="k", lw=0.5)
ax[0].set_xlabel("ATAA vs control, meta log2FC (random effects)"); ax[0].set_ylabel("TAAD vs control, meta log2FC (random effects)")
ax[0].set_title(f"Gene-level concordance (rho={r_all[0]:.2f}, n={len(cmp)})", fontsize=10); ax[0].legend(fontsize=6.5, frameon=False, markerscale=2)
ax[1].scatter(cc.diff_lfc, cc.lfc_direct, s=4, c="lightgrey", rasterized=True)
ax[1].scatter(sig.diff_lfc, sig.lfc_direct, s=8, c="#C44E52", rasterized=True, label=f"FDR_diff<0.05 (n={len(sig)})")
ax[1].axhline(0, c="k", lw=0.5); ax[1].axvline(0, c="k", lw=0.5)
ax[1].set_xlabel("Predicted difference: log2FC(TAAD) - log2FC(ATAA), meta"); ax[1].set_ylabel("Observed TAAD vs ATAA log2FC\n(GSE318877 aortic media, direct)")
ax[1].set_title(f"Direct contrast cohort: rho={rho_all[0]:.2f} (all), {rho_sig[0]:.2f} (sig.)", fontsize=10); ax[1].legend(fontsize=7, frameon=False)
vc = cmp["class"].value_counts().drop("none", errors="ignore")
ax[2].barh(vc.index[::-1], vc.values[::-1], color=[cols[c] for c in vc.index[::-1]])
for i, v in enumerate(vc.values[::-1]): ax[2].text(v, i, f" {v}", va="center", fontsize=8)
ax[2].set_xlabel("number of genes"); ax[2].set_title("Classification of consensus DEGs", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig2_gene_level_comparison.png"), dpi=200); fig.savefig(os.path.join(FIG, "Fig2_gene_level_comparison.pdf")); plt.close(fig)
for cl in ["TAAD_specific", "ATAA_specific", "shared_concordant", "discordant"]:
    s = cmp[cmp["class"] == cl].copy()
    s["abs_z"] = (s.z_TAAD.abs() if "TAAD" in cl else s.z_ATAA.abs()) if cl not in ("shared_concordant", "discordant") else (s.z_TAAD.abs() + s.z_ATAA.abs())
    s.sort_values("abs_z", ascending=False).head(50).to_csv(os.path.join(META, f"top_{cl}.csv"))
