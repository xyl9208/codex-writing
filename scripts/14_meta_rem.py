"""Meta-analysis v2: random-effects (DerSimonian-Laird) inference as the PRIMARY test, with
Hartung-Knapp-Sidik-Jonkman (HKSJ) adjustment and prediction intervals as small-k sensitivity analyses.

Also:
  * retention of the v1 (Stouffer) consensus sets (set intersection, Jaccard) -- not a ratio of totals
  * gene classification v2: cross-disease difference test (FDR_diff) is the primary criterion;
    TOST equivalence (|log2FC| < 0.5, 90 % CI) adds a higher evidence tier for 'near-zero in the other disease';
    'only in one disease, difference not established' is reported as undetermined, not as specific
  * LOCO robustness under the v2 rule
  * GSE190635 inclusion sensitivity (TAAD, 7 cohorts)
  * TAAD 4-of-6 cohort subsampling to separate 'number of cohorts' from 'true cross-cohort consistency'
  * GSEA (Hallmark/KEGG/Reactome) on the v2 (z_rem) rankings and comparison with v1 NES
"""
import os, json, itertools, warnings, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DE = os.path.join(ROOT, "results", "de"); META1 = os.path.join(ROOT, "results", "meta"); OUT = os.path.join(ROOT, "results", "meta_v2"); os.makedirs(OUT, exist_ok=True)
REG = json.load(open(os.path.join(ROOT, "results", "processed", "cohorts.json"))); DESIGN = json.load(open(os.path.join(META1, "design.json")))
EQUIV_BOUND = 0.5  # log2 units (~1.41-fold), pre-specified equivalence margin

def load(name):
    d = pd.read_csv(os.path.join(DE, f"{name}_de.csv"), index_col=0)
    d = d[~d.index.duplicated()].dropna(subset=["pvalue", "log2FC"])
    d["z"] = np.sign(d["log2FC"]) * stats.norm.isf(d["pvalue"].clip(1e-300) / 2)
    d["se"] = d["se"].abs().replace(0, np.nan); d["n"] = sum(REG[name]["groups"].values())
    return d

def rem_meta(dfs, min_cohorts=None):
    k = len(dfs); min_cohorts = min_cohorts or max(2, int(np.ceil(0.75 * k)))
    genes = pd.Index(sorted(set().union(*[set(d.index) for d in dfs])))
    names = [d["cohort"].iloc[0] for d in dfs]
    Y = pd.DataFrame({n: d["log2FC"].reindex(genes) for n, d in zip(names, dfs)}); V = pd.DataFrame({n: (d["se"] ** 2).reindex(genes) for n, d in zip(names, dfs)})
    P = pd.DataFrame({n: d["pvalue"].reindex(genes) for n, d in zip(names, dfs)}); Zc = pd.DataFrame({n: d["z"].reindex(genes) for n, d in zip(names, dfs)})
    W = pd.Series({n: np.sqrt(d["n"].iloc[0]) for n, d in zip(names, dfs)})
    present = Y.notna() & V.notna() & np.isfinite(V); nk = present.sum(axis=1); keep = nk >= min_cohorts
    Y, V, P, Zc, present, nk = Y[keep], V[keep], P[keep], Zc[keep], present[keep], nk[keep]
    V = V.where(present); w = 1 / V
    yfix = (w * Y).sum(axis=1) / w.sum(axis=1)
    Q = (w * (Y.sub(yfix, axis=0)) ** 2).sum(axis=1); dfq = nk - 1
    C = w.sum(axis=1) - (w ** 2).sum(axis=1) / w.sum(axis=1); tau2 = ((Q - dfq) / C).clip(lower=0)
    wr = (1 / V.add(tau2, axis=0)).where(present); yrem = (wr * Y).sum(axis=1) / wr.sum(axis=1); se_rem = np.sqrt(1 / wr.sum(axis=1))
    z_rem = yrem / se_rem; p_rem = 2 * stats.norm.sf(np.abs(z_rem))
    # HKSJ: t-distribution with k-1 df and variance scaled by q
    q = (wr * (Y.sub(yrem, axis=0)) ** 2).sum(axis=1) / (nk - 1)
    se_hk = np.sqrt(q.clip(lower=1e-12)) * se_rem  # standard HKSJ (unmodified)
    se_hk_mod = np.sqrt(q.clip(lower=1.0)) * se_rem  # modified HKSJ (q >= 1), the conservative variant
    p_hk = 2 * stats.t.sf(np.abs(yrem / se_hk), df=nk - 1); p_hk_mod = 2 * stats.t.sf(np.abs(yrem / se_hk_mod), df=nk - 1)
    # prediction interval (k >= 3): yrem +/- t_{k-2} * sqrt(se_rem^2 + tau2)
    tcrit = np.where(nk >= 3, stats.t.ppf(0.975, df=np.clip(nk - 2, 1, None)), np.nan)
    pi_half = tcrit * np.sqrt(se_rem ** 2 + tau2)
    # Stouffer (secondary)
    Wm = present * W; zst = (Zc.fillna(0) * Wm).sum(axis=1) / np.sqrt((Wm ** 2).sum(axis=1))
    sign = np.sign(yrem)
    same = (np.sign(Y).eq(sign, axis=0) & present).sum(axis=1) / nk
    nominal = ((P < 0.05) & np.sign(Y).eq(sign, axis=0) & present).sum(axis=1) / nk
    out = pd.DataFrame({"lfc_rem": yrem, "se_rem": se_rem, "z_rem": z_rem, "p_rem": p_rem, "fdr_rem": multipletests(p_rem, method="fdr_bh")[1],
                        "ci95_low": yrem - 1.96 * se_rem, "ci95_high": yrem + 1.96 * se_rem, "pi95_low": yrem - pi_half, "pi95_high": yrem + pi_half,
                        "p_hksj": p_hk, "fdr_hksj": multipletests(p_hk, method="fdr_bh")[1], "p_hksj_mod": p_hk_mod, "fdr_hksj_mod": multipletests(p_hk_mod, method="fdr_bh")[1],
                        "tau2": tau2, "I2": ((Q - dfq) / Q).clip(0, 1) * 100, "n_cohorts": nk, "frac_same_direction": same, "frac_nominal_same_dir": nominal,
                        "z_stouffer": zst, "fdr_stouffer": multipletests(2 * stats.norm.sf(np.abs(zst)), method="fdr_bh")[1]})
    out["consensus_v2"] = (out.fdr_rem < 0.05) & (out.frac_same_direction >= 0.75) & (out.frac_nominal_same_dir >= 0.5)
    out["consensus_hksj"] = (out.fdr_hksj < 0.05) & (out.frac_same_direction >= 0.75) & (out.frac_nominal_same_dir >= 0.5)
    out["direction"] = np.where(out.lfc_rem > 0, "up", "down")
    # TOST equivalence to zero (90 % CI within +/- EQUIV_BOUND)
    out["equiv_zero"] = ((yrem + EQUIV_BOUND) / se_rem > stats.norm.ppf(0.95)) & ((yrem - EQUIV_BOUND) / se_rem < -stats.norm.ppf(0.95))
    for n in names: out[f"lfc_{n}"] = Y[n]; out[f"p_{n}"] = P[n]
    return out

# ------------------------------------------------------------------ primary meta v2 + retention + LOCO
res = {}; summ = {}
for dis in ["ATAA", "TAAD"]:
    dfs = [load(c) for c in DESIGN[dis]["cohorts"]]
    m = rem_meta(dfs).sort_values("p_rem"); m.to_csv(os.path.join(OUT, f"{dis}_meta_v2.csv")); res[dis] = m
    v1 = pd.read_csv(os.path.join(META1, f"{dis}_meta.csv"), index_col=0); s1 = set(v1.index[v1.consensus_DEG]); s2 = set(m.index[m.consensus_v2]); s3 = set(m.index[m.consensus_hksj])
    deg = m[m.consensus_v2]
    summ[dis] = dict(n_genes=int(len(m)), consensus_v1_stouffer=len(s1), consensus_v2_rem=len(s2), intersection_v1_v2=len(s1 & s2), retention_of_v1=len(s1 & s2) / len(s1), jaccard_v1_v2=len(s1 & s2) / len(s1 | s2),
                     consensus_hksj=len(s3), intersection_v2_hksj=len(s2 & s3), consensus_hksj_mod=int(((m.fdr_hksj_mod < 0.05) & (m.frac_same_direction >= 0.75) & (m.frac_nominal_same_dir >= 0.5)).sum()),
                     n_up=int((deg.direction == "up").sum()), n_down=int((deg.direction == "down").sum()), median_I2_deg=float(deg.I2.median()), median_tau2_deg=float(deg.tau2.median()),
                     frac_deg_PI_excludes_zero=float(((deg.pi95_low > 0) | (deg.pi95_high < 0)).mean()))
    print(f"[{dis}] v2 REM consensus={len(s2)} (up {summ[dis]['n_up']}, down {summ[dis]['n_down']}); v1 Stouffer={len(s1)}; intersection={len(s1&s2)} (retention {len(s1&s2)/len(s1)*100:.0f}%, Jaccard {summ[dis]['jaccard_v1_v2']:.2f}); "
          f"HKSJ consensus={len(s3)} (modified HKSJ {summ[dis]['consensus_hksj_mod']}); DEGs whose 95% prediction interval excludes 0: {summ[dis]['frac_deg_PI_excludes_zero']*100:.0f}%")
    rows = []
    for drop in DESIGN[dis]["cohorts"]:
        sub = [d for d in dfs if d["cohort"].iloc[0] != drop]; ml = rem_meta(sub)
        kept = ml.consensus_v2.reindex(deg.index).fillna(False)
        rows.append(dict(disease=dis, dropped=drop, n_consensus_without=int(ml.consensus_v2.sum()), frac_consensus_retained=float(kept.mean()),
                         frac_direction_retained=float((np.sign(ml.lfc_rem.reindex(deg.index)) == np.sign(deg.lfc_rem)).mean()), rho_lfc=float(stats.spearmanr(ml.lfc_rem.reindex(deg.index).fillna(0), deg.lfc_rem)[0])))
    loco = pd.DataFrame(rows); loco.to_csv(os.path.join(OUT, f"LOCO_v2_{dis}.csv"), index=False); print(loco.round(3).to_string(index=False))
    summ[dis]["loco_retention_min"] = float(loco.frac_consensus_retained.min()); summ[dis]["loco_retention_median"] = float(loco.frac_consensus_retained.median())

# ------------------------------------------------------------------ classification v2
A, T = res["ATAA"], res["TAAD"]; common = A.index.intersection(T.index)
c = pd.DataFrame({"lfc_ATAA": A.loc[common, "lfc_rem"], "se_ATAA": A.loc[common, "se_rem"], "fdr_ATAA": A.loc[common, "fdr_rem"], "DEG_ATAA": A.loc[common, "consensus_v2"], "equiv_ATAA": A.loc[common, "equiv_zero"], "I2_ATAA": A.loc[common, "I2"], "pi_low_ATAA": A.loc[common, "pi95_low"], "pi_high_ATAA": A.loc[common, "pi95_high"],
                  "lfc_TAAD": T.loc[common, "lfc_rem"], "se_TAAD": T.loc[common, "se_rem"], "fdr_TAAD": T.loc[common, "fdr_rem"], "DEG_TAAD": T.loc[common, "consensus_v2"], "equiv_TAAD": T.loc[common, "equiv_zero"], "I2_TAAD": T.loc[common, "I2"], "pi_low_TAAD": T.loc[common, "pi95_low"], "pi_high_TAAD": T.loc[common, "pi95_high"]})
c["diff_lfc"] = c.lfc_TAAD - c.lfc_ATAA; c["se_diff"] = np.sqrt(c.se_TAAD ** 2 + c.se_ATAA ** 2); c["z_diff"] = c.diff_lfc / c.se_diff
c["p_diff"] = 2 * stats.norm.sf(c.z_diff.abs()); c["fdr_diff"] = multipletests(c.p_diff, method="fdr_bh")[1]
def classify(r):
    a, t, d = r.DEG_ATAA, r.DEG_TAAD, r.fdr_diff < 0.05
    if a and t:
        if np.sign(r.lfc_ATAA) != np.sign(r.lfc_TAAD): return "discordant"
        return "shared_magnitude_differs" if d else "shared_concordant"
    if t and not a:
        if not d: return "TAAD_only_undetermined"
        return "TAAD_enriched_equivATAA" if r.equiv_ATAA else "TAAD_enriched"
    if a and not t:
        if not d: return "ATAA_only_undetermined"
        return "ATAA_enriched_equivTAAD" if r.equiv_TAAD else "ATAA_enriched"
    return "none"
c["class_v2"] = c.apply(classify, axis=1)
# for continuity: map v1 class
v1c = pd.read_csv(os.path.join(META1, "ATAA_vs_TAAD_gene_comparison.csv"), index_col=0)["class"]; c["class_v1"] = v1c.reindex(c.index)
c = c.sort_values("p_diff"); c.to_csv(os.path.join(OUT, "gene_classes_v2.csv"))
vc = c.class_v2.value_counts(); print("\nGene classes v2:"); print(vc.to_string())
rA, rT = stats.spearmanr(c.lfc_ATAA, c.lfc_TAAD)[0], None
nA, nT, nAT = int(c.DEG_ATAA.sum()), int(c.DEG_TAAD.sum()), int((c.DEG_ATAA & c.DEG_TAAD).sum())
print(f"log2FC correlation ATAA vs TAAD (v2): rho={rA:.3f}; DEG overlap {nAT} of ATAA {nA} / TAAD {nT} (expected {nA*nT/len(c):.1f}; hypergeom P={stats.hypergeom.sf(nAT-1, len(c), nA, nT):.1e}); Jaccard={nAT/(nA+nT-nAT):.3f}")
summ["classes_v2"] = vc.to_dict(); summ["rho_lfc_v2"] = float(rA); summ["overlap"] = dict(ATAA=nA, TAAD=nT, overlap=nAT, expected=nA * nT / len(c))

# ------------------------------------------------------------------ GSE190635 inclusion sensitivity
dfs7 = [load(cc) for cc in DESIGN["TAAD"]["cohorts"] + DESIGN["excluded"]]
m7 = rem_meta(dfs7); s6 = set(T.index[T.consensus_v2]); s7 = set(m7.index[m7.consensus_v2])
cc7 = pd.DataFrame({"lfc_ATAA": A.loc[common, "lfc_rem"], "se_ATAA": A.loc[common, "se_rem"], "DEG_ATAA": A.loc[common, "consensus_v2"], "equiv_ATAA": A.loc[common, "equiv_zero"],
                    "lfc_TAAD": m7.lfc_rem.reindex(common), "se_TAAD": m7.se_rem.reindex(common), "DEG_TAAD": m7.consensus_v2.reindex(common).fillna(False), "equiv_TAAD": m7.equiv_zero.reindex(common).fillna(False)}).dropna(subset=["lfc_TAAD"])
cc7["z_diff"] = (cc7.lfc_TAAD - cc7.lfc_ATAA) / np.sqrt(cc7.se_TAAD ** 2 + cc7.se_ATAA ** 2); cc7["fdr_diff"] = multipletests(2 * stats.norm.sf(cc7.z_diff.abs()), method="fdr_bh")[1]
cc7["class_v2"] = cc7.apply(classify, axis=1)
sens = dict(TAAD_consensus_6=len(s6), TAAD_consensus_7_with_GSE190635=len(s7), intersection=len(s6 & s7), retention_of_6=len(s6 & s7) / len(s6), rho_lfc_6_vs_7=float(stats.spearmanr(T.lfc_rem.reindex(common), m7.lfc_rem.reindex(common), nan_policy="omit")[0]),
            classes_with_GSE190635=cc7.class_v2.value_counts().to_dict(), median_I2_deg_7=float(m7[m7.consensus_v2].I2.median()))
json.dump(sens, open(os.path.join(OUT, "GSE190635_inclusion_sensitivity.json"), "w"), indent=1); m7.to_csv(os.path.join(OUT, "TAAD_meta_v2_with_GSE190635.csv"))
print("\nGSE190635 inclusion sensitivity:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in sens.items()})

# ------------------------------------------------------------------ TAAD 4-of-6 subsampling vs ATAA (same rule, shared gene universe)
rho = pd.read_csv(os.path.join(META1, "cohort_concordance_rho.csv"), index_col=0)
ataa_genes = set(A.index); rows = []
for combo in itertools.combinations(DESIGN["TAAD"]["cohorts"], 4):
    dfs4 = [load(cc) for cc in combo]; m4 = rem_meta(dfs4)
    uni = ataa_genes & set(m4.index)
    n4 = int(m4.consensus_v2.reindex(uni).sum()); nA4 = int(A.consensus_v2.reindex(uni).sum())
    sub = rho.loc[list(combo), list(combo)].values; iu = np.triu_indices(4, 1)
    # LOCO retention within the 4-cohort subset
    d4 = m4[m4.consensus_v2]; rets = []
    for drop in combo:
        ml = rem_meta([d for d in dfs4 if d["cohort"].iloc[0] != drop], min_cohorts=3); rets.append(float(ml.consensus_v2.reindex(d4.index).fillna(False).mean()))
    rows.append(dict(cohorts="+".join(combo), n_total=sum(sum(REG[cc]["groups"].values()) for cc in combo), median_rho=float(np.median(sub[iu])), n_consensus_shared_universe=n4, n_consensus_ATAA_same_universe=nA4,
                     loco_retention_min=min(rets), loco_retention_median=float(np.median(rets)), median_I2_deg=float(d4.I2.median()) if len(d4) else np.nan))
sub4 = pd.DataFrame(rows); sub4.to_csv(os.path.join(OUT, "TAAD_4of6_subsets.csv"), index=False)
ia = rho.loc[DESIGN["ATAA"]["cohorts"], DESIGN["ATAA"]["cohorts"]].values[np.triu_indices(4, 1)]
la = pd.read_csv(os.path.join(OUT, "LOCO_v2_ATAA.csv"))
print("\nTAAD 4-of-6 subsets (15): median concordance rho range %.2f-%.2f (median %.2f) vs ATAA 4 cohorts %.2f; consensus DEGs (shared universe) range %d-%d (median %d) vs ATAA %d-%d; LOCO min retention median %.2f (range %.2f-%.2f) vs ATAA %.2f" % (
    sub4.median_rho.min(), sub4.median_rho.max(), sub4.median_rho.median(), np.median(ia), sub4.n_consensus_shared_universe.min(), sub4.n_consensus_shared_universe.max(), sub4.n_consensus_shared_universe.median(),
    sub4.n_consensus_ATAA_same_universe.min(), sub4.n_consensus_ATAA_same_universe.max(), sub4.loco_retention_min.median(), sub4.loco_retention_min.min(), sub4.loco_retention_min.max(), la.frac_consensus_retained.min()))
summ["TAAD_4of6"] = dict(median_rho_range=[float(sub4.median_rho.min()), float(sub4.median_rho.max())], consensus_range=[int(sub4.n_consensus_shared_universe.min()), int(sub4.n_consensus_shared_universe.max())],
                         ATAA_consensus_same_universe_range=[int(sub4.n_consensus_ATAA_same_universe.min()), int(sub4.n_consensus_ATAA_same_universe.max())], loco_min_retention_range=[float(sub4.loco_retention_min.min()), float(sub4.loco_retention_min.max())],
                         ATAA_median_rho=float(np.median(ia)), ATAA_loco_min_retention=float(la.frac_consensus_retained.min()))
json.dump(summ, open(os.path.join(OUT, "summary_v2.json"), "w"), indent=1, default=float)

# ------------------------------------------------------------------ GSEA on v2 rankings (z_rem) and comparison with v1
import gseapy as gp
GS = os.path.join(ROOT, "data", "genesets"); P2 = os.path.join(ROOT, "results", "pathways_v2"); os.makedirs(P2, exist_ok=True)
GMT = {"Hallmark": "h.all.v2024.1.Hs.symbols.gmt", "KEGG": "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt", "Reactome": "c2.cp.reactome.v2024.1.Hs.symbols.gmt"}
rank = {"meta_ATAA_v2": A["z_rem"].dropna().sort_values(ascending=False), "meta_TAAD_v2": T["z_rem"].dropna().sort_values(ascending=False), "meta_diff_v2": c["z_diff"].dropna().sort_values(ascending=False)}
allr = []
for tag, rnk in rank.items():
    rnk = rnk[~rnk.index.duplicated()]
    for lib, f in GMT.items():
        r = gp.prerank(rnk=rnk, gene_sets=os.path.join(GS, f), permutation_num=1000, min_size=15, max_size=500, seed=7, threads=4, outdir=None, verbose=False).res2d
        r = r.rename(columns={"Term": "term", "NES": "nes", "FDR q-val": "fdr"}); r["library"] = lib; r["ranking"] = tag; allr.append(r[["term", "library", "nes", "fdr", "ranking"]])
g2 = pd.concat(allr); g2["nes"] = g2.nes.astype(float); g2["fdr"] = g2.fdr.astype(float); g2.to_csv(os.path.join(P2, "gsea_meta_v2.csv"), index=False)
nes1 = pd.read_csv(os.path.join(ROOT, "results", "pathways", "NES_matrix.csv"), index_col=[0, 1])
piv = g2.pivot_table(index=["library", "term"], columns="ranking", values="nes"); fdr2 = g2.pivot_table(index=["library", "term"], columns="ranking", values="fdr")
for a_, b_ in [("meta_ATAA_v2", "meta_ATAA"), ("meta_TAAD_v2", "meta_TAAD"), ("meta_diff_v2", "meta_diff_TAADminusATAA")]:
    j = piv[a_].to_frame("v2").join(nes1[b_].rename("v1")).dropna(); print(f"NES {a_} vs {b_}: rho={stats.spearmanr(j.v2, j.v1)[0]:.3f} (n={len(j)})")
hm = piv.loc["Hallmark"].join(fdr2.loc["Hallmark"], rsuffix="_fdr")
hm.to_csv(os.path.join(P2, "hallmark_NES_v2.csv")); print(hm.sort_values("meta_diff_v2").round(2).to_string())
