"""F1 (bulk): Tier 1 programme scores C, I, R per patient in each cohort (label-free standardisation within cohort),
per-cohort Hedges' g (disease vs control) with bootstrap CI, DL random-effects pooling per disease, and the
cohort-level ATAA-vs-TAAD comparison (3 tests: C, I, R; Holm).  Exploratory: every Tier 1 module (BH).
GSE318877: patient-level (mean of punches) C / I / R, dissection (subacute/chronic) vs dilatation, estimation only."""
import os, sys, json, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from v21_lib import *
PROC = os.path.join(ROOT, "results", "processed"); OUT = os.path.join(ROOT, "results", "v21_bulk"); os.makedirs(OUT, exist_ok=True)
REG = json.load(open(os.path.join(PROC, "cohorts.json"))); DESIGN = json.load(open(os.path.join(ROOT, "results", "meta", "design.json")))
MOD = load_tier1(); cohorts = {c: "ATAA" for c in DESIGN["ATAA"]["cohorts"]}; cohorts.update({c: "TAAD" for c in DESIGN["TAAD"]["cohorts"]})

def rem(g, se):
    g, se = np.asarray(g, float), np.asarray(se, float); w = 1 / se ** 2; yf = np.sum(w * g) / np.sum(w); Q = np.sum(w * (g - yf) ** 2); df = len(g) - 1
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w); tau2 = max(0, (Q - df) / C) if C > 0 else 0; wr = 1 / (se ** 2 + tau2); y = np.sum(wr * g) / np.sum(wr); s = np.sqrt(1 / np.sum(wr))
    return y, s, 2 * stats.norm.sf(abs(y / s)), max(0, (Q - df) / Q * 100) if Q > 0 else 0
def g_se(g, nx_, ny): return np.sqrt((nx_ + ny) / (nx_ * ny) + g ** 2 / (2 * (nx_ + ny)))

allscores = []; eff = []
for c, dis in cohorts.items():
    e = pd.read_csv(os.path.join(PROC, f"{c}_expr.csv"), index_col=0); m = pd.read_csv(os.path.join(PROC, f"{c}_meta.csv"), index_col=0); e.columns = e.columns.astype(str); m.index = m.index.astype(str)
    ms, nmod = module_scores(e, MOD); pr = program_scores(ms); s = pd.concat([pr, ms], axis=1); s["group"] = m.loc[s.index, "group"]; s["cohort"] = c; s["disease"] = dis; allscores.append(s)
    for k in ["C", "I", "R"] + list(ms.columns):
        x, y = s.loc[s.group == dis, k].dropna(), s.loc[s.group == "Control", k].dropna()
        g, lo, hi = hedges_g(x.values, y.values); eff.append(dict(cohort=c, disease=dis, measure=k, g=g, ci_low=lo, ci_high=hi, se=g_se(g, len(x), len(y)), p_mwu=stats.mannwhitneyu(x, y).pvalue, n_case=len(x), n_ctrl=len(y), n_genes=int(nmod.get(k, np.nan)) if k in nmod else np.nan))
S = pd.concat(allscores); S.to_csv(os.path.join(OUT, "bulk_patient_scores.csv")); E = pd.DataFrame(eff); E.to_csv(os.path.join(OUT, "bulk_effects_per_cohort.csv"), index=False)
rows = []
for k in E.measure.unique():
    r = {"measure": k}
    for dis in ["ATAA", "TAAD"]:
        s = E[(E.measure == k) & (E.disease == dis)]; y, se, p, i2 = rem(s.g, s.se); r.update({f"g_{dis}": y, f"se_{dis}": se, f"ci_low_{dis}": y - 1.96 * se, f"ci_high_{dis}": y + 1.96 * se, f"p_{dis}": p, f"I2_{dis}": i2, f"n_cohorts_{dis}": len(s), f"n_cohorts_same_dir_{dis}": int((np.sign(s.g) == np.sign(y)).sum())})
    r["diff_TAAD_minus_ATAA"] = r["g_TAAD"] - r["g_ATAA"]; r["se_diff"] = np.sqrt(r["se_TAAD"] ** 2 + r["se_ATAA"] ** 2); r["z_diff"] = r["diff_TAAD_minus_ATAA"] / r["se_diff"]; r["p_diff"] = 2 * stats.norm.sf(abs(r["z_diff"]))
    rows.append(r)
M = pd.DataFrame(rows).set_index("measure")
prim = M.loc[["C", "I", "R"]].copy(); prim["p_diff_holm"] = multipletests(prim.p_diff, method="holm")[1]; prim["family"] = "F1_primary"
expl = M.drop(index=["C", "I", "R"]).copy(); expl["p_diff_bh"] = multipletests(expl.p_diff, method="fdr_bh")[1]; expl["family"] = "exploratory"
M2 = pd.concat([prim, expl]); M2.to_csv(os.path.join(OUT, "bulk_pooled_and_disease_difference.csv"))
pd.set_option("display.width", 250)
print("F1 primary (cohort-level ATAA vs TAAD; Holm over C, I, R):"); print(prim[["g_ATAA", "ci_low_ATAA", "ci_high_ATAA", "g_TAAD", "ci_low_TAAD", "ci_high_TAAD", "diff_TAAD_minus_ATAA", "p_diff", "p_diff_holm"]].round(3).to_string())
print("\nExploratory Tier 1 modules (pooled g; BH on disease difference):"); print(expl[["g_ATAA", "p_ATAA", "n_cohorts_same_dir_ATAA", "g_TAAD", "p_TAAD", "n_cohorts_same_dir_TAAD", "diff_TAAD_minus_ATAA", "p_diff_bh"]].round(3).to_string())

# GSE318877 patient level (estimation only)
er = pd.read_csv(os.path.join(PROC, "GSE318877_regions_log2tpm.csv"), index_col=0); mr = pd.read_csv(os.path.join(PROC, "GSE318877_regions_meta.csv"), index_col=0)
er = er[~er.index.duplicated()]; er = er[(er > 1).mean(axis=1) > 0.25]
# patient-level expression = mean log2TPM over punches, then label-free module scores across the 7 patients
pe = er.T.groupby(mr.loc[er.columns, "subject"].astype(str)).mean().T; ms7, n7 = module_scores(pe, MOD); pr7 = program_scores(ms7)
pr7["group"] = mr.drop_duplicates("subject").set_index("subject")["group"].reindex(pr7.index.astype(int) if pr7.index.dtype != object else pr7.index).values if False else [mr.loc[mr.subject.astype(str) == s, "group"].iloc[0] for s in pr7.index]
pr7 = pd.concat([pr7, ms7], axis=1); pr7.to_csv(os.path.join(OUT, "GSE318877_patient_scores.csv"))
rows = []
for k in ["C", "I", "R"] + list(ms7.columns):
    x, y = pr7.loc[pr7.group == "TAAD", k], pr7.loc[pr7.group == "ATAA", k]; g, lo, hi = hedges_g(x.values, y.values); obs, p, na = exact_perm_test(x.values, y.values)
    rows.append(dict(measure=k, mean_dissection=x.mean(), mean_dilatation=y.mean(), diff=obs, g=g, ci_low=lo, ci_high=hi, p_exact=p, n_assignments=na, note="estimation only (3 vs 4; min attainable two-sided P=0.057)"))
D = pd.DataFrame(rows).set_index("measure"); D.to_csv(os.path.join(OUT, "GSE318877_patient_estimates.csv"))
print("\nGSE318877 patient level (subacute/chronic dissection n=3 vs dilatation n=4; estimation only):"); print(D.loc[["C", "I", "R", "SMC_contractile", "Calcium_handling", "Hypoxia", "Glycolysis", "Oxidative_stress_NFE2L2", "MYC_ribosome", "p53_DNA_damage", "NFkB_IL6_inflammation", "IFN_alpha_response", "MHCII_antigen_presentation", "ECM_collagen"]][["mean_dissection", "mean_dilatation", "g", "ci_low", "ci_high", "p_exact"]].round(3).to_string())
print("\nraw patient scores:"); print(pr7[["group", "C", "I", "R"]].round(2).to_string())
