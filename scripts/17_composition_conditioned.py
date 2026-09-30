"""A7: composition-signature-conditioned SENSITIVITY analysis (not a decomposition).

For each bulk cohort, the disease effect on each Tier 1 programme score (C, I, R and the individual modules) is
re-estimated in a linear model that conditions on three pre-specified cell-signature scores
(sig_SMC_contractile, sig_Inflammatory_myeloid, sig_T_cell; label-free scores from control-cell markers).
Rules: (i) genes shared between the outcome module and any of the three covariate marker sets are removed from the
outcome module before scoring; (ii) if fewer than 5 genes remain the result is recorded as 'not estimable';
(iii) unconditioned and conditioned effects (standardised coefficients with 95 % CI) are reported side by side;
(iv) with <= 3 covariates and n between 8 and 53 per cohort, degrees of freedom are reported.
"""
import os, sys, json, numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from v21_lib import *
PROC = os.path.join(ROOT, "results", "processed"); OUT = os.path.join(ROOT, "results", "v21_conditioned"); os.makedirs(OUT, exist_ok=True)
REG = json.load(open(os.path.join(PROC, "cohorts.json"))); DESIGN = json.load(open(os.path.join(ROOT, "results", "meta", "design.json")))
MOD = load_tier1(); COV = ["sig_SMC_contractile", "sig_Inflammatory_myeloid", "sig_T_cell"]; covgenes = set().union(*[set(MOD[c]) for c in COV])
cohorts = {c: "ATAA" for c in DESIGN["ATAA"]["cohorts"]}; cohorts.update({c: "TAAD" for c in DESIGN["TAAD"]["cohorts"]})
rows = []
for c, dis in cohorts.items():
    e = pd.read_csv(os.path.join(PROC, f"{c}_expr.csv"), index_col=0); m = pd.read_csv(os.path.join(PROC, f"{c}_meta.csv"), index_col=0); e.columns = e.columns.astype(str); m.index = m.index.astype(str)
    cov_scores, _ = module_scores(e, {k: MOD[k] for k in COV}); cov_scores = (cov_scores - cov_scores.mean()) / cov_scores.std()
    # outcome modules with covariate-marker genes removed
    mod_clean = {k: [g for g in v if g not in covgenes] for k, v in MOD.items() if not k.startswith("sig_")}
    ms_full, n_full = module_scores(e, {k: MOD[k] for k in mod_clean}); ms_clean, n_clean = module_scores(e, mod_clean)
    pr_full = program_scores(ms_full); pr_clean = program_scores(ms_clean)
    y_group = (m.loc[e.columns, "group"] == dis).astype(float)
    for k in ["C", "I", "R"] + list(mod_clean):
        yf = pr_full[k] if k in ("C", "I", "R") else (ms_full[k] if k in ms_full else None); yc = pr_clean[k] if k in ("C", "I", "R") else (ms_clean[k] if k in ms_clean else None)
        if yf is None or yc is None or yc.isna().all():
            rows.append(dict(cohort=c, disease=dis, measure=k, status="not estimable (<5 genes after removing covariate-marker genes)")); continue
        n_removed = (n_full.get(k, np.nan) - n_clean.get(k, np.nan)) if k in n_full else np.nan
        yz = (yc - yc.mean()) / yc.std()
        X0 = sm.add_constant(y_group.values); r0 = sm.OLS(yz.values, X0).fit()
        X1 = sm.add_constant(np.column_stack([y_group.values, cov_scores.loc[yz.index].values])); r1 = sm.OLS(yz.values, X1).fit()
        rows.append(dict(cohort=c, disease=dis, measure=k, status="ok", n=len(yz), df_resid_conditioned=int(r1.df_resid), n_genes_used=int(n_clean.get(k, 0)) if k in n_clean else np.nan, n_genes_removed=n_removed,
                         beta_unconditioned=r0.params[1], ci_low_unc=r0.conf_int()[1][0], ci_high_unc=r0.conf_int()[1][1], beta_conditioned=r1.params[1], ci_low_cond=r1.conf_int()[1][0], ci_high_cond=r1.conf_int()[1][1],
                         p_unconditioned=r0.pvalues[1], p_conditioned=r1.pvalues[1], max_vif_proxy=float(np.max(np.abs(np.corrcoef(X1[:, 1:].T)[0, 1:])))))
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "conditioned_sensitivity_per_cohort.csv"), index=False)
ok = R[R.status == "ok"]
def rem(g, se):
    g, se = np.asarray(g, float), np.asarray(se, float); w = 1 / se ** 2; yf = np.sum(w * g) / np.sum(w); Q = np.sum(w * (g - yf) ** 2); df = len(g) - 1; C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0, (Q - df) / C) if C > 0 else 0; wr = 1 / (se ** 2 + tau2); y = np.sum(wr * g) / np.sum(wr); s = np.sqrt(1 / np.sum(wr)); return y, s
srows = []
for (dis, k), d in ok.groupby(["disease", "measure"]):
    se_u = (d.ci_high_unc - d.ci_low_unc) / 3.92; se_c = (d.ci_high_cond - d.ci_low_cond) / 3.92
    yu, su = rem(d.beta_unconditioned, se_u); yc_, sc_ = rem(d.beta_conditioned, se_c)
    srows.append(dict(disease=dis, measure=k, n_cohorts=len(d), pooled_beta_unconditioned=yu, ci_unc=f"{yu-1.96*su:.2f},{yu+1.96*su:.2f}", pooled_beta_conditioned=yc_, ci_cond=f"{yc_-1.96*sc_:.2f},{yc_+1.96*sc_:.2f}", retained=("yes" if np.sign(yc_) == np.sign(yu) and abs(yc_ / sc_) > 1.96 else "attenuated/uncertain")))
S = pd.DataFrame(srows); S.to_csv(os.path.join(OUT, "conditioned_sensitivity_pooled.csv"), index=False)
pd.set_option("display.width", 250); print(S[S.measure.isin(["C", "I", "R", "SMC_contractile", "Calcium_handling", "Glycolysis", "Oxidative_stress_NFE2L2", "MYC_ribosome", "p53_DNA_damage", "NFkB_IL6_inflammation", "IFN_alpha_response", "MHCII_antigen_presentation", "ECM_collagen"])].round(3).to_string(index=False))
print("\nnot estimable:", R[R.status != "ok"][["cohort", "measure"]].values.tolist()[:10])
