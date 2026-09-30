"""Table 2: evidence matrix (estimation / inference / estimability separated), assembled from the v2.1 result files.

Rows = pre-specified comparisons; each row records direction, effect (Hedges' g or standardised difference) with 95 % CI,
effective n, whether the module is independent of module construction (Tier 1 = yes), the inference column (only for the
primary tests of the plan appendix), the estimation column (CI excludes / includes 0) and the estimability column.
Also applies BH to the exploratory ordered-trend results of GSE26155 (not done in 16_within_cohort.py).
"""
import os, sys, numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(ROOT, "results"); MS = os.path.join(ROOT, "manuscript")
def est(lo, hi):
    if np.isnan(lo) or np.isnan(hi): return "not estimated"
    return "CI excludes 0" if (lo > 0 or hi < 0) else "CI includes 0"
def fmt(g, lo, hi): return f"{g:+.2f} [{lo:.2f}, {hi:.2f}]" if not np.isnan(g) else "NA"
def direction(g): return "NA" if np.isnan(g) else ("higher" if g > 0 else "lower")
rows = []
# ---------------------------------------------------------------- F1 bulk (cohort-level comparison of diseases)
B = pd.read_csv(os.path.join(R, "v21_bulk", "bulk_pooled_and_disease_difference.csv"))
for _, r in B.iterrows():
    prim = r.family == "F1_primary"
    for dis in ["ATAA", "TAAD"]:
        rows.append(dict(block="1. Bulk cohorts (disease vs control; REM over cohorts)", dataset=f"{r[f'n_cohorts_{dis}']} {dis} cohorts", unit="patient within cohort; cohort-level pooling", comparison=f"{dis} vs control", measure=r.measure,
                         direction=direction(r[f"g_{dis}"]), effect=fmt(r[f"g_{dis}"], r[f"ci_low_{dis}"], r[f"ci_high_{dis}"]), n=f"{int(r[f'n_cohorts_{dis}'])} cohorts ({int(r[f'n_cohorts_same_dir_{dis}'])} same direction)",
                         tier1="yes", inference="", estimation=est(r[f"ci_low_{dis}"], r[f"ci_high_{dis}"]), estimability="credible", note=f"I2 = {r[f'I2_{dis}']:.0f} %"))
    lo, hi = r.diff_TAAD_minus_ATAA - 1.96 * r.se_diff, r.diff_TAAD_minus_ATAA + 1.96 * r.se_diff
    rows.append(dict(block="1. Bulk cohorts (disease vs control; REM over cohorts)", dataset="10 cohorts", unit="cohort-level difference of pooled g", comparison="TAAD effect minus ATAA effect", measure=r.measure,
                     direction=direction(r.diff_TAAD_minus_ATAA), effect=fmt(r.diff_TAAD_minus_ATAA, lo, hi), n="4 vs 6 cohorts", tier1="yes",
                     inference=(("passed" if r.p_diff_holm < 0.05 else "not passed") + f" (Holm P = {r.p_diff_holm:.3f})") if prim else f"exploratory (BH P = {r.p_diff_bh:.3f})",
                     estimation=est(lo, hi), estimability="credible", note="F1 primary family (C, I, R)" if prim else ""))
# ---------------------------------------------------------------- F2 single-cell patient level
S = pd.read_csv(os.path.join(R, "v21_sc", "sc_patient_level_tests.csv"))
for ds, dis in [("GSE155468", "ATAA"), ("GSE213740", "TAAD")]:
    d = S[(S.dataset == ds) & (S.celltype == "SMC_contractile")]
    for k in ["C", "I", "R", "Hypoxia", "Glycolysis", "Oxidative_stress_NFE2L2", "MYC_ribosome", "p53_DNA_damage", "NFkB_IL6_inflammation", "IFN_alpha_response", "MHCII_antigen_presentation", "ECM_collagen"]:
        r = d[d.measure == k].iloc[0]
        rows.append(dict(block="2. Single-cell, contractile SMC pseudo-bulk (patient level)", dataset=ds, unit="patient (>=20 SMC)", comparison=f"{dis} vs control", measure=k, direction=direction(r.g), effect=fmt(r.g, r.ci_low, r.ci_high),
                         n=f"{int(r.n_case)} vs {int(r.n_ctrl)}", tier1="yes", inference=(("passed" if r.p_exact < 0.05 else "not passed") + f" (exact P = {r.p_exact:.3f}, {int(r.n_assignments)} assignments)") if r.primary else f"exploratory (exact P = {r.p_exact:.3f}; BH = {r.p_bh_exploratory:.2f})",
                         estimation=est(r.ci_low, r.ci_high), estimability="limited precision (n small)", note="F2 primary" if r.primary else ""))
    for ct in ["Fibroblast", "Endothelial", "Macrophage", "T_cell"]:
        for k in ["C", "I", "R"]:
            q = S[(S.dataset == ds) & (S.celltype == ct) & (S.measure == k)]
            if len(q) == 0: continue
            r = q.iloc[0]
            rows.append(dict(block="2. Single-cell, other cell types (patient level)", dataset=ds, unit=f"patient (>=20 {ct})", comparison=f"{dis} vs control", measure=k, direction=direction(r.g), effect=fmt(r.g, r.ci_low, r.ci_high),
                             n=f"{int(r.n_case)} vs {int(r.n_ctrl)}", tier1="yes", inference=f"exploratory (exact P = {r.p_exact:.3f})", estimation=est(r.ci_low, r.ci_high), estimability="limited precision (n small)", note=""))
# GSE189795 (acute external) if available
f = os.path.join(R, "v21_stage", "GSE189795_tests.csv")
if os.path.exists(f):
    E2 = pd.read_csv(f); d = E2[E2.celltype == "SMC_contractile"]
    for k in ["C", "I", "R", "Hypoxia", "Glycolysis", "Oxidative_stress_NFE2L2", "MYC_ribosome", "p53_DNA_damage", "NFkB_IL6_inflammation"]:
        q = d[d.measure == k]
        if len(q) == 0: continue
        r = q.iloc[0]
        rows.append(dict(block="3. Stage: acute dissection (external replicate)", dataset="GSE189795", unit="patient (>=20 SMC)", comparison="acute TAAD vs control", measure=k, direction=direction(r.g), effect=fmt(r.g, r.ci_low, r.ci_high),
                         n=f"{int(r.n_acute)} vs {int(r.n_ctrl)}", tier1="yes", inference=(("passed" if r.p_exact < 0.05 else "not passed") + f" (exact P = {r.p_exact:.3f}, {int(r.n_assignments)} assignments)") if r.primary else f"exploratory (exact P = {r.p_exact:.3f}; BH = {r.p_bh_exploratory:.2f})",
                         estimation=est(r.ci_low, r.ci_high), estimability="limited precision (n small)", note="F2 primary" if r.primary else ""))
f = os.path.join(R, "v21_stage", "GSE222318_stage_estimates.csv")
if os.path.exists(f):
    E1 = pd.read_csv(f); d = E1[E1.celltype == "SMC_contractile"]
    for con in ["acute_vs_control", "subacute_vs_control", "chronic_vs_control", "acute_vs_nonacute", "acute_vs_subacute_samePosition"]:
        for k in ["C", "I", "R"]:
            q = d[(d.contrast == con) & (d.measure == k)]
            if len(q) == 0: continue
            r = q.iloc[0]
            rows.append(dict(block="3. Stage: within-cohort stages (estimation only)", dataset="GSE222318", unit="patient (>=20 SMC)", comparison=con.replace("_", " "), measure=k, direction=direction(r.g), effect=fmt(r.g, r.ci_low, r.ci_high),
                             n=f"{int(r.n_A)} vs {int(r.n_B)} (positions {r.positions_A} vs {r.positions_B})", tier1="yes", inference="none (Position confounded with stage)", estimation=est(r.ci_low, r.ci_high),
                             estimability="limited precision; sampling position differs between stages" if con not in ("acute_vs_subacute_samePosition",) else "limited precision (same sampling position)", note=""))
# ---------------------------------------------------------------- F3 within-patient regions
W = pd.read_csv(os.path.join(R, "v21_within", "GSE140947_paired_tests.csv"))
for con in ["aneurysm_belly_minus_neck", "donor_mid_minus_distal"]:
    for k in ["C", "I", "R"]:
        r = W[(W.contrast == con) & (W.measure == k)].iloc[0]
        rows.append(dict(block="4. Within-patient region (paired)", dataset="GSE140947", unit="within-patient difference", comparison=con.replace("_", " "), measure=k, direction=direction(r.mean_delta), effect=fmt(r.mean_delta, r.ci_low, r.ci_high),
                         n=f"{int(r.n_pairs)} pairs", tier1="yes", inference=(("passed" if r.p_exact_signflip < 0.05 else "not passed") + f" (exact sign-flip P = {r.p_exact_signflip:.3f}, 64 assignments)") if r.primary else f"exploratory (exact sign-flip P = {r.p_exact_signflip:.3f})",
                         estimation=est(r.ci_low, r.ci_high), estimability="limited precision (6 pairs)", note="F3 primary" if r.primary else ""))
# ---------------------------------------------------------------- GSE318877 patient level and layer gradients (estimation only)
P = pd.read_csv(os.path.join(R, "v21_bulk", "GSE318877_patient_estimates.csv"))
for k in ["C", "I", "R"]:
    r = P[P.measure == k].iloc[0]
    rows.append(dict(block="5. Direct dissection-vs-dilatation cohort (estimation only)", dataset="GSE318877", unit="patient (mean of medial punches)", comparison="dissection vs dilatation", measure=k, direction=direction(r.g), effect=fmt(r.g, r.ci_low, r.ci_high),
                     n="3 vs 4", tier1="yes", inference="none (3 vs 4; minimum attainable two-sided P = 0.057)", estimation=est(r.ci_low, r.ci_high), estimability="limited precision (n = 7)", note=f"descriptive exact P = {r.p_exact:.2f}"))
L = pd.read_csv(os.path.join(R, "v21_within", "GSE318877_layer_gradient_estimates.csv"))
for k in ["C", "I", "R"]:
    r = L[L.measure == k].iloc[0]
    rows.append(dict(block="5. Direct dissection-vs-dilatation cohort (estimation only)", dataset="GSE318877", unit="within-patient layer gradient (adventitial minus luminal)", comparison="gradient: dissection vs dilatation", measure=k, direction=direction(r["diff"]), effect=fmt(r.g, r.ci_low, r.ci_high),
                     n="3 vs 4", tier1="yes", inference="none", estimation=est(r.ci_low, r.ci_high), estimability="limited precision (n = 7)", note=f"mean gradient dissection {r.mean_gradient_dissection:+.2f}, dilatation {r.mean_gradient_dilatation:+.2f}"))
# ---------------------------------------------------------------- GSE26155 ordered dilatation (exploratory; BH)
O = pd.read_csv(os.path.join(R, "v21_within", "GSE26155_ordered_trend.csv")); O["p_bh"] = multipletests(O.p_perm, method="fdr_bh")[1]; O.to_csv(os.path.join(R, "v21_within", "GSE26155_ordered_trend_bh.csv"), index=False)
for k in ["C", "I", "R", "Calcium_handling", "ECM_collagen", "Matrix_degradation", "IFN_alpha_response", "IFNG_specific", "MHCII_antigen_presentation", "NFkB_IL6_inflammation", "sig_T_cell", "sig_Macrophage", "sig_B_cell"]:
    r = O[O.measure == k].iloc[0]
    rows.append(dict(block="6. Ordered dilatation within one cohort (exploratory)", dataset="GSE26155", unit="patient (intima-media, TAV)", comparison="No < Borderline < Yes dilatation (trend)", measure=k, direction=("increasing" if r.rho_spearman > 0 else "decreasing"),
                     effect=f"Spearman rho = {r.rho_spearman:+.2f}; means {r.mean_No:+.2f} / {r.mean_Borderline:+.2f} / {r.mean_Yes:+.2f}", n=f"{int(r.n_No)} / {int(r.n_Borderline)} / {int(r.n_Yes)}", tier1="yes",
                     inference=f"exploratory (JT permutation P = {r.p_perm:.3f}; BH = {r.p_bh:.3f})", estimation="", estimability="credible", note=""))
# ---------------------------------------------------------------- conditioned sensitivity
Cn = pd.read_csv(os.path.join(R, "v21_conditioned", "conditioned_sensitivity_pooled.csv"))
for dis in ["ATAA", "TAAD"]:
    for k in ["C", "I", "R", "Glycolysis", "Oxidative_stress_NFE2L2", "MYC_ribosome", "p53_DNA_damage", "NFkB_IL6_inflammation", "IFN_alpha_response", "MHCII_antigen_presentation", "ECM_collagen"]:
        r = Cn[(Cn.disease == dis) & (Cn.measure == k)].iloc[0]
        rows.append(dict(block="7. Composition-signature-conditioned sensitivity (bulk)", dataset=f"{int(r.n_cohorts)} {dis} cohorts", unit="patient within cohort; REM over cohorts", comparison=f"{dis} vs control, conditioned on 3 cell-signature scores", measure=k,
                         direction=direction(r.pooled_beta_conditioned), effect=f"beta_unc {r.pooled_beta_unconditioned:+.2f} [{r.ci_unc}] -> beta_cond {r.pooled_beta_conditioned:+.2f} [{r.ci_cond}]", n=f"{int(r.n_cohorts)} cohorts", tier1="yes (covariate-marker genes removed)",
                         inference="sensitivity", estimation=("retained" if r.retained == "yes" else "attenuated / uncertain"), estimability="credible", note=""))
T = pd.DataFrame(rows); T.to_csv(os.path.join(MS, "Table2_evidence_matrix.csv"), index=False)
# short main-text version (C, I, R only)
short = T[T.measure.isin(["C", "I", "R"]) & ~T.block.str.startswith("7.")].copy()
short.to_csv(os.path.join(MS, "Table2_evidence_matrix_main.csv"), index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_colwidth", 60)
print(short[["block", "dataset", "comparison", "measure", "effect", "n", "inference", "estimation", "estimability"]].to_string(index=False))
print(f"\nfull table: {len(T)} rows; main table: {len(short)} rows")
