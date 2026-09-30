"""Implementation check of the v2 random-effects / HKSJ code (14_meta_rem.py) against statsmodels combine_effects, with
IDENTICAL inputs (per-cohort log2FC and SE, same missing handling) for ALL genes in {ATAA,TAAD}_meta_v2.csv.

Convention difference: statsmodels leaves a negative DerSimonian-Laird tau2 untruncated, whereas 14_meta_rem.py truncates it
at 0 (tau2.clip(lower=0), the standard DL convention). Two comparisons are therefore reported for every gene:
  raw   = statsmodels output as is (untruncated tau2; can give negative/inf weights when v_i + tau2 <= 0)
  trunc = statsmodels DL tau2 truncated at 0, then RE weights 1/(v_i+tau2), pooled effect, SE and the Knapp-Hartung variance
          q = (1/(k-1)) * sum(w_i (y_i - yhat)^2) / sum(w_i) recomputed from statsmodels' own eff/var_eff inputs, i.e. exactly
          combine_effects' formulas (mean_effect_re, sd_eff_w_re, var_hksj_re) with tau2 replaced by max(0, tau2).
Genes where truncation was applied (tau2_sm < 0) are flagged. Ours: p_rem = 2*norm.sf(|lfc/se|); p_hksj = 2*t.sf(|lfc/se_hk|, k-1)
with se_hk = sqrt(q)*se_rem (unmodified HKSJ; the 'modified' variant p_hksj_mod is not checked here). Main results are not touched."""
import os, json, warnings, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.meta_analysis import combine_effects
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DE = os.path.join(ROOT, "results", "de"); M2 = os.path.join(ROOT, "results", "meta_v2")
DESIGN = json.load(open(os.path.join(ROOT, "results", "meta", "design.json")))

def load(name):  # identical filtering to 14_meta_rem.load()
    d = pd.read_csv(os.path.join(DE, f"{name}_de.csv"), index_col=0)
    d = d[~d.index.duplicated()].dropna(subset=["pvalue", "log2FC"]); d["se"] = d["se"].abs().replace(0, np.nan); return d

def p_norm(eff, se): return 2 * stats.norm.sf(np.abs(eff / se))
def p_t(eff, var, k): return 2 * stats.t.sf(np.abs(eff / np.sqrt(var)), df=k - 1)

rows = []
for dis in ["ATAA", "TAAD"]:
    m = pd.read_csv(os.path.join(M2, f"{dis}_meta_v2.csv"), index_col=0); des = [load(c) for c in DESIGN[dis]["cohorts"]]
    Y = np.column_stack([d["log2FC"].reindex(m.index).values for d in des]); V = np.column_stack([(d["se"] ** 2).reindex(m.index).values for d in des])
    present = np.isfinite(Y) & np.isfinite(V)  # same as Y.notna() & V.notna() & isfinite(V) in rem_meta()
    for i, g in enumerate(m.index):
        y, v = Y[i, present[i]], V[i, present[i]]; k = len(y)
        r = combine_effects(y, v, method_re="dl", use_t=True)
        # (b) raw statsmodels values, untruncated tau2
        lfc_raw, se_raw, var_kh_raw = float(r.mean_effect_re), float(r.sd_eff_w_re), float(r.var_hksj_re)
        # (a) truncation-consistent recomputation with statsmodels' own inputs (r.eff, r.var_eff) and its formulas
        tau2_t = max(0.0, float(r.tau2)); w = 1 / (r.var_eff + tau2_t); wsum = w.sum(); lfc_t = (w * r.eff).sum() / wsum; se_t = np.sqrt(1 / wsum)
        var_kh_t = (w * (r.eff - lfc_t) ** 2).sum() / ((k - 1) * wsum)
        rows.append(dict(disease=dis, gene=g, k=k, tau2_ours=m.at[g, "tau2"], tau2_sm=float(r.tau2), truncated=bool(r.tau2 < 0),
                         lfc_ours=m.at[g, "lfc_rem"], se_ours=m.at[g, "se_rem"], p_rem_ours=m.at[g, "p_rem"], p_hksj_ours=m.at[g, "p_hksj"],
                         lfc_sm_trunc=lfc_t, se_sm_trunc=se_t, p_rem_sm_trunc=p_norm(lfc_t, se_t), p_hksj_sm_trunc=p_t(lfc_t, var_kh_t, k),
                         lfc_sm_raw=lfc_raw, se_sm_raw=se_raw, p_rem_sm_raw=p_norm(lfc_raw, se_raw), p_hksj_sm_raw=p_t(lfc_raw, var_kh_raw, k) if var_kh_raw > 0 else np.nan,
                         var_hksj_sm_raw=var_kh_raw))
df = pd.DataFrame(rows)
for tag in ["trunc", "raw"]:
    for q in ["lfc", "se", "p_rem", "p_hksj"]: df[f"d_{q}_{tag}"] = (df[f"{q}_ours"] - df[f"{q}_sm_{tag}"]).abs()
    df[f"ratio_p_hksj_{tag}"] = df.p_hksj_ours / df[f"p_hksj_sm_{tag}"]; df[f"ratio_p_rem_{tag}"] = df.p_rem_ours / df[f"p_rem_sm_{tag}"]
df["d_tau2_trunc"] = (df.tau2_ours - df.tau2_sm.clip(lower=0)).abs(); df["d_tau2_raw"] = (df.tau2_ours - df.tau2_sm).abs()
df.to_csv(os.path.join(M2, "hksj_check.csv"), index=False)

def fold(r): return np.maximum(r, 1 / r)  # max fold discrepancy of a P ratio, either direction
def block(d, tag):
    ours = d[["lfc_ours", "se_ours", "p_rem_ours", "p_hksj_ours"]].values; sm = d[[f"lfc_sm_{tag}", f"se_sm_{tag}", f"p_rem_sm_{tag}", f"p_hksj_sm_{tag}"]].values
    ok = np.isfinite(ours).all(1) & np.isfinite(sm).all(1)
    fr, fh = fold(d[f"ratio_p_rem_{tag}"]), fold(d[f"ratio_p_hksj_{tag}"]); fr_ok, fh_ok = np.isfinite(fr), np.isfinite(fh)
    return dict(n_genes=int(len(d)), n_truncation_applied=int(d.truncated.sum()), n_nonfinite_ours=int((~np.isfinite(ours).all(1)).sum()), n_nonfinite_statsmodels=int((~np.isfinite(sm).all(1)).sum()),
                n_statsmodels_se_nonfinite=int((~np.isfinite(d[f"se_sm_{tag}"])).sum()), n_statsmodels_p_hksj_nonfinite=int((~np.isfinite(d[f"p_hksj_sm_{tag}"])).sum()),
                n_compared_both_finite=int(ok.sum()), max_abs_d_tau2=float(np.nanmax(d[f"d_tau2_{tag}"])), max_abs_d_lfc=float(np.nanmax(d[f"d_lfc_{tag}"])), max_abs_d_se=float(np.nanmax(d[f"d_se_{tag}"])),
                max_abs_d_p_rem=float(np.nanmax(d[f"d_p_rem_{tag}"])), max_abs_d_p_hksj=float(np.nanmax(d[f"d_p_hksj_{tag}"])),
                max_p_rem_fold_ratio_finite=float(fr[fr_ok].max()), n_p_rem_ratio_nonfinite=int((~fr_ok).sum()), max_p_hksj_fold_ratio_finite=float(fh[fh_ok].max()), n_p_hksj_ratio_nonfinite=int((~fh_ok).sum()),
                p_hksj_ratio_min=float(np.nanmin(d[f"ratio_p_hksj_{tag}"])), p_hksj_ratio_max=float(np.nanmax(d[f"ratio_p_hksj_{tag}"])),
                n_genes_all_abs_diff_below_1e_minus_8=int((ok & (d[[f"d_lfc_{tag}", f"d_se_{tag}", f"d_p_rem_{tag}", f"d_p_hksj_{tag}"]].values < 1e-8).all(1)).sum()),
                n_genes_p_hksj_fold_ratio_above_1p01=int((fh > 1.01).sum()))
summ = {}
for dis in ["ATAA", "TAAD", "ALL"]:
    d = df if dis == "ALL" else df[df.disease == dis]
    summ[dis] = dict(truncation_consistent=block(d, "trunc"), raw_statsmodels_untruncated=block(d, "raw"))
    summ[dis]["truncation_consistent"]["n_truncated_genes_with_p_hksj_fold_ratio_above_1p01"] = int((fold(d.ratio_p_hksj_trunc)[d.truncated] > 1.01).sum())
    summ[dis]["raw_statsmodels_untruncated"]["n_truncated_genes_with_p_hksj_fold_ratio_above_1p01"] = int((fold(d.ratio_p_hksj_raw)[d.truncated] > 1.01).sum())
summ["note"] = ("trunc: statsmodels DL tau2 truncated at 0 and RE effect/SE/Knapp-Hartung variance recomputed with combine_effects' formulas from its own inputs; "
                "raw: statsmodels output with untruncated tau2. ours = 14_meta_rem.py (tau2.clip(lower=0)); p_hksj is the unmodified HKSJ variant. All genes of {dis}_meta_v2.csv checked.")
json.dump(summ, open(os.path.join(M2, "hksj_check_summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
a, b = summ["ALL"]["truncation_consistent"], summ["ALL"]["raw_statsmodels_untruncated"]
print(f"\nALL: n={a['n_genes']} genes checked; tau2_sm<0 (truncation applied) in {a['n_truncation_applied']}. Truncation-consistent: max|d_lfc|={a['max_abs_d_lfc']:.1e}, max|d_se|={a['max_abs_d_se']:.1e}, "
      f"max|d_p_rem|={a['max_abs_d_p_rem']:.1e}, max|d_p_hksj|={a['max_abs_d_p_hksj']:.1e}, max P_hksj fold ratio={a['max_p_hksj_fold_ratio_finite']:.12f}, non-finite ours/sm={a['n_nonfinite_ours']}/{a['n_nonfinite_statsmodels']}. "
      f"Raw untruncated: max|d_lfc|={b['max_abs_d_lfc']:.3g}, max|d_se|={b['max_abs_d_se']:.3g}, max|d_p_rem|={b['max_abs_d_p_rem']:.3g}, max|d_p_hksj|={b['max_abs_d_p_hksj']:.3g}, "
      f"max finite P_hksj fold ratio={b['max_p_hksj_fold_ratio_finite']:.3g}, non-finite ours/sm={b['n_nonfinite_ours']}/{b['n_nonfinite_statsmodels']}, genes with P_hksj fold ratio>1.01: {b['n_genes_p_hksj_fold_ratio_above_1p01']} ({b['n_truncated_genes_with_p_hksj_fold_ratio_above_1p01']} of them with tau2_sm<0).")
ex = df[df.truncated].sort_values("ratio_p_hksj_raw").head(5)
print("\nExamples of genes with negative statsmodels tau2 (truncation applied), most discrepant raw HKSJ P first:")
print(ex[["disease", "gene", "k", "tau2_sm", "tau2_ours", "p_hksj_ours", "p_hksj_sm_trunc", "p_hksj_sm_raw", "ratio_p_hksj_trunc", "ratio_p_hksj_raw"]].to_string(index=False))
