"""Strict implementation check of the v2 random-effects / HKSJ code against statsmodels with IDENTICAL inputs
(per-cohort log2FC and SE as used by 14_meta_rem.py, same missing handling)."""
import os, json, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.meta_analysis import combine_effects
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DE = os.path.join(ROOT, "results", "de"); M2 = os.path.join(ROOT, "results", "meta_v2")
REG = json.load(open(os.path.join(ROOT, "results", "processed", "cohorts.json"))); DESIGN = json.load(open(os.path.join(ROOT, "results", "meta", "design.json")))
rng = np.random.default_rng(1); rows = []
for dis in ["ATAA", "TAAD"]:
    m = pd.read_csv(os.path.join(M2, f"{dis}_meta_v2.csv"), index_col=0)
    des = {c: pd.read_csv(os.path.join(DE, f"{c}_de.csv"), index_col=0) for c in DESIGN[dis]["cohorts"]}
    des = {c: d[~d.index.duplicated()].dropna(subset=["pvalue", "log2FC"]) for c, d in des.items()}  # identical filtering to 14_meta_rem.load()
    genes = list(m.index[m.consensus_v2][:40]) + list(rng.choice(m.index, 60, replace=False))
    for g in genes:
        y, v = [], []
        for c, d in des.items():
            if g in d.index and np.isfinite(d.loc[g, "se"]) and d.loc[g, "se"] != 0 and np.isfinite(d.loc[g, "log2FC"]):
                y.append(d.loc[g, "log2FC"]); v.append(abs(d.loc[g, "se"]) ** 2)
        if len(y) < 2: continue
        r = combine_effects(np.array(y), np.array(v), method_re="dl", use_t=True); s = r.summary_frame(); re = s.loc["random effect"]
        k = len(y); p_re = 2 * stats.norm.sf(abs(re["eff"] / re["sd_eff"]))
        # KH: statsmodels applies Knapp-Hartung when use_t=True to the random-effect row's sd
        sd_hk = float(np.sqrt(np.atleast_1d(r.var_hksj_re)[0])) if hasattr(r, 'var_hksj_re') else float(getattr(r, 'sd_eff_w_re_hksj')); p_kh = 2 * stats.t.sf(abs(re["eff"] / sd_hk), df=k - 1)
        rows.append(dict(disease=dis, gene=g, k=k, lfc_ours=m.loc[g, "lfc_rem"], lfc_sm=re["eff"], se_ours=m.loc[g, "se_rem"], sd_sm=re["sd_eff"], tau2_ours=m.loc[g, "tau2"], tau2_sm=r.tau2,
                         p_rem_ours=m.loc[g, "p_rem"], p_hksj_ours=m.loc[g, "p_hksj"], p_kh_sm=p_kh))
df = pd.DataFrame(rows); df["d_lfc"] = (df.lfc_ours - df.lfc_sm).abs(); df["d_tau2"] = (df.tau2_ours - df.tau2_sm).abs(); df["ratio_p_hksj"] = df.p_hksj_ours / df.p_kh_sm
df.to_csv(os.path.join(M2, "hksj_check.csv"), index=False)
print(f"n genes checked={len(df)}; max |lfc diff|={df.d_lfc.max():.2e}; max |tau2 diff|={df.d_tau2.max():.2e}; HKSJ p ratio ours/statsmodels: min={df.ratio_p_hksj.min():.3f} median={df.ratio_p_hksj.median():.3f} max={df.ratio_p_hksj.max():.3f}")
print(df.sort_values("ratio_p_hksj").head(3).round(4).to_string(index=False)); print(df.sort_values("ratio_p_hksj").tail(3).round(4).to_string(index=False))
