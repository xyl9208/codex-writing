"""Within-cohort structure (v2.1).

B1  GSE140947: within-patient aneurysm belly minus neck (6 pairs) for Tier 1 programme scores C, I, R.
    PRIMARY (F3): delta R, exact sign-flip permutation (64 assignments).  Donor mid minus distal (6 pairs) as background.
B4  GSE318877: per patient, per medial layer (I lumen-side, L, M, A adventitia-side; replicates averaged) -> label-free
    module scores across the 28 patient-layer profiles -> within-patient gradient (A - I) of C, I, R; dissection vs
    dilatation difference of gradients: ESTIMATION ONLY (3 vs 4).
B2  GSE26155: ordered dilatation groups (No < Borderline < Yes; intima-media, TAV) -> Jonckheere-Terpstra trend of
    C, I, R (permutation P), estimation of group means.
"""
import os, sys, io, json, gzip, warnings, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from v21_lib import *
warnings.filterwarnings("ignore")
PROC = os.path.join(ROOT, "results", "processed"); GEO = os.path.join(ROOT, "data", "geo"); OUT = os.path.join(ROOT, "results", "v21_within"); os.makedirs(OUT, exist_ok=True)
MOD = load_tier1(); pd.set_option("display.width", 250)

# ---------------------------------------------------------------- B1 GSE140947 paired regions
cnt = pd.read_csv(os.path.join(PROC, "GSE140947_regions_counts.csv"), index_col=0); meta = pd.read_csv(os.path.join(PROC, "GSE140947_regions_meta.csv"), index_col=0)
cnt = cnt[~cnt.index.duplicated()]; lcpm = np.log2(cnt / cnt.sum(axis=0) * 1e6 + 1); lcpm = lcpm[(cnt >= 10).sum(axis=1) >= 6]
ms, n = module_scores(lcpm, MOD); pr = program_scores(ms); s = pd.concat([pr, ms], axis=1).join(meta[["group", "region", "subject"]]); s.to_csv(os.path.join(OUT, "GSE140947_region_scores.csv"))
rows = []
for grp, r1, r2, label in [("ATAA", "AN_BELLY", "AN_NECK", "aneurysm_belly_minus_neck"), ("Control", "MID_ASC", "DISTAL_ASC", "donor_mid_minus_distal")]:
    a = s[(s.group == grp) & (s.region == r1)].set_index("subject"); b = s[(s.group == grp) & (s.region == r2)].set_index("subject"); common = a.index.intersection(b.index)
    for k in ["C", "I", "R"] + list(ms.columns):
        d = (a.loc[common, k] - b.loc[common, k]).dropna().values
        if len(d) < 3: continue
        obs, p, na = exact_sign_perm_paired(d); ci = stats.t.interval(0.95, len(d) - 1, loc=d.mean(), scale=stats.sem(d)) if len(d) > 1 else (np.nan, np.nan)
        rows.append(dict(contrast=label, measure=k, n_pairs=len(d), mean_delta=obs, ci_low=ci[0], ci_high=ci[1], p_exact_signflip=p, n_assignments=na, primary=(label == "aneurysm_belly_minus_neck" and k == "R"), deltas=";".join(f"{v:.2f}" for v in d)))
B1 = pd.DataFrame(rows); B1.to_csv(os.path.join(OUT, "GSE140947_paired_tests.csv"), index=False)
print("B1 GSE140947 within-patient belly - neck (PRIMARY: delta R):"); print(B1[B1.measure.isin(["C", "I", "R", "SMC_contractile", "Hypoxia", "Glycolysis", "Oxidative_stress_NFE2L2", "MYC_ribosome", "p53_DNA_damage", "NFkB_IL6_inflammation"])][["contrast", "measure", "n_pairs", "mean_delta", "ci_low", "ci_high", "p_exact_signflip", "primary"]].round(3).to_string(index=False))

# ---------------------------------------------------------------- B4 GSE318877 layers within patient
er = pd.read_csv(os.path.join(PROC, "GSE318877_regions_log2tpm.csv"), index_col=0); mr = pd.read_csv(os.path.join(PROC, "GSE318877_regions_meta.csv"), index_col=0)
er = er[~er.index.duplicated()]; er = er[(er > 1).mean(axis=1) > 0.25]
key = mr["subject"].astype(str) + "|" + mr["region"]; pl = er.T.groupby(key.loc[er.columns]).mean().T  # patient-layer profiles (replicates averaged)
ms4, n4 = module_scores(pl, MOD); pr4 = program_scores(ms4); s4 = pd.concat([pr4, ms4], axis=1)
s4["subject"] = [c.split("|")[0] for c in s4.index]; s4["layer"] = [c.split("|")[1] for c in s4.index]; s4["group"] = [mr.loc[mr.subject.astype(str) == sb, "group"].iloc[0] for sb in s4.subject]
s4.to_csv(os.path.join(OUT, "GSE318877_patient_layer_scores.csv"))
order = {"I": 0, "L": 1, "M": 2, "A": 3}; rows = []
for k in ["C", "I", "R"] + list(ms4.columns):
    grad = {}
    for sb, d in s4.groupby("subject"):
        d = d.set_index("layer")
        if "A" in d.index and "I" in d.index: grad[sb] = d.loc["A", k] - d.loc["I", k]
    g_ = pd.Series(grad); grp = s4.drop_duplicates("subject").set_index("subject")["group"].reindex(g_.index)
    x, y = g_[grp == "TAAD"].values, g_[grp == "ATAA"].values
    gg, lo, hi = hedges_g(x, y); obs, p, na = exact_perm_test(x, y)
    rows.append(dict(measure=k, n_dissection=len(x), n_dilatation=len(y), mean_gradient_dissection=x.mean(), mean_gradient_dilatation=y.mean(), diff=obs, g=gg, ci_low=lo, ci_high=hi, p_exact_descriptive=p, gradients=";".join(f"{s}:{v:.2f}" for s, v in g_.items())))
B4 = pd.DataFrame(rows); B4.to_csv(os.path.join(OUT, "GSE318877_layer_gradient_estimates.csv"), index=False)
print("\nB4 GSE318877 within-patient adventitia-side minus lumen-side gradient (estimation only):"); print(B4[B4.measure.isin(["C", "I", "R", "SMC_contractile", "Glycolysis", "Oxidative_stress_NFE2L2", "MYC_ribosome", "NFkB_IL6_inflammation", "ECM_collagen"])][["measure", "mean_gradient_dissection", "mean_gradient_dilatation", "diff", "g", "ci_low", "ci_high", "p_exact_descriptive"]].round(3).to_string(index=False))
print(s4.pivot_table(index="subject", columns="layer", values="R").round(2).join(s4.drop_duplicates("subject").set_index("subject")["group"]).to_string())

# ---------------------------------------------------------------- B2 GSE26155 ordered dilatation (re-read matrix incl. borderline)
def read_series_matrix(path):
    meta = {}; rows_ = []; in_table = False; header = None
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"): in_table = True; continue
            if line.startswith("!series_matrix_table_end"): break
            if in_table:
                parts = [p.strip('"') for p in line.split("\t")]
                if header is None: header = parts
                else: rows_.append(parts)
            elif line.startswith("!Sample_"):
                k_, *vals = line.split("\t"); vals = [v.strip('"') for v in vals]; k_ = k_[1:]
                if k_ == "Sample_characteristics_ch1":
                    for i, v in enumerate(vals):
                        if ":" in v: kk, vv = v.split(":", 1); meta.setdefault("char:" + kk.strip(), [None] * len(vals))[i] = vv.strip()
                else: meta[k_] = vals
    e = pd.DataFrame(rows_, columns=header).set_index(header[0]).apply(pd.to_numeric, errors="coerce"); m = pd.DataFrame(meta); m.index = m["Sample_geo_accession"]; return m, e
m26, e26 = read_series_matrix(os.path.join(GEO, "GSE26155_series_matrix.txt.gz"))
keep = (m26["Sample_source_name_ch1"] == "Aorta intima-media"); m26 = m26[keep]; e26 = e26[m26.index]
lines = open(os.path.join(GEO, "GPL5175_full.txt"), errors="replace").read().split("\n"); start = next(i for i, l in enumerate(lines) if l.startswith("ID\t"))
ann = pd.read_csv(io.StringIO("\n".join(lines[start:])), sep="\t", dtype=str, quoting=3, on_bad_lines="skip"); ann = ann[ann["category"] == "main"]
def first_sym(ga):
    if not isinstance(ga, str) or ga == "---": return None
    s_ = [c.split(" // ")[1].strip() for c in ga.split(" /// ") if len(c.split(" // ")) > 1]; s_ = [x for x in s_ if x and x != "---"]; return s_[0] if s_ else None
ann["symbol"] = ann["gene_assignment"].map(first_sym); ann = ann.dropna(subset=["symbol"]).set_index("ID"); e26 = e26.loc[e26.index.intersection(ann.index)]
e26["sym"] = ann.loc[e26.index, "symbol"].values; e26["mean"] = e26.drop(columns="sym").mean(axis=1); e26 = e26.sort_values("mean", ascending=False).drop_duplicates("sym").set_index("sym").drop(columns="mean")
ms2, n2 = module_scores(e26, MOD); pr2 = program_scores(ms2); s2 = pd.concat([pr2, ms2], axis=1); s2["dilation"] = m26.loc[s2.index, "char:aorta dilation"].values
s2["ordinal"] = s2.dilation.map({"No": 0, "Borderline": 1, "Yes": 2}); s2.to_csv(os.path.join(OUT, "GSE26155_ordered_scores.csv"))
def jt_test(values, groups, n_perm=5000, seed=0):
    """Jonckheere-Terpstra statistic with permutation P (two-sided), groups ordinal 0<1<2."""
    v = np.asarray(values, float); g = np.asarray(groups); levels = sorted(set(g))
    def jt(vv):
        s_ = 0
        for i in range(len(levels)):
            for j in range(i + 1, len(levels)):
                a, b = vv[g == levels[i]], vv[g == levels[j]]; s_ += np.sum(a[:, None] < b[None, :]) + 0.5 * np.sum(a[:, None] == b[None, :])
        return s_
    obs = jt(v); rng_ = np.random.default_rng(seed); perm = np.array([jt(rng_.permutation(v)) for _ in range(n_perm)])
    p = min(1, 2 * min(np.mean(perm >= obs), np.mean(perm <= obs))); return obs, p
rows = []
for k in ["C", "I", "R"] + list(ms2.columns):
    d = s2[[k, "ordinal"]].dropna(); obs, p = jt_test(d[k].values, d["ordinal"].values)
    means = d.groupby("ordinal")[k].mean(); rows.append(dict(measure=k, n_No=int((d.ordinal == 0).sum()), n_Borderline=int((d.ordinal == 1).sum()), n_Yes=int((d.ordinal == 2).sum()), mean_No=means.get(0, np.nan), mean_Borderline=means.get(1, np.nan), mean_Yes=means.get(2, np.nan), rho_spearman=stats.spearmanr(d[k], d.ordinal)[0], JT=obs, p_perm=p))
B2 = pd.DataFrame(rows); B2.to_csv(os.path.join(OUT, "GSE26155_ordered_trend.csv"), index=False)
print("\nB2 GSE26155 ordered dilatation (No < Borderline < Yes; intima-media TAV):"); print(B2[B2.measure.isin(["C", "I", "R", "SMC_contractile", "Calcium_handling", "IFN_alpha_response", "MHCII_antigen_presentation", "ECM_collagen", "sig_T_cell", "Glycolysis", "NFkB_IL6_inflammation"])].round(3).to_string(index=False))
