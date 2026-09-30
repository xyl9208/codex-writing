"""F2 (single-cell, patient level): Tier 1 programme scores C, I, R within each cell type, computed from patient x cell-type
pseudo-bulk (mean log-normalised expression, >= 20 cells), label-free standardisation across ALL patients of the dataset
within the cell type.  Primary inference (one per dataset): R in SMC_contractile, disease vs control, exact enumeration
(GSE155468 8 vs 3 -> 165 assignments; GSE213740 6 vs 3 -> 84).  Everything else is estimation / exploratory (BH).
Also: control-referenced version with full-pipeline permutation as a check."""
import os, sys, itertools, warnings, numpy as np, pandas as pd, scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from v21_lib import *
warnings.filterwarnings("ignore")
SCD = os.path.join(ROOT, "results", "sc"); OUT = os.path.join(ROOT, "results", "v21_sc"); os.makedirs(OUT, exist_ok=True)
MOD = load_tier1()
a = sc.read_h5ad(os.path.join(SCD, "aorta_sc_integrated.h5ad"))
pb, ncells = celltype_pseudobulk(a, min_cells=20)
cols = pd.DataFrame([c.split("|") for c in pb.columns], columns=["sample", "celltype"], index=pb.columns)
smeta = a.obs.drop_duplicates("sample").set_index("sample")[["dataset", "group"]]; cols = cols.join(smeta, on="sample"); cols["n_cells"] = ncells
pb.to_csv(os.path.join(OUT, "pseudobulk_mean_lognorm.csv")); cols.to_csv(os.path.join(OUT, "pseudobulk_index.csv"))
CTS = ["SMC_contractile", "SMC_modulated", "Fibromyocyte", "Fibroblast", "Endothelial", "Macrophage", "Inflammatory_myeloid", "T_cell", "NK", "B_cell", "Plasma", "Mast"]

def scores_for(ds, ct):
    cc = cols.index[(cols.dataset == ds) & (cols.celltype == ct)]
    if len(cc) < 4: return None
    x = pb[cc]; x = x[(x > 0).mean(axis=1) >= 0.2]  # label-free detectability filter within cell type
    ms, n = module_scores(x, MOD, detect_frac=None); pr = program_scores(ms); s = pd.concat([pr, ms], axis=1); s["group"] = cols.loc[cc, "group"].values; s["sample"] = cols.loc[cc, "sample"].values; s["n_cells"] = cols.loc[cc, "n_cells"].values
    return s, n

rows = []; scores_all = []
for ds, dis in [("GSE155468", "ATAA"), ("GSE213740", "TAAD")]:
    for ct in CTS:
        r = scores_for(ds, ct)
        if r is None: continue
        s, n = r; s["dataset"] = ds; s["celltype"] = ct; scores_all.append(s)
        for k in ["C", "I", "R"] + [m for m in MOD if m in s.columns]:
            x, y = s.loc[s.group == dis, k].dropna(), s.loc[s.group == "Control", k].dropna()
            if len(x) < 2 or len(y) < 2: continue
            g, lo, hi = hedges_g(x.values, y.values); obs, p, na = exact_perm_test(x.values, y.values)
            rows.append(dict(dataset=ds, disease=dis, celltype=ct, measure=k, n_case=len(x), n_ctrl=len(y), n_genes=int(n.get(k, 0)) if k in n else np.nan, mean_diff=obs, g=g, ci_low=lo, ci_high=hi, p_exact=p, n_assignments=na,
                             primary=(ct == "SMC_contractile" and k == "R")))
R = pd.DataFrame(rows); R["p_bh_exploratory"] = np.nan
mask = ~R.primary; R.loc[mask, "p_bh_exploratory"] = multipletests(R.loc[mask, "p_exact"], method="fdr_bh")[1]
R.to_csv(os.path.join(OUT, "sc_patient_level_tests.csv"), index=False); pd.concat(scores_all).to_csv(os.path.join(OUT, "sc_patient_celltype_scores.csv"), index=False)
pd.set_option("display.width", 250)
print("PRIMARY (one per dataset): R in contractile SMC, exact enumeration"); print(R[R.primary][["dataset", "n_case", "n_ctrl", "mean_diff", "g", "ci_low", "ci_high", "p_exact", "n_assignments"]].round(3).to_string(index=False))
for ds in ["GSE155468", "GSE213740"]:
    sub = R[(R.dataset == ds) & (R.celltype == "SMC_contractile") & (R.measure.isin(["C", "I", "R"] + INJURY_SIX + ["Calcium_handling", "IFN_alpha_response", "MHCII_antigen_presentation", "ECM_collagen"]))]
    print(f"\n{ds} contractile SMC (estimation; BH exploratory):"); print(sub[["measure", "n_genes", "mean_diff", "g", "ci_low", "ci_high", "p_exact", "p_bh_exploratory"]].round(3).to_string(index=False))
piv = R[R.measure.isin(["C", "I", "R"])].pivot_table(index="celltype", columns=["dataset", "measure"], values="g").round(2); print("\nHedges' g by cell type (C / I / R):"); print(piv.to_string())
# raw patient scores for the primary comparison
sa = pd.concat(scores_all); prim = sa[sa.celltype == "SMC_contractile"][["dataset", "sample", "group", "n_cells", "C", "I", "R"]].round(2); print("\nraw patient scores (contractile SMC):"); print(prim.to_string(index=False))

# check: control-referenced z with full-pipeline permutation (statistic: mean R case - mean R control), SMC only
def ctrl_ref_R(x, groups, ctrl_mask):
    mu = x.loc[:, ctrl_mask].mean(axis=1); sd = x.loc[:, ctrl_mask].std(axis=1).clip(lower=0.05); z = x.sub(mu, axis=0).div(sd, axis=0).clip(-10, 10)
    ms = pd.DataFrame({k: z.loc[[g for g in v if g in z.index]].mean(axis=0) for k, v in MOD.items() if sum(g in z.index for g in v) >= 5}); pr = program_scores(ms); return pr["R"]
chk = []
for ds, dis in [("GSE155468", "ATAA"), ("GSE213740", "TAAD")]:
    cc = cols.index[(cols.dataset == ds) & (cols.celltype == "SMC_contractile")]; x = pb[cc]; x = x[(x > 0).mean(axis=1) >= 0.2]; grp = cols.loc[cc, "group"].values
    ctrl = grp == "Control"; obs = ctrl_ref_R(x, grp, ctrl); stat_obs = obs[~ctrl].mean() - obs[ctrl].mean(); n1 = ctrl.sum(); vals = []
    for idx in itertools.combinations(range(len(cc)), n1):
        cm = np.zeros(len(cc), bool); cm[list(idx)] = True; r_ = ctrl_ref_R(x, grp, cm); vals.append(r_[~cm].mean() - r_[cm].mean())
    vals = np.array(vals); p = min(1, 2 * min(np.mean(vals >= stat_obs - 1e-12), np.mean(vals <= stat_obs + 1e-12)))
    chk.append(dict(dataset=ds, stat_obs=stat_obs, p_fullpipeline_exact=p, n_assignments=len(vals)))
print("\ncheck (control-referenced z, full-pipeline exact permutation, R in contractile SMC):"); print(pd.DataFrame(chk).round(4).to_string(index=False)); pd.DataFrame(chk).to_csv(os.path.join(OUT, "sc_primary_fullpipeline_check.csv"), index=False)
