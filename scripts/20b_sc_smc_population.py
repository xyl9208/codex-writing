"""EXPLORATORY / SENSITIVITY (single-cell, patient level; not a primary test).  Reviewer point: the primary single-cell tests use
only cells still classified as contractile SMC, so a disease change that consists of cells LEAVING the contractile state could be
missed.  Per dataset (GSE155468 ATAA vs Control; GSE213740 TAAD vs Control), with the patient as the unit:
 (i)  programme scores C / I / R (and the six injury modules, MHCII, IFN-alpha) in the WHOLE SMC population: patient pseudo-bulk
      (mean log-normalised expression over >= 20 cells) of SMC_all = SMC_contractile + SMC_modulated (primary definition) and of
      SMC_all_plus_fibromyocyte (broader mural definition, sensitivity); label-free standardisation across all patients of the dataset;
 (ii) SMC-state composition per patient: fraction SMC_modulated among SMC_all; fraction Fibromyocyte among mural cells
      (captured-cell proportions: subject to dissociation bias, NOT tissue cell counts);
 (iii) within-state scores in SMC_modulated (>= 20 cells); the contractile-state values are reused from sc_patient_level_tests.csv.
Effects: Hedges' g with bootstrap CI, exact enumeration permutation P, BH across all tests computed here (one exploratory family).
Memory: only the rows of SMC-lineage cells are read from the h5ad (h5py block reads of the on-disk CSR), no full load."""
import os, sys, resource, warnings, h5py, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from v21_lib import *
warnings.filterwarnings("ignore")
H5 = os.path.join(ROOT, "results", "sc", "aorta_sc_integrated.h5ad"); OUT = os.path.join(ROOT, "results", "v21_sc"); os.makedirs(OUT, exist_ok=True)
MOD = load_tier1(); MIN_CELLS = 20; DS = [("GSE155468", "ATAA"), ("GSE213740", "TAAD")]
DEFS = {"SMC_all": ["SMC_contractile", "SMC_modulated"], "SMC_all_plus_fibromyocyte": ["SMC_contractile", "SMC_modulated", "Fibromyocyte"]}
MEAS = ["C", "I", "R"] + INJURY_SIX + ["MHCII_antigen_presentation", "IFN_alpha_response"]

def read_obs_cat(f, k): g = f["obs"][k]; return pd.Categorical.from_codes(g["codes"][:], g["categories"][:].astype(str))
def read_rows_csr(f, key, rows, ncol, block=4000):
    """Selected rows of an on-disk CSR matrix, read block by block (only blocks that contain wanted rows are touched)."""
    ip = f[f"{key}/indptr"][:]; rows = np.sort(rows); parts = []
    for s in range(0, len(ip) - 1, block):
        e = min(s + block, len(ip) - 1); r = rows[(rows >= s) & (rows < e)]
        if len(r) == 0: continue
        a0, b0 = ip[s], ip[e]; m = sp.csr_matrix((f[f"{key}/data"][a0:b0], f[f"{key}/indices"][a0:b0], ip[s:e + 1] - a0), shape=(e - s, ncol)); parts.append(m[r - s])
    return sp.vstack(parts).tocsr()

with h5py.File(H5, "r") as f:
    obs = pd.DataFrame({k: read_obs_cat(f, k) for k in ["sample", "dataset", "group", "celltype"]}, index=f["obs/_index"][:].astype(str)); genes = f["var/_index"][:].astype(str)
    keep = np.where(obs.celltype.isin(DEFS["SMC_all_plus_fibromyocyte"]).values)[0]; X = read_rows_csr(f, "X", keep, len(genes))
a = ad.AnnData(X=X, obs=obs.iloc[keep].copy(), var=pd.DataFrame(index=genes)); a.obs["state"] = a.obs["celltype"].astype(str)
smeta = obs.drop_duplicates("sample").set_index("sample")[["dataset", "group"]].astype(str); smeta.index = smeta.index.astype(str)
print(f"SMC-lineage cells read: {a.n_obs} of {len(obs)} (genes {a.n_vars}); X range {X.data.min():.3f}-{X.data.max():.3f} (log-normalised)")

# (ii) composition per patient (all patients; every patient has >= 20 SMC_all cells)
cnt = pd.crosstab(obs["sample"].astype(str), obs["celltype"].astype(str)).reindex(columns=["SMC_contractile", "SMC_modulated", "Fibromyocyte"], fill_value=0)
comp = smeta.join(cnt).reset_index().rename(columns={"index": "sample", "SMC_contractile": "n_SMC_contractile", "SMC_modulated": "n_SMC_modulated", "Fibromyocyte": "n_Fibromyocyte"})
comp["frac_modulated_of_SMC"] = comp.n_SMC_modulated / (comp.n_SMC_contractile + comp.n_SMC_modulated)
comp["frac_fibromyocyte_of_mural"] = comp.n_Fibromyocyte / (comp.n_SMC_contractile + comp.n_SMC_modulated + comp.n_Fibromyocyte)
comp = comp[["dataset", "sample", "group", "n_SMC_contractile", "n_SMC_modulated", "n_Fibromyocyte", "frac_modulated_of_SMC", "frac_fibromyocyte_of_mural"]].sort_values(["dataset", "group", "sample"])
comp.to_csv(os.path.join(OUT, "smc_state_composition.csv"), index=False)

# (i) / (iii) pseudo-bulk per scope -> label-free module / programme scores across the patients of each dataset
def scope_scores(label, states):
    b = a[a.obs.state.isin(states).values].copy(); b.obs["celltype"] = label; pb, n = celltype_pseudobulk(b, min_cells=MIN_CELLS)
    cols = pd.DataFrame([c.split("|") for c in pb.columns], columns=["sample", "celltype"], index=pb.columns).join(smeta, on="sample"); cols["n_cells"] = n; out = {}
    for ds, _ in DS:
        cc = cols.index[cols.dataset == ds]
        if len(cc) < 4: continue
        x = pb[cc]; x = x[(x > 0).mean(axis=1) >= 0.2]  # label-free detectability filter within the scope (as in 20_sc_patient_tests.py)
        ms, ng = module_scores(x, MOD); s = pd.concat([program_scores(ms), ms], axis=1); s["group"] = cols.loc[cc, "group"].values; s["sample"] = cols.loc[cc, "sample"].values; s["n_cells"] = cols.loc[cc, "n_cells"].values; out[ds] = (s, ng)
    return pb, out

rows, pat = [], []
def add_tests(s, ds, dis, definition, scope, measures, atype, ng=None):
    for k in measures:
        x, y = s.loc[s.group == dis, k].dropna(), s.loc[s.group == "Control", k].dropna()
        if len(x) < 2 or len(y) < 2: continue
        g, lo, hi = hedges_g(x.values, y.values); d, p, na = exact_perm_test(x.values, y.values)
        rows.append(dict(dataset=ds, definition=definition, celltype_scope=scope, measure=k, n_case=len(x), n_ctrl=len(y), n_genes=(int(ng[k]) if ng is not None and k in ng else np.nan), mean_diff=d, g=g, ci_low=lo, ci_high=hi, p_exact=p, n_assignments=na, analysis_type=atype))

scopes = [("SMC_all", "whole_SMC", DEFS["SMC_all"], "exploratory_whole_SMC_pseudobulk"), ("SMC_all_plus_fibromyocyte", "whole_mural", DEFS["SMC_all_plus_fibromyocyte"], "exploratory_whole_mural_pseudobulk"),
          ("SMC_modulated", "within_state", ["SMC_modulated"], "exploratory_within_state"), ("SMC_contractile", "within_state", ["SMC_contractile"], "check_only")]
pb_check = None
for definition, scope, states, atype in scopes:
    pb, out = scope_scores(definition, states)
    if definition == "SMC_contractile": pb_check = (pb, out); continue  # recomputed only to verify the row-extraction path against the stored primary pipeline
    for ds, dis in DS:
        if ds not in out: continue
        s, ng = out[ds]; s = s.assign(dataset=ds, definition=definition, celltype_scope=scope); pat.append(s)
        add_tests(s, ds, dis, definition, scope, MEAS if scope != "within_state" else ["C", "I", "R"], atype, ng)
for ds, dis in DS:
    c = comp[comp.dataset == ds]
    add_tests(c, ds, dis, "SMC_all", "composition", ["frac_modulated_of_SMC"], "exploratory_composition"); add_tests(c, ds, dis, "SMC_all_plus_fibromyocyte", "composition", ["frac_fibromyocyte_of_mural"], "exploratory_composition")
R = pd.DataFrame(rows); R["p_bh_exploratory"] = multipletests(R.p_exact, method="fdr_bh")[1]  # one exploratory family: every test computed in this script

# reuse the contractile-state values from the primary pipeline (not recomputed; R there is the primary test)
prim = pd.read_csv(os.path.join(OUT, "sc_patient_level_tests.csv")); prim = prim[(prim.celltype == "SMC_contractile") & prim.measure.isin(["C", "I", "R"])].copy()
prim["definition"] = "SMC_contractile"; prim["celltype_scope"] = "within_state"; prim["analysis_type"] = np.where(prim.primary, "reused_primary_sc_patient_level_tests", "reused_exploratory_sc_patient_level_tests")
R = pd.concat([R, prim[R.columns]], ignore_index=True)
R = R[["dataset", "definition", "celltype_scope", "measure", "n_case", "n_ctrl", "mean_diff", "g", "ci_low", "ci_high", "p_exact", "n_assignments", "p_bh_exploratory", "analysis_type", "n_genes"]]; R.to_csv(os.path.join(OUT, "smc_population_tests.csv"), index=False)
P = pd.concat(pat, ignore_index=True); lead = ["dataset", "definition", "celltype_scope", "sample", "group", "n_cells", "C", "I", "R", "n_injury_modules"]; P = P[lead + [c for c in P.columns if c not in lead]]
P.to_csv(os.path.join(OUT, "smc_population_patient_scores.csv"), index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
# check: the contractile-only pseudo-bulk / scores recomputed through the h5py row-extraction path must match the stored primary pipeline
pbs = pd.read_csv(os.path.join(OUT, "pseudobulk_mean_lognorm.csv"), index_col=0); pbc, outc = pb_check; common = [c for c in pbc.columns if c in pbs.columns]
sc_old = pd.read_csv(os.path.join(OUT, "sc_patient_celltype_scores.csv")); sc_old = sc_old[sc_old.celltype == "SMC_contractile"].set_index(["dataset", "sample"])
d_scores = max(np.abs(outc[ds][0].set_index("sample")[["C", "I", "R"]].values - sc_old.loc[ds].loc[outc[ds][0]["sample"], ["C", "I", "R"]].values).max() for ds, _ in DS)
print(f"check vs stored primary pipeline (SMC_contractile): {len(common)} patient columns, max |diff| pseudo-bulk = {np.abs(pbc[common].values - pbs[common].loc[pbc.index].values).max():.2e}, max |diff| C/I/R = {d_scores:.2e}")
sm = pd.read_csv(os.path.join(OUT, "sc_patient_level_tests.csv")); sm = sm[(sm.celltype == "SMC_modulated") & sm.measure.isin(["C", "I", "R"])].set_index(["dataset", "measure"])
rm = R[(R.definition == "SMC_modulated")].set_index(["dataset", "measure"]); print(f"check within-state SMC_modulated C/I/R vs stored table: max |diff g| = {np.abs(rm.g - sm.loc[rm.index, 'g']).max():.2e}")

K = ["dataset", "definition", "measure", "n_case", "n_ctrl", "mean_diff", "g", "ci_low", "ci_high", "p_exact", "n_assignments", "p_bh_exploratory"]
print("\n(i) EXPLORATORY: C / I / R in the whole SMC population (patient pseudo-bulk; exact enumeration; BH over all tests in this script)")
print(R[R.celltype_scope.isin(["whole_SMC", "whole_mural"]) & R.measure.isin(["C", "I", "R"])][K].round(3).to_string(index=False))
print("\n(i) EXPLORATORY: injury modules, MHCII, IFN-alpha in SMC_all")
print(R[(R.definition == "SMC_all") & R.measure.isin(INJURY_SIX + ["MHCII_antigen_presentation", "IFN_alpha_response"])][["dataset", "measure", "n_genes", "mean_diff", "g", "ci_low", "ci_high", "p_exact", "p_bh_exploratory"]].round(3).to_string(index=False))
print("\n(ii) SMC-state composition per patient (captured-cell proportions: subject to dissociation bias, NOT tissue cell counts)")
print(comp.round(3).to_string(index=False)); print(R[R.celltype_scope == "composition"][K].round(3).to_string(index=False))
print("\n(iii) comparison across scopes: Hedges' g [p_exact] for C / I / R (SMC_contractile rows reused from sc_patient_level_tests.csv; R there is the primary test)")
cmpt = R[R.measure.isin(["C", "I", "R"]) & (R.celltype_scope != "composition")].copy(); cmpt["cell"] = cmpt.apply(lambda r: f"{r.g:.2f} [{r.p_exact:.3f}] n={r.n_case}v{r.n_ctrl}", axis=1)
print(cmpt.pivot_table(index=["dataset", "definition"], columns="measure", values="cell", aggfunc="first").reindex(columns=["C", "I", "R"]).to_string())
print("\nraw patient scores (SMC_all):"); print(P[P.definition == "SMC_all"][["dataset", "sample", "group", "n_cells", "C", "I", "R"]].round(2).to_string(index=False))
print(f"\npeak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.2f} GB")
