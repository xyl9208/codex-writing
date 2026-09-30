"""Pre-specified gene modules and module-level (patient-level) evidence in bulk cohorts.

Modules are defined ONCE from the v2 meta-analysis (disease-associated consensus genes intersected with curated
pathway / cell-type marker sets) and are then tested unchanged in every subsequent data structure
(cell-type-resolved scRNA-seq, disease stage, within-patient regions, medial layers).

Bulk evidence per module: per-sample control-referenced module score (mean z over module genes, sign-oriented so that
the disease-associated direction is positive) -> Hedges' g (disease vs control) per cohort with 95 % CI ->
DerSimonian-Laird pooled g per disease -> cross-disease difference; GSE318877 (subacute/chronic dissection vs dilatation)
scored at the PATIENT level (mean over punches).
"""
import os, json, warnings, numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
M2 = os.path.join(ROOT, "results", "meta_v2"); GS = os.path.join(ROOT, "data", "genesets"); SCD = os.path.join(ROOT, "results", "sc"); PROC = os.path.join(ROOT, "results", "processed")
OUT = os.path.join(ROOT, "results", "modules"); os.makedirs(OUT, exist_ok=True)
REG = json.load(open(os.path.join(PROC, "cohorts.json"))); DESIGN = json.load(open(os.path.join(ROOT, "results", "meta", "design.json")))
A = pd.read_csv(os.path.join(M2, "ATAA_meta_v2.csv"), index_col=0); T = pd.read_csv(os.path.join(M2, "TAAD_meta_v2.csv"), index_col=0)

def gmt(path):
    d = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t"); d[f[0]] = set(f[2:])
    return d
sets = {}
for f in ["h.all.v2024.1.Hs.symbols.gmt", "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt", "c2.cp.reactome.v2024.1.Hs.symbols.gmt"]: sets.update(gmt(os.path.join(GS, f)))
markers = gmt(os.path.join(SCD, "aorta_celltype_markers.gmt"))
def U(*names): return set().union(*[sets.get(n, markers.get(n, set())) for n in names])

def pool(m, direction, relaxed=False):
    """disease-associated genes: v2 consensus (primary) or relaxed (fdr_rem<0.05 & direction consistent) if the consensus set is too small."""
    base = m[m.consensus_v2] if not relaxed else m[(m.fdr_rem < 0.05) & (m.frac_same_direction >= 0.75)]
    return base[base.direction == direction]

DEFS = {
 # TAAD-associated programmes
 "T_SMC_contractile_loss":      ("TAAD", "down", U("KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION", "REACTOME_SMOOTH_MUSCLE_CONTRACTION", "HALLMARK_MYOGENESIS", "SMC_contractile")),
 "T_Hypoxia_glycolysis":        ("TAAD", "up",   U("HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS", "KEGG_GLYCOLYSIS_GLUCONEOGENESIS")),
 "T_Oxidative_metal_stress":    ("TAAD", "up",   U("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "REACTOME_RESPONSE_TO_METAL_IONS", "REACTOME_METALLOTHIONEINS_BIND_METALS", "REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2", "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES")),
 "T_MYC_ribosome_biogenesis":   ("TAAD", "up",   U("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "REACTOME_RRNA_PROCESSING", "REACTOME_MAJOR_PATHWAY_OF_RRNA_PROCESSING_IN_THE_NUCLEOLUS_AND_CYTOSOL")),
 "T_p53_DNA_damage":            ("TAAD", "up",   U("HALLMARK_P53_PATHWAY", "HALLMARK_DNA_REPAIR", "REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS", "REACTOME_STABILIZATION_OF_P53")),
 "T_Innate_inflammation":       ("TAAD", "up",   U("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "Inflammatory_myeloid")),
 "T_Endothelial_loss":          ("TAAD", "down", U("Endothelial")),
 "T_Adventitial_fibroblast_loss": ("TAAD", "down", U("Fibroblast")),
 # ATAA-associated programmes
 "A_IFN_I_response":            ("ATAA", "up",   U("HALLMARK_INTERFERON_ALPHA_RESPONSE", "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE")),
 "A_MHCII_antigen_presentation": ("ATAA", "up",  U("KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION", "KEGG_ALLOGRAFT_REJECTION", "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION", "KEGG_CELL_ADHESION_MOLECULES_CAMS")),
 "A_ECM_collagen":              ("ATAA", "up",   U("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", "REACTOME_COLLAGEN_FORMATION", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "REACTOME_COLLAGEN_CHAIN_TRIMERIZATION")),
 "A_Lymphoid_infiltration":     ("ATAA", "up",   U("T_cell", "NK", "B_cell", "Plasma", "KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION")),
}
rows = []; modules = {}
for name, (dis, direction, gs) in DEFS.items():
    m = A if dis == "ATAA" else T; relaxed = False
    p = pool(m, direction); g = p.index.intersection(gs)
    if len(g) < 8: p = pool(m, direction, relaxed=True); g = p.index.intersection(gs); relaxed = True
    g = p.loc[g].reindex(p.loc[g].z_rem.abs().sort_values(ascending=False).index).index[:60].tolist()
    modules[name] = g; rows.append(dict(module=name, disease=dis, direction=direction, n_genes=len(g), relaxed_pool=relaxed, genes=";".join(g)))
# shared programmes: fdr_rem<0.05 in both diseases with the same direction
common = A.index.intersection(T.index); both = (A.loc[common, "fdr_rem"] < 0.05) & (T.loc[common, "fdr_rem"] < 0.05) & (np.sign(A.loc[common, "lfc_rem"]) == np.sign(T.loc[common, "lfc_rem"]))
for name, direction in [("S_Shared_up", "up"), ("S_Shared_down", "down")]:
    g = common[both & ((A.loc[common, "lfc_rem"] > 0) if direction == "up" else (A.loc[common, "lfc_rem"] < 0))]
    g = (A.loc[g, "z_rem"].abs() + T.loc[g, "z_rem"].abs()).sort_values(ascending=False).index[:60].tolist()
    modules[name] = g; rows.append(dict(module=name, disease="both", direction=direction, n_genes=len(g), relaxed_pool=False, genes=";".join(g)))
tab = pd.DataFrame(rows); tab.to_csv(os.path.join(OUT, "modules_table.csv"), index=False)
with open(os.path.join(OUT, "modules.gmt"), "w") as fh:
    for k, v in modules.items(): fh.write("\t".join([k, "meta_v2"] + v) + "\n")
print(tab[["module", "disease", "direction", "n_genes", "relaxed_pool"]].to_string(index=False))
for k, v in modules.items(): print(k, ":", ", ".join(v[:12]))
# overlap between modules and the composition marker sets used later (documented for the conditioned sensitivity analysis)
ov = {k: {ct: len(set(v) & markers[ct]) for ct in ["SMC_contractile", "Inflammatory_myeloid", "T_cell", "Endothelial", "Fibroblast"]} for k, v in modules.items()}
pd.DataFrame(ov).T.to_csv(os.path.join(OUT, "module_marker_overlap.csv")); print("\nmodule x marker-set overlap:\n", pd.DataFrame(ov).T.to_string())

# ------------------------------------------------------------------ bulk module evidence (patient level)
def hedges(x, y):
    nx_, ny = len(x), len(y); sp = np.sqrt(((nx_ - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / (nx_ + ny - 2))
    g = (np.mean(x) - np.mean(y)) / sp * (1 - 3 / (4 * (nx_ + ny) - 9)); se = np.sqrt((nx_ + ny) / (nx_ * ny) + g ** 2 / (2 * (nx_ + ny)))
    return g, se
def rem(g, se):
    w = 1 / se ** 2; yf = np.sum(w * g) / np.sum(w); Q = np.sum(w * (g - yf) ** 2); df = len(g) - 1; C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0, (Q - df) / C) if C > 0 else 0; wr = 1 / (se ** 2 + tau2); y = np.sum(wr * g) / np.sum(wr); s = np.sqrt(1 / np.sum(wr))
    return y, s, 2 * stats.norm.sf(abs(y / s)), max(0, (Q - df) / Q * 100) if Q > 0 else 0
def module_scores(expr, meta, ref_group="Control"):
    ref = meta.index[meta.group == ref_group]
    mu = expr[ref].mean(axis=1); sd = expr[ref].std(axis=1).clip(lower=0.1); z = expr.sub(mu, axis=0).div(sd, axis=0).clip(-10, 10)
    out = {}
    for k, v in modules.items():
        g = [x for x in v if x in z.index]; sgn = 1 if tab.set_index("module").loc[k, "direction"] == "up" else -1
        out[k] = sgn * z.loc[g].mean(axis=0) if len(g) >= 5 else np.nan
    return pd.DataFrame(out)
cohorts = {c: "ATAA" for c in DESIGN["ATAA"]["cohorts"]}; cohorts.update({c: "TAAD" for c in DESIGN["TAAD"]["cohorts"]})
eff = []; allscores = []
for c, dis in cohorts.items():
    e = pd.read_csv(os.path.join(PROC, f"{c}_expr.csv"), index_col=0); m = pd.read_csv(os.path.join(PROC, f"{c}_meta.csv"), index_col=0); e.columns = e.columns.astype(str); m.index = m.index.astype(str); e = e[~e.index.duplicated()]
    s = module_scores(e, m.loc[e.columns]); s["group"] = m.loc[s.index, "group"]; s["cohort"] = c; allscores.append(s)
    for k in modules:
        x, y = s.loc[s.group == dis, k].dropna(), s.loc[s.group == "Control", k].dropna()
        if len(x) < 2 or len(y) < 2: continue
        g, se = hedges(x.values, y.values); eff.append(dict(cohort=c, disease=dis, module=k, g=g, se=se, ci_low=g - 1.96 * se, ci_high=g + 1.96 * se, p_mwu=stats.mannwhitneyu(x, y).pvalue, n_case=len(x), n_ctrl=len(y)))
pd.concat(allscores).to_csv(os.path.join(OUT, "bulk_module_scores.csv")); eff = pd.DataFrame(eff); eff.to_csv(os.path.join(OUT, "bulk_module_effects_per_cohort.csv"), index=False)
mrows = []
for k in modules:
    r = {"module": k}
    for dis in ["ATAA", "TAAD"]:
        s = eff[(eff.module == k) & (eff.disease == dis)]; y, se, p, i2 = rem(s.g.values, s.se.values)
        r.update({f"g_{dis}": y, f"se_{dis}": se, f"p_{dis}": p, f"I2_{dis}": i2, f"n_cohorts_{dis}": len(s), f"frac_same_dir_{dis}": float((np.sign(s.g) == np.sign(y)).mean()), f"n_cohorts_p05_{dis}": int((s.p_mwu < 0.05).sum())})
    r["diff_TAAD_minus_ATAA"] = r["g_TAAD"] - r["g_ATAA"]; r["z_diff"] = r["diff_TAAD_minus_ATAA"] / np.sqrt(r["se_TAAD"] ** 2 + r["se_ATAA"] ** 2); r["p_diff"] = 2 * stats.norm.sf(abs(r["z_diff"]))
    mrows.append(r)
mt = pd.DataFrame(mrows).set_index("module"); mt.to_csv(os.path.join(OUT, "bulk_module_meta.csv"))
pd.set_option("display.width", 250); print("\nBulk module-level evidence (pooled Hedges' g, disease vs control):"); print(mt[["g_ATAA", "p_ATAA", "n_cohorts_p05_ATAA", "g_TAAD", "p_TAAD", "n_cohorts_p05_TAAD", "diff_TAAD_minus_ATAA", "p_diff"]].round(3).to_string())

# GSE318877 (subacute/chronic dissection vs dilatation), patient-level scores from punch-level log2TPM (within-cohort z)
er = pd.read_csv(os.path.join(PROC, "GSE318877_regions_log2tpm.csv"), index_col=0); mr = pd.read_csv(os.path.join(PROC, "GSE318877_regions_meta.csv"), index_col=0)
er = er[~er.index.duplicated()]; er = er[(er > 1).mean(axis=1) > 0.25]; z = er.sub(er.mean(axis=1), axis=0).div(er.std(axis=1).clip(lower=0.1), axis=0)
sc_ = pd.DataFrame({k: (1 if tab.set_index("module").loc[k, "direction"] == "up" else -1) * z.loc[[x for x in v if x in z.index]].mean(axis=0) for k, v in modules.items()})
sc_ = sc_.join(mr[["subject", "group", "region"]]); sc_.to_csv(os.path.join(OUT, "GSE318877_punch_module_scores.csv"))
pt = sc_.groupby("subject").agg({**{k: "mean" for k in modules}, "group": "first"}); pt.to_csv(os.path.join(OUT, "GSE318877_patient_module_scores.csv"))
rows = []
for k in modules:
    x, y = pt.loc[pt.group == "TAAD", k], pt.loc[pt.group == "ATAA", k]; g, se = hedges(x.values, y.values)
    rows.append(dict(module=k, mean_dissection=x.mean(), mean_dilatation=y.mean(), g=g, ci_low=g - 1.96 * se, ci_high=g + 1.96 * se, p_exact_mwu=stats.mannwhitneyu(x, y, method="exact").pvalue, n=f"{len(x)} vs {len(y)}"))
d318 = pd.DataFrame(rows).set_index("module"); d318.to_csv(os.path.join(OUT, "GSE318877_patient_module_tests.csv")); print("\nGSE318877 patient-level (3 subacute/chronic dissection vs 4 dilatation):"); print(d318.round(3).to_string())
