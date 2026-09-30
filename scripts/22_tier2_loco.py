"""Tier 2 disease signatures: fold-internal rebuilding + leave-one-cohort-out evaluation.

For each held-out cohort, the disease signatures (REM consensus genes of the remaining cohorts of that disease intersected
with the same pathway themes as in 15_modules.py, top 60 by |z_rem|) are rebuilt WITHOUT the held-out cohort, then scored
in the held-out cohort (label-free z, sign-oriented) and the disease-vs-control Hedges' g is recorded.
This evaluates the signature-construction procedure, not a fixed gene list."""
import os, sys, json, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from v21_lib import *
PROC = os.path.join(ROOT, "results", "processed"); DE = os.path.join(ROOT, "results", "de"); GS = os.path.join(ROOT, "data", "genesets"); SCD = os.path.join(ROOT, "results", "sc")
OUT = os.path.join(ROOT, "results", "v21_tier2_loco"); os.makedirs(OUT, exist_ok=True)
REG = json.load(open(os.path.join(PROC, "cohorts.json"))); DESIGN = json.load(open(os.path.join(ROOT, "results", "meta", "design.json")))
def gmt(path): return {l.split("\t")[0]: set(l.rstrip("\n").split("\t")[2:]) for l in open(path)}
S = {}
for f in ["h.all.v2024.1.Hs.symbols.gmt", "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt", "c2.cp.reactome.v2024.1.Hs.symbols.gmt"]: S.update(gmt(os.path.join(GS, f)))
MK = gmt(os.path.join(SCD, "aorta_celltype_markers.gmt"))
def U(*n): return set().union(*[S.get(x, MK.get(x, set())) for x in n])
THEMES = {"T_SMC_contraction_related": ("TAAD", "down", U("KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION", "REACTOME_SMOOTH_MUSCLE_CONTRACTION", "HALLMARK_MYOGENESIS", "SMC_contractile")),
          "T_hypoxia_glycolysis": ("TAAD", "up", U("HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS", "KEGG_GLYCOLYSIS_GLUCONEOGENESIS")),
          "T_oxidative_metal_stress": ("TAAD", "up", U("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "REACTOME_RESPONSE_TO_METAL_IONS", "REACTOME_METALLOTHIONEINS_BIND_METALS", "REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2", "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES")),
          "T_MYC_ribosome": ("TAAD", "up", U("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "REACTOME_RRNA_PROCESSING")),
          "T_p53_DNA_damage": ("TAAD", "up", U("HALLMARK_P53_PATHWAY", "HALLMARK_DNA_REPAIR", "REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS", "REACTOME_STABILIZATION_OF_P53")),
          "T_NFkB_inflammation": ("TAAD", "up", U("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "Inflammatory_myeloid")),
          "T_endothelial_related": ("TAAD", "down", U("Endothelial")), "T_fibroblast_related": ("TAAD", "down", U("Fibroblast")),
          "A_interferon": ("ATAA", "up", U("HALLMARK_INTERFERON_ALPHA_RESPONSE", "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE")),
          "A_antigen_presentation_adhesion": ("ATAA", "up", U("KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION", "KEGG_ALLOGRAFT_REJECTION", "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION", "KEGG_CELL_ADHESION_MOLECULES_CAMS")),
          "A_ECM_collagen": ("ATAA", "up", U("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", "REACTOME_COLLAGEN_FORMATION", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "REACTOME_COLLAGEN_CHAIN_TRIMERIZATION")),
          "A_lymphocyte_related": ("ATAA", "up", U("T_cell", "NK", "B_cell", "Plasma", "KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION"))}
def load_de(c):
    d = pd.read_csv(os.path.join(DE, f"{c}_de.csv"), index_col=0); d = d[~d.index.duplicated()].dropna(subset=["pvalue", "log2FC"]); d["se"] = d["se"].abs().replace(0, np.nan); return d
def rem_consensus(names):
    dfs = [load_de(c) for c in names]; k = len(dfs); min_c = max(2, int(np.ceil(0.75 * k))); genes = pd.Index(sorted(set().union(*[set(d.index) for d in dfs])))
    Y = pd.DataFrame({c: d["log2FC"].reindex(genes) for c, d in zip(names, dfs)}); V = pd.DataFrame({c: (d["se"] ** 2).reindex(genes) for c, d in zip(names, dfs)}); P = pd.DataFrame({c: d["pvalue"].reindex(genes) for c, d in zip(names, dfs)})
    present = Y.notna() & V.notna() & np.isfinite(V); nk = present.sum(axis=1); keep = nk >= min_c; Y, V, P, present, nk = Y[keep], V[keep], P[keep], present[keep], nk[keep]
    V = V.where(present); w = 1 / V; yf = (w * Y).sum(axis=1) / w.sum(axis=1); Q = (w * (Y.sub(yf, axis=0)) ** 2).sum(axis=1); C = w.sum(axis=1) - (w ** 2).sum(axis=1) / w.sum(axis=1)
    tau2 = ((Q - (nk - 1)) / C).clip(lower=0); wr = (1 / V.add(tau2, axis=0)).where(present); yr = (wr * Y).sum(axis=1) / wr.sum(axis=1); se = np.sqrt(1 / wr.sum(axis=1)); z = yr / se
    fdr = multipletests(2 * stats.norm.sf(np.abs(z)), method="fdr_bh")[1]; sign = np.sign(yr); same = (np.sign(Y).eq(sign, axis=0) & present).sum(axis=1) / nk; nominal = ((P < 0.05) & np.sign(Y).eq(sign, axis=0) & present).sum(axis=1) / nk
    out = pd.DataFrame({"z": z, "lfc": yr}); out["consensus"] = (fdr < 0.05) & (same >= 0.75) & (nominal >= 0.5); out["relaxed"] = (fdr < 0.05) & (same >= 0.75); return out
def build(names):
    sigs = {}
    for th, (dis, direction, gs) in THEMES.items():
        if not any(c in DESIGN[dis]["cohorts"] for c in names): continue
        m = rem_consensus([c for c in names if c in DESIGN[dis]["cohorts"]])
        pool = m[m.consensus & ((m.lfc > 0) if direction == "up" else (m.lfc < 0))]; g = pool.index.intersection(gs)
        if len(g) < 8: pool = m[m.relaxed & ((m.lfc > 0) if direction == "up" else (m.lfc < 0))]; g = pool.index.intersection(gs)
        g = pool.loc[g].z.abs().sort_values(ascending=False).index[:60].tolist()
        if len(g) >= 8: sigs[th] = (g, direction)
    return sigs
cohorts = {c: "ATAA" for c in DESIGN["ATAA"]["cohorts"]}; cohorts.update({c: "TAAD" for c in DESIGN["TAAD"]["cohorts"]})
rows = []
for held, dis in cohorts.items():
    train = [c for c in cohorts if c != held]; sigs = build(train)
    e = pd.read_csv(os.path.join(PROC, f"{held}_expr.csv"), index_col=0); m = pd.read_csv(os.path.join(PROC, f"{held}_meta.csv"), index_col=0); e.columns = e.columns.astype(str); m.index = m.index.astype(str)
    ms, n = module_scores(e, {k: v[0] for k, v in sigs.items()})
    for th, (g, direction) in sigs.items():
        if th not in ms: continue
        s = ms[th] * (1 if direction == "up" else -1); x, y = s[m.loc[s.index, "group"] == dis], s[m.loc[s.index, "group"] == "Control"]
        gg, lo, hi = hedges_g(x.values, y.values); rows.append(dict(held_out=held, held_disease=dis, signature=th, signature_disease=THEMES[th][0], n_genes=int(n[th]), g=gg, ci_low=lo, ci_high=hi, p_mwu=stats.mannwhitneyu(x, y).pvalue, same_disease=(THEMES[th][0] == dis)))
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "tier2_fold_internal_loco.csv"), index=False)
pd.set_option("display.width", 250)
summ = R.groupby(["signature", "held_disease"]).agg(n_folds=("g", "size"), median_g=("g", "median"), min_g=("g", "min"), max_g=("g", "max"), n_p05=("p_mwu", lambda p: int((p < 0.05).sum()))).round(2)
print("Tier 2 signatures rebuilt inside each fold and scored in the held-out cohort (g = disease vs control, oriented to the signature direction):"); print(summ.to_string()); summ.to_csv(os.path.join(OUT, "tier2_fold_internal_summary.csv"))
