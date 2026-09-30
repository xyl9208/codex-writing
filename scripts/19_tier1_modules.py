"""Tier 1 pathway modules: externally defined gene sets (MSigDB v2024.1 + explicitly listed core genes + control-cell
marker sets).  Membership does NOT depend on any differential-expression result of this study.

Outputs results/modules_v21/tier1_modules.gmt (full sets), tier1_modules_noRP.gmt (ribosomal-protein genes removed;
sensitivity), tier1_table.csv (membership + source), tier2_rename.csv.
"""
import os, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); GS = os.path.join(ROOT, "data", "genesets"); SCD = os.path.join(ROOT, "results", "sc")
OUT = os.path.join(ROOT, "results", "modules_v21"); os.makedirs(OUT, exist_ok=True)
def gmt(path):
    return {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(path)}
S = {}
for f in ["h.all.v2024.1.Hs.symbols.gmt", "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt", "c2.cp.reactome.v2024.1.Hs.symbols.gmt"]: S.update(gmt(os.path.join(GS, f)))
CTRL = gmt(os.path.join(SCD, "aorta_celltype_markers_control.gmt"))
ALL = set().union(*[set(v) for v in S.values()])
def U(*n): return sorted(set().union(*[set(S[x]) for x in n]))
def I(a, b): return sorted(set(S[a]) & set(S[b]))
CORE_SMC = ["MYH11", "ACTA2", "CNN1", "TAGLN", "MYOCD", "LMOD1", "SMTN", "MYLK", "DES", "CALD1", "TPM2", "MYL9", "ACTG2", "SYNPO2", "PLN", "KCNMB1", "SORBS1", "CASQ2"]
CORE_CA = ["RYR2", "CASQ1", "CASQ2", "ATP1A2", "PLN", "ATP2A2", "TRPC1", "CACNA1C", "CACNB2", "SLC8A1", "ITPR1", "KCNJ8", "ABCC9"]
CORE_GLY = ["HK1", "HK2", "PFKFB3", "PFKP", "ALDOA", "TPI1", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA", "SLC2A1", "SLC16A3", "PDK1"]
MT = sorted({g for g in ALL if g.startswith(("MT1", "MT2A", "MT3", "MT4"))} | {"MT1X", "MT1E", "MT1F", "MT1G", "MT1M", "MT2A"})
MMP = sorted({g for g in ALL if g.startswith(("MMP", "TIMP", "ADAMTS"))})
MHCII = ["HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DQA2", "HLA-DQB2", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "CD74", "CIITA", "RFX5", "RFXAP", "RFXANK"]
ifna = set(S["HALLMARK_INTERFERON_ALPHA_RESPONSE"]); ifng_spec = sorted(set(S["HALLMARK_INTERFERON_GAMMA_RESPONSE"]) - ifna)
MODULES = {
 "SMC_contractile":        (sorted(set(I("KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION", "REACTOME_SMOOTH_MUSCLE_CONTRACTION")) | set(CORE_SMC)), "KEGG_VSMC_CONTRACTION ∩ REACTOME_SMOOTH_MUSCLE_CONTRACTION ∪ listed core genes", "functional_phenotype"),
 "Calcium_handling":       (sorted(set(CORE_CA) | set(S.get("REACTOME_ION_HOMEOSTASIS", []))), "REACTOME_ION_HOMEOSTASIS ∪ listed core calcium-handling genes", "functional_phenotype"),
 "Cell_matrix_adhesion":   (I("KEGG_FOCAL_ADHESION", "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS"), "KEGG_FOCAL_ADHESION ∩ REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS", "functional_phenotype"),
 "ECM_collagen":           (U("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", "REACTOME_COLLAGEN_FORMATION"), "REACTOME ECM organization ∪ collagen formation", "remodelling"),
 "Matrix_degradation":     (sorted(set(S["REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"]) | set(MMP)), "REACTOME_DEGRADATION_OF_THE_ECM ∪ MMP/TIMP/ADAMTS family (MSigDB symbols)", "remodelling"),
 "Hypoxia":                (S["HALLMARK_HYPOXIA"], "HALLMARK_HYPOXIA", "injury_response"),
 "Glycolysis":             (sorted(set(I("HALLMARK_GLYCOLYSIS", "KEGG_GLYCOLYSIS_GLUCONEOGENESIS")) | set(CORE_GLY)), "HALLMARK_GLYCOLYSIS ∩ KEGG_GLYCOLYSIS ∪ listed core enzymes", "injury_response"),
 "Oxidative_stress_NFE2L2": (sorted(set(S["HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]) | set(S.get("REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2", [])) | set(MT)), "HALLMARK_ROS ∪ REACTOME_NFE2L2 events ∪ metallothioneins", "injury_response"),
 "MYC_ribosome":           (U("HALLMARK_MYC_TARGETS_V1", "REACTOME_RRNA_PROCESSING"), "HALLMARK_MYC_TARGETS_V1 ∪ REACTOME_RRNA_PROCESSING", "injury_response"),
 "p53_DNA_damage":         (U("HALLMARK_P53_PATHWAY", "REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS"), "HALLMARK_P53_PATHWAY ∪ REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS", "injury_response"),
 "NFkB_IL6_inflammation":  (U("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING"), "HALLMARK_TNFA_SIGNALING_VIA_NFKB ∪ HALLMARK_IL6_JAK_STAT3_SIGNALING", "injury_response"),
 "IFN_alpha_response":     (S["HALLMARK_INTERFERON_ALPHA_RESPONSE"], "HALLMARK_INTERFERON_ALPHA_RESPONSE", "injury_response"),
 "IFNG_specific":          (ifng_spec, "HALLMARK_INTERFERON_GAMMA_RESPONSE minus genes shared with HALLMARK_INTERFERON_ALPHA_RESPONSE (custom subset)", "injury_response"),
 "MHCII_antigen_presentation": (MHCII, "explicitly listed MHC class II genes + CD74 + CIITA/RFX transactivators", "immune"),
}
for ct in ["Endothelial", "Fibroblast", "Inflammatory_myeloid", "Macrophage", "T_cell", "NK", "B_cell", "Plasma", "SMC_contractile"]:
    if ct in CTRL: MODULES[f"sig_{ct}"] = (CTRL[ct], "control-cell one-vs-rest markers (GSE155468 + GSE213740 controls)", "cell_signature")
INJURY_SIX = ["Hypoxia", "Glycolysis", "Oxidative_stress_NFE2L2", "MYC_ribosome", "p53_DNA_damage", "NFkB_IL6_inflammation"]
rows = []
with open(os.path.join(OUT, "tier1_modules.gmt"), "w") as f1, open(os.path.join(OUT, "tier1_modules_noRP.gmt"), "w") as f2:
    for k, (g, src, cat) in MODULES.items():
        g = sorted(set(g)); g2 = [x for x in g if not x.startswith(("RPL", "RPS", "MRPL", "MRPS"))]
        f1.write("\t".join([k, src] + g) + "\n"); f2.write("\t".join([k, src + " (ribosomal proteins removed)"] + g2) + "\n")
        rows.append(dict(module=k, category=cat, n_genes=len(g), n_genes_noRP=len(g2), source=src, in_injury_composite=k in INJURY_SIX, genes=";".join(g)))
t = pd.DataFrame(rows); t.to_csv(os.path.join(OUT, "tier1_table.csv"), index=False)
pd.DataFrame([("T_SMC_contractile_loss", "TAAD signature: SMC-contraction-related (down)"), ("T_Endothelial_loss", "TAAD signature: endothelial-related (down)"), ("T_Adventitial_fibroblast_loss", "TAAD signature: fibroblast-related (down)"),
              ("T_Hypoxia_glycolysis", "TAAD signature: hypoxia-glycolysis (up)"), ("T_Oxidative_metal_stress", "TAAD signature: oxidative/metal stress (up)"), ("T_MYC_ribosome_biogenesis", "TAAD signature: MYC-ribosome (up)"),
              ("T_p53_DNA_damage", "TAAD signature: p53-DNA damage (up)"), ("T_Innate_inflammation", "TAAD signature: NF-kB inflammation (up)"), ("A_IFN_I_response", "ATAA signature: interferon response (alpha/beta/gamma sets; up)"),
              ("A_MHCII_antigen_presentation", "ATAA signature: antigen presentation and immune adhesion-related (MHC-I/II + adhesion; up)"), ("A_ECM_collagen", "ATAA signature: ECM-collagen (up)"), ("A_Lymphoid_infiltration", "ATAA signature: lymphocyte-related (up)"),
              ("S_Shared_up", "concordant candidates, both diseases REM-FDR<0.05 (up; n=96 pool)"), ("S_Shared_down", "concordant candidates (down)")], columns=["old", "new"]).to_csv(os.path.join(OUT, "tier2_rename.csv"), index=False)
print(t[["module", "category", "n_genes", "n_genes_noRP", "in_injury_composite"]].to_string(index=False))
print("SMC_contractile:", ", ".join(MODULES["SMC_contractile"][0][:25]))
