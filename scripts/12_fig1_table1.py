"""Figure 1 (study design schematic) and Table 1 (cohort summary)."""
import os, json, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REG = json.load(open(os.path.join(ROOT, "results", "processed", "cohorts.json"))); DESIGN = json.load(open(os.path.join(ROOT, "results", "meta", "design.json")))
de = pd.read_csv(os.path.join(ROOT, "results", "de", "de_summary.csv")).set_index("cohort")
FIG = os.path.join(ROOT, "figures"); TAB = os.path.join(ROOT, "manuscript"); os.makedirs(TAB, exist_ok=True)

rows = []
extra = {"GSE26155": ("Folkersen et al. 2011 (ASAP)", "Ascending aorta intima-media; dilated (>45 mm) vs non-dilated (<40 mm) TAV patients", "Affymetrix Human Exon 1.0 ST"),
         "GSE140947": ("Chen et al. 2020", "Ascending aortic media (aneurysm neck + belly per patient) vs organ-donor media", "RNA-seq (NextSeq 500)"),
         "GSE235161": ("Ma et al. 2024", "Thoracic aortic aneurysm tissue (surgical) vs non-aneurysmal aorta (CABG/donor)", "RNA-seq (FPKM)"),
         "GSE202267": ("2022", "Thoracic aortic aneurysm vs non-aneurysmal aorta (surgical controls)", "RNA-seq (counts)"),
         "GSE52093": ("Pan et al. 2014", "Acute Stanford type A dissection, ascending aorta vs normal ascending aorta", "Illumina HumanHT-12 V4"),
         "GSE153434": ("Zhou et al. 2020", "Acute type A dissection ascending aorta vs normal (heart-transplant donors)", "RNA-seq (HiSeq X Ten)"),
         "GSE98770": ("Kimura et al. 2017", "Acute type A dissection intima-media vs transplant-donor ascending aorta", "Agilent SurePrint G3 8x60K"),
         "GSE294606": ("2025", "Aortic dissection tissue vs control aorta", "RNA-seq (counts)"),
         "GSE267434": ("2024", "Dissected ascending aorta vs normal ascending aorta", "RNA-seq (FPKM)"),
         "GSE147026": ("Li et al. 2020", "Acute dissection aortic media vs control media", "RNA-seq"),
         "GSE190635": ("2021", "Aortic dissection vs healthy control aorta (excluded: label inconsistency)", "Affymetrix HG-U133 Plus 2.0"),
         "GSE318877": ("Kimura et al. 2026", "Micro-regional (100 um) punches of ascending aortic media; dissection (3 pts) vs dilatation (4 pts); 70 punches", "Micro-regional RNA-seq")}
for c in list(DESIGN["ATAA"]["cohorts"]) + list(DESIGN["TAAD"]["cohorts"]) + DESIGN["excluded"] + DESIGN["direct"]:
    g = REG[c]["groups"]; case = [k for k in g if k != "Control"]
    role = "ATAA meta-analysis" if c in DESIGN["ATAA"]["cohorts"] else "TAAD meta-analysis" if c in DESIGN["TAAD"]["cohorts"] else "excluded (QC)" if c in DESIGN["excluded"] else "direct TAAD-vs-ATAA contrast"
    rows.append(dict(Accession=c, Reference=extra[c][0], Disease=REG[c]["disease"], Design=extra[c][1], Platform=extra[c][2], Cases="; ".join(f"{k}={v}" for k, v in g.items() if k != "Control") if c != "GSE318877" else "TAAD=3 subjects (34 punches)", Controls=g.get("Control", "ATAA=4 subjects (36 punches)"),
                     DE_method=de.loc[c, "method"], DEGs_FDR05=int(de.loc[c, "n_FDR05"]), Role=role))
rows.append(dict(Accession="GSE155468", Reference="Li et al. 2020", Disease="ATAA", Design="scRNA-seq of ascending aorta: 8 ATAA vs 3 controls", Platform="10x Genomics", Cases="ATAA=8", Controls=3, DE_method="pseudo-bulk limma", DEGs_FDR05="-", Role="single-cell validation"))
rows.append(dict(Accession="GSE213740", Reference="2022", Disease="TAAD", Design="scRNA-seq of ascending aortic wall: 6 sporadic type A dissection vs 3 donors", Platform="10x Genomics", Cases="TAAD=6", Controls=3, DE_method="pseudo-bulk limma", DEGs_FDR05="-", Role="single-cell validation"))
t1 = pd.DataFrame(rows); t1.to_csv(os.path.join(TAB, "Table1_cohorts.csv"), index=False)
print(t1[["Accession", "Disease", "Cases", "Controls", "DE_method", "DEGs_FDR05", "Role"]].to_string(index=False))

# ---------------------------------------------------------------- Figure 1: schematic
fig, ax = plt.subplots(figsize=(13, 8.2)); ax.set_xlim(0, 13); ax.set_ylim(0, 8.2); ax.axis("off")
def box(x, y, w, h, text, fc, fs=8.5, weight="normal"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12", fc=fc, ec="k", lw=0.8))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, weight=weight)
def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12, lw=0.9, color="k"))
ax.text(6.5, 7.9, "GEO systematic search (Sep 2026): human ascending/thoracic aortic aneurysm and dissection tissue transcriptomes", ha="center", fontsize=10, weight="bold")
box(0.3, 6.2, 4.0, 1.3, "ATAA vs non-diseased aorta\n4 bulk cohorts, n = 115\nGSE26155 (53) . GSE140947 (12)\nGSE235161 (42) . GSE202267 (8)", "#D8F3EF", weight="bold")
box(4.6, 6.2, 4.0, 1.3, "TAAD vs non-diseased aorta\n6 bulk cohorts, n = 76\nGSE52093 (12) . GSE153434 (20) . GSE98770 (11)\nGSE294606 (13) . GSE267434 (12) . GSE147026 (8)", "#F8DCDC", weight="bold")
box(8.9, 6.2, 3.9, 1.3, "Independent validation resources\nDirect contrast: GSE318877 (TAAD vs ATAA media, 70 punches)\nscRNA-seq: GSE155468 (ATAA) + GSE213740 (TAAD), 73,755 cells\nExcluded after QC: GSE190635 (label inconsistency)", "#EEEEEE", fs=7.5)
box(0.3, 4.6, 8.3, 1.2, "Harmonised per-cohort processing (symbol harmonisation, re-processing from raw where needed) -> per-cohort DE\n(limma-trend for arrays/FPKM; DESeq2 for counts) -> QC: PCA, cross-cohort concordance of DE profiles", "#FFFFFF", fs=8)
arrow(2.3, 6.2, 2.3, 5.8); arrow(6.6, 6.2, 6.6, 5.8)
box(0.3, 3.1, 4.0, 1.1, "Random-effects meta-analysis (Stouffer + DL)\nATAA consensus DEGs; LOCO robustness", "#D8F3EF", fs=8)
box(4.6, 3.1, 4.0, 1.1, "Random-effects meta-analysis (Stouffer + DL)\nTAAD consensus DEGs; LOCO robustness", "#F8DCDC", fs=8)
arrow(2.3, 4.6, 2.3, 4.2); arrow(6.6, 4.6, 6.6, 4.2)
box(0.3, 1.5, 8.3, 1.2, "ATAA vs TAAD comparison\nGene classes (shared / TAAD-specific / ATAA-specific / discordant; z_diff) . GSEA on meta-rankings (Hallmark, KEGG, Reactome) and ORA\nssGSEA cell-composition inference (aorta scRNA-derived markers) . STRING PPI hubs . control-referenced signature scores + LOCO classifier", "#FFF4D6", fs=7.8)
arrow(2.3, 3.1, 3.5, 2.7); arrow(6.6, 3.1, 5.4, 2.7)
box(8.9, 1.5, 3.9, 2.7, "Validation\n\n- Direct TAAD-vs-ATAA contrast (GSE318877):\n  gene-level, pathway-level, signature scores\n- scRNA-seq: cell-type-resolved pseudo-bulk DE,\n  cellular source of gene classes, composition\n- LOCO classifier AUC across 10 cohorts", "#EEEEEE", fs=7.8)
arrow(8.6, 2.1, 8.9, 2.1); arrow(10.85, 6.2, 10.85, 4.2)
box(0.3, 0.2, 12.5, 0.9, "Output: distinct and shared molecular programmes of ascending aortic aneurysm and acute type A dissection", "#E9E3F5", fs=9, weight="bold")
arrow(4.45, 1.5, 4.45, 1.1)
fig.savefig(os.path.join(FIG, "Fig1_study_design.png"), dpi=200, bbox_inches="tight"); fig.savefig(os.path.join(FIG, "Fig1_study_design.pdf"), bbox_inches="tight"); plt.close(fig)
