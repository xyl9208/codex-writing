"""Main figures v2 (six) and supplementary figures for the v2.1 manuscript.

Fig 1  evidence structure and cohort comparability
Fig 2  Tier 1 programmes in bulk cohorts (per-cohort and pooled effects; module heatmap; composition-conditioned sensitivity)
Fig 3  patient-level programme scores within cell types (scRNA-seq)
Fig 4  stage (GSE222318 within-cohort stages, estimation only; GSE189795 acute external replicate)
Fig 5  within-patient region (GSE140947), medial layer (GSE318877) and ordered dilatation (GSE26155)
Fig 6  programme organisation: (a) Hedges' g of C vs I per dataset (no R inference: different SDs), (b) within-dataset mean differences on the score scale with delta-R iso-lines, (c) state comparison table (no arrows)
FigS   Tier 2 fold-internal evaluation; composition-conditioned sensitivity (all modules)
"""
import os, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
from matplotlib.patches import FancyBboxPatch
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures"); os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False, "axes.spines.right": False})
COL = {"ATAA": "#1f77b4", "TAAD": "#d62728", "Control": "#7f7f7f", "Acute": "#d62728", "Subacute": "#ff7f0e", "Chronic": "#9467bd"}
DESIGN = json.load(open(os.path.join(R, "meta", "design.json"))); ATAA_C = DESIGN["ATAA"]["cohorts"]; TAAD_C = DESIGN["TAAD"]["cohorts"]
LAB = {"C": "C (SMC contractile)", "I": "I (injury composite)", "R": "R = I + C (relative state)"}
def save(fig, name): fig.savefig(os.path.join(FIG, name + ".png"), dpi=220, bbox_inches="tight"); fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight"); plt.close(fig); print("saved", name)
def panel(ax, letter, dy=1.06): ax.text(-0.1, dy, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")

# ====================================================================== Figure 1
def fig1():
    rho = pd.read_csv(os.path.join(R, "meta", "cohort_concordance_rho.csv"), index_col=0); S = json.load(open(os.path.join(R, "meta_v2", "summary_v2.json")))
    fig = plt.figure(figsize=(11, 7.8)); gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[1.2, 1], hspace=0.4, wspace=0.3)
    # (a) evidence structure
    ax = fig.add_subplot(gs[0, :]); ax.set_xlim(0, 10); ax.set_ylim(-0.7, 6.6); ax.axis("off"); panel(ax, "a")
    def box(x, y, w, h, text, fc, fs=7.2, ec="#444444"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=fc, ec=ec, lw=0.8)); ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, wrap=True)
    rows = [("Bulk transcriptomes\n10 cohorts: 4 ATAA (63 vs 52), 6 TAAD (41 vs 35)\n+ GSE190635 (flagged, sensitivity only)", "Unit: patient within cohort;\ncohorts compared at cohort level", "F1 primary: C, I, R (Holm)\nTAAD effect vs ATAA effect", "#dbe9f6"),
            ("scRNA-seq, patient x cell-type pseudo-bulk\nGSE155468 (ATAA 8 vs 3; 7 vs 3 with >= 20 contractile SMC), GSE213740 (TAAD 6 vs 3),\nGSE189795 (acute TAAD 5 vs 4; no contractile-SMC cluster recovered)", "Unit: patient (>= 20 cells of the type)\nlabel-free scores within cell type", "F2 primary: R in contractile SMC (exact enumeration;\none test per dataset; not estimable if the cell type is missing)", "#fde4d8"),
            ("Within-patient regions\nGSE140947: aneurysm belly vs neck (6 pairs)", "Unit: within-patient difference", "F3 primary: delta R\n(exact sign-flip, 64 assignments)", "#e6f2e0"),
            ("Stage and wall layer\nGSE222318 (acute/subacute/chronic; sampling position differs)\nGSE318877 (3 dissection vs 4 dilatation; 4 medial layers)", "Unit: patient; patient x layer", "Estimation only\n(no inferential test; CI and raw scores)", "#f3e6f7"),
            ("Ordered dilatation within one cohort\nGSE26155: No (31) < Borderline (6) < Yes (22)", "Unit: patient (intima-media, TAV)", "Exploratory trend (Jonckheere-Terpstra,\npermutation P; BH)", "#f7f3d9")]
    y0 = 5.1
    for i, (a_, b_, c_, fc) in enumerate(rows):
        y = y0 - i * 1.12; box(0.05, y, 4.3, 1.0, a_, fc, fs=6.5); box(4.5, y, 2.6, 1.0, b_, "#ffffff", fs=6.5); box(7.25, y, 2.7, 1.0, c_, "#ffffff", fs=6.5)
    ax.text(2.2, 6.2, "Data structure", ha="center", fontsize=8, fontweight="bold"); ax.text(5.8, 6.2, "Unit of analysis", ha="center", fontsize=8, fontweight="bold"); ax.text(8.6, 6.2, "Inference status (plan lock after exploration)", ha="center", fontsize=8, fontweight="bold")
    ax.text(5, -0.45, "Programme scores (Tier 1, external definitions): C = SMC contractile module; I = mean of six injury-response modules (hypoxia, glycolysis, oxidative stress/NFE2L2, MYC-ribosome, p53-DNA damage, NF-kB/IL-6); R = I + C", ha="center", fontsize=7, style="italic")
    # (b) concordance heatmap
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "b"); order = ATAA_C + TAAD_C + ["GSE190635"]; M = rho.loc[order, order]
    sns.heatmap(M, cmap="RdBu_r", vmin=-0.7, vmax=0.7, center=0, annot=True, fmt=".2f", annot_kws={"size": 5.5}, ax=ax, cbar_kws=dict(label="Spearman rho of log2FC profiles", shrink=0.8), square=True)
    ax.set_xticklabels([f"{c}" for c in order], rotation=90, fontsize=6.5); ax.set_yticklabels(order, rotation=0, fontsize=6.5)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels(): lbl.set_color(COL["ATAA"] if lbl.get_text() in ATAA_C else (COL["TAAD"] if lbl.get_text() in TAAD_C else "#555555"))
    ax.axhline(4, color="k", lw=1); ax.axvline(4, color="k", lw=1); ax.axhline(10, color="k", lw=0.6, ls="--"); ax.axvline(10, color="k", lw=0.6, ls="--"); ax.set_title("Cross-cohort concordance of case-vs-control profiles")
    # (c) consensus and robustness
    ax = fig.add_subplot(gs[1, 1]); panel(ax, "c"); ax.axis("off")
    sub = pd.read_csv(os.path.join(R, "meta_v2", "TAAD_4of6_subsets.csv"))
    cells = [["REM consensus DEGs", f"{S['ATAA']['consensus_v2_rem']}", f"{S['TAAD']['consensus_v2_rem']}", f"{int(sub.n_consensus_shared_universe.min())}-{int(sub.n_consensus_shared_universe.max())}"],
             ["LOCO retention of consensus\n(median)", f"{S['ATAA']['loco_retention_median']:.2f}", f"{S['TAAD']['loco_retention_median']:.2f}", f"{sub.loco_retention_min.median():.2f} (min)"],
             ["Median pairwise rho", f"{S['TAAD_4of6']['ATAA_median_rho']:.2f}", f"{np.median(sub.median_rho):.2f}", f"{sub.median_rho.min():.2f}-{sub.median_rho.max():.2f}"],
             ["Consensus under HKSJ", f"{S['ATAA']['consensus_hksj']}", f"{S['TAAD']['consensus_hksj']}", "-"],
             ["Disease-difference genes\n(enriched in ATAA / TAAD)", f"{S['classes_v2']['ATAA_enriched'] + S['classes_v2']['ATAA_enriched_equivTAAD']}", f"{S['classes_v2']['TAAD_enriched'] + S['classes_v2']['TAAD_enriched_equivATAA']}", f"{S['classes_v2']['shared_concordant']} shared\n{S['classes_v2']['discordant']} discordant"]]
    tab = ax.table(cellText=cells, colLabels=["", "ATAA\n(4 cohorts)", "TAAD\n(6 cohorts)", "TAAD 4-of-6\nsubsets"], colWidths=[0.4, 0.2, 0.2, 0.2], loc="upper center", cellLoc="left")
    tab.auto_set_font_size(False); tab.set_fontsize(6.5)
    for (r_, c_), cell in tab.get_celld().items():
        cell.set_edgecolor("#dddddd"); cell.set_height(0.15)
        if r_ == 0: cell.set_text_props(fontweight="bold", color=(COL["ATAA"] if c_ == 1 else COL["TAAD"] if c_ in (2, 3) else "k"))
    ax.set_title("Reproducibility of gene-level findings (REM primary)", fontsize=8.5, loc="left")
    ax.text(0, -0.02, "Gene-level asymmetry is a statement about statistical reproducibility across\nthe included cohorts, not about biological stability.", transform=ax.transAxes, fontsize=6.5, style="italic", va="top")
    save(fig, "Fig1_evidence_structure")

# ====================================================================== Figure 2
def fig2():
    E = pd.read_csv(os.path.join(R, "v21_bulk", "bulk_effects_per_cohort.csv")); P = pd.read_csv(os.path.join(R, "v21_bulk", "bulk_pooled_and_disease_difference.csv")).set_index("measure")
    Cn = pd.read_csv(os.path.join(R, "v21_conditioned", "conditioned_sensitivity_pooled.csv"))
    fig = plt.figure(figsize=(12, 8.2)); gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], hspace=0.5, wspace=0.35)
    order = ATAA_C + TAAD_C
    for j, k in enumerate(["C", "I", "R"]):
        ax = fig.add_subplot(gs[0, j]); panel(ax, "abc"[j])
        d = E[E.measure == k].set_index("cohort").loc[order]; y = np.arange(len(order))[::-1]
        for yi, (c, r) in zip(y, d.iterrows()):
            col = COL[r.disease]; ax.errorbar(r.g, yi, xerr=[[r.g - r.ci_low], [r.ci_high - r.g]], fmt="s", color=col, ms=4, capsize=2, lw=0.9)
        for dis, yy in [("ATAA", -1.0), ("TAAD", -2.0)]:
            g, lo, hi = P.loc[k, f"g_{dis}"], P.loc[k, f"ci_low_{dis}"], P.loc[k, f"ci_high_{dis}"]
            ax.plot([lo, g, hi, g, lo], [yy, yy + 0.25, yy, yy - 0.25, yy], color=COL[dis], lw=1); ax.fill([lo, g, hi, g], [yy, yy + 0.25, yy, yy - 0.25], color=COL[dis], alpha=0.5)
        ax.axvline(0, color="k", lw=0.6); ax.axhline(len(TAAD_C) - 0.5, color="#aaaaaa", lw=0.6, ls=":")
        ax.set_yticks(list(y) + [-1, -2]); ax.set_yticklabels(order + ["ATAA pooled (REM)", "TAAD pooled (REM)"])
        for lbl in ax.get_yticklabels(): lbl.set_color(COL["ATAA"] if ("ATAA" in lbl.get_text() or lbl.get_text() in ATAA_C) else COL["TAAD"])
        ax.set_xlabel("Hedges' g, disease vs control (95 % CI)"); ax.set_title(LAB[k])
        holm = P.loc[k, "p_diff_holm"]; ax.text(0.02, 0.02, f"TAAD - ATAA difference: {P.loc[k, 'diff_TAAD_minus_ATAA']:+.2f} [{P.loc[k, 'diff_TAAD_minus_ATAA'] - 1.96 * P.loc[k, 'se_diff']:.2f}, {P.loc[k, 'diff_TAAD_minus_ATAA'] + 1.96 * P.loc[k, 'se_diff']:.2f}]\nHolm P = {holm:.3f}" if holm >= 0.001 else f"TAAD - ATAA difference: {P.loc[k, 'diff_TAAD_minus_ATAA']:+.2f}\nHolm P < 0.001", transform=ax.transAxes, fontsize=6.8, va="bottom", bbox=dict(fc="white", ec="#cccccc", lw=0.5))
        ax.set_ylim(-2.8, len(order) - 0.3)
    # (d) module heatmap
    ax = fig.add_subplot(gs[1, 0:2]); panel(ax, "d")
    mods = ["SMC_contractile", "Calcium_handling", "Cell_matrix_adhesion", "ECM_collagen", "Matrix_degradation", "Hypoxia", "Glycolysis", "Oxidative_stress_NFE2L2", "MYC_ribosome", "p53_DNA_damage", "NFkB_IL6_inflammation", "IFN_alpha_response", "IFNG_specific", "MHCII_antigen_presentation", "sig_SMC_contractile", "sig_Endothelial", "sig_Fibroblast", "sig_Inflammatory_myeloid", "sig_Macrophage", "sig_T_cell", "sig_NK", "sig_B_cell", "sig_Plasma"]
    H = P.loc[mods, ["g_ATAA", "g_TAAD", "diff_TAAD_minus_ATAA"]].T; H.index = ["ATAA vs control", "TAAD vs control", "TAAD - ATAA"]
    ann = H.round(2).astype(str)
    for m in mods:
        for i, dis in enumerate(["ATAA", "TAAD"]):
            if P.loc[m, f"ci_low_{dis}"] > 0 or P.loc[m, f"ci_high_{dis}"] < 0: ann.loc[H.index[i], m] += "*"
        q = P.loc[m, "p_diff_bh"]
        if q < 0.05: ann.loc["TAAD - ATAA", m] += "#"
    sns.heatmap(H, cmap="RdBu_r", vmin=-2.3, vmax=2.3, center=0, annot=ann.values, fmt="", annot_kws={"size": 5.6}, ax=ax, cbar_kws=dict(label="Hedges' g (REM pooled)", shrink=0.6, pad=0.01))
    ax.set_xticklabels([m.replace("_", " ").replace("sig ", "signature: ") for m in mods], rotation=60, ha="right", fontsize=6.5); ax.set_yticklabels(H.index, rotation=0)
    ax.set_title("Tier 1 modules and control-cell signatures: pooled effects per disease and their difference (* CI excludes 0; # BH < 0.05 for the difference)", fontsize=8)
    # (e) conditioned sensitivity
    ax = fig.add_subplot(gs[1, 2]); panel(ax, "e")
    show = ["C", "I", "R", "Glycolysis", "Oxidative_stress_NFE2L2", "MYC_ribosome", "p53_DNA_damage", "NFkB_IL6_inflammation", "IFN_alpha_response", "MHCII_antigen_presentation", "ECM_collagen"]
    y = np.arange(len(show))[::-1]
    for dis, off, mk in [("ATAA", 0.18, "o"), ("TAAD", -0.18, "s")]:
        d = Cn[Cn.disease == dis].set_index("measure").loc[show]
        lo_u = d.ci_unc.str.split(",").str[0].astype(float); hi_u = d.ci_unc.str.split(",").str[1].astype(float); lo_c = d.ci_cond.str.split(",").str[0].astype(float); hi_c = d.ci_cond.str.split(",").str[1].astype(float)
        ax.errorbar(d.pooled_beta_unconditioned, y + off, xerr=[d.pooled_beta_unconditioned - lo_u, hi_u - d.pooled_beta_unconditioned], fmt=mk, color=COL[dis], mfc="white", ms=4, capsize=1.5, lw=0.7, label=f"{dis} unconditioned")
        ax.errorbar(d.pooled_beta_conditioned, y + off - 0.0, xerr=[d.pooled_beta_conditioned - lo_c, hi_c - d.pooled_beta_conditioned], fmt=mk, color=COL[dis], ms=4, capsize=1.5, lw=0.7, alpha=0.9, label=f"{dis} conditioned on 3 cell signatures")
    ax.set_yticks(y); ax.set_yticklabels([s.replace("_", " ") for s in show]); ax.axvline(0, color="k", lw=0.6); ax.set_xlabel("standardised disease coefficient (REM pooled, 95 % CI)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=5.8, frameon=False, ncol=2); ax.set_title("Composition-signature-conditioned sensitivity", fontsize=8)
    save(fig, "Fig2_bulk_programmes")

# ====================================================================== Figure 3
def fig3():
    S = pd.read_csv(os.path.join(R, "v21_sc", "sc_patient_celltype_scores.csv")); T = pd.read_csv(os.path.join(R, "v21_sc", "sc_patient_level_tests.csv"))
    fig = plt.figure(figsize=(12.5, 8)); gs = fig.add_gridspec(2, 4, hspace=0.6, wspace=0.5)
    for j, (ds, dis) in enumerate([("GSE155468", "ATAA"), ("GSE213740", "TAAD")]):
        d = S[(S.dataset == ds) & (S.celltype == "SMC_contractile")]
        ax = fig.add_subplot(gs[0, j * 2]); panel(ax, "ac"[j], 1.12)
        for grp, dd in d.groupby("group"):
            ax.scatter(dd.C, dd.I, s=np.clip(dd.n_cells / 25, 12, 120), color=COL.get(grp, COL[dis]), alpha=0.85, edgecolor="k", lw=0.4, label=f"{grp} (n = {len(dd)})")
        lim = max(abs(d.C).max(), abs(d.I).max()) * 1.15; ax.plot([-lim, lim], [lim, -lim], color="#999999", lw=0.7, ls="--"); ax.text(lim * 0.55, -lim * 0.95, "R = 0", fontsize=6.5, color="#777777")
        ax.axhline(0, color="#dddddd", lw=0.6); ax.axvline(0, color="#dddddd", lw=0.6); ax.set_xlabel("C (contractile)"); ax.set_ylabel("I (injury composite)"); ax.legend(frameon=False, loc="upper right"); ax.set_title(f"{ds} ({dis}): per-patient C vs I, contractile SMC\n(point size = number of cells)", fontsize=7.5, pad=8)
        ax = fig.add_subplot(gs[0, j * 2 + 1]); panel(ax, "bd"[j], 1.12)
        long = d.melt(id_vars=["group"], value_vars=["C", "I", "R"], var_name="measure", value_name="score"); sns.stripplot(data=long, x="measure", y="score", hue="group", dodge=True, palette={dis: COL[dis], "Control": COL["Control"]}, size=4.5, ax=ax, jitter=0.12, edgecolor="k", linewidth=0.3)
        sns.pointplot(data=long, x="measure", y="score", hue="group", dodge=0.4, palette={dis: COL[dis], "Control": COL["Control"]}, errorbar=None, markers="_", markersize=14, linestyle="none", ax=ax, legend=False)
        ax.axhline(0, color="#dddddd", lw=0.6); ax.legend(frameon=False, fontsize=6.5); ax.set_xlabel(""); ax.set_ylabel("label-free score")
        t = T[(T.dataset == ds) & (T.celltype == "SMC_contractile") & (T.measure.isin(["C", "I", "R"]))].set_index("measure")
        txt = "\n".join([f"{k}: g = {t.loc[k, 'g']:+.2f} [{t.loc[k, 'ci_low']:.2f}, {t.loc[k, 'ci_high']:.2f}]" + (f", exact P = {t.loc[k, 'p_exact']:.3f} (primary)" if k == "R" else "") for k in ["C", "I", "R"]])
        ax.text(0.02, 0.98, txt, transform=ax.transAxes, fontsize=6.3, va="top", bbox=dict(fc="white", ec="#cccccc", lw=0.5)); ax.set_title(f"{ds}: {dis} vs control, contractile SMC", fontsize=7.5, pad=8)
    # (e) forest by cell type
    cts = ["SMC_contractile", "SMC_modulated", "Fibromyocyte", "Fibroblast", "Endothelial", "Macrophage", "Inflammatory_myeloid", "T_cell", "NK"]
    for j, k in enumerate(["C", "I", "R"]):
        ax = fig.add_subplot(gs[1, j]); panel(ax, "efg"[j]); y = np.arange(len(cts))[::-1]
        for ds, dis, off in [("GSE155468", "ATAA", 0.17), ("GSE213740", "TAAD", -0.17)]:
            d = T[(T.dataset == ds) & (T.measure == k)].set_index("celltype").reindex(cts)
            ax.errorbar(d.g, y + off, xerr=[np.nan_to_num(d.g - d.ci_low), np.nan_to_num(d.ci_high - d.g)], fmt="o" if dis == "ATAA" else "s", color=COL[dis], ms=3.8, capsize=1.5, lw=0.7, label=f"{ds} ({dis} vs control)")
            for yi, (ct, r) in zip(y + off, d.iterrows()):
                if not np.isnan(r.g): ax.text(3.7, yi, f"{int(r.n_case)}v{int(r.n_ctrl)}", fontsize=5.2, va="center", color=COL[dis])
        ax.axvline(0, color="k", lw=0.6); ax.set_yticks(y); ax.set_yticklabels([c.replace("_", " ") for c in cts]); ax.set_xlim(-5, 5); ax.set_xlabel(f"Hedges' g, {k} (95 % CI, clipped)"); ax.set_title(LAB[k], fontsize=8)
        if j == 0: ax.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1)
    ax = fig.add_subplot(gs[1, 3]); ax.axis("off"); panel(ax, "h")
    chk = pd.read_csv(os.path.join(R, "v21_sc", "sc_primary_fullpipeline_check.csv"))
    ax.text(0, 1, "Primary tests (F2), contractile SMC, R = I + C\n\n" + "\n".join([f"{r.dataset}: {int(r.n_case)} vs {int(r.n_ctrl)} patients\n   g = {r.g:+.2f} [{r.ci_low:.2f}, {r.ci_high:.2f}]\n   exact P = {r.p_exact:.3f} ({int(r.n_assignments)} assignments)" for _, r in T[T.primary].iterrows()]) +
            "\n\nFull-pipeline check (control-referenced z,\npermutation of labels through scoring):\n" + "\n".join([f"{r.dataset}: P = {r.p_fullpipeline_exact:.3f}" for _, r in chk.iterrows()]) + "\n\nAll cell-type estimates: limited precision\n(3 controls per dataset; CI from bootstrap).", transform=ax.transAxes, fontsize=6.8, va="top", family="monospace")
    save(fig, "Fig3_sc_patient_level")

# ====================================================================== Figure 4 (stage)
def fig4():
    f1 = os.path.join(R, "v21_stage", "GSE222318_patient_celltype_scores.csv"); f2 = os.path.join(R, "v21_stage", "GSE189795_patient_celltype_scores.csv")
    if not (os.path.exists(f1) and os.path.exists(f2) and os.path.exists(os.path.join(R, "v21_stage", "GSE189795_tests.csv"))): print("Fig 4 skipped: stage results not yet available"); return
    S1 = pd.read_csv(f1); E1 = pd.read_csv(os.path.join(R, "v21_stage", "GSE222318_stage_estimates.csv")); S2 = pd.read_csv(f2); E2 = pd.read_csv(os.path.join(R, "v21_stage", "GSE189795_tests.csv"))
    fig = plt.figure(figsize=(12, 7.4)); gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.4)
    d = S1[S1.celltype == "SMC_contractile"].copy(); order = ["Control", "Acute", "Subacute", "Chronic"]; d["group"] = pd.Categorical(d.group, order)
    mk = {p: m for p, m in zip(sorted(d.Position.unique()), ["o", "s", "^", "D"])}
    for j, k in enumerate(["C", "I", "R"]):
        ax = fig.add_subplot(gs[0, j]); panel(ax, "abc"[j])
        for i, grp in enumerate(order):
            dd = d[d.group == grp]
            for _, r in dd.iterrows(): ax.scatter(i + np.random.default_rng(abs(hash(r["sample"])) % 1000).uniform(-0.18, 0.18), r[k], marker=mk[r.Position], color=COL[grp], s=28, edgecolor="k", lw=0.4)
            if len(dd): ax.hlines(dd[k].mean(), i - 0.3, i + 0.3, color=COL[grp], lw=1.5)
        ax.set_xticks(range(4)); ax.set_xticklabels([f"{g}\n(n = {int((d.group == g).sum())})" for g in order]); ax.axhline(0, color="#dddddd", lw=0.6); ax.set_ylabel("label-free score, contractile SMC"); ax.set_title(f"GSE222318 within-cohort stages: {LAB[k]}", fontsize=8)
        if j == 0:
            from matplotlib.lines import Line2D
            ax.legend(handles=[Line2D([], [], marker=m, color="w", markerfacecolor="#888888", markeredgecolor="k", markersize=6, label=f"position: {p}") for p, m in mk.items()], frameon=False, fontsize=6, loc="lower left")
    # (d) contrasts
    ax = fig.add_subplot(gs[1, 0:2]); panel(ax, "d")
    cons = ["acute_vs_control", "subacute_vs_control", "chronic_vs_control", "acute_vs_nonacute", "acute_vs_subacute_samePosition"]; lab = {"acute_vs_control": "acute vs control", "subacute_vs_control": "subacute vs control", "chronic_vs_control": "chronic vs control", "acute_vs_nonacute": "acute vs non-acute (subacute + chronic)", "acute_vs_subacute_samePosition": "acute vs subacute (same sampling position)"}
    e = E1[(E1.celltype == "SMC_contractile") & (E1.measure.isin(["C", "I", "R"]))]; y = np.arange(len(cons))[::-1]
    for k, off, col in [("C", 0.22, "#2ca02c"), ("I", 0.0, "#ff7f0e"), ("R", -0.22, "#1f1f1f")]:
        dd = e[e.measure == k].set_index("contrast").reindex(cons); ax.errorbar(dd.g, y + off, xerr=[np.nan_to_num(dd.g - dd.ci_low), np.nan_to_num(dd.ci_high - dd.g)], fmt="o", color=col, ms=4, capsize=1.5, lw=0.8, label=LAB[k])
        for yi, (c, r) in zip(y + off, dd.iterrows()):
            if not np.isnan(r.g): ax.text(9.7, yi, f"{int(r.n_A)} vs {int(r.n_B)}; {r.positions_A} vs {r.positions_B}".replace("Intima-medium", "IM").replace("Whole", "W"), fontsize=5.3, va="center", color=col)
    ax.axvline(0, color="k", lw=0.6); ax.set_yticks(y); ax.set_yticklabels([lab[c] for c in cons]); ax.set_xlim(-6, 13.5); ax.set_xlabel("Hedges' g (95 % CI, clipped); estimation only, no inferential test; IM = intima-media, W = whole wall"); ax.legend(frameon=False, fontsize=6.5, loc="lower left"); ax.set_title("GSE222318 contractile SMC: stage contrasts (sampling position confounded with stage)", fontsize=8)
    # (e) GSE189795: contractile SMC if recovered in >= 4 patients, otherwise modulated SMC as an exploratory substitute (primary test not estimable)
    ax = fig.add_subplot(gs[1, 2]); panel(ax, "e", 1.16)
    ct = "SMC_contractile" if (S2.celltype == "SMC_contractile").sum() >= 4 else "SMC_modulated"; d2 = S2[S2.celltype == ct]
    long = d2.melt(id_vars=["group"], value_vars=["C", "I", "R"], var_name="measure", value_name="score"); sns.stripplot(data=long, x="measure", y="score", hue="group", dodge=True, palette={"Acute": COL["TAAD"], "Control": COL["Control"]}, size=4.5, ax=ax, jitter=0.12, edgecolor="k", linewidth=0.3)
    sns.pointplot(data=long, x="measure", y="score", hue="group", dodge=0.4, palette={"Acute": COL["TAAD"], "Control": COL["Control"]}, errorbar=None, markers="_", markersize=14, linestyle="none", ax=ax, legend=False)
    t = E2[(E2.celltype == ct) & (E2.measure.isin(["C", "I", "R"]))].set_index("measure")
    tag = "primary" if ct == "SMC_contractile" else "exploratory substitute"
    txt = "\n".join([f"{k}: g = {t.loc[k, 'g']:+.2f} [{t.loc[k, 'ci_low']:.2f}, {t.loc[k, 'ci_high']:.2f}], exact P = {t.loc[k, 'p_exact']:.3f}" + (f" ({tag})" if k == "R" else "") for k in ["C", "I", "R"]])
    ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=5.8, va="bottom", bbox=dict(fc="white", ec="#cccccc", lw=0.5)); ax.axhline(0, color="#dddddd", lw=0.6); ax.legend(frameon=False, fontsize=6.5, loc="upper right"); ax.set_xlabel(""); ax.set_ylabel(f"label-free score, {ct.replace('_', ' ')}")
    if ct == "SMC_contractile": ax.set_title(f"GSE189795 acute TAAD ({int(t.loc['R', 'n_acute'])}) vs control ({int(t.loc['R', 'n_ctrl'])}), contractile SMC", fontsize=8)
    else: ax.set_title(f"GSE189795 acute TAAD ({int(t.loc['R', 'n_acute'])}) vs control ({int(t.loc['R', 'n_ctrl'])})\nno contractile-SMC cluster; modulated SMC (exploratory)", fontsize=7.5, pad=10)
    save(fig, "Fig4_stage")

# ====================================================================== Figure 5 (within patient)
def fig5():
    reg = pd.read_csv(os.path.join(R, "v21_within", "GSE140947_region_scores.csv"), index_col=0); B1 = pd.read_csv(os.path.join(R, "v21_within", "GSE140947_paired_tests.csv"))
    lay = pd.read_csv(os.path.join(R, "v21_within", "GSE318877_patient_layer_scores.csv"), index_col=0); O = pd.read_csv(os.path.join(R, "v21_within", "GSE26155_ordered_scores.csv"), index_col=0)
    Ot = pd.read_csv(os.path.join(R, "v21_within", "GSE26155_ordered_trend_bh.csv")).set_index("measure")
    fig = plt.figure(figsize=(12, 8)); gs = fig.add_gridspec(2, 6, hspace=0.55, wspace=0.7)
    # (a-c) GSE140947 paired
    for j, k in enumerate(["C", "I", "R"]):
        ax = fig.add_subplot(gs[0, j * 2:j * 2 + 2]); panel(ax, "abc"[j], 1.14)
        for grp, r1, r2, col, x0 in [("ATAA", "AN_NECK", "AN_BELLY", COL["ATAA"], 0), ("Control", "DISTAL_ASC", "MID_ASC", COL["Control"], 3)]:
            a = reg[(reg.group == grp) & (reg.region == r1)].set_index("subject")[k]; b = reg[(reg.group == grp) & (reg.region == r2)].set_index("subject")[k]; common = a.index.intersection(b.index)
            for s in common: ax.plot([x0, x0 + 1], [a[s], b[s]], color=col, lw=0.9, marker="o", ms=3.5, alpha=0.85)
        ax.set_xticks([0, 1, 3, 4]); ax.set_xticklabels(["neck", "belly", "distal", "mid"]); ax.text(0.2, 1.01, "aneurysm (6 patients)", transform=ax.transAxes, ha="center", va="bottom", fontsize=7, color=COL["ATAA"]); ax.text(0.8, 1.01, "donor (6)", transform=ax.transAxes, ha="center", va="bottom", fontsize=7, color=COL["Control"])
        r = B1[(B1.contrast == "aneurysm_belly_minus_neck") & (B1.measure == k)].iloc[0]
        ax.text(0.02, 0.02, f"belly - neck: {r.mean_delta:+.2f} [{r.ci_low:.2f}, {r.ci_high:.2f}]\nexact sign-flip P = {r.p_exact_signflip:.3f}" + (" (primary, F3)" if k == "R" else ""), transform=ax.transAxes, fontsize=6.5, va="bottom", bbox=dict(fc="white", ec="#cccccc", lw=0.5))
        ax.set_ylabel("label-free score"); ax.set_title(f"GSE140947 within-patient regions: {LAB[k]}", fontsize=8, pad=16)
    # (d) GSE318877 layers
    ax = fig.add_subplot(gs[1, 0:2]); panel(ax, "d"); order = ["I", "L", "M", "A"]
    for k, ls in [("C", "-"), ("I", "--")]:
        for sb, d in lay.groupby("subject"):
            d = d.set_index("layer").reindex(order); col = COL["TAAD"] if d.group.dropna().iloc[0] == "TAAD" else COL["ATAA"]; ax.plot(range(4), d[k], color=col, ls=ls, lw=0.9, marker="o" if k == "C" else "^", ms=3, alpha=0.8)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color=COL["TAAD"], label="dissection (3)"), Line2D([], [], color=COL["ATAA"], label="dilatation (4)"), Line2D([], [], color="k", ls="-", marker="o", ms=3, label="C"), Line2D([], [], color="k", ls="--", marker="^", ms=3, label="I")], frameon=False, fontsize=6.3, ncol=2)
    ax.set_xticks(range(4)); ax.set_xticklabels(["luminal (I)", "L", "M", "adventitial (A)"]); ax.axhline(0, color="#dddddd", lw=0.6); ax.set_ylabel("label-free score (28 patient-layer profiles)"); ax.set_title("GSE318877 medial layers within patient (estimation only)", fontsize=8)
    # (e) GSE26155 ordered
    show = ["C", "I", "IFN_alpha_response", "MHCII_antigen_presentation", "sig_T_cell", "ECM_collagen"]
    ax = fig.add_subplot(gs[1, 2:6]); panel(ax, "e")
    long = O.melt(id_vars=["dilation"], value_vars=show, var_name="measure", value_name="score"); long["dilation"] = pd.Categorical(long.dilation, ["No", "Borderline", "Yes"])
    sns.boxplot(data=long, x="measure", y="score", hue="dilation", palette=["#c7d7e8", "#8fb3d9", COL["ATAA"]], ax=ax, fliersize=0, linewidth=0.6, width=0.7); sns.stripplot(data=long, x="measure", y="score", hue="dilation", dodge=True, palette=["#c7d7e8", "#8fb3d9", COL["ATAA"]], size=2.2, ax=ax, edgecolor="k", linewidth=0.2, legend=False)
    ax.set_xticks(range(len(show))); ax.set_xticklabels(["C", "I", "IFN-alpha response", "MHC-II antigen\npresentation", "T-cell signature", "ECM-collagen"], fontsize=6.8); ax.set_xlabel(""); ax.set_ylabel("label-free score"); ax.axhline(0, color="#dddddd", lw=0.6)
    for i, m in enumerate(show): ax.text(i, ax.get_ylim()[1] * 0.98, f"rho {Ot.loc[m, 'rho_spearman']:+.2f}\nP {Ot.loc[m, 'p_perm']:.3f}\nBH {Ot.loc[m, 'p_bh']:.3f}", ha="center", va="top", fontsize=5.8)
    ax.legend(title="dilatation (No 31 / Borderline 6 / Yes 22)", frameon=False, fontsize=6.3, title_fontsize=6.3, loc="lower left"); ax.set_title("GSE26155 ordered dilatation (intima-media, TAV): exploratory trend", fontsize=8)
    save(fig, "Fig5_within_patient")

# ====================================================================== Figure 6 (programme organisation)
def fig6():
    from matplotlib.lines import Line2D
    E = pd.read_csv(os.path.join(R, "v21_bulk", "bulk_effects_per_cohort.csv")); T = pd.read_csv(os.path.join(R, "v21_sc", "sc_patient_level_tests.csv")); P = pd.read_csv(os.path.join(R, "v21_bulk", "GSE318877_patient_estimates.csv")).set_index("measure")
    pooled = pd.read_csv(os.path.join(R, "v21_bulk", "bulk_pooled_and_disease_difference.csv")).set_index("measure")
    B = pd.read_csv(os.path.join(R, "v21_bulk", "bulk_patient_scores.csv"), index_col=0); S = pd.read_csv(os.path.join(R, "v21_sc", "sc_patient_celltype_scores.csv")); G = pd.read_csv(os.path.join(R, "v21_bulk", "GSE318877_patient_scores.csv"))
    W = pd.read_csv(os.path.join(R, "v21_within", "GSE140947_paired_tests.csv")); Cn = pd.read_csv(os.path.join(R, "v21_conditioned", "conditioned_sensitivity_pooled.csv"))
    BR = "#8c564b"; f189 = os.path.join(R, "v21_stage", "GSE189795_tests.csv"); s189 = os.path.join(R, "v21_stage", "GSE189795_patient_celltype_scores.csv")
    t189 = pd.read_csv(f189) if os.path.exists(f189) else pd.DataFrame(columns=["celltype", "measure", "g"]); t189 = t189[t189.celltype == "SMC_contractile"].set_index("measure")
    n_diss, n_dil = int((G.group == "TAAD").sum()), int((G.group == "ATAA").sum())
    fig = plt.figure(figsize=(13.5, 11)); gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.1], height_ratios=[1.15, 1], hspace=0.45, wspace=0.28)
    # (a) standardised effects (Hedges' g), disease vs control; GSE318877 (reference = dilatation) in a separate inset
    ax = fig.add_subplot(gs[0, 0]); panel(ax, "a")
    for c in ATAA_C + TAAD_C:
        d = E[E.cohort == c].set_index("measure"); dis = d.disease.iloc[0]; ax.errorbar(d.loc["C", "g"], d.loc["I", "g"], xerr=[[d.loc["C", "g"] - d.loc["C", "ci_low"]], [d.loc["C", "ci_high"] - d.loc["C", "g"]]], yerr=[[d.loc["I", "g"] - d.loc["I", "ci_low"]], [d.loc["I", "ci_high"] - d.loc["I", "g"]]], fmt="o", color=COL[dis], ms=5, capsize=1.5, lw=0.6, alpha=0.8)
        ax.annotate(c.replace("GSE", ""), (d.loc["C", "g"], d.loc["I", "g"]), fontsize=5.5, xytext=(3, 3), textcoords="offset points", color=COL[dis])
    for dis in ["ATAA", "TAAD"]: ax.scatter(pooled.loc["C", f"g_{dis}"], pooled.loc["I", f"g_{dis}"], marker="D", s=90, color=COL[dis], edgecolor="k", lw=0.8, zorder=5)
    for ds, dis in [("GSE155468", "ATAA"), ("GSE213740", "TAAD")]:
        t = T[(T.dataset == ds) & (T.celltype == "SMC_contractile")].set_index("measure"); ax.scatter(t.loc["C", "g"], t.loc["I", "g"], marker="*", s=110, color=COL[dis], edgecolor="k", lw=0.6, zorder=5); ax.annotate(f"{ds} SMC", (t.loc["C", "g"], t.loc["I", "g"]), fontsize=5.5, xytext=(-40 if dis == "TAAD" else 6, 8 if dis == "TAAD" else -12), textcoords="offset points", color=COL[dis])
    if "C" in t189.index: ax.scatter(t189.loc["C", "g"], t189.loc["I", "g"], marker="*", s=110, color=COL["TAAD"], edgecolor="k", lw=0.6, zorder=5); ax.annotate("GSE189795 SMC (acute)", (t189.loc["C", "g"], t189.loc["I", "g"]), fontsize=5.5, xytext=(4, 4), textcoords="offset points", color=COL["TAAD"])
    else: ax.text(0.02, 0.02, "GSE189795 (acute TAAD, scRNA): no contractile-SMC cluster annotated; not shown", transform=ax.transAxes, fontsize=5.5, color="#555555", va="bottom")
    lim = 4.2; ax.axhline(0, color="#dddddd", lw=0.6); ax.axvline(0, color="#dddddd", lw=0.6); ax.set_xlim(-lim, 2); ax.set_ylim(-1.5, lim); ax.set_xlabel("effect on C (contractile), Hedges' g"); ax.set_ylabel("effect on I (injury composite), Hedges' g")
    ins = ax.inset_axes([0.62, 0.60, 0.36, 0.33]); ins.errorbar(P.loc["C", "g"], P.loc["I", "g"], xerr=[[P.loc["C", "g"] - P.loc["C", "ci_low"]], [P.loc["C", "ci_high"] - P.loc["C", "g"]]], yerr=[[P.loc["I", "g"] - P.loc["I", "ci_low"]], [P.loc["I", "ci_high"] - P.loc["I", "g"]]], fmt="X", color=BR, ms=7, mec="k", mew=0.5, capsize=2, lw=0.8)
    ins.axhline(0, color="#dddddd", lw=0.6); ins.axvline(0, color="#dddddd", lw=0.6); ins.set_xlim(-5, 2.5); ins.set_ylim(-4.5, 2.5); ins.tick_params(labelsize=5, length=2, pad=1); ins.set_xlabel("g, C", fontsize=5.5, labelpad=1); ins.set_ylabel("g, I", fontsize=5.5, labelpad=1)
    ins.set_title(f"GSE318877 dissection ({n_diss}) vs dilatation ({n_dil}), 95 % CI\nreference = dilatation, not control", fontsize=5.6, color=BR, pad=2)
    H = [Line2D([], [], marker="o", ls="", color=COL["ATAA"], ms=5, label="bulk cohort, ATAA vs control (95 % CI)"), Line2D([], [], marker="o", ls="", color=COL["TAAD"], ms=5, label="bulk cohort, TAAD vs control (95 % CI)"),
         Line2D([], [], marker="D", ls="", mfc="#bbbbbb", mec="k", ms=7, label="REM pooled over bulk cohorts (per disease)"), Line2D([], [], marker="*", ls="", mfc="#bbbbbb", mec="k", ms=10, label="patient-level contractile SMC (scRNA), disease vs control"),
         Line2D([], [], marker="X", ls="", mfc=BR, mec="k", ms=7, label="GSE318877 dissection vs dilatation (inset): reference = dilatation, not control")]
    ax.legend(handles=H, frameon=False, fontsize=6.2, loc="upper left", bbox_to_anchor=(0.0, -0.1), ncol=1)
    ax.set_title("Standardised effects on C and I (Hedges' g) across datasets\n(g_C and g_I are scaled by different SDs, so delta R cannot be read from this plane; see b)", fontsize=7.5)
    # (b) within-dataset mean differences on the score scale (delta R = delta I + delta C is valid here)
    rows = []
    def add(ds, comp, ref_name, case, ref, kind, dis): rows.append(dict(dataset=ds, comparison=comp, reference=ref_name, delta_C=case.C.mean() - ref.C.mean(), delta_I=case.I.mean() - ref.I.mean(), n_case=len(case), n_ref=len(ref), kind=kind, dis=dis))
    for c in ATAA_C + TAAD_C:
        d = B[B.cohort == c]; dis = d.disease.iloc[0]; add(c, f"{dis} vs control (bulk, patient level)", "control", d[d.group == dis], d[d.group == "Control"], "bulk", dis)
    for ds, dis in [("GSE155468", "ATAA"), ("GSE213740", "TAAD")]:
        d = S[(S.dataset == ds) & (S.celltype == "SMC_contractile")]; add(ds, f"{dis} vs control (contractile SMC, patient level)", "control", d[d.group == dis], d[d.group == "Control"], "smc", dis)
    if os.path.exists(s189):
        d = pd.read_csv(s189); d = d[d.celltype == "SMC_contractile"]
        if len(d): add("GSE189795", "acute TAAD vs control (contractile SMC, patient level)", "control", d[d.group == "Acute"], d[d.group == "Control"], "smc", "TAAD")
    add("GSE318877", "dissection vs dilatation (patient level)", "dilatation", G[G.group == "TAAD"], G[G.group == "ATAA"], "318877", "TAAD")
    w = W[W.contrast == "aneurysm_belly_minus_neck"].set_index("measure"); rows.append(dict(dataset="GSE140947", comparison="aneurysm belly minus neck (within patient)", reference="neck (same patient)", delta_C=w.loc["C", "mean_delta"], delta_I=w.loc["I", "mean_delta"], n_case=int(w.loc["C", "n_pairs"]), n_ref=int(w.loc["C", "n_pairs"]), kind="within", dis="ATAA"))
    D = pd.DataFrame(rows); D["delta_R"] = D.delta_C + D.delta_I
    D[["dataset", "comparison", "reference", "delta_C", "delta_I", "delta_R", "n_case", "n_ref"]].to_csv(os.path.join(R, "v21_bulk", "programme_mean_differences.csv"), index=False)
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "b"); x0, x1, y0, y1 = -1.85, 0.6, -1.25, 2.05
    for k in np.arange(-2.5, 2.6, 0.5):
        z = abs(k) < 1e-9; ax.plot([x0, x1], [k - x0, k - x1], color="#555555" if z else "#c8c8c8", lw=0.9 if z else 0.5, ls="--" if z else ":", zorder=0)
        lab = "delta R = 0 (delta I = -delta C)" if z else f"delta R = {k:+.1f}"
        if k - x1 >= y0: ax.text(x1 + 0.03, k - x1, lab, fontsize=6 if z else 5.2, color="#555555" if z else "#999999", ha="left", va="center")   # line leaves through the right edge: label outside
        elif k - y0 >= x0: ax.text(k - y0, y0, lab + " ", fontsize=5.2, color="#999999", rotation=-45, transform_rotates_text=True, rotation_mode="anchor", ha="right", va="bottom")   # leaves through the bottom edge: label inside, along the line
    MK = {"bulk": ("o", 50), "smc": ("*", 130), "318877": ("X", 70), "within": ("^", 60)}; SUF = {"smc": " SMC", "within": " belly-neck", "318877": " diss-dil", "bulk": ""}
    POS = {"GSE147026": (-1.42, 1.9, "right"), "GSE153434": (-1.33, 1.72, "right"), "GSE52093": (-1.5, 1.3, "right"), "GSE98770": (-1.5, 0.98, "right"), "GSE267434": (-0.74, 1.5, "left"), "GSE202267": (-0.74, 1.02, "left")}
    OFF = {("GSE294606", "bulk"): (-6, -3, "right"), ("GSE213740", "smc"): (-7, 2, "right"), ("GSE155468", "smc"): (6, -8, "left"), ("GSE318877", "318877"): (-7, -8, "right"), ("GSE140947", "within"): (6, -8, "left"), ("GSE26155", "bulk"): (5, 1, "left")}
    for _, r in D.iterrows():
        col = BR if r.kind == "318877" else COL[r.dis]; m, s = MK[r.kind]; ax.scatter(r.delta_C, r.delta_I, marker=m, s=s, color=col, edgecolor="k", lw=0.5, zorder=5); lab = f"{r.dataset.replace('GSE', '')}{SUF[r.kind]}: {r.delta_R:+.2f}"
        if r.kind == "bulk" and r.dataset in POS: px, py, ha = POS[r.dataset]; ax.annotate(lab, (r.delta_C, r.delta_I), xytext=(px, py), textcoords="data", fontsize=5.5, ha=ha, va="center", color=col, zorder=6, arrowprops=dict(arrowstyle="-", color=col, lw=0.4, shrinkA=0, shrinkB=2))
        else: dx, dy, ha = OFF.get((r.dataset, r.kind), (5, 2, "left")); ax.annotate(lab, (r.delta_C, r.delta_I), xytext=(dx, dy), textcoords="offset points", fontsize=5.5, ha=ha, color=col, zorder=6)
    ax.axhline(0, color="#dddddd", lw=0.6); ax.axvline(0, color="#dddddd", lw=0.6); ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_xlabel("delta C = mean C (case) - mean C (reference), score units"); ax.set_ylabel("delta I = mean I (case) - mean I (reference), score units")
    H = [Line2D([], [], marker="o", ls="", color=COL["ATAA"], ms=5, label="bulk cohort, ATAA vs control"), Line2D([], [], marker="o", ls="", color=COL["TAAD"], ms=5, label="bulk cohort, TAAD vs control"), Line2D([], [], marker="*", ls="", mfc="#bbbbbb", mec="k", ms=10, label="contractile SMC (scRNA, patient level), disease vs control"),
         Line2D([], [], marker="^", ls="", color=COL["ATAA"], mec="k", ms=6, label="GSE140947 aneurysm belly minus neck (within patient, 6 pairs)"), Line2D([], [], marker="X", ls="", mfc=BR, mec="k", ms=7, label="GSE318877 dissection vs dilatation: reference = dilatation, not control")]
    ax.legend(handles=H, frameon=False, fontsize=6.2, loc="upper left", bbox_to_anchor=(0.0, -0.1), ncol=1)
    ax.set_title("Within-dataset mean differences on the label-free score scale (case minus reference)\niso-lines: constant delta R = delta I + delta C (dashed: delta R = 0); point labels give dataset: delta R\ncomparable within a dataset, not across datasets; GSE318877 uses dilatation, not control, as reference", fontsize=7.5)
    # (c) state comparison table; composition-conditioned sensitivity read at render time
    ax = fig.add_subplot(gs[1, :]); ax.axis("off"); ax.text(0.0, 1.0, "c", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")
    def g_(m, dis): return f"{pooled.loc[m, f'g_{dis}']:+.2f}" + ("*" if (pooled.loc[m, f'ci_low_{dis}'] > 0 or pooled.loc[m, f'ci_high_{dis}'] < 0) else "")
    def ret(m, dis): v = Cn[(Cn.disease == dis) & (Cn.measure == m)].retained; return "n/a" if len(v) == 0 else ("yes" if v.iloc[0] == "yes" else "uncertain")
    def cond(ms): return "retained after conditioning on 3 cell signatures (ATAA: " + " / ".join(ret(m, "ATAA") for m in ms) + "; TAAD: " + " / ".join(ret(m, "TAAD") for m in ms) + ")"
    met = ["Glycolysis", "Oxidative_stress_NFE2L2", "MYC_ribosome"]
    rows = [["C  SMC contractile", g_("C", "ATAA"), g_("C", "TAAD"), "yes (Holm P<0.001)", "decreases with dilatation, GSE26155 (BH 0.02); point estimates lower in SMC of both scRNA sets (exact P 0.17; n small);\n" + cond(["C"])],
            ["I  injury composite", g_("I", "ATAA"), g_("I", "TAAD"), f"borderline\n(Holm P={pooled.loc['I', 'p_diff_holm']:.2f})", "no trend with dilatation (GSE26155); SMC level: CI includes 0 in both scRNA sets;\n" + cond(["I"])],
            ["R  relative state (I + C)", g_("R", "ATAA"), g_("R", "TAAD"), f"no (P = {pooled.loc['R', 'p_diff']:.2f})", "F2 and F3 primary tests not passed; GSE318877 estimation only;\n" + cond(["R"])],
            ["Glycolysis / oxidative /\nMYC-ribosome", " / ".join([f"{pooled.loc[m, 'g_ATAA']:+.1f}" for m in met]), " / ".join([f"{pooled.loc[m, 'g_TAAD']:+.1f}" for m in met]), "yes (BH 0.01-0.02)", "no trend with dilatation;\n" + cond(met)],
            ["IFN-alpha response", g_("IFN_alpha_response", "ATAA"), g_("IFN_alpha_response", "TAAD"), f"no (BH {pooled.loc['IFN_alpha_response', 'p_diff_bh']:.2f})", "increases with dilatation (BH 0.01)"],
            ["MHC-II antigen\npresentation", g_("MHCII_antigen_presentation", "ATAA"), g_("MHCII_antigen_presentation", "TAAD"), f"borderline\n(BH {pooled.loc['MHCII_antigen_presentation', 'p_diff_bh']:.2f})", "increases with dilatation (BH 0.01);\nhigher in ATAA SMC (exact P 0.07; n small)"],
            ["T-cell signature", g_("sig_T_cell", "ATAA"), g_("sig_T_cell", "TAAD"), f"borderline\n(BH {pooled.loc['sig_T_cell', 'p_diff_bh']:.2f})", "increases with dilatation (BH 0.01)"],
            ["ECM-collagen", g_("ECM_collagen", "ATAA"), g_("ECM_collagen", "TAAD"), f"no (BH {pooled.loc['ECM_collagen', 'p_diff_bh']:.2f})", "increases with dilatation (BH 0.03)"],
            ["Endothelial / fibroblast\nsignatures", " / ".join([f"{pooled.loc[m, 'g_ATAA']:+.1f}" for m in ["sig_Endothelial", "sig_Fibroblast"]]), " / ".join([f"{pooled.loc[m, 'g_TAAD']:+.1f}" for m in ["sig_Endothelial", "sig_Fibroblast"]]), "yes\n(BH 0.002 / 0.007)", "composition level; not a within-cell finding"]]
    tab = ax.table(cellText=rows, colLabels=["programme (Tier 1)", "ATAA vs\ncontrol", "TAAD vs\ncontrol", "differs between\ndiseases?", "within-cohort / patient-level evidence; composition-signature-conditioned sensitivity (Fig. 2e)"], colWidths=[0.15, 0.09, 0.09, 0.12, 0.55], loc="upper center", cellLoc="left", bbox=[0.02, 0.15, 0.96, 0.85])
    tab.auto_set_font_size(False); tab.set_fontsize(6.3)
    for (r_, c_), cell in tab.get_celld().items():
        cell.set_edgecolor("#dddddd"); cell.PAD = 0.02
        if r_ == 0: cell.set_text_props(fontweight="bold")
    ax.text(0.02, 0.11, "* pooled 95 % CI excludes 0. Tissue-level state comparison; no temporal or causal ordering is implied. Supported: cohort-level differences in C and in some metabolic/stress modules.\nNot determined: whether the relative organisation of the two programmes differs (R and patient-level tests inconclusive).", transform=ax.transAxes, fontsize=6.3, va="top")
    save(fig, "Fig6_programme_organisation")

# ====================================================================== Supplementary
def figS():
    T2 = pd.read_csv(os.path.join(R, "v21_tier2_loco", "tier2_fold_internal_loco.csv"))
    fig, ax = plt.subplots(figsize=(9, 4.2)); order = sorted(T2.signature.unique(), key=lambda s: (s[0] != "T", s))
    sns.stripplot(data=T2, x="signature", y="g", hue="held_disease", order=order, palette=COL, dodge=True, size=4, ax=ax, edgecolor="k", linewidth=0.3)
    sns.pointplot(data=T2, x="signature", y="g", hue="held_disease", order=order, palette=COL, dodge=0.4, errorbar=None, markers="_", markersize=12, linestyle="none", ax=ax, legend=False)
    ax.axhline(0, color="k", lw=0.6); ax.set_xticklabels([s.replace("_", " ") for s in order], rotation=60, ha="right", fontsize=6.5); ax.set_ylabel("Hedges' g in held-out cohort (disease vs control,\noriented to signature direction)"); ax.set_xlabel("")
    ax.legend(title="held-out cohort disease", frameon=False, fontsize=6.5, title_fontsize=6.5); ax.set_title("Tier 2 signatures rebuilt inside each leave-one-cohort-out fold and scored in the held-out cohort", fontsize=8)
    save(fig, "FigS_tier2_fold_internal")
    Cn = pd.read_csv(os.path.join(R, "v21_conditioned", "conditioned_sensitivity_pooled.csv")); mods = [m for m in Cn.measure.unique() if m not in ("C", "I", "R")]
    fig, axes = plt.subplots(1, 2, figsize=(10, 6), sharey=True)
    for ax, dis in zip(axes, ["ATAA", "TAAD"]):
        d = Cn[Cn.disease == dis].set_index("measure").loc[mods]; y = np.arange(len(mods))[::-1]
        lo_u = d.ci_unc.str.split(",").str[0].astype(float); hi_u = d.ci_unc.str.split(",").str[1].astype(float); lo_c = d.ci_cond.str.split(",").str[0].astype(float); hi_c = d.ci_cond.str.split(",").str[1].astype(float)
        ax.errorbar(d.pooled_beta_unconditioned, y + 0.15, xerr=[d.pooled_beta_unconditioned - lo_u, hi_u - d.pooled_beta_unconditioned], fmt="o", mfc="white", color=COL[dis], ms=4, capsize=1.5, lw=0.7, label="unconditioned")
        ax.errorbar(d.pooled_beta_conditioned, y - 0.15, xerr=[d.pooled_beta_conditioned - lo_c, hi_c - d.pooled_beta_conditioned], fmt="o", color=COL[dis], ms=4, capsize=1.5, lw=0.7, label="conditioned on sig_SMC_contractile, sig_Inflammatory_myeloid, sig_T_cell")
        ax.set_yticks(y); ax.set_yticklabels([m.replace("_", " ") for m in mods]); ax.axvline(0, color="k", lw=0.6); ax.set_title(f"{dis} vs control (REM over cohorts)"); ax.set_xlabel("standardised coefficient (95 % CI)"); ax.legend(frameon=False, fontsize=6, loc="lower right")
    fig.tight_layout(); save(fig, "FigS_conditioned_sensitivity")

if __name__ == "__main__":
    import sys
    which = sys.argv[1:] or ["1", "2", "3", "4", "5", "6", "S"]
    for w in which: {"1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5, "6": fig6, "S": figS}[w]()
