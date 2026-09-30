"""Control-referenced signature scoring and cross-cohort machine learning.

Every disease sample is expressed relative to the non-diseased controls of its own cohort
(z-score per gene: (x - mean_ctrl) / sd_ctrl), which removes platform/batch offsets and makes
samples from different cohorts comparable.

1. Disease-specific signature scores (TAAD-specific / ATAA-specific / shared gene sets from the meta-analysis)
   -> score = mean z(up genes) - mean z(down genes); compared between TAAD samples and ATAA samples across cohorts.
2. Leave-one-cohort-out (LOCO) classifier TAAD vs ATAA (elastic-net logistic regression) on control-referenced z profiles;
   AUC on each held-out cohort; consensus gene weights; compact panel.
3. Independent tests: (i) GSE318877 direct cohort (within-cohort z, no controls) - TAAD vs ATAA subjects and micro-regions;
   (ii) scRNA-seq pseudo-bulk (GSE155468 ATAA vs control; GSE213740 TAAD vs control) control-referenced z.
"""
import os, json, warnings, numpy as np, pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "results", "processed"); META = os.path.join(ROOT, "results", "meta"); OUT = os.path.join(ROOT, "results", "ml"); os.makedirs(OUT, exist_ok=True)
FIG = os.path.join(ROOT, "figures"); SCD = os.path.join(ROOT, "results", "sc")
REG = json.load(open(os.path.join(PROC, "cohorts.json"))); DESIGN = json.load(open(os.path.join(META, "design.json")))
cmp = pd.read_csv(os.path.join(META, "ATAA_vs_TAAD_gene_comparison.csv"), index_col=0)
COHORTS = {c: "ATAA" for c in DESIGN["ATAA"]["cohorts"]}; COHORTS.update({c: "TAAD" for c in DESIGN["TAAD"]["cohorts"]})

def load_expr(name):
    e = pd.read_csv(os.path.join(PROC, f"{name}_expr.csv"), index_col=0); m = pd.read_csv(os.path.join(PROC, f"{name}_meta.csv"), index_col=0)
    e.columns = e.columns.astype(str); m.index = m.index.astype(str); e = e[~e.index.duplicated()]
    return e, m.loc[e.columns]

def ctrl_z(e, m, min_sd=0.1):
    ctrl = m.index[m.group == "Control"]
    mu = e[ctrl].mean(axis=1); sd = e[ctrl].std(axis=1).clip(lower=min_sd)
    # genes not expressed/variable in the cohort get NaN
    z = e.sub(mu, axis=0).div(sd, axis=0)
    z = z[(e[ctrl].std(axis=1) > 0)]
    return z.clip(-10, 10)

# ---------------------------------------------------------------- build the control-referenced sample x gene matrix (disease samples only)
Zs = []; info = []
for c, dis in COHORTS.items():
    e, m = load_expr(c); z = ctrl_z(e, m)
    dz = z[m.index[m.group != "Control"]]
    Zs.append(dz); info += [dict(sample=s, cohort=c, disease=dis) for s in dz.columns]
info = pd.DataFrame(info).set_index("sample")
Z = pd.concat(Zs, axis=1, join="outer")
genes_ok = Z.notna().mean(axis=1) >= 0.9
Z = Z[genes_ok].fillna(0.0); print("control-referenced matrix:", Z.shape, info.disease.value_counts().to_dict())
Z.T.to_csv(os.path.join(OUT, "control_referenced_z_disease_samples.csv")); info.to_csv(os.path.join(OUT, "control_referenced_z_sample_info.csv"))

# ---------------------------------------------------------------- 1. signature scores
def sig_score(zmat, up, down):
    up = [g for g in up if g in zmat.index]; down = [g for g in down if g in zmat.index]
    s = zmat.loc[up].mean(axis=0) if up else 0
    return s - (zmat.loc[down].mean(axis=0) if down else 0)
sets = {}
for cl in ["TAAD_specific", "ATAA_specific", "shared_concordant"]:
    s = cmp[cmp["class"] == cl]; ref = s.lfc_TAAD if "TAAD" in cl else s.lfc_ATAA
    sets[cl] = (s.index[ref > 0].tolist(), s.index[ref < 0].tolist())
scores = pd.DataFrame({cl: sig_score(Z, *ud) for cl, ud in sets.items()}).join(info)
scores.to_csv(os.path.join(OUT, "signature_scores_disease_samples.csv"))
print("\nSignature scores (control-referenced), median by disease:")
print(scores.groupby("disease")[list(sets)].median().round(2).to_string())
for cl in sets:
    t, a = scores.loc[scores.disease == "TAAD", cl], scores.loc[scores.disease == "ATAA", cl]
    print(f"  {cl}: TAAD vs ATAA Mann-Whitney p={stats.mannwhitneyu(t, a).pvalue:.2e}; AUC(TAAD vs ATAA)={roc_auc_score((scores.disease=='TAAD').astype(int), scores[cl]):.3f}")
# per-cohort effect (Hedges g of disease vs control z=0 baseline is implicit): report per cohort medians
print(scores.groupby(["disease", "cohort"])[list(sets)].median().round(2).to_string())

# ---------------------------------------------------------------- 2. LOCO classifier TAAD vs ATAA
X = Z.T.loc[info.index]; y = (info.disease == "TAAD").astype(int).values; groups = info.cohort.values
feat = cmp.index[(cmp["class"] != "none")].intersection(X.columns)  # consensus DEGs of either disease
Xf = X[feat]
rows = []; preds = pd.Series(index=X.index, dtype=float)
for held in np.unique(groups):
    tr, te = groups != held, groups == held
    clf = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=0.05, max_iter=5000, class_weight="balanced")
    clf.fit(Xf[tr], y[tr]); p = clf.predict_proba(Xf[te])[:, 1]; preds[Xf.index[te]] = p
    truth = y[te][0]
    rows.append(dict(held_out=held, disease=COHORTS[held], n=int(te.sum()), mean_prob_TAAD=float(p.mean()), frac_called_TAAD=float((p > 0.5).mean()), correct=float(((p > 0.5).astype(int) == truth).mean())))
loco = pd.DataFrame(rows); print("\nLOCO classifier (TAAD vs ATAA):"); print(loco.round(3).to_string(index=False))
auc_pooled = roc_auc_score(y, preds.loc[X.index]); print(f"pooled LOCO AUC = {auc_pooled:.3f}; overall accuracy = {((preds.loc[X.index] > 0.5).astype(int) == y).mean():.3f}")
loco.to_csv(os.path.join(OUT, "LOCO_classifier_results.csv"), index=False)
# final model on all cohorts -> weights
clf = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=0.05, max_iter=5000, class_weight="balanced").fit(Xf, y)
w = pd.Series(clf.coef_[0], index=feat); w = w[w != 0].sort_values()
w.to_frame("weight").join(cmp[["class", "lfc_ATAA", "lfc_TAAD", "z_diff", "fdr_diff"]]).to_csv(os.path.join(OUT, "classifier_gene_weights.csv"))
print(f"non-zero weights: {len(w)}; top TAAD-favouring: {', '.join(w.index[-15:][::-1])}; top ATAA-favouring: {', '.join(w.index[:15])}")

# ---------------------------------------------------------------- 3a. GSE318877 direct cohort: apply signature scores (within-cohort z) at subject and micro-region level
e, m = load_expr("GSE318877")
zz = e.sub(e.mean(axis=1), axis=0).div(e.std(axis=1).clip(lower=0.1), axis=0)
d_sc = pd.DataFrame({cl: sig_score(zz, *ud) for cl, ud in sets.items()}).join(m[["group"]])
print("\nGSE318877 subject-level signature scores:"); print(d_sc.round(2).to_string())
er = pd.read_csv(os.path.join(PROC, "GSE318877_regions_log2tpm.csv"), index_col=0); mr = pd.read_csv(os.path.join(PROC, "GSE318877_regions_meta.csv"), index_col=0)
er = er[~er.index.duplicated()]; er = er[(er > 1).mean(axis=1) > 0.25]
zr = er.sub(er.mean(axis=1), axis=0).div(er.std(axis=1).clip(lower=0.1), axis=0)
r_sc = pd.DataFrame({cl: sig_score(zr, *ud) for cl, ud in sets.items()}).join(mr[["group", "subject", "region"]])
r_sc.to_csv(os.path.join(OUT, "GSE318877_region_signature_scores.csv"))
# subject-level test (mean over punches per subject) and mixed-model-like permutation at subject level
subj = r_sc.groupby("subject").agg({**{cl: "mean" for cl in sets}, "group": "first"})
for cl in sets:
    t, a = subj.loc[subj.group == "TAAD", cl], subj.loc[subj.group == "ATAA", cl]
    print(f"  GSE318877 {cl}: TAAD subjects mean={t.mean():.2f}, ATAA subjects mean={a.mean():.2f}; Mann-Whitney (subject level, n=3 vs 4) p={stats.mannwhitneyu(t, a).pvalue:.3f}; region-level AUC={roc_auc_score((r_sc.group=='TAAD').astype(int), r_sc[cl]):.3f}")
# classifier probabilities in GSE318877 (within-cohort z, features filled 0 if absent)
Xd = zz.reindex(feat).fillna(0).T
pd_ = pd.Series(clf.predict_proba(Xd)[:, 1], index=Xd.index).to_frame("prob_TAAD").join(m[["group"]]); print(pd_.round(3).to_string())
pd_.to_csv(os.path.join(OUT, "GSE318877_classifier_probs.csv"))

# ---------------------------------------------------------------- 3b. scRNA pseudo-bulk (sample level) control-referenced z
pbc = pd.read_csv(os.path.join(SCD, "pseudobulk_counts_celltype_sample.csv"), index_col=0)
samples = sorted(set(c.split("|")[0] for c in pbc.columns))
pb = pd.DataFrame({s: pbc[[c for c in pbc.columns if c.startswith(s + "|")]].sum(axis=1) for s in samples})
pb = pb[(pb >= 10).sum(axis=1) >= 3]; lcpm = np.log2(pb / pb.sum(axis=0) * 1e6 + 1)
cm = pd.read_csv(os.path.join(SCD, "cell_metadata.csv"), index_col=0).drop_duplicates("sample").set_index("sample")[["dataset", "group"]]
sc_rows = []
for ds in ["GSE155468", "GSE213740"]:
    ss = cm.index[cm.dataset == ds]; mm = cm.loc[ss]; z = ctrl_z(lcpm[ss], mm)
    dz = z[mm.index[mm.group != "Control"]]
    s = pd.DataFrame({cl: sig_score(dz, *ud) for cl, ud in sets.items()}); s["dataset"] = ds; s["disease"] = mm.loc[dz.columns, "group"].values
    Xs = dz.reindex(feat).fillna(0).T; s["prob_TAAD"] = clf.predict_proba(Xs)[:, 1]
    sc_rows.append(s)
scv = pd.concat(sc_rows); scv.to_csv(os.path.join(OUT, "scRNA_pseudobulk_signature_scores.csv"))
print("\nscRNA pseudo-bulk (control-referenced) signature scores and classifier probabilities:"); print(scv.round(2).to_string())
print("  AUC TAAD vs ATAA (pseudo-bulk, classifier prob):", round(roc_auc_score((scv.disease == "TAAD").astype(int), scv.prob_TAAD), 3),
      "| TAAD-specific score AUC:", round(roc_auc_score((scv.disease == "TAAD").astype(int), scv.TAAD_specific), 3),
      "| ATAA-specific score AUC (ATAA vs TAAD):", round(roc_auc_score((scv.disease == "ATAA").astype(int), scv.ATAA_specific), 3))

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 4, figsize=(17, 4.4))
order = list(DESIGN["ATAA"]["cohorts"]) + list(DESIGN["TAAD"]["cohorts"])
pal = {"ATAA": "#2A9D8F", "TAAD": "#C44E52"}
sns.boxplot(data=scores, x="cohort", y="TAAD_specific", hue="disease", order=order, palette=pal, ax=ax[0], fliersize=0, dodge=False)
sns.stripplot(data=scores, x="cohort", y="TAAD_specific", order=order, color="k", size=2.5, ax=ax[0])
ax[0].axhline(0, c="grey", lw=0.6, ls="--"); ax[0].set_ylabel("TAAD-specific signature score\n(control-referenced z)"); ax[0].tick_params(axis="x", rotation=90, labelsize=7); ax[0].legend(fontsize=7, frameon=False); ax[0].set_title("Bulk cohorts", fontsize=10)
sns.boxplot(data=scores, x="cohort", y="ATAA_specific", hue="disease", order=order, palette=pal, ax=ax[1], fliersize=0, dodge=False)
sns.stripplot(data=scores, x="cohort", y="ATAA_specific", order=order, color="k", size=2.5, ax=ax[1])
ax[1].axhline(0, c="grey", lw=0.6, ls="--"); ax[1].set_ylabel("ATAA-specific signature score"); ax[1].tick_params(axis="x", rotation=90, labelsize=7); ax[1].get_legend().remove(); ax[1].set_title("Bulk cohorts", fontsize=10)
fpr, tpr, _ = roc_curve(y, preds.loc[X.index]); ax[2].plot(fpr, tpr, c="#C44E52", lw=2, label=f"LOCO pooled AUC={auc_pooled:.2f}")
fpr2, tpr2, _ = roc_curve((scv.disease == "TAAD").astype(int), scv.prob_TAAD); ax[2].plot(fpr2, tpr2, c="#264653", lw=2, ls="--", label=f"scRNA pseudo-bulk AUC={roc_auc_score((scv.disease == 'TAAD').astype(int), scv.prob_TAAD):.2f}")
ax[2].plot([0, 1], [0, 1], "k:", lw=0.8); ax[2].set_xlabel("1 - specificity"); ax[2].set_ylabel("sensitivity"); ax[2].legend(fontsize=7, frameon=False); ax[2].set_title("TAAD vs ATAA classifier", fontsize=10)
sns.boxplot(data=r_sc, x="group", y="TAAD_specific", palette=pal, ax=ax[3], fliersize=0); sns.stripplot(data=r_sc, x="group", y="TAAD_specific", hue="subject", size=3, ax=ax[3], palette="tab10")
ax[3].set_title("GSE318877 micro-regions (direct)", fontsize=10); ax[3].set_ylabel("TAAD-specific score (within-cohort z)"); ax[3].legend(fontsize=6, frameon=False, title="subject", title_fontsize=6)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig8_signature_ml.png"), dpi=200); fig.savefig(os.path.join(FIG, "Fig8_signature_ml.pdf")); plt.close(fig)
