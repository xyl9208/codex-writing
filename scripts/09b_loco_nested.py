"""Leakage-free leave-one-cohort-out classifier (TAAD vs ATAA) -- v2.1: feature selection uses the REM consensus rule inside each training fold.

For every held-out cohort, the meta-analysis and the feature selection (consensus DEGs of either disease)
are re-computed using ONLY the remaining cohorts; the elastic-net model is then fitted on the remaining
cohorts and applied to the held-out cohort.  A second scheme holds out one ATAA and one TAAD cohort at a
time (24 pairs) so that every test fold contains both classes.
"""
import os, json, warnings, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "results", "processed"); DE = os.path.join(ROOT, "results", "de"); META = os.path.join(ROOT, "results", "meta"); OUT = os.path.join(ROOT, "results", "ml")
REG = json.load(open(os.path.join(PROC, "cohorts.json"))); DESIGN = json.load(open(os.path.join(META, "design.json")))
COHORTS = {c: "ATAA" for c in DESIGN["ATAA"]["cohorts"]}; COHORTS.update({c: "TAAD" for c in DESIGN["TAAD"]["cohorts"]})

def load_de(name):
    d = pd.read_csv(os.path.join(DE, f"{name}_de.csv"), index_col=0)
    d = d[~d.index.duplicated()].dropna(subset=["pvalue", "log2FC"])
    d["z"] = np.sign(d["log2FC"]) * stats.norm.isf(d["pvalue"].clip(1e-300) / 2); d["n"] = sum(REG[name]["groups"].values())
    return d

def consensus(names):
    """v2.1 rule: DerSimonian-Laird random-effects consensus DEGs (REM-FDR<0.05, direction >=75 %, nominal >=50 %) restricted to the given cohorts."""
    dfs = [load_de(c) for c in names]; k = len(dfs); min_c = max(2, int(np.ceil(0.75 * k)))
    genes = pd.Index(sorted(set().union(*[set(d.index) for d in dfs]))); nm = [d["cohort"].iloc[0] for d in dfs]
    Y = pd.DataFrame({n: d["log2FC"].reindex(genes) for n, d in zip(nm, dfs)}); V = pd.DataFrame({n: (d["se"].abs().replace(0, np.nan) ** 2).reindex(genes) for n, d in zip(nm, dfs)}); P = pd.DataFrame({n: d["pvalue"].reindex(genes) for n, d in zip(nm, dfs)})
    present = Y.notna() & V.notna() & np.isfinite(V); nk = present.sum(axis=1); keep = nk >= min_c; Y, V, P, present, nk = Y[keep], V[keep], P[keep], present[keep], nk[keep]
    V = V.where(present); w = 1 / V; yf = (w * Y).sum(axis=1) / w.sum(axis=1); Q = (w * (Y.sub(yf, axis=0)) ** 2).sum(axis=1); C = w.sum(axis=1) - (w ** 2).sum(axis=1) / w.sum(axis=1)
    tau2 = ((Q - (nk - 1)) / C).clip(lower=0); wr = (1 / V.add(tau2, axis=0)).where(present); yr = (wr * Y).sum(axis=1) / wr.sum(axis=1); se = np.sqrt(1 / wr.sum(axis=1))
    fdr = multipletests(2 * stats.norm.sf(np.abs(yr / se)), method="fdr_bh")[1]; sign = np.sign(yr)
    same = (np.sign(Y).eq(sign, axis=0) & present).sum(axis=1) / nk; nominal = ((P < 0.05) & np.sign(Y).eq(sign, axis=0) & present).sum(axis=1) / nk
    return set(yr.index[(fdr < 0.05) & (same >= 0.75) & (nominal >= 0.5)])

Z = pd.read_csv(os.path.join(OUT, "control_referenced_z_disease_samples.csv"), index_col=0)  # samples x genes
info = pd.read_csv(os.path.join(OUT, "control_referenced_z_sample_info.csv"), index_col=0).loc[Z.index]
y = (info.disease == "TAAD").astype(int).values; groups = info.cohort.values

def fit_predict(train_cohorts, test_cohorts):
    feat = set()
    for dis in ["ATAA", "TAAD"]:
        names = [c for c in DESIGN[dis]["cohorts"] if c in train_cohorts]
        if len(names) >= 2: feat |= consensus(names)
    feat = sorted(feat & set(Z.columns))
    tr = np.isin(groups, train_cohorts); te = np.isin(groups, test_cohorts)
    clf = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=0.05, max_iter=5000, class_weight="balanced").fit(Z.loc[tr, feat], y[tr])
    return pd.Series(clf.predict_proba(Z.loc[te, feat])[:, 1], index=Z.index[te]), len(feat)

# ---- scheme 1: leave one cohort out (nested feature selection)
rows = []; preds = pd.Series(index=Z.index, dtype=float)
for held in sorted(set(groups)):
    p, nfeat = fit_predict([c for c in COHORTS if c != held], [held]); preds[p.index] = p
    truth = y[groups == held][0]
    rows.append(dict(held_out=held, disease=COHORTS[held], n=len(p), n_features=nfeat, mean_prob_TAAD=p.mean(), correct=((p > 0.5).astype(int) == truth).mean()))
r1 = pd.DataFrame(rows); auc1 = roc_auc_score(y, preds.loc[Z.index]); acc1 = ((preds.loc[Z.index] > 0.5).astype(int) == y).mean()
print("Scheme 1 - LOCO with nested (training-only) feature selection:"); print(r1.round(3).to_string(index=False))
print(f"pooled AUC = {auc1:.3f}; accuracy = {acc1:.3f}")

# ---- scheme 2: leave one ATAA + one TAAD cohort out (both classes in every test fold)
rows2 = []; aucs = []
for a in DESIGN["ATAA"]["cohorts"]:
    for t in DESIGN["TAAD"]["cohorts"]:
        p, nfeat = fit_predict([c for c in COHORTS if c not in (a, t)], [a, t]); yt = y[np.isin(groups, [a, t])]
        auc = roc_auc_score((info.loc[p.index, "disease"] == "TAAD").astype(int), p); aucs.append(auc)
        rows2.append(dict(held_ATAA=a, held_TAAD=t, n=len(p), n_features=nfeat, AUC=auc, accuracy=((p > 0.5).astype(int) == (info.loc[p.index, "disease"] == "TAAD").astype(int)).mean()))
r2 = pd.DataFrame(rows2)
print("\nScheme 2 - leave-one-ATAA-and-one-TAAD-cohort-out (24 folds): median fold AUC = %.3f (IQR %.3f-%.3f), mean = %.3f" % (np.median(aucs), np.percentile(aucs, 25), np.percentile(aucs, 75), np.mean(aucs)))
print(r2.round(3).to_string(index=False))
r1.to_csv(os.path.join(OUT, "LOCO_nested_results.csv"), index=False); r2.to_csv(os.path.join(OUT, "LOCO_pairs_nested_results.csv"), index=False)
# 4 x 6 fold matrix figure (AUC and accuracy), no CI computed from the 24 non-independent folds
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
M = r2.pivot(index="held_ATAA", columns="held_TAAD", values="AUC"); Acc = r2.pivot(index="held_ATAA", columns="held_TAAD", values="accuracy")
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
sns.heatmap(M, annot=True, fmt=".2f", vmin=0.5, vmax=1, cmap="Blues", ax=ax[0], cbar_kws=dict(label="AUC")); ax[0].set_title("Nested LOCO (one ATAA + one TAAD cohort held out): AUC", fontsize=9)
sns.heatmap(Acc, annot=True, fmt=".2f", vmin=0.5, vmax=1, cmap="Greens", ax=ax[1], cbar_kws=dict(label="accuracy")); ax[1].set_title("accuracy (threshold 0.5)", fontsize=9)
for a_ in ax: a_.set_xlabel("held-out TAAD cohort"); a_.set_ylabel("held-out ATAA cohort"); a_.tick_params(labelsize=7)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures", "FigS_nested_LOCO_matrix.png"), dpi=200); plt.close(fig)
json.dump(dict(pooled_auc_nested=auc1, accuracy_nested=acc1, pair_auc_median=float(np.median(aucs)), pair_auc_mean=float(np.mean(aucs))), open(os.path.join(OUT, "LOCO_nested_summary.json"), "w"), indent=1)
