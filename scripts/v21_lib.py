"""Shared v2.1 functions: Tier 1 module loading, label-free standardised module scores, C / I / R programme scores,
exact enumeration permutation tests, Hedges' g with bootstrap CI."""
import os, itertools, numpy as np, pandas as pd
from scipy import stats
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T1 = os.path.join(ROOT, "results", "modules_v21", "tier1_modules.gmt")
INJURY_SIX = ["Hypoxia", "Glycolysis", "Oxidative_stress_NFE2L2", "MYC_ribosome", "p53_DNA_damage", "NFkB_IL6_inflammation"]

def load_tier1(path=T1):
    return {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(path)}

def module_scores(expr, modules, min_genes=5, detect_frac=None):
    """expr: genes x samples (log scale).  Label-free: z per gene across ALL samples; module score = mean z of members.
    detect_frac: optional label-free detectability filter (gene kept if expressed > 0 in >= detect_frac of samples)."""
    x = expr[~expr.index.duplicated()].dropna(how="all")
    if detect_frac is not None: x = x[(x > 0).mean(axis=1) >= detect_frac]
    sd = x.std(axis=1); x = x[sd > 0]; z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1), axis=0)
    out, n = {}, {}
    for k, g in modules.items():
        gg = [y for y in g if y in z.index]
        if len(gg) >= min_genes: out[k] = z.loc[gg].mean(axis=0); n[k] = len(gg)
    return pd.DataFrame(out), pd.Series(n)

def program_scores(ms):
    """C = SMC_contractile score; I = equal-weight mean of the six injury modules (each standardised across samples);
    R = I + C (relative expression state: injury response relative to contractile phenotype)."""
    d = pd.DataFrame(index=ms.index)
    d["C"] = ms["SMC_contractile"] if "SMC_contractile" in ms else np.nan
    avail = [k for k in INJURY_SIX if k in ms]
    if avail:
        zz = (ms[avail] - ms[avail].mean()) / ms[avail].std(); d["I"] = zz.mean(axis=1); d["n_injury_modules"] = len(avail)
    else: d["I"] = np.nan; d["n_injury_modules"] = 0
    d["R"] = d["I"] + d["C"]
    return d

def hedges_g(x, y, n_boot=2000, seed=0):
    x, y = np.asarray(x, float), np.asarray(y, float); nx_, ny = len(x), len(y)
    def g_of(a, b):
        sp = np.sqrt(((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) / (len(a) + len(b) - 2))
        return (np.mean(a) - np.mean(b)) / sp * (1 - 3 / (4 * (len(a) + len(b)) - 9)) if sp > 0 else np.nan
    g = g_of(x, y); rng = np.random.default_rng(seed)
    if nx_ < 2 or ny < 2: return g, np.nan, np.nan
    bs = [g_of(rng.choice(x, nx_, replace=True), rng.choice(y, ny, replace=True)) for _ in range(n_boot)]
    return g, float(np.nanpercentile(bs, 2.5)), float(np.nanpercentile(bs, 97.5))

def exact_perm_test(x, y):
    """Exact enumeration two-sided permutation test of the mean difference (SciPy convention: 2*min tail, capped at 1).
    Returns (observed diff, p, n_assignments)."""
    x, y = np.asarray(x, float), np.asarray(y, float); allv = np.concatenate([x, y]); n1 = len(x); N = len(allv)
    obs = x.mean() - y.mean(); diffs = []
    for idx in itertools.combinations(range(N), n1):
        m = np.zeros(N, bool); m[list(idx)] = True; diffs.append(allv[m].mean() - allv[~m].mean())
    diffs = np.array(diffs); eps = 1e-12
    p_ge = np.mean(diffs >= obs - eps); p_le = np.mean(diffs <= obs + eps)
    return float(obs), float(min(1.0, 2 * min(p_ge, p_le))), len(diffs)

def exact_sign_perm_paired(d):
    """Exact sign-flip permutation test (two-sided) for paired differences d (2^n assignments)."""
    d = np.asarray(d, float); n = len(d); obs = d.mean(); vals = []
    for signs in itertools.product([1, -1], repeat=n): vals.append(np.mean(d * np.array(signs)))
    vals = np.array(vals); eps = 1e-12
    return float(obs), float(min(1.0, 2 * min(np.mean(vals >= obs - eps), np.mean(vals <= obs + eps)))), len(vals)

def celltype_pseudobulk(adata, min_cells=20, layer=None):
    """Mean log-normalised expression per sample x cell type (>= min_cells). Returns (genes x 'sample|celltype' DataFrame, cell counts)."""
    X = (adata.layers[layer] if layer else adata.X).tocsr(); rows = []; counts = {}
    for (s, ct), idx in adata.obs.groupby(["sample", "celltype"]).indices.items():
        if len(idx) < min_cells: continue
        rows.append(pd.Series(np.asarray(X[idx].mean(axis=0)).ravel(), index=adata.var_names, name=f"{s}|{ct}")); counts[f"{s}|{ct}"] = len(idx)
    return pd.DataFrame(rows).T, pd.Series(counts)
