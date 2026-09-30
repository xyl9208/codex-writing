"""Cell-type marker sets derived from CONTROL cells only (both scRNA datasets), for Tier 1 cell-related signatures."""
import os, warnings, pandas as pd, scanpy as sc
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "results", "sc")
a = sc.read_h5ad(os.path.join(OUT, "aorta_sc_integrated.h5ad")); a = a[a.obs.group == "Control"].copy()
order = ["SMC_contractile", "SMC_modulated", "Fibromyocyte", "Fibroblast", "Endothelial", "Macrophage", "Inflammatory_myeloid", "T_cell", "NK", "B_cell", "Plasma", "Mast"]
vc = a.obs.celltype.value_counts(); keep = [c for c in order if vc.get(c, 0) >= 30]; a = a[a.obs.celltype.isin(keep)].copy()
sc.tl.rank_genes_groups(a, "celltype", method="wilcoxon", pts=True)
rows = []
with open(os.path.join(OUT, "aorta_celltype_markers_control.gmt"), "w") as fh:
    for ct in keep:
        df = sc.get.rank_genes_groups_df(a, ct)
        df = df[(df.logfoldchanges > 1.5) & (df.pvals_adj < 0.01) & (df.pct_nz_group > 0.3) & (df.pct_nz_reference < 0.25)]
        df = df[~df.names.str.startswith(("MT-", "RPL", "RPS", "MALAT1"))].head(50)
        if len(df) >= 8: fh.write("\t".join([ct, "scRNA_control_cells_GSE155468_GSE213740"] + df.names.tolist()) + "\n")
        rows.append(dict(celltype=ct, n_control_cells=int(vc[ct]), n_markers=len(df), top=", ".join(df.names[:8])))
print(pd.DataFrame(rows).to_string(index=False))
