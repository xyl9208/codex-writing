"""Harmonised preprocessing of all bulk transcriptomic cohorts.

Outputs (results/processed/):
  <cohort>_expr.csv    gene x sample, log2-scale expression (arrays / FPKM / log2-CPM)
  <cohort>_counts.csv  gene x sample raw counts (RNA-seq cohorts only)
  <cohort>_meta.csv    sample metadata: group (Control / ATAA / TAAD), subject, covariates
  cohorts.json         cohort registry
"""
import gzip, re, json, os, sys, glob, tarfile, io
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEO = os.path.join(ROOT, "data", "geo")
OUT = os.path.join(ROOT, "results", "processed")
os.makedirs(OUT, exist_ok=True)

# ----------------------------------------------------------------------------------------------
# Gene symbol harmonisation (NCBI gene_info): alias -> official symbol, Entrez -> symbol
# ----------------------------------------------------------------------------------------------
gi = pd.read_csv(os.path.join(GEO, "Homo_sapiens.gene_info.gz"), sep="\t",
                 usecols=["GeneID", "Symbol", "Synonyms", "type_of_gene", "dbXrefs"], dtype=str)
gi = gi[gi.type_of_gene != "unknown"]
OFFICIAL = set(gi.Symbol)
ENTREZ2SYM = dict(zip(gi.GeneID, gi.Symbol))
ENSG2SYM = {}
for xr, s in zip(gi.dbXrefs, gi.Symbol):
    m = re.search(r"Ensembl:(ENSG\d+)", xr or "")
    if m: ENSG2SYM.setdefault(m.group(1), s)
alias_counts = {}
for syns, sym in zip(gi.Synonyms, gi.Symbol):
    if syns == "-": continue
    for a in syns.split("|"):
        alias_counts.setdefault(a, set()).add(sym)
ALIAS2SYM = {a: list(v)[0] for a, v in alias_counts.items() if len(v) == 1 and a not in OFFICIAL}
PROTEIN_CODING = set(gi.loc[gi.type_of_gene == "protein-coding", "Symbol"])

def harmonise(symbols):
    out = []
    for s in symbols:
        if not isinstance(s, str) or s in ("", "---", "nan", "NA"): out.append(None); continue
        s = s.strip()
        out.append(s if s in OFFICIAL else ALIAS2SYM.get(s, s))
    return out

def collapse_max_mean(df, symbols):
    """Collapse probes -> genes keeping the probe with the highest mean (log-scale matrices)."""
    df = df.copy(); df["__sym"] = harmonise(symbols)
    df = df.dropna(subset=["__sym"])
    df["__mean"] = df.drop(columns="__sym").mean(axis=1)
    df = df.sort_values("__mean", ascending=False).drop_duplicates("__sym")
    return df.set_index("__sym").drop(columns="__mean").sort_index()

def collapse_sum(df, symbols):
    """Collapse to genes by summing (counts / FPKM)."""
    df = df.copy(); df["__sym"] = harmonise(symbols)
    df = df.dropna(subset=["__sym"])
    return df.groupby("__sym").sum(numeric_only=True).sort_index()

def quantile_normalize(df):
    rank_mean = df.stack().groupby(df.rank(method="first").stack().astype(int)).mean()
    return df.rank(method="min").stack().astype(int).map(rank_mean).unstack()

def log2cpm(counts):
    lib = counts.sum(axis=0)
    return np.log2((counts / lib * 1e6) + 1)

def read_series_matrix(path):
    """Return (meta DataFrame indexed by GSM, expression DataFrame)."""
    meta = {}; rows = []; in_table = False; header = None
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"): in_table = True; continue
            if line.startswith("!series_matrix_table_end"): break
            if in_table:
                parts = [p.strip('"') for p in line.split("\t")]
                if header is None: header = parts
                else: rows.append(parts)
            elif line.startswith("!Sample_"):
                key, *vals = line.split("\t")
                vals = [v.strip('"') for v in vals]
                key = key[1:]
                if key == "Sample_characteristics_ch1":
                    for i, v in enumerate(vals):
                        if ":" in v:
                            k, val = v.split(":", 1)
                            meta.setdefault("char:" + k.strip(), [None] * len(vals))[i] = val.strip()
                else:
                    meta[key] = vals
    expr = pd.DataFrame(rows, columns=header).set_index(header[0])
    expr = expr.apply(pd.to_numeric, errors="coerce")
    m = pd.DataFrame(meta)
    m.index = m["Sample_geo_accession"]
    return m, expr

def save(name, expr, meta, counts=None, info=None):
    expr.to_csv(os.path.join(OUT, f"{name}_expr.csv"))
    meta.to_csv(os.path.join(OUT, f"{name}_meta.csv"))
    if counts is not None: counts.to_csv(os.path.join(OUT, f"{name}_counts.csv"))
    reg = json.load(open(os.path.join(OUT, "cohorts.json"))) if os.path.exists(os.path.join(OUT, "cohorts.json")) else {}
    reg[name] = dict(n_genes=int(expr.shape[0]), n_samples=int(expr.shape[1]),
                     groups=meta["group"].value_counts().to_dict(), counts=counts is not None, **(info or {}))
    json.dump(reg, open(os.path.join(OUT, "cohorts.json"), "w"), indent=1)
    print(f"[{name}] genes={expr.shape[0]} samples={expr.shape[1]} groups={meta['group'].value_counts().to_dict()}")

# ==============================================================================================
# ATAA cohorts
# ==============================================================================================
def gse26155():
    m, e = read_series_matrix(os.path.join(GEO, "GSE26155_series_matrix.txt.gz"))
    keep = (m["Sample_source_name_ch1"] == "Aorta intima-media") & (m["char:aorta dilation"].isin(["Yes", "No"]))
    m = m[keep]
    meta = pd.DataFrame({"group": np.where(m["char:aorta dilation"] == "Yes", "ATAA", "Control"),
                         "subject": m["char:patient id"] if "char:patient id" in m else m["Sample_title"],
                         "title": m["Sample_title"]}, index=m.index)
    # GPL5175 transcript-cluster annotation
    lines = open(os.path.join(GEO, "GPL5175_full.txt"), errors="replace").read().split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("ID\t"))
    ann = pd.read_csv(io.StringIO("\n".join(lines[start:])), sep="\t", dtype=str, quoting=3, on_bad_lines="skip")
    ann = ann[ann["category"] == "main"]
    def first_sym(ga):
        if not isinstance(ga, str) or ga == "---": return None
        syms = [chunk.split(" // ")[1].strip() for chunk in ga.split(" /// ") if len(chunk.split(" // ")) > 1]
        syms = [s for s in syms if s and s != "---"]
        return syms[0] if syms else None
    ann["symbol"] = ann["gene_assignment"].map(first_sym)
    ann = ann.dropna(subset=["symbol"]).set_index("ID")
    e = e.loc[e.index.intersection(ann.index), meta.index]
    expr = collapse_max_mean(e, ann.loc[e.index, "symbol"].values)
    save("GSE26155", expr, meta, info=dict(disease="ATAA", platform="Affymetrix Human Exon 1.0 ST (RMA)", tissue="ascending aorta intima-media", role="discovery"))

def gse140947():
    df = pd.read_csv(os.path.join(GEO, "GSE140947_bulkrna_24_samp_genecounts.csv.gz"))
    sym = df["gene_id"].str.split("_", n=1).str[1]
    counts24 = df.drop(columns="gene_id"); counts24.index = df["gene_id"]
    counts24 = collapse_sum(counts24, sym.values)
    # collapse regions to subject level (aneurysm belly+neck; donor mid+distal)
    subj = {}
    for c in counts24.columns:
        m_ = re.match(r"(AN)_(BELLY|NECK)_(\d+)", c) or re.match(r"(MID|DISTAL)_ASC_(\d+)", c)
        if c.startswith("AN"): subj[c] = f"ATAA_{c.split('_')[-1]}"
        else: subj[c] = f"Control_{c.split('_')[-1]}"
    counts = counts24.T.groupby(pd.Series(subj)).sum().T
    meta = pd.DataFrame({"group": [c.split("_")[0] for c in counts.columns], "subject": counts.columns}, index=counts.columns)
    save("GSE140947", log2cpm(counts), meta, counts=counts, info=dict(disease="ATAA", platform="RNA-seq (Illumina NextSeq 500), media; regions summed per subject", tissue="ascending aorta media", role="discovery"))
    # also keep region-level matrix for supplementary use
    meta24 = pd.DataFrame({"group": ["ATAA" if c.startswith("AN") else "Control" for c in counts24.columns],
                           "region": [c.rsplit("_", 1)[0] for c in counts24.columns],
                           "subject": [subj[c] for c in counts24.columns]}, index=counts24.columns)
    counts24.to_csv(os.path.join(OUT, "GSE140947_regions_counts.csv")); meta24.to_csv(os.path.join(OUT, "GSE140947_regions_meta.csv"))

def gse235161():
    df = pd.read_csv(os.path.join(GEO, "GSE235161_expr_fpkm_42sample.csv.gz"), index_col=0)
    df.index = df.index.astype(str)
    sym = [ENTREZ2SYM.get(i) for i in df.index]
    fpkm = collapse_sum(df, sym)
    # total-RNA library: rRNA / 7SL / snRNA dominate FPKM -> restrict to protein-coding genes and re-scale to TPM-like values
    fpkm = fpkm.loc[fpkm.index.isin(PROTEIN_CODING)]
    tpm = fpkm / fpkm.sum(axis=0) * 1e6
    x = np.log2(tpm + 1)
    x = quantile_normalize(x)
    code = pd.Series(x.columns, index=x.columns).str[0]
    meta = pd.DataFrame({"group": np.where(code == "N", "Control", "ATAA"),
                         "subtype_code": code.values, "subject": x.columns}, index=x.columns)
    save("GSE235161", x, meta, info=dict(disease="ATAA", platform="RNA-seq FPKM (protein-coding, TPM-rescaled, quantile-normalised)", tissue="thoracic (ascending) aorta", role="validation"))

def gse202267():
    df = pd.read_csv(os.path.join(GEO, "GSE202267_GEO_TAA_AAA_readcounts.csv.gz"), index_col=0)
    df = df[[c for c in df.columns if not c.startswith("AAA")]]
    counts = collapse_sum(df, df.index.values)
    meta = pd.DataFrame({"group": ["ATAA" if c.startswith("TAA") else "Control" for c in counts.columns], "subject": counts.columns}, index=counts.columns)
    save("GSE202267", log2cpm(counts), meta, counts=counts, info=dict(disease="ATAA", platform="RNA-seq counts", tissue="thoracic aorta", role="validation"))

# ==============================================================================================
# TAAD cohorts
# ==============================================================================================
def gse52093():
    m, e = read_series_matrix(os.path.join(GEO, "GSE52093_series_matrix.txt.gz"))
    meta = pd.DataFrame({"group": np.where(m["char:disease state"].str.contains("dissection"), "TAAD", "Control"),
                         "sex": m["char:gender"], "subject": m["Sample_title"]}, index=m.index)
    # raw (non-normalised) intensities with detection p-values -> match columns to GSMs by rank correlation
    raw = pd.read_csv(os.path.join(GEO, "GSE52093_non-normalized.txt.gz"), sep="\t", index_col=0)
    sig = raw[[c for c in raw.columns if c.startswith("SAMPLE")]]
    det = raw[[c for c in raw.columns if c.startswith("Detection")]]; det.columns = sig.columns
    common = sig.index.intersection(e.index)
    corr = pd.DataFrame({g: [sig.loc[common, s].corr(e.loc[common, g], method="spearman") for s in sig.columns] for g in e.columns}, index=sig.columns)
    mapping = corr.idxmax(axis=0)  # GSM -> SAMPLE column
    assert mapping.is_unique and (corr.max(axis=0) > 0.95).all(), corr.max(axis=0)
    sig = sig[mapping.values]; sig.columns = mapping.index; det = det[mapping.values]; det.columns = mapping.index
    expressed = (det < 0.05).sum(axis=1) >= 3
    x = np.log2(sig.clip(lower=1))
    x = quantile_normalize(x)
    x = x.loc[expressed[expressed].index]
    # the annot file has a header block starting with '^' lines; find header row
    with gzip.open(os.path.join(GEO, "GPL10558.annot.gz"), "rt", errors="replace") as fh:
        lines = fh.read().split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("ID\t"))
    end = next((i for i, l in enumerate(lines) if l.startswith("!platform_table_end")), len(lines))
    ann = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", dtype=str, quoting=3).set_index("ID")
    sym = ann.reindex(x.index)["Gene symbol"]
    expr = collapse_max_mean(x[meta.index], sym.values)
    save("GSE52093", expr, meta, info=dict(disease="TAAD", platform="Illumina HumanHT-12 V4 (log2, quantile)", tissue="ascending aorta", role="discovery"))

def gse153434():
    df = pd.read_csv(os.path.join(GEO, "GSE153434_all.counts.txt.gz"), sep="\t")
    df = df[df["type_of_gene"].isin(["protein-coding", "ncRNA", "lncRNA"])]
    cols = [c for c in df.columns if c.startswith("Normal") or c.startswith("Patient")]
    counts = collapse_sum(df[cols].set_index(df["Symbol"]), df["Symbol"].values)
    meta = pd.DataFrame({"group": ["TAAD" if c.startswith("Patient") else "Control" for c in counts.columns], "subject": counts.columns}, index=counts.columns)
    save("GSE153434", log2cpm(counts), meta, counts=counts, info=dict(disease="TAAD", platform="RNA-seq counts (Illumina HiSeq X Ten)", tissue="ascending aorta", role="discovery"))

def gse98770():
    """Re-processed from raw Agilent Feature Extraction files: gProcessedSignal -> log2 -> quantile normalisation;
    probes retained if gIsWellAboveBG in >= 3 arrays; slide/batch recorded."""
    m, _ = read_series_matrix(os.path.join(GEO, "GSE98770-GPL14550_series_matrix.txt.gz"))
    meta = pd.DataFrame({"group": np.where(m["char:tissue"].str.contains("dissected"), "TAAD", "Control"),
                         "sex": m["char:gender"], "age": m["char:age"].str.replace("y", ""), "subject": m.index}, index=m.index)
    sig = {}; wab = {}
    for f in sorted(glob.glob(os.path.join(GEO, "GSE98770", "GSM26115*_GE1_*.txt.gz"))):
        gsm = os.path.basename(f).split("_")[0]
        meta.loc[gsm, "slide"] = os.path.basename(f).split("_")[2]
        with gzip.open(f, "rt", errors="replace") as fh:
            lines = fh.read().split("\n")
        start = next(i for i, l in enumerate(lines) if l.startswith("FEATURES"))
        d = pd.read_csv(io.StringIO("\n".join(lines[start:])), sep="\t", usecols=["ControlType", "ProbeName", "gProcessedSignal", "gIsWellAboveBG"], low_memory=False)
        d = d[d.ControlType == 0].groupby("ProbeName").agg({"gProcessedSignal": "mean", "gIsWellAboveBG": "max"})
        sig[gsm] = d["gProcessedSignal"]; wab[gsm] = d["gIsWellAboveBG"]
    sig = pd.DataFrame(sig); wab = pd.DataFrame(wab)
    keep = (wab == 1).sum(axis=1) >= 3
    x = quantile_normalize(np.log2(sig.clip(lower=1)))[keep[keep].index if False else sig.columns].loc[keep[keep].index]
    lines = open(os.path.join(GEO, "GPL14550_full.txt"), errors="replace").read().split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("ID\t"))
    ann = pd.read_csv(io.StringIO("\n".join(lines[start:])), sep="\t", dtype=str, quoting=3, on_bad_lines="skip").set_index("ID")
    x = x.loc[x.index.intersection(ann.index)]
    sym = ann.loc[x.index, "GENE_SYMBOL"]
    expr = collapse_max_mean(x[meta.index], sym.values)
    save("GSE98770", expr, meta, info=dict(disease="TAAD", platform="Agilent SurePrint G3 8x60K v2 (re-processed from raw: log2, quantile)", tissue="ascending aorta intima-media", role="validation"))

def gse190635():
    m, e = read_series_matrix(os.path.join(GEO, "GSE190635_series_matrix.txt.gz"))
    meta = pd.DataFrame({"group": np.where(m["char:disease state"].str.contains("dissection"), "TAAD", "Control"), "subject": m.index}, index=m.index)
    with gzip.open(os.path.join(GEO, "GPL570.annot.gz"), "rt", errors="replace") as fh:
        lines = fh.read().split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("ID\t"))
    end = next((i for i, l in enumerate(lines) if l.startswith("!platform_table_end")), len(lines))
    ann = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", dtype=str, quoting=3).set_index("ID")
    sym = ann.reindex(e.index)["Gene symbol"].str.split("///").str[0].str.strip()
    expr = collapse_max_mean(e[meta.index], sym.values)
    save("GSE190635", expr, meta, info=dict(disease="TAAD", platform="Affymetrix HG-U133 Plus 2.0 (MAS5 log2)", tissue="aorta", role="validation"))

def gse294606():
    df = pd.read_csv(os.path.join(GEO, "GSE294606_proceeded_data.txt.gz"), encoding="latin-1",
                     usecols=lambda c: c in ("gene_name",) or c.startswith("count."))
    cnt = df[[c for c in df.columns if c.startswith("count.")]]; cnt.columns = [c.replace("count.", "") for c in cnt.columns]
    counts = collapse_sum(cnt, df["gene_name"].values)
    meta = pd.DataFrame({"group": ["TAAD" if c.startswith("A") else "Control" for c in counts.columns], "subject": counts.columns}, index=counts.columns)
    save("GSE294606", log2cpm(counts), meta, counts=counts, info=dict(disease="TAAD", platform="RNA-seq counts", tissue="aorta", role="validation"))

def gse267434():
    # transcript -> gene via Ensembl GTF
    t2g = {}
    with gzip.open(os.path.join(GEO, "Homo_sapiens.GRCh38.110.gtf.gz"), "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.split("\t")
            if f[2] != "transcript": continue
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            gn = re.search(r'gene_name "([^"]+)"', f[8])
            gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            t2g[tid] = gn.group(1) if gn else ENSG2SYM.get(gid)
    df = pd.read_csv(os.path.join(GEO, "GSE267434_transcript_sample_FPKM.txt.gz"), sep="\t", index_col=0)
    sym = [t2g.get(t.split(".")[0]) for t in df.index]
    fpkm = collapse_sum(df, sym)
    ctrl = {"S1", "S2", "S3", "S4", "S5", "S13"}
    meta = pd.DataFrame({"group": ["Control" if c in ctrl else "TAAD" for c in fpkm.columns], "subject": fpkm.columns}, index=fpkm.columns)
    save("GSE267434", np.log2(fpkm + 1), meta, info=dict(disease="TAAD", platform="RNA-seq FPKM (transcript-level summed)", tissue="ascending aorta", role="validation"))

def gse147026():
    df = pd.read_csv(os.path.join(GEO, "GSE147026_mRNA-ADVSCK.All.anno.txt.gz"), sep="\t")
    ad = ["A00710", "A00711", "A00717", "A00720"]; ck = ["GA1", "Normol3", "X262708", "X262982"]
    x = df[ad + ck]; x.index = df["ensembl_gene_id"]
    # sanity: sign of provided log2FC agrees with AD-CK difference
    d = np.log2(x[ad] + 1).mean(axis=1) - np.log2(x[ck] + 1).mean(axis=1)
    assert pd.Series(d.values).corr(df["log2FC"].reset_index(drop=True), method="spearman") > 0.8
    sym = [s if isinstance(s, str) else ENSG2SYM.get(g) for s, g in zip(df["external_gene_name"], df["ensembl_gene_id"])]
    fpkm = collapse_sum(x, sym)
    meta = pd.DataFrame({"group": ["TAAD"] * 4 + ["Control"] * 4, "subject": ad + ck}, index=ad + ck)
    save("GSE147026", np.log2(fpkm + 1), meta, info=dict(disease="TAAD", platform="RNA-seq (normalised expression)", tissue="aortic media", role="validation"))

# ==============================================================================================
# Direct ATAA-vs-TAAD cohort: GSE318877 micro-regional RNA-seq of aortic media
# ==============================================================================================
def gse318877():
    files = sorted(glob.glob(os.path.join(GEO, "GSE318877", "*.csv.gz")))
    # sample annotation from GEO (subject id, micro-region, disease)
    import requests
    cache = os.path.join(GEO, "GSE318877_gsm_meta.tsv")
    if not os.path.exists(cache):
        txt = requests.get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318877&targ=gsm&form=text&view=brief", timeout=120).text
        rows = []; cur = {}
        for l in txt.splitlines():
            if l.startswith("^SAMPLE"):
                if cur: rows.append(cur)
                cur = {"gsm": l.split(" = ")[1]}
            elif l.startswith("!Sample_title"): cur["title"] = l.split(" = ", 1)[1]
            elif l.startswith("!Sample_characteristics_ch1"):
                k, v = l.split(" = ", 1)[1].split(": ", 1); cur[k] = v
        rows.append(cur)
        pd.DataFrame(rows).to_csv(cache, sep="\t", index=False)
    gm = pd.read_csv(cache, sep="\t", dtype=str).set_index("gsm")
    cnt = {}; tpm = {}; meta = []
    for f in files:
        gsm = os.path.basename(f).split("_", 1)[0]
        d = pd.read_csv(f, usecols=["Name", "Total gene reads", "TPM"])
        d = d.groupby("Name").sum()
        cnt[gsm] = d["Total gene reads"]; tpm[gsm] = d["TPM"]
        reg = gm.loc[gsm, "group"]
        meta.append(dict(sample=gsm, subject=gm.loc[gsm, "subject id"], region=reg[0], region_rep=reg,
                         group="TAAD" if gm.loc[gsm, "disease"] == "Dissection" else "ATAA"))
    counts = pd.DataFrame(cnt); tpmdf = pd.DataFrame(tpm)
    meta = pd.DataFrame(meta).set_index("sample")
    counts = collapse_sum(counts, counts.index.values); tpmdf = collapse_sum(tpmdf, tpmdf.index.values)
    counts.to_csv(os.path.join(OUT, "GSE318877_regions_counts.csv")); meta.to_csv(os.path.join(OUT, "GSE318877_regions_meta.csv"))
    np.log2(tpmdf + 1).to_csv(os.path.join(OUT, "GSE318877_regions_log2tpm.csv"))
    # subject-level pseudo-bulk (sum of all punches)
    subj_counts = counts.T.groupby(meta["subject"]).sum().T
    subj_meta = meta.drop_duplicates("subject").set_index("subject")[["group"]]
    subj_meta["subject"] = subj_meta.index
    subj_counts = subj_counts[subj_meta.index]
    save("GSE318877", log2cpm(subj_counts), subj_meta, counts=subj_counts,
         info=dict(disease="ATAA_vs_TAAD", platform="micro-regional RNA-seq (aortic media punches, summed per subject)", tissue="ascending aorta media", role="direct comparison"))

if __name__ == "__main__":
    which = sys.argv[1:] or ["gse26155", "gse140947", "gse235161", "gse202267", "gse52093", "gse153434", "gse98770", "gse190635", "gse294606", "gse267434", "gse147026", "gse318877"]
    for w in which:
        print("==>", w, flush=True)
        globals()[w]()
