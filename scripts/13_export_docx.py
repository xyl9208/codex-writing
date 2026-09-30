"""Export manuscript.md to a Word document (JTM submission style) with embedded figures and Table 1."""
import os, re, pandas as pd
from docx import Document
from docx.shared import Pt, Inches, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MS = os.path.join(ROOT, "manuscript"); FIG = os.path.join(ROOT, "figures")
md = open(os.path.join(MS, "manuscript.md"), encoding="utf-8").read().split("\n")
doc = Document()
st = doc.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(12)
for s in doc.sections: s.left_margin = s.right_margin = Cm(2.5); s.top_margin = s.bottom_margin = Cm(2.5)

def add_runs(par, text):
    # minimal markdown inline: **bold**, *italic*, ^sup^
    tokens = re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*|\^[^^]+\^)", text)
    for t in tokens:
        if not t: continue
        if t.startswith("**") and t.endswith("**"): r = par.add_run(t[2:-2]); r.bold = True
        elif t.startswith("*") and t.endswith("*"): r = par.add_run(t[1:-1]); r.italic = True
        elif t.startswith("^") and t.endswith("^"): r = par.add_run(t[1:-1]); r.font.superscript = True
        else: par.add_run(t)

FIGMAP = {"Fig. 1": ["Fig1_study_design.png", "Fig1b_cohort_concordance.png"], "Fig. 2": ["Fig2_gene_level_comparison.png"], "Fig. 3": ["Fig3_hallmark_NES_heatmap.png", "Fig3b_pathway_scatter.png"],
          "Fig. 4": ["Fig4_cell_composition.png"], "Fig. 5": ["Fig5_signature_ml.png"], "Fig. 6": ["Fig6_PPI_hubs.png"], "Fig. 7": ["Fig7_scRNA_overview.png"], "Fig. 8": ["Fig8_sc_validation.png"]}
in_legends = False; first = True
for line in md:
    line = line.rstrip()
    if not line or line == "---": continue
    if line.startswith("# ") and first:
        p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER; r = p.add_run(line[2:]); r.bold = True; r.font.size = Pt(15); first = False; continue
    if line.startswith("## "):
        h = doc.add_heading(line[3:], level=1); in_legends = line[3:].startswith("Figure legends"); continue
    if line.startswith("### "): doc.add_heading(line[4:], level=2); continue
    if line.startswith("**Table 1**"):
        p = doc.add_paragraph(); add_runs(p, line)
        t1 = pd.read_csv(os.path.join(MS, "Table1_cohorts.csv"))
        cols = ["Accession", "Disease", "Design", "Platform", "Cases", "Controls", "DE_method", "DEGs_FDR05", "Role"]
        table = doc.add_table(rows=1, cols=len(cols)); table.style = "Table Grid"
        for i, c in enumerate(cols): table.rows[0].cells[i].text = c.replace("_", " ")
        for _, r in t1.iterrows():
            cells = table.add_row().cells
            for i, c in enumerate(cols): cells[i].text = str(r[c])
        for row in table.rows:
            for c in row.cells:
                for p in c.paragraphs:
                    for run in p.runs: run.font.size = Pt(8)
        continue
    if in_legends and line.startswith("**Fig."):
        key = re.match(r"\*\*(Fig\. \d)", line).group(1)
        for f in FIGMAP.get(key, []):
            doc.add_picture(os.path.join(FIG, f), width=Inches(6.3))
        p = doc.add_paragraph(); add_runs(p, line); continue
    if re.match(r"^\d+\. ", line):
        p = doc.add_paragraph(); p.paragraph_format.left_indent = Cm(0.8); p.paragraph_format.first_line_indent = Cm(-0.8); add_runs(p, line); continue
    p = doc.add_paragraph(); add_runs(p, line)
    if line.startswith("**Keywords"): p.paragraph_format.space_after = Pt(12)
out = os.path.join(MS, "manuscript_JTM.docx"); doc.save(out); print("saved", out)
