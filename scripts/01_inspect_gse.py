"""Inspect GEO series: sample titles/characteristics + supplementary files."""
import requests, re, sys, time
def brief(acc, targ="self"):
    url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ={targ}&form=text&view=brief"
    r=requests.get(url,timeout=120); r.raise_for_status(); return r.text
def suppl_list(acc):
    stub=acc[:-3]+"nnn"
    out=[]
    for sub in ["suppl","matrix"]:
        url=f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{acc}/{sub}/"
        r=requests.get(url,timeout=120)
        if r.ok:
            out+= [f"{sub}/{m}" for m in re.findall(r'href="([^"]+)"',r.text) if not m.startswith("/") and not m.startswith("?")]
    return out
for acc in sys.argv[1:]:
    print("="*100); print(acc)
    s=brief(acc,"self")
    for k in ["!Series_title","!Series_summary","!Series_overall_design","!Series_platform_id","!Series_sample_id","!Series_type"]:
        vals=[l.split(" = ",1)[1] for l in s.splitlines() if l.startswith(k)]
        if k=="!Series_sample_id": print(k, len(vals)); continue
        for v in vals: print(k, v[:600])
    g=brief(acc,"gsm")
    cur=None
    for l in g.splitlines():
        if l.startswith("^SAMPLE"): cur=l.split(" = ")[1]; title=None; ch=[]; src=None
        elif l.startswith("!Sample_title"): title=l.split(" = ",1)[1]
        elif l.startswith("!Sample_source_name_ch1"): src=l.split(" = ",1)[1]
        elif l.startswith("!Sample_characteristics_ch1"): ch.append(l.split(" = ",1)[1])
        elif l.startswith("!Sample_data_processing") and cur:
            pass
        elif l.startswith("!Sample_platform_id") and cur:
            print(f"  {cur}\t{title}\t{src}\t{' | '.join(ch)[:200]}"); cur=None
    print("  FILES:", suppl_list(acc))
    time.sleep(0.5)
