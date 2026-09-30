"""Search GEO DataSets (gds) for human aortic dissection / ascending aortic aneurysm expression series."""
import requests, time, json, sys, re
import xml.etree.ElementTree as ET

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def esearch(term, retmax=300):
    r = requests.get(BASE+"esearch.fcgi", params=dict(db="gds", term=term, retmax=retmax, retmode="json"), timeout=60)
    r.raise_for_status()
    return r.json()["esearchresult"]["idlist"]

def esummary(ids):
    out = []
    for i in range(0, len(ids), 100):
        chunk = ids[i:i+100]
        r = requests.get(BASE+"esummary.fcgi", params=dict(db="gds", id=",".join(chunk), retmode="json"), timeout=120)
        r.raise_for_status()
        res = r.json()["result"]
        for uid in chunk:
            d = res.get(uid)
            if d: out.append(d)
        time.sleep(0.4)
    return out

queries = {
 "dissection": '("aortic dissection"[Title] OR "aortic dissection"[Description]) AND "Homo sapiens"[Organism] AND gse[Entry Type] AND ("expression profiling by array"[DataSet Type] OR "expression profiling by high throughput sequencing"[DataSet Type])',
 "aneurysm": '(("ascending aortic aneurysm"[Description]) OR ("thoracic aortic aneurysm"[Description]) OR ("ascending aorta"[Description] AND aneurysm[Description])) AND "Homo sapiens"[Organism] AND gse[Entry Type] AND ("expression profiling by array"[DataSet Type] OR "expression profiling by high throughput sequencing"[DataSet Type])',
}
allrows = {}
for k, q in queries.items():
    ids = esearch(q)
    print(k, len(ids), file=sys.stderr)
    for d in esummary(ids):
        allrows.setdefault(d["accession"], {})
        allrows[d["accession"]].update(dict(acc=d["accession"], title=d["title"], n=d.get("n_samples"), gpl=d.get("gpl"), type=d.get("gdstype"), date=d.get("pdat"), summary=d.get("summary","")[:400]))
        allrows[d["accession"]].setdefault("hits", set()).add(k)

rows = sorted(allrows.values(), key=lambda x: (-int(x["n"] or 0)))
for r in rows:
    print(f"{r['acc']}\t{r['n']}\tGPL{r['gpl']}\t{r['date']}\t{'/'.join(sorted(r['hits']))}\t{r['type']}\n   {r['title']}\n   {r['summary'][:300]}")
json.dump([{k:(list(v) if isinstance(v,set) else v) for k,v in r.items()} for r in rows], open("data/geo_search_results.json","w"), indent=1)
