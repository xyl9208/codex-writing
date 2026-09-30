"""STRING protein-protein interaction networks and hub genes for TAAD-specific and ATAA-specific consensus DEGs."""
import os, json, time, math, requests, numpy as np, pandas as pd, networkx as nx
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
META = os.path.join(ROOT, "results", "meta"); OUT = os.path.join(ROOT, "results", "ppi"); os.makedirs(OUT, exist_ok=True); FIG = os.path.join(ROOT, "figures")
cmp = pd.read_csv(os.path.join(META, "ATAA_vs_TAAD_gene_comparison.csv"), index_col=0)
API = "https://string-db.org/api"

def string_network(genes, score=700):
    ids = []
    for i in range(0, len(genes), 500):
        r = requests.post(f"{API}/tsv/get_string_ids", data={"identifiers": "\r".join(genes[i:i+500]), "species": 9606, "limit": 1, "echo_query": 1, "caller_identity": "aorta_meta"}, timeout=180)
        r.raise_for_status(); ids.append(pd.read_csv(pd.io.common.StringIO(r.text), sep="\t")); time.sleep(1)
    ids = pd.concat(ids).drop_duplicates("queryItem")
    mapping = dict(zip(ids.stringId, ids.queryItem))
    edges = []
    sids = ids.stringId.tolist()
    r = requests.post(f"{API}/tsv/network", data={"identifiers": "\r".join(sids), "species": 9606, "required_score": score, "caller_identity": "aorta_meta"}, timeout=300)
    r.raise_for_status(); net = pd.read_csv(pd.io.common.StringIO(r.text), sep="\t")
    net["a"] = net.stringId_A.map(mapping); net["b"] = net.stringId_B.map(mapping)
    return ids, net[["a", "b", "score"]].drop_duplicates()

def analyse(cl, top_n=1200):
    s = cmp[cmp["class"] == cl].copy()
    s["absz"] = s.z_TAAD.abs() if "TAAD" in cl else s.z_ATAA.abs()
    genes = s.sort_values("absz", ascending=False).head(top_n).index.tolist()
    ids, net = string_network(genes)
    G = nx.Graph(); G.add_weighted_edges_from(net[["a", "b", "score"]].values.tolist())
    deg = pd.Series(dict(G.degree())).sort_values(ascending=False)
    # Maximal Clique Centrality (cytoHubba MCC)
    mcc = {}
    for n in G.nodes:
        mcc[n] = sum(math.factorial(len(c) - 1) if len(c) > 2 else 1 for c in nx.find_cliques(G) if n in c) if G.degree(n) < 200 else np.nan
    hubs = pd.DataFrame({"degree": deg, "MCC": pd.Series(mcc), "betweenness": pd.Series(nx.betweenness_centrality(G, k=min(300, G.number_of_nodes()), seed=1))})
    hubs = hubs.join(s[["lfc_ATAA", "lfc_TAAD", "z_ATAA", "z_TAAD", "z_diff", "fdr_diff"]]).sort_values("degree", ascending=False)
    hubs.to_csv(os.path.join(OUT, f"hubs_{cl}.csv")); net.to_csv(os.path.join(OUT, f"edges_{cl}.csv"), index=False)
    print(f"[{cl}] genes submitted={len(genes)} mapped={len(ids)} nodes={G.number_of_nodes()} edges={G.number_of_edges()}")
    print("  top hubs (degree):", ", ".join(f"{g}({int(d)})" for g, d in deg.head(25).items()))
    return G, hubs, s

fig, axes = plt.subplots(1, 2, figsize=(16, 8))
for ax, cl in zip(axes, ["TAAD_specific", "ATAA_specific"]):
    G, hubs, s = analyse(cl)
    top = hubs.index[:40]; H = G.subgraph(top)
    pos = nx.spring_layout(H, seed=3, k=0.9)
    ref = s.lfc_TAAD if "TAAD" in cl else s.lfc_ATAA
    colors = ["#C44E52" if ref.get(n, 0) > 0 else "#4C72B0" for n in H.nodes]
    sizes = [60 + 12 * hubs.loc[n, "degree"] for n in H.nodes]
    nx.draw_networkx_edges(H, pos, ax=ax, alpha=0.25, width=0.6); nx.draw_networkx_nodes(H, pos, ax=ax, node_color=colors, node_size=sizes, edgecolors="k", linewidths=0.4)
    nx.draw_networkx_labels(H, pos, ax=ax, font_size=7)
    ax.set_title(f"{cl.replace('_', '-')} DEGs: top-40 STRING hubs (red up / blue down in disease)", fontsize=10); ax.axis("off")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig6_PPI_hubs.png"), dpi=200); fig.savefig(os.path.join(FIG, "Fig6_PPI_hubs.pdf")); plt.close(fig)
