# Evidence for using all QC genes as the gene universe (README, "Data"):
# how many IFN-responding genes and Hallmark genes survive each HVG cutoff, and which responders are missed.
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

adata = ad.read_h5ad("data/processed/kang.h5ad")
train = adata[adata.obs["split"] == "train"].copy()

# same symbol matching as 04_build_masks.py
hgnc_file = sorted(Path("data/raw").glob("hgnc_complete_set_*.txt"))[-1]
hgnc = pd.read_csv(hgnc_file, sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
hgnc = hgnc.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["symbol"]
symbol = pd.Series(train.var_names.map(hgnc), index=train.var_names).fillna(train.var["gene_symbols"])

hallmark = {}
for line in open("data/raw/h.all.v2024.1.Hs.symbols.gmt"):
    name, _, *members = line.rstrip("\n").split("\t")
    hallmark[name] = set(symbol.index[symbol.isin(members)])

# seurat_v3 rank for every gene, on training cells, per donor (the first plan)
sc.pp.highly_variable_genes(train, flavor="seurat_v3", layer="counts", n_top_genes=train.n_vars, batch_key="donor")
rank = train.var["highly_variable_rank"].rank(method="first")
detected = pd.Series(np.asarray((train.layers["counts"] > 0).mean(axis=0)).ravel(), index=train.var_names)

# responding genes: |stim - ctrl| > 0.5 in mean log-normalised expression, in at least one cell type
sc.pp.normalize_total(train)
sc.pp.log1p(train)
responders = {}
for cell_type in train.obs["cell_type"].cat.categories:
    cells = train[train.obs["cell_type"] == cell_type]
    if cells.n_obs < 100:
        continue
    stim = (cells.obs["condition"] == "stim").values
    diff = pd.Series(np.asarray(cells.X[stim].mean(axis=0) - cells.X[~stim].mean(axis=0)).ravel(), index=train.var_names)
    for gene in diff.index[diff.abs() > 0.5]:
        if abs(diff[gene]) > abs(responders.get(gene, (0, ""))[0]):
            responders[gene] = (diff[gene], cell_type)
print("genes responding to IFN-beta in at least one cell type:", len(responders))

ifn = hallmark["HALLMARK_INTERFERON_ALPHA_RESPONSE"]
in_hallmark = set().union(*hallmark.values())
rows = []
for n in [2000, 3000, 5000, 8000, train.n_vars]:
    kept = set(rank.index[rank <= n])
    sizes = [len(s & kept) for s in hallmark.values()]
    rows.append({
        "genes": n,
        "responders kept": len(kept & set(responders)),
        "IFN-alpha genes kept": len(kept & ifn),
        "Hallmark genes kept": len(kept & in_hallmark),
        "sets with >= 12 genes": sum(size >= 12 for size in sizes),
        "median set size": int(np.median(sizes)),
    })
print(pd.DataFrame(rows).to_string(index=False))

print("\nresponders outside the top 5,000 HVGs:")
for gene in sorted(responders, key=lambda g: -abs(responders[g][0])):
    if rank[gene] > 5000:
        d, cell_type = responders[gene]
        print(f"  {symbol[gene]:<10} rank {int(rank[gene]):>5}  {d:+.2f} in {cell_type:<18} detected in {detected[gene]:.0%} of cells")

print("\nmedian share of cells a gene is detected in, by HVG rank:")
for lo, hi in [(1, 500), (501, 2000), (2001, 5000), (5001, 8000), (8001, train.n_vars)]:
    print(f"  rank {lo}-{hi}: {detected[rank.between(lo, hi)].median():.1%}")
