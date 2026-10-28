# Weaker-signal benchmark: remove the n genes with the largest stimulation effect from the data. The genes are ranked by
# |t| of stimulated vs control log expression in training cells (megakaryocytes excluded), so the choice comes from the data
# and never from the gene sets. Writes data/processed/kang_drop<n>.h5ad; the other scripts read it when SC_TAG=drop<n>.
# Usage: python src/extras/drop250/drop_genes.py 250
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

n = int(sys.argv[1])
adata = ad.read_h5ad("data/processed/kang.h5ad")
obs = adata.obs
cell_type = obs["cell_type"].astype(str).values
rows = np.flatnonzero((obs["split"] == "train").values & (cell_type != "Megakaryocytes"))
stim = (obs["condition"] == "stim").values[rows]

counts = adata.layers["counts"].tocsr().astype(np.float32)
total = np.asarray(counts.sum(axis=1)).ravel()
x = counts.multiply(1 / total[:, None]).tocsr()
x.data = np.log1p(x.data * np.median(total[(obs["split"] == "train").values]))  # the encoder's shifted log
xs = x[rows]
mean1, mean0 = np.asarray(xs[stim].mean(axis=0)).ravel(), np.asarray(xs[~stim].mean(axis=0)).ravel()
square = lambda m: np.asarray(m.multiply(m).mean(axis=0)).ravel()
var1, var0 = square(xs[stim]) - mean1**2, square(xs[~stim]) - mean0**2
t = (mean1 - mean0) / np.sqrt(var1 / stim.sum() + var0 / (~stim).sum() + 1e-12)
drop = np.argsort(-np.abs(t))[:n]

Path(f"reports/tables/kang_drop{n}").mkdir(parents=True, exist_ok=True)
pd.DataFrame({"gene_id": adata.var_names[drop], "symbol": adata.var["gene_symbols"].values[drop], "t": t[drop]}) \
    .to_csv(f"reports/tables/kang_drop{n}/removed_genes_drop{n}.csv", index=False, float_format="%.5g")
keep = np.setdiff1d(np.arange(adata.n_vars), drop)
adata[:, keep].copy().write_h5ad(f"data/processed/kang_drop{n}.h5ad")
print(f"removed {n} genes; {len(keep)} left; strongest: {', '.join(adata.var['gene_symbols'].values[drop[:10]])}")
