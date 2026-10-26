# Shared setup for the scoring and check scripts (imported by 08_evaluate.py and 10_checks.py, not run on its own).
# The cells, the splits and the AUROC rule are defined here once. SC_TAG selects a derived benchmark (preprint, weaker-signal
# benchmark): SC_TAG=drop250 reads data/processed/*_drop250.* and writes results/kang_drop250/ and reports/tables/kang_drop250/.
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from gene_sets import alpha_name, data, gamma_name, main, run, tag  # noqa: F401  (re-exported to the scoring scripts)

results_dir = Path(f"results/{run}")
tables_dir = Path(f"reports/tables/{run}")

min_variance = 0.01  # a latent whose posterior means vary less across validation cells is inactive

obs = ad.read_h5ad(f"data/processed/{data}.h5ad", backed="r").obs  # the cells are the same in every benchmark
stim = (obs["condition"] == "stim").values
cell_type = obs["cell_type"].astype(str).values
donor = obs["donor"].astype(str).values
split = obs["split"].values
if data == "kang":
    types = sorted(set(cell_type) - {"Megakaryocytes"})  # no nucleus, so no interferon response to find
else:  # a cell type is scored when the test split has at least 20 cells of both conditions
    test = split == "test"
    n = pd.crosstab(cell_type[test], stim[test])
    types = sorted(n.index[(n >= 20).all(axis=1)])
rows = {s: np.flatnonzero((split == s) & np.isin(cell_type, types)) for s in ["val", "test"]}
val_all = np.flatnonzero(split == "val")
full = int((split == "train").sum())


def auroc_by_type(score, which, keep=None):
    """AUROC of stimulated vs control cells within each cell type, so a latent cannot score by separating cell types.
    `keep` (a boolean array over all cells) restricts the cells, for example to one donor."""
    r = rows[which] if keep is None else rows[which][keep[rows[which]]]
    values = {}
    for t in types:
        k = r[cell_type[r] == t]
        if 0 < stim[k].sum() < len(k):
            values[t] = roc_auc_score(stim[k], score[k])
    return pd.Series(values)


def auroc(score, which, keep=None):
    """Mean over the cell types of auroc_by_type."""
    return auroc_by_type(score, which, keep).mean()
