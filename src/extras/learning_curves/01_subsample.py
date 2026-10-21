# Training subsamples for the learning curves, with the shuffled masks and co-expression modules rebuilt on each, as in
# src/04_build_masks.py and src/05_build_coexpression_null.py. One nested series of subsamples per seed; validation
# and test cells stay fixed. Choices are explained in the README.
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sizes = [500, 1000, 2000, 4000, 8000]
seeds = range(10)

masks = np.load("data/processed/masks_hallmark.npz")
real = masks["real"]
set_sizes = real.sum(axis=0)
adata = ad.read_h5ad("data/processed/kang.h5ad")
assert (adata.var_names == masks["genes"]).all()
counts = adata.layers["counts"].tocsr()
train = np.flatnonzero((adata.obs["split"] == "train").values)
group = adata.obs.iloc[train].groupby(["donor", "condition", "cell_type"], observed=True).ngroup().values

out = {}
for s in seeds:
    # Stratified and nested: cells are put in random order inside each donor x condition x cell type group and
    # keyed by their relative position in it; the n smallest keys give every group its share of n cells, and each
    # subsample contains the smaller ones.
    rng = np.random.default_rng(s)
    key = np.zeros(len(train))
    for g in np.unique(group):
        members = np.flatnonzero(group == g)
        key[members] = (rng.permutation(len(members)) + 0.5) / len(members)
    order = np.lexsort((rng.random(len(train)), key))

    for n in sizes:
        rows = np.sort(train[order[:n]])
        x = counts[rows]
        rng = np.random.default_rng([s, n])

        # shuffled mask: gene labels permuted inside 25 bins of mean expression in the subsample
        mean = np.asarray(x.mean(axis=0)).ravel()
        bins = pd.qcut(pd.Series(mean).rank(method="first"), 25, labels=False).values
        perm = np.arange(len(mean))
        for b in range(25):
            idx = np.flatnonzero(bins == b)
            perm[idx] = rng.permutation(idx)

        # co-expression modules: a random seed gene detected in >= 1% of the subsample and its most correlated genes
        norm = ad.AnnData(x.astype(np.float32))
        sc.pp.normalize_total(norm)
        sc.pp.log1p(norm)
        z = norm.X.toarray()
        z -= z.mean(axis=0)
        sd = z.std(axis=0)
        sd[sd == 0] = 1
        z /= sd
        detected = np.asarray((x > 0).mean(axis=0)).ravel()
        seed_genes = rng.choice(np.flatnonzero(detected >= 0.01), real.shape[1], replace=False)
        corr = z[:, seed_genes].T @ z / len(rows)  # only the rows of the seed genes are needed
        modules = np.zeros(real.shape, dtype=bool)
        for k, size in enumerate(set_sizes):
            modules[np.argsort(-corr[k])[:size], k] = True

        out.setdefault(f"rows_{n}", []).append(rows)
        out.setdefault(f"perms_{n}", []).append(perm)
        out.setdefault(f"modules_{n}", []).append(modules)
        print(f"seed {s}, {n} cells: {len(np.unique(group[np.isin(train, rows)]))} of {group.max() + 1} groups, "
              f"{(detected >= 0.01).sum()} genes eligible as module seeds, modules cover {modules.any(axis=1).sum()} genes")

# nesting check: every subsample contains the smaller ones of its seed
for small, large in zip(sizes, sizes[1:]):
    assert all(np.isin(a, b).all() for a, b in zip(out[f"rows_{small}"], out[f"rows_{large}"]))
np.savez_compressed("data/processed/subsamples_hallmark.npz", sizes=sizes, **{k: np.array(v) for k, v in out.items()})
