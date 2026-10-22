# Response-matched random mask. Run after 04_build_masks.py for the same benchmark (SC_DATA, SC_SETS):
#   python src/04b_build_responsematch.py
# Like the within-set random mask, only genes in at least one set are permuted and only among themselves, so the genes in no set (and
# so the genes the free latents reach) stay as in the gene set mask. The permutation stays inside bins that cross
#   - 5 expression-level bins (the 25 bins of 04_build_masks.py, mean count in training cells, taken five at a time), and
#   - 5 bins of stimulation effect: the rank of |t| of stimulated against control training cells on the shifted-log expression
#     (cells of the scored cell types; as in 10_checks.py, F), cut at ranks 250, 1,000, 3,000 and 6,000.
# Every set therefore keeps its size, its overlaps, its expression level and its number of genes in each effect bin, exactly; in particular the
# number of its genes among the 250 strongest. Quintile bins did not do this (the interferon-alpha set kept 25 of its 63 strongest genes), and finer
# expression bins left too few genes to exchange. The mask uses the stimulation labels of the training cells, so it is a stricter null than the
# other random masks. Writes
# data/processed/masks_<collection>_responsematch<suffix>.npz (perms: seed x genes, with its own random seeds 2000 + seed).
import os

import anndata as ad
import numpy as np
import pandas as pd

from gene_sets import alpha_name, collection, data, file_tag, tag
from scoring import cell_type, split, stim, types

n_seeds = int(os.environ.get("SC_NSEEDS", 10))
n_effect = 5
masks = np.load(f"data/processed/masks_{collection}{file_tag}.npz")
real, bins, sets = masks["real"], masks["bins"], list(masks["sets"])

adata = ad.read_h5ad(f"data/processed/{data}{tag}.h5ad")
assert (masks["genes"] == adata.var_names).all()
counts = adata.layers["counts"].tocsr().astype(np.float32)
total = np.asarray(counts.sum(axis=1)).ravel()
train = np.flatnonzero(split == "train")
logx = counts.multiply(1 / total[:, None]).tocsr()
logx.data = np.log1p(logx.data * np.median(total[train]))

rows = np.flatnonzero((split == "train") & np.isin(cell_type, types))
x, s = logx[rows], stim[rows]
mean_stim = np.asarray(x[s].mean(axis=0)).ravel()
mean_ctrl = np.asarray(x[~s].mean(axis=0)).ravel()
square = lambda m: np.asarray(m.multiply(m).mean(axis=0)).ravel()
var_stim, var_ctrl = square(x[s]) - mean_stim**2, square(x[~s]) - mean_ctrl**2
t_stat = (mean_stim - mean_ctrl) / np.sqrt(var_stim / s.sum() + var_ctrl / (~s).sum() + 1e-12)
rank = np.argsort(np.argsort(-np.abs(t_stat)))  # 0 for the strongest gene
effect_bin = np.searchsorted([250, 1000, 3000, 6000], rank, side="right")  # 0: the 250 strongest ... 4: the rest

in_set = real.any(axis=1)
expression_bin = bins // 5
cell = expression_bin * n_effect + effect_bin  # one bin per expression level x effect level
perms = np.zeros((n_seeds, real.shape[0]), dtype=np.int32)
for seed in range(n_seeds):
    rng = np.random.default_rng(2000 + seed)
    perm = np.arange(real.shape[0])
    for b in np.unique(cell[in_set]):
        idx = np.flatnonzero((cell == b) & in_set)
        perm[idx] = rng.permutation(idx)
    assert (real[perm].any(axis=1) == in_set).all()
    perms[seed] = perm

# what the control keeps: for every set, the number of its genes in each expression bin and in each effect bin
for seed in range(n_seeds):
    shuffled = real[perms[seed]]
    for kind, labels, n in [("expression", expression_bin, 5), ("effect", effect_bin, n_effect)]:
        for b in range(n):
            assert (shuffled[labels == b].sum(axis=0) == real[labels == b].sum(axis=0)).all(), (kind, b)
kept = [(real & real[p]).sum(axis=0) / real.sum(axis=0) for p in perms]
target = sets.index(alpha_name)
strong = effect_bin == 0  # the 250 strongest genes
print(f"{collection}{file_tag}: {int(in_set.sum())} genes in sets, {len(np.unique(cell[in_set]))} bins; {1 - np.mean(kept):.0%} of set members changed")
setshuffle = np.load(f"data/processed/masks_{collection}_setshuffle{file_tag}.npz")["perms"]
for name, ps in [("responsematch", perms), ("within-set random", setshuffle)]:
    print(f"known-answer set {alpha_name}, {name}: {np.mean([(real[p][:, target] & strong).sum() for p in ps]):.1f} of the {int((real[:, target] & strong).sum())} "
          f"strongest genes it holds in the gene set mask; {np.mean([(real[:, target] & real[p][:, target]).sum() for p in ps]):.1f} of its {int(real[:, target].sum())} genes retained")
np.savez_compressed(f"data/processed/masks_{collection}_responsematch{file_tag}.npz", perms=perms, effect_bin=effect_bin, t_stat=t_stat)
