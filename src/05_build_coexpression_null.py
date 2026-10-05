# Co-expression null for the Hallmark masks: for every Hallmark set, a module of the same size built only from
# gene-gene correlation in training cells (a random seed gene and its most correlated genes). It asks whether
# curated knowledge adds anything beyond groups of genes that are simply expressed together.
import anndata as ad
import numpy as np
import scanpy as sc

masks = np.load("data/processed/masks_hallmark.npz")
real = masks["real"]
sizes = real.sum(axis=0)

adata = ad.read_h5ad("data/processed/kang.h5ad")
assert (adata.var_names == masks["genes"]).all()
train = adata[adata.obs["split"] == "train"]

x = ad.AnnData(train.layers["counts"])
sc.pp.normalize_total(x)
sc.pp.log1p(x)
stim = (train.obs["condition"] == "stim").values

# Pearson correlation between all genes, training cells only (~0.6 GB for 12,034 genes)
z = x.X.toarray()
z -= z.mean(axis=0)
sd = z.std(axis=0)
sd[sd == 0] = 1  # genes never seen in training cells get correlation 0
z /= sd
corr = z.T @ z / z.shape[0]
del z

# Seeds only among genes detected in at least 1% of training cells: the nearest neighbours of a gene seen in a
# handful of cells are noise, and the module would not be co-expressed at all.
detected = np.asarray((train.layers["counts"] > 0).mean(axis=0)).ravel()
candidates = np.flatnonzero(detected >= 0.01)
print("genes eligible as seeds:", len(candidates), "of", len(detected))

# Module k has the size of Hallmark set k and keeps its column position, but carries no biological label.
seeds = range(10)
modules = np.zeros((len(seeds), *real.shape), dtype=bool)
seed_genes = np.zeros((len(seeds), real.shape[1]), dtype=np.int32)
for s in seeds:
    rng = np.random.default_rng(s)
    seed_genes[s] = rng.choice(candidates, real.shape[1], replace=False)
    for k, (seed, size) in enumerate(zip(seed_genes[s], sizes)):
        modules[s, np.argsort(-corr[seed])[:size], k] = True  # the seed itself comes first


def coherence(mask):
    # mean correlation between the genes of each set
    out = []
    for k in range(mask.shape[1]):
        idx = np.flatnonzero(mask[:, k])
        c = corr[np.ix_(idx, idx)]
        out.append((c.sum() - len(idx)) / (len(idx) * (len(idx) - 1)))
    return np.array(out)


# A fair null should be at least as co-expressed as the real sets (Hallmark sets were built to be coherent in bulk data).
print(f"mean within-set correlation, median over sets: real {np.median(coherence(real)):.3f}, "
      f"co-expression modules {np.median([np.median(coherence(modules[s])) for s in seeds]):.3f}")

# Known-answer check: the module most raised by IFN-beta in each seed, and the share of its genes that belong to the
# Hallmark interferon sets (share, because modules and sets differ in size)
sets = list(masks["sets"])
ifn = real[:, sets.index("HALLMARK_INTERFERON_ALPHA_RESPONSE")] | real[:, sets.index("HALLMARK_INTERFERON_GAMMA_RESPONSE")]
for s in seeds:
    m = modules[s].astype(np.float32)
    score = (x.X @ m) / m.sum(axis=0)
    diff = np.asarray(score[stim].mean(axis=0) - score[~stim].mean(axis=0)).ravel()
    best = diff.argmax()
    members = modules[s][:, best]
    print(f"  seed {s}: most raised module +{diff[best]:.3f}, {members.sum()} genes, "
          f"{ifn[members].mean():.0%} in Hallmark IFN-alpha or IFN-gamma")

np.savez_compressed(
    "data/processed/coexpression_hallmark.npz",
    genes=masks["genes"], sets=masks["sets"], modules=modules, seed_genes=seed_genes,
)
