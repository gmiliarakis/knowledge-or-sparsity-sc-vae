# Raw counts -> cell QC -> train/val/test split. Choices are explained in the README.
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

raw = "data/raw/"

# symbols repeat in this annotation, so genes are indexed by Ensembl ID
genes = pd.read_csv(raw + "GSE96583_batch2.genes.tsv.gz", sep="\t", header=None, names=["gene_ids", "gene_symbols"])
print("genes:", len(genes), "| symbols used more than once:", genes["gene_symbols"].duplicated().sum())
genes = genes.set_index("gene_ids")

samples = {}
for condition, gsm, mtx in [("ctrl", "GSM2560248", "2.1"), ("stim", "GSM2560249", "2.2")]:
    a = sc.read_mtx(raw + f"{gsm}_{mtx}.mtx.gz").T  # files are genes x cells
    a.obs_names = pd.read_csv(raw + f"{gsm}_barcodes.tsv.gz", header=None)[0].values
    a.var = genes.copy()
    samples[condition] = a

# some barcodes occur in both runs; the authors' metadata names the stim copies "-11"
ctrl, stim = samples["ctrl"], samples["stim"]
clash = stim.obs_names.isin(ctrl.obs_names)
print("barcodes found in both samples:", clash.sum())
stim.obs_names = [b[:-2] + "-11" if c else b for b, c in zip(stim.obs_names, clash)]

adata = ad.concat(samples, label="sample", merge="same")

meta = pd.read_csv(raw + "GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t", index_col=0)
assert set(meta.index) == set(adata.obs_names)
meta = meta.loc[adata.obs_names]
assert (meta["stim"].values == adata.obs["sample"].values).all()

adata.obs["condition"] = adata.obs.pop("sample")
adata.obs["donor"] = meta["ind"].astype(str).values
adata.obs["cell_type"] = meta["cell"].astype("category").values
adata.obs["demuxlet"] = meta["multiplets"].values
adata.layers["counts"] = adata.X.copy()

print(adata)
print(pd.crosstab(adata.obs["condition"], adata.obs["demuxlet"]))

# cell QC
print("cells without a cell type:", adata.obs["cell_type"].isna().sum())
adata = adata[(adata.obs["demuxlet"] == "singlet") & adata.obs["cell_type"].notna()].copy()
print("singlets with a cell type:", adata.n_obs)

# no mitochondrial filter: these genes have no counts in the deposited matrices
mt = adata.var["gene_symbols"].str.startswith("MT-").values
print("MT- genes:", mt.sum(), "| their total counts:", adata.X[:, mt].sum())

sc.pp.calculate_qc_metrics(adata, percent_top=[20], log1p=True, inplace=True)


def is_outlier(x, k):
    return (x - x.median()).abs() > k * median_abs_deviation(x)


# 5 MADs, thresholds per condition x cell type
flags = []
for _, obs in adata.obs.groupby(["condition", "cell_type"], observed=True):
    flags.append(pd.DataFrame({
        "counts": is_outlier(obs["log1p_total_counts"], 5),
        "genes": is_outlier(obs["log1p_n_genes_by_counts"], 5),
        "top20": is_outlier(obs["pct_counts_in_top_20_genes"], 5),
    }))
flags = pd.concat(flags).loc[adata.obs_names]
outlier = flags.any(axis=1)

print("QC flags per rule, cell type and condition:")
print(flags.assign(removed=outlier).groupby([adata.obs["cell_type"], adata.obs["condition"]], observed=True).sum())
metrics = ["total_counts", "n_genes_by_counts", "pct_counts_in_top_20_genes"]
print(adata.obs.groupby("cell_type", observed=True)[metrics].median().round(1))
adata = adata[~outlier.values].copy()

sc.pp.filter_genes(adata, min_cells=20)
print("after QC:", adata.n_obs, "cells x", adata.n_vars, "genes")

# split inside every donor x condition x cell type group; at least one val and one test cell each
val_frac, test_frac = 0.15, 0.15
rng = np.random.default_rng(0)
split = np.full(adata.n_obs, "train", dtype=object)
groups = adata.obs.groupby(["donor", "condition", "cell_type"], observed=True).indices
print("smallest group:", min(len(idx) for idx in groups.values()), "cells")
for idx in groups.values():
    idx = rng.permutation(idx)
    n_val = max(1, round(val_frac * len(idx)))
    n_test = max(1, round(test_frac * len(idx)))
    split[idx[:n_val]] = "val"
    split[idx[n_val : n_val + n_test]] = "test"
adata.obs["split"] = pd.Categorical(split, categories=["train", "val", "test"])
print(pd.crosstab([adata.obs["cell_type"], adata.obs["condition"]], adata.obs["split"]))

Path("data/processed").mkdir(exist_ok=True)
adata.write_h5ad("data/processed/kang.h5ad")
print(adata)
