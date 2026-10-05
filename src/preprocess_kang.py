from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

raw = "data/raw/"

genes = pd.read_csv(raw + "GSE96583_batch2.genes.tsv.gz", sep="\t", header=None, names=["gene_ids", "gene_symbols"])
genes = genes.set_index("gene_ids")  # symbols are not unique in this hg19 annotation

samples = {}
for condition, gsm, mtx in [("ctrl", "GSM2560248", "2.1"), ("stim", "GSM2560249", "2.2")]:
    a = sc.read_mtx(raw + f"{gsm}_{mtx}.mtx.gz").T
    a.obs_names = pd.read_csv(raw + f"{gsm}_barcodes.tsv.gz", header=None)[0].values
    a.var = genes.copy()
    samples[condition] = a

# 313 barcodes occur in both samples; the authors' metadata calls the stim copies "-11" instead of "-1"
ctrl, stim = samples["ctrl"], samples["stim"]
clash = stim.obs_names.isin(ctrl.obs_names)
stim.obs_names = [b[:-2] + "-11" if c else b for b, c in zip(stim.obs_names, clash)]

adata = ad.concat(samples, label="sample", merge="same")

meta = pd.read_csv(raw + "GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t", index_col=0)
assert set(meta.index) == set(adata.obs_names)
meta = meta.loc[adata.obs_names]
assert (meta["stim"].values == adata.obs["sample"].values).all()

adata.obs["condition"] = adata.obs.pop("sample")
adata.obs["donor"] = meta["ind"].astype(str).values
adata.obs["cell_type"] = meta["cell"].astype("category").values  # 9 cells unlabelled
adata.obs["demuxlet"] = meta["multiplets"].values
adata.layers["counts"] = adata.X.copy()

print(adata)
print(pd.crosstab(adata.obs["condition"], adata.obs["demuxlet"]))

# QC, doublet removal and gene filtering follow here (step 3)

Path("data/processed").mkdir(exist_ok=True)
adata.write_h5ad("data/processed/kang.h5ad")
