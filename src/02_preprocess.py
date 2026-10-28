# Raw counts -> cell QC -> train/val/test split, for the Kang et al. data or for one cytokine contrast of the Human Cytokine Dictionary
# (a cytokine against PBS). Choices are explained in the preprint.
#   python src/02_preprocess.py kang           reads data/raw/ (01_download.py), writes data/processed/kang.h5ad
#   python src/02_preprocess.py IFN-beta        reads data/raw/cytokine_ifn_beta_cells.h5ad, writes data/processed/cytokine_ifn_beta.h5ad
#   python src/02_preprocess.py --choose DEGs.csv   the second cytokine, from the table of differentially expressed genes
import argparse
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

from sources import CYTOKINE_URL as URL


def preprocess_kang():
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


def choose_cytokine(choose_path):
    candidates = ["TNF-alpha", "IL-6", "IL-2"]
    reference = "IFN-beta"  # always analysed; its count is printed for comparison
    min_abs_log_fc, max_adj_p, min_degs_per_type = 0.5, 0.05, 50

    parts = []
    columns = ["gene", "log_fc", "adj_p_value", "celltype", "cytokine", "num_donors", "well_biased"]
    for chunk in pd.read_csv(choose_path, usecols=columns, chunksize=1_000_000):
        parts.append(chunk[chunk["cytokine"].isin(candidates + [reference])])
    deg = pd.concat(parts)
    missing = set(candidates + [reference]) - set(deg["cytokine"])
    assert not missing, f"not in the table: {missing}"

    biased = deg["well_biased"].astype(str).str.lower().isin(["true", "1"])
    deg["is_deg"] = (deg["adj_p_value"] < max_adj_p) & (deg["log_fc"].abs() >= min_abs_log_fc) & ~biased
    print("rows per cytokine:", deg.groupby("cytokine").size().to_dict())
    print("rows excluded as well-biased:", int(biased.sum()), "of", len(deg))
    print("donors in the table:", deg.groupby("cytokine")["num_donors"].max().to_dict())

    per_type = deg.groupby(["cytokine", "celltype"])["is_deg"].sum().unstack(0).fillna(0).astype(int)
    print("\nDEGs per cell type:")
    print(per_type.to_string())
    total = per_type.sum()
    print("\nDEGs summed over cell types:", total.to_dict())
    chosen = total[candidates].idxmax()
    print(f"\nchosen by the rule: {chosen}")
    print(f"responder cell types of {chosen} (>= {min_degs_per_type} DEGs; the scoring rule of >= 20 test cells is applied later):")
    print(sorted(per_type.index[per_type[chosen] >= min_degs_per_type]))


def preprocess_cytokine(cytokine, control, cells_per_group, seed, all_genes, out):
    slug = re.sub("[^a-z0-9]+", "_", cytokine.lower()).strip("_")
    out_file = Path(out or f"data/processed/cytokine_{slug}.h5ad")
    adata = ad.read_h5ad(f"data/raw/cytokine_{slug}_cells.h5ad")
    # anndata saves strings as categories in natural order (Donor1, Donor2, ..., Donor10). The subsampling and the split below visit the groups in
    # the order of their keys and draw random numbers group by group, so the columns get the types they had when the counts were read:
    # strings (in text order, Donor1, Donor10, ...) and, for the cell type, categories in text order.
    for column in ["cytokine", "condition", "donor", "sample"]:
        adata.obs[column] = adata.obs[column].astype(str).astype(object).values
    adata.obs["cell_type"] = adata.obs["cell_type"].astype(str).astype("category").values
    # cell QC, the rules of the Kang data: 5 MADs, thresholds per condition x cell type. The counts include mitochondrial
    # genes, but no mitochondrial filter is applied, as for Kang et al.
    sc.pp.calculate_qc_metrics(adata, percent_top=[20], log1p=True, inplace=True)


    def is_outlier(x, k):
        return (x - x.median()).abs() > k * median_abs_deviation(x)


    flags = []
    for _, o in adata.obs.groupby(["condition", "cell_type"], observed=True):
        flags.append(pd.DataFrame({
            "counts": is_outlier(o["log1p_total_counts"], 5),
            "genes": is_outlier(o["log1p_n_genes_by_counts"], 5),
            "top20": is_outlier(o["pct_counts_in_top_20_genes"], 5),
        }))
    flags = pd.concat(flags).loc[adata.obs_names]
    outlier = flags.any(axis=1)
    print(f"QC removes {outlier.sum()} of {adata.n_obs} cells ({outlier.mean():.1%})")
    print(flags.assign(removed=outlier).groupby([adata.obs["cell_type"], adata.obs["condition"]], observed=True).sum())
    adata = adata[~outlier.values].copy()

    # the same number of cells in every donor x condition group (at most --cells-per-group)
    keep = []
    for (donor, condition), idx in adata.obs.groupby(["donor", "condition"], observed=True).indices.items():
        rng = np.random.default_rng([seed, int(re.sub(r"\D", "", donor)), condition == "stim", 1])
        keep.append(rng.permutation(idx)[: cells_per_group])
        if len(idx) < cells_per_group:
            print(f"  {donor} {condition}: {len(idx)} cells after QC, fewer than {cells_per_group}")
    adata = adata[np.sort(np.concatenate(keep))].copy()
    sc.pp.filter_genes(adata, min_cells=20)
    if not all_genes:
        # the genes of the Kang data (current HGNC symbols), so that the gene sets, the masks and the Reactome selection chosen there
        # keep their sizes and sparsity; no response information is used
        kang = ad.read_h5ad("data/processed/kang.h5ad", backed="r")
        hgnc = pd.read_csv(sorted(Path("data/raw").glob("hgnc_complete_set_*.txt"))[-1], sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
        hgnc = hgnc.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["symbol"]
        kang_symbols = set(pd.Series(kang.var_names.map(hgnc), index=kang.var_names).fillna(kang.var["gene_symbols"]))
        shared = adata.var_names.isin(kang_symbols)
        print(f"genes with >= 20 cells: {adata.n_vars}; shared with the {len(kang_symbols)} Kang genes: {shared.sum()}")
        adata = adata[:, shared].copy()
    print("after QC and subsampling:", adata.n_obs, "cells x", adata.n_vars, "genes")

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

    adata.uns["source"] = {"dataset": "Human Cytokine Dictionary, 10M PBMCs, 12 donors, 90 cytokines (Parse Biosciences; CC BY-NC 4.0)", "url": URL,
                           "contrast": f"{cytokine} vs {control}", "seed": seed}
    out_file.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_file)
    print(adata, "->", out_file)



parser = argparse.ArgumentParser()
parser.add_argument("dataset", nargs="?", help="kang, or the label of a cytokine in the cytokine column, e.g. IFN-beta")
parser.add_argument("--choose", metavar="DEGS_CSV", help="choose the second cytokine from the DEG table of the atlas (DEGs.csv), then stop")
parser.add_argument("--control", default="PBS")
parser.add_argument("--cells-per-group", type=int, default=1200, help="cells kept per donor x condition, after QC")
parser.add_argument("--seed", type=int, default=0)
parser.add_argument("--out", help="default data/processed/cytokine_<cytokine>.h5ad")
parser.add_argument("--all-genes", action="store_true", help="keep every gene with >= 20 cells, not only the genes shared with the Kang data")
args = parser.parse_args()

if args.choose:
    choose_cytokine(args.choose)
elif args.dataset in (None, "kang"):
    preprocess_kang()
else:
    preprocess_cytokine(args.dataset, args.control, args.cells_per_group, args.seed, args.all_genes, args.out)
