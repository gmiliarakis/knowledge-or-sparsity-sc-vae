# Oesinghaus data (10M PBMCs, the Human Cytokine Dictionary) -> one stimulation contrast (a cytokine against PBS) -> cell QC -> train/val/test split. The counts are
# read from the published zip in place (class RangeFile below), so the 41 GB file is never downloaded. Choices are explained in the README.
#   python src/03c_preprocess_oesinghaus.py IFN-beta --dry-run     plan only: the cells and the size of the read
#   python src/03c_preprocess_oesinghaus.py --choose DEGs.csv      the second cytokine, from the table of differentially expressed genes
#   python src/03c_preprocess_oesinghaus.py IFN-beta               writes data/processed/oesinghaus_ifn_beta.h5ad
import argparse
import io
import re
import struct
from collections import OrderedDict
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.stats import median_abs_deviation


URL = "https://parse-wget.s3.us-west-2.amazonaws.com/10m/10M_PBMC_12donor_90cytokines_h5ad_20260203.zip"
CHUNK = 143667  # elements per stored chunk of X/data, for the size estimate


# Read an HDF5 file (an h5ad) that is stored, uncompressed, inside a zip on a web server, without downloading the zip.
# The zip entry sits at a fixed byte offset, so range requests on the zip are range reads on the HDF5 file. h5py opens the
# object like a local file; blocks are cached because HDF5 reads many small pieces.
class RangeFile(io.RawIOBase):
    def __init__(self, url, block=16 * 2**20, cache_blocks=32):
        self.url, self.block, self.cache_blocks = url, block, cache_blocks
        self.session = requests.Session()
        head = self._fetch(0, 4095)
        signature, _, flags, method, *_ = struct.unpack("<IHHHHHIIIHH", head[:30])
        name_length, extra_length = struct.unpack("<HH", head[26:30])
        assert signature == 0x04034B50, "not a zip file"
        assert method == 0, "the zip entry is compressed, so it cannot be read in place"
        extra = head[30 + name_length:30 + name_length + extra_length]
        self.size = struct.unpack("<I", head[22:26])[0]  # uncompressed size, or 0xFFFFFFFF when zip64 holds it
        while extra:  # zip64 extra field: uncompressed size first
            field, length = struct.unpack("<HH", extra[:4])
            if field == 1 and self.size == 0xFFFFFFFF:
                self.size = struct.unpack("<Q", extra[4:12])[0]
            extra = extra[4 + length:]
        self.name = head[30:30 + name_length].decode()
        self.offset = 30 + name_length + extra_length
        self.position, self.cache = 0, OrderedDict()

    def _fetch(self, start, end):
        for attempt in range(8):
            try:
                r = self.session.get(self.url, headers={"Range": f"bytes={start}-{end}"}, timeout=60)
                r.raise_for_status()
                return r.content
            except requests.RequestException:
                if attempt == 7:
                    raise
        return b""

    def _block(self, i):
        if i in self.cache:
            self.cache.move_to_end(i)
            return self.cache[i]
        start = i * self.block
        data = self._fetch(self.offset + start, self.offset + min(start + self.block, self.size) - 1)
        self.cache[i] = data
        if len(self.cache) > self.cache_blocks:
            self.cache.popitem(last=False)
        return data

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=0):
        self.position = offset if whence == 0 else self.position + offset if whence == 1 else self.size + offset
        return self.position

    def readinto(self, buffer):
        want = min(len(buffer), self.size - self.position)
        done = 0
        while done < want:
            i, within = divmod(self.position, self.block)
            piece = self._block(i)[within:within + want - done]
            buffer[done:done + len(piece)] = piece
            done += len(piece)
            self.position += len(piece)
        return done

    def read(self, n=-1):
        n = self.size - self.position if n < 0 else n
        buffer = bytearray(min(n, self.size - self.position))
        self.readinto(buffer)
        return bytes(buffer)

parser = argparse.ArgumentParser()
parser.add_argument("cytokine", nargs="?", help="label in the cytokine column, e.g. IFN-beta")
parser.add_argument("--choose", metavar="DEGS_CSV", help="choose the second cytokine from the DEG table of the atlas (DEGs.csv), then stop")
parser.add_argument("--control", default="PBS")
parser.add_argument("--read-per-group", type=int, default=1800, help="cells read per donor x condition, before QC")
parser.add_argument("--cells-per-group", type=int, default=1200, help="cells kept per donor x condition, after QC")
parser.add_argument("--seed", type=int, default=0)
parser.add_argument("--out", help="default data/processed/oesinghaus_<cytokine>.h5ad")
parser.add_argument("--all-genes", action="store_true", help="keep every gene with >= 20 cells, not only the genes shared with the Kang data")
parser.add_argument("--dry-run", action="store_true")
args = parser.parse_args()

# The cytokine besides IFN-beta, by the rule in the README (Data): the candidate with the most DEGs summed over cell types, and its
# responder cell types. Reads only the DEG table of the atlas, never model results.
if args.choose:
    candidates = ["TNF-alpha", "IL-6", "IL-2"]
    reference = "IFN-beta"  # always analysed; its count is printed for comparison
    min_abs_log_fc, max_adj_p, min_degs_per_type = 0.5, 0.05, 50

    parts = []
    columns = ["gene", "log_fc", "adj_p_value", "celltype", "cytokine", "num_donors", "well_biased"]
    for chunk in pd.read_csv(args.choose, usecols=columns, chunksize=1_000_000):
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
    raise SystemExit
assert args.cytokine, "give a cytokine label, or --choose DEGS_CSV"
slug = re.sub("[^a-z0-9]+", "_", args.cytokine.lower()).strip("_")
out_file = Path(args.out or f"data/processed/oesinghaus_{slug}.h5ad")

# cell annotations, read without the counts
rf = RangeFile(URL, block=2**21, cache_blocks=256)
h = h5py.File(rf, "r", rdcc_nbytes=256 * 2**20, rdcc_nslots=100003)


def column(name):
    g = h["obs"][name]
    return pd.Categorical.from_codes(g["codes"][:], g["categories"][:].astype(str))


obs = pd.DataFrame({k: column(k) for k in ["cytokine", "donor", "cell_type", "sample"]})
codes = obs["sample"].cat.codes.values
starts = np.r_[0, np.flatnonzero(np.diff(codes)) + 1]  # rows come in runs of consecutive cells from one sample
ends = np.r_[starts[1:], len(codes)]
runs = pd.DataFrame({"start": starts, "end": ends, "donor": obs["donor"].values[starts].astype(str),
                     "cytokine": obs["cytokine"].values[starts].astype(str)})
runs["n"] = runs["end"] - runs["start"]
runs = runs[runs["cytokine"].isin([args.cytokine, args.control])]
runs["condition"] = np.where(runs["cytokine"] == args.cytokine, "stim", "ctrl")
assert (runs["condition"] == "stim").any(), f"no cells for {args.cytokine}; labels: {sorted(obs['cytokine'].cat.categories)}"
print(f"{args.cytokine}: {runs.loc[runs.condition == 'stim', 'n'].sum()} cells in {(runs.condition == 'stim').sum()} runs; "
      f"{args.control}: {runs.loc[runs.condition == 'ctrl', 'n'].sum()} cells in {(runs.condition == 'ctrl').sum()} runs")

# whole runs per donor x condition until enough cells are read: scattered single cells would read far more of the file.
# The random state depends on donor and condition only, so the control cells are the same for every cytokine.
chosen = []
for (donor, condition), g in runs.groupby(["donor", "condition"]):
    rng = np.random.default_rng([args.seed, int(re.sub(r"\D", "", donor)), condition == "stim"])
    order = rng.permutation(len(g))
    cum = np.cumsum(g["n"].values[order])
    keep = order[: int(np.searchsorted(cum, args.read_per_group)) + 1]
    if cum[-1] < args.read_per_group:
        print(f"  {donor} {condition}: only {cum[-1]} cells available")
    chosen.append(g.iloc[keep])
chosen = pd.concat(chosen).sort_values("start")
indptr = h["X/indptr"][:]
nnz = int((indptr[chosen["end"].values] - indptr[chosen["start"].values]).sum())
print(f"plan: {chosen['n'].sum()} cells in {len(chosen)} runs, {nnz} counts, about {nnz * 12 / 2**30:.1f} GiB read "
      f"(about {len(chosen) * 2 * 0.5:.0f} MB to {nnz * 12 / 2**20:.0f} MB transferred)")
if args.dry_run:
    print(chosen.groupby(["donor", "condition"])["n"].sum().unstack().to_string())
    raise SystemExit

# counts of the chosen runs
data, indices, row_nnz, rows = [], [], [], []
for i, (s, e) in enumerate(zip(chosen["start"], chosen["end"])):
    a, b = int(indptr[s]), int(indptr[e])
    data.append(h["X/data"][a:b])
    indices.append(h["X/indices"][a:b].astype(np.int32))
    row_nnz.append(np.diff(indptr[s:e + 1]))
    rows.append(np.arange(s, e))
    if i % 100 == 0:
        print(f"  read {i}/{len(chosen)} runs", flush=True)
rows = np.concatenate(rows)
genes = h["var/_index"][:].astype(str)
counts = csr_matrix((np.concatenate(data), np.concatenate(indices), np.r_[0, np.cumsum(np.concatenate(row_nnz))]),
                    shape=(len(rows), len(genes)))
assert (counts.data == np.round(counts.data)).all(), "the counts are not integers"

adata = ad.AnnData(counts)
adata.obs_names = [f"row{r}" for r in rows]  # row in the atlas file
adata.var_names = genes
print("duplicated gene symbols:", adata.var_names.duplicated().sum())
adata.var["gene_symbols"] = genes
adata.var_names_make_unique()
sel = obs.iloc[rows]
adata.obs["cytokine"] = sel["cytokine"].astype(str).values
adata.obs["condition"] = np.where(adata.obs["cytokine"] == args.cytokine, "stim", "ctrl")
adata.obs["donor"] = sel["donor"].astype(str).values
adata.obs["cell_type"] = sel["cell_type"].astype(str).astype("category").values
adata.obs["sample"] = sel["sample"].astype(str).values
adata.layers["counts"] = adata.X.copy()
print(adata)
print(pd.crosstab(adata.obs["cell_type"], adata.obs["condition"]))

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
    rng = np.random.default_rng([args.seed, int(re.sub(r"\D", "", donor)), condition == "stim", 1])
    keep.append(rng.permutation(idx)[: args.cells_per_group])
    if len(idx) < args.cells_per_group:
        print(f"  {donor} {condition}: {len(idx)} cells after QC, fewer than {args.cells_per_group}")
adata = adata[np.sort(np.concatenate(keep))].copy()
sc.pp.filter_genes(adata, min_cells=20)
if not args.all_genes:
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
                       "contrast": f"{args.cytokine} vs {args.control}", "seed": args.seed}
out_file.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(out_file)
print(adata, "->", out_file)
