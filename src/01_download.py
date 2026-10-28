# Downloads: the Kang et al. files (GEO), the gene sets and gene tables, and the cells of the Human Cytokine Dictionary that the study uses.
# The cytokine archive (a 41 GB zip) is never downloaded: the counts of the chosen cells are read from it in place by range requests
# (class RangeFile) and saved as a small file in data/raw/, which 02_preprocess.py then reads.
#   python src/01_download.py                      all three, with IFN-beta as the cytokine
#   python src/01_download.py --cytokine IL-2       the cells of another cytokine too
#   python src/01_download.py --dry-run             plan only, for the cytokine: the cells and the size of the read
# Kang et al. 2018 (GSE96583), batch 2 only: PBMCs from 8 donors, ctrl vs 6 h IFN-beta.
import argparse
import io
import re
import struct
import urllib.request
from collections import OrderedDict
from datetime import date
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import requests
from scipy.sparse import csr_matrix

from sources import CYTOKINE_URL as URL


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


def download_kang():
    geo = "https://ftp.ncbi.nlm.nih.gov/geo"
    files = [
        "samples/GSM2560nnn/GSM2560248/suppl/GSM2560248_2.1.mtx.gz",
        "samples/GSM2560nnn/GSM2560248/suppl/GSM2560248_barcodes.tsv.gz",
        "samples/GSM2560nnn/GSM2560249/suppl/GSM2560249_2.2.mtx.gz",
        "samples/GSM2560nnn/GSM2560249/suppl/GSM2560249_barcodes.tsv.gz",
        "series/GSE96nnn/GSE96583/suppl/GSE96583_batch2.genes.tsv.gz",
        "series/GSE96nnn/GSE96583/suppl/GSE96583_batch2.total.tsne.df.tsv.gz",
    ]

    # NCBI returns 403 for urllib's default user agent
    opener = urllib.request.build_opener()
    opener.addheaders = [("User-Agent", "Mozilla/5.0")]
    urllib.request.install_opener(opener)

    raw = Path("data/raw")
    raw.mkdir(parents=True, exist_ok=True)
    for f in files:
        out = raw / Path(f).name
        if not out.exists():
            print("downloading", out.name)
            urllib.request.urlretrieve(f"{geo}/{f}", out)


def download_genesets():
    # reactome.org refuses Python's default User-Agent
    urllib.request.install_opener(urllib.request.build_opener())
    urllib.request._opener.addheaders = [("User-Agent", "Mozilla/5.0")]

    msigdb = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs"
    files = {
        "h.all.v2024.1.Hs.symbols.gmt": f"{msigdb}/h.all.v2024.1.Hs.symbols.gmt",
        "c2.cp.reactome.v2024.1.Hs.symbols.gmt": f"{msigdb}/c2.cp.reactome.v2024.1.Hs.symbols.gmt",
    }

    raw = Path("data/raw")
    raw.mkdir(parents=True, exist_ok=True)

    # HGNC updates this table regularly, so keep a dated copy and reuse it once downloaded
    if not list(raw.glob("hgnc_complete_set_*.txt")):
        hgnc = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
        files[f"hgnc_complete_set_{date.today()}.txt"] = hgnc

    # Reactome's pathway names and hierarchy, used to leave out the Disease branch of the Reactome sets. They change with every
    # release, so keep dated copies and reuse them once downloaded
    reactome = "https://reactome.org/download/current"
    for name in ["ReactomePathways", "ReactomePathwaysRelation"]:
        if not list(raw.glob(f"{name}_*.txt")):
            files[f"{name}_{date.today()}.txt"] = f"{reactome}/{name}.txt"

    for name, url in files.items():
        out = raw / name
        if not out.exists():
            print("downloading", name)
            urllib.request.urlretrieve(url, out)


def fetch_cytokine(cytokine, control, read_per_group, seed, dry_run):
    """The counts of the cells of one cytokine and of the control, read in place from the atlas, saved to data/raw/cytokine_<slug>_cells.h5ad."""
    slug = re.sub("[^a-z0-9]+", "_", cytokine.lower()).strip("_")
    extract_file = Path(f"data/raw/cytokine_{slug}_cells.h5ad")
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
    runs = runs[runs["cytokine"].isin([cytokine, control])]
    runs["condition"] = np.where(runs["cytokine"] == cytokine, "stim", "ctrl")
    assert (runs["condition"] == "stim").any(), f"no cells for {cytokine}; labels: {sorted(obs['cytokine'].cat.categories)}"
    print(f"{cytokine}: {runs.loc[runs.condition == 'stim', 'n'].sum()} cells in {(runs.condition == 'stim').sum()} runs; "
          f"{control}: {runs.loc[runs.condition == 'ctrl', 'n'].sum()} cells in {(runs.condition == 'ctrl').sum()} runs")

    # whole runs per donor x condition until enough cells are read: scattered single cells would read far more of the file.
    # The random state depends on donor and condition only, so the control cells are the same for every cytokine.
    chosen = []
    for (donor, condition), g in runs.groupby(["donor", "condition"]):
        rng = np.random.default_rng([seed, int(re.sub(r"\D", "", donor)), condition == "stim"])
        order = rng.permutation(len(g))
        cum = np.cumsum(g["n"].values[order])
        keep = order[: int(np.searchsorted(cum, read_per_group)) + 1]
        if cum[-1] < read_per_group:
            print(f"  {donor} {condition}: only {cum[-1]} cells available")
        chosen.append(g.iloc[keep])
    chosen = pd.concat(chosen).sort_values("start")
    indptr = h["X/indptr"][:]
    nnz = int((indptr[chosen["end"].values] - indptr[chosen["start"].values]).sum())
    print(f"plan: {chosen['n'].sum()} cells in {len(chosen)} runs, {nnz} counts, about {nnz * 12 / 2**30:.1f} GiB read "
          f"(about {len(chosen) * 2 * 0.5:.0f} MB to {nnz * 12 / 2**20:.0f} MB transferred)")
    if dry_run:
        print(chosen.groupby(["donor", "condition"])["n"].sum().unstack().to_string())
        h.close()  # close the HDF5 file and its reader explicitly: leaving them to the interpreter at exit can crash it
        rf.close()
        return

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
    h.close()  # the counts are read; see the note at the dry run
    rf.close()
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
    adata.obs["condition"] = np.where(adata.obs["cytokine"] == cytokine, "stim", "ctrl")
    adata.obs["donor"] = sel["donor"].astype(str).values
    adata.obs["cell_type"] = sel["cell_type"].astype(str).astype("category").values
    adata.obs["sample"] = sel["sample"].astype(str).values
    adata.layers["counts"] = adata.X.copy()
    print(adata)
    print(pd.crosstab(adata.obs["cell_type"], adata.obs["condition"]))

    extract_file.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(extract_file)
    print(adata, "->", extract_file)


parser = argparse.ArgumentParser()
parser.add_argument("--cytokine", action="append", help="label in the cytokine column, e.g. IFN-beta (repeat for several; default IFN-beta)")
parser.add_argument("--control", default="PBS")
parser.add_argument("--read-per-group", type=int, default=1800, help="cells read per donor x condition, before QC")
parser.add_argument("--seed", type=int, default=0)
parser.add_argument("--dry-run", action="store_true", help="plan the read of the cytokine cells and stop")
args = parser.parse_args()

if not args.dry_run:
    download_kang()
    download_genesets()
for cytokine in args.cytokine or ["IFN-beta"]:
    fetch_cytokine(cytokine, args.control, args.read_per_group, args.seed, args.dry_run)
