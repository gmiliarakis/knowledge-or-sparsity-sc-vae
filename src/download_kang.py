# Kang et al. 2018 (GSE96583), batch 2 only: PBMCs from 8 donors, ctrl vs 6 h IFN-beta.
from pathlib import Path
import urllib.request

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
