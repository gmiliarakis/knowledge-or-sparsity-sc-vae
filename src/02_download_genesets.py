# MSigDB v2024.1 gene sets (Hallmark, Reactome) and the HGNC gene table used to match their symbols.
from datetime import date
from pathlib import Path
import urllib.request

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
