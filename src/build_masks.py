# Gene-set masks (genes x sets, True = gene in set) and their shuffled versions. Choices are explained in the README.
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

adata = ad.read_h5ad("data/processed/kang.h5ad")
genes = adata.var_names  # Ensembl IDs; all QC genes are the gene universe

# today's official symbol for each Ensembl ID (HGNC); genes HGNC does not know keep their 2017 symbol
hgnc_file = sorted(Path("data/raw").glob("hgnc_complete_set_*.txt"))[-1]
hgnc = pd.read_csv(hgnc_file, sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
hgnc = hgnc.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["symbol"]
symbol = pd.Series(genes.map(hgnc), index=genes).fillna(adata.var["gene_symbols"])
print("HGNC table:", hgnc_file.name)
print("genes whose symbol changed since 2017:", (symbol != adata.var["gene_symbols"]).sum())

# 25 bins of equal size by mean expression in training cells (as in scanpy's score_genes),
# so a shuffled gene is replaced by one with similar expression
train = (adata.obs["split"] == "train").values
mean = np.asarray(adata.layers["counts"][train].mean(axis=0)).ravel()
bins = pd.qcut(pd.Series(mean).rank(method="first"), 25, labels=False).values

# log-normalised training cells for the known-answer check
x = ad.AnnData(adata.layers["counts"][train])
sc.pp.normalize_total(x)
sc.pp.log1p(x)
stim = (adata.obs["condition"][train] == "stim").values


# known-answer check: mean expression of each set's genes per cell, stim minus ctrl
def stim_minus_ctrl(m):
    score = (x.X @ m) / m.sum(axis=0)
    return np.asarray(score[stim].mean(axis=0) - score[~stim].mean(axis=0)).ravel()


levels = [0, 0.25, 0.5, 0.75, 1]
seeds = range(10)

for name, gmt in [("hallmark", "h.all.v2024.1.Hs.symbols.gmt"), ("reactome", "c2.cp.reactome.v2024.1.Hs.symbols.gmt")]:
    sets = {}
    for line in open("data/raw/" + gmt):
        set_name, _, *members = line.rstrip("\n").split("\t")
        sets[set_name] = set(members)

    mask = pd.DataFrame({s: symbol.isin(m).values for s, m in sets.items()}, index=genes)
    all_members = set().union(*sets.values())
    print(f"\n{name}: {len(sets)} sets, {len(all_members)} genes, {len(all_members - set(symbol))} not found in our data")

    mask = mask.loc[:, mask.sum() >= 12]
    real = mask.values
    print("sets with >= 12 genes:", real.shape[1], "| genes in at least one set:", real.any(axis=1).sum(), "of", len(genes))

    # Shuffle by permuting gene labels inside each expression bin: every set keeps its size and its overlaps
    # with other sets, only which genes are in it changes. A gene at position i takes the memberships of
    # gene perm[i]. Levels are nested: for one seed, the genes permuted at 25% are among those at 50%.
    perms = np.zeros((len(levels), len(seeds), len(genes)), dtype=np.int32)
    for s in seeds:
        rng = np.random.default_rng(s)
        order = [rng.permutation(np.flatnonzero(bins == b)) for b in range(25)]
        for i, level in enumerate(levels):
            perm = np.arange(len(genes))
            for idx in order:
                chosen = idx[: round(level * len(idx))]
                perm[chosen] = rng.permutation(chosen)
            perms[i, s] = perm

    # Share of set members that actually changed; a permuted gene can land on a member of its own set
    for i, level in enumerate(levels):
        kept = [(real & real[perms[i, s]]).sum(axis=0) / real.sum(axis=0) for s in seeds]
        print(f"  {level:.0%} of genes permuted: {1 - np.mean(kept):.0%} of set members changed")

    # Known-answer check
    diff = pd.Series(stim_minus_ctrl(real.astype(np.float32)), index=mask.columns).sort_values(ascending=False)
    print("  sets most up in stim (real mask):")
    print(diff.head(5).round(3).to_string())
    if name == "hallmark":
        ifn = mask.columns.get_loc("HALLMARK_INTERFERON_ALPHA_RESPONSE")
        shuffled = [stim_minus_ctrl(real[perms[-1, s]].astype(np.float32))[ifn] for s in seeds]
        print(f"  IFN-alpha set, stim minus ctrl: real {diff['HALLMARK_INTERFERON_ALPHA_RESPONSE']:.3f}, "
              f"100% shuffled {np.mean(shuffled):.3f} (range {min(shuffled):.3f} to {max(shuffled):.3f})")

    np.savez_compressed(
        f"data/processed/masks_{name}.npz",
        genes=np.array(genes, dtype=str), sets=np.array(mask.columns, dtype=str), real=real, perms=perms, levels=levels, bins=bins,
    )
