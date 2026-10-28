# Masks of one benchmark (SC_DATA, SC_TAG), for both collections. Three stages, each a function below, run for the Hallmark and then the Reactome sets:
#   gene_set_masks        the gene set masks (genes x sets, True = gene in set) and their shuffled versions (all-gene shuffle and sets-only shuffle)
#   expression_matched_masks  the expression-matched random mask
#   coexpression_groups   the co-expression null (groups of co-expressed genes of the sizes of the sets)
# Choices are explained in the preprint. The Reactome selection of the Kang data needs the Hallmark masks, so Hallmark runs first.
# Details of each stage:

# --- gene_set_masks
# Gene set masks and their shuffled versions: the all-gene shuffle (fully random mask) and the sets-only shuffle (within-set random mask).

# --- expression_matched_masks
# Expression-matched (response-matched) random mask. Like the within-set random mask, only genes in at least one set are permuted and only among themselves, so the genes in no set (and
# so the genes the free latents reach) stay as in the gene set mask. The permutation stays inside bins that cross
#   - 5 expression-level bins (the 25 bins of 03_build_masks.py, mean count in training cells, taken five at a time), and
#   - 5 bins of stimulation effect: the rank of |t| of stimulated against control training cells on the shifted-log expression
#     (cells of the scored cell types; as in 07_checks.py, F), cut at ranks 250, 1,000, 3,000 and 6,000.
# Every set therefore keeps its size, its overlaps, its expression level and its number of genes in each effect bin, exactly; in particular the
# number of its genes among the 250 strongest. Quintile bins did not do this (the interferon-alpha set kept 25 of its 63 strongest genes), and finer
# expression bins left too few genes to exchange. The mask uses the stimulation labels of the training cells, so it is a stricter null than the
# other random masks. Writes
# data/processed/masks_<collection>_responsematch<suffix>.npz (perms: seed x genes, with its own random seeds 2000 + seed).

# --- coexpression_groups
# Co-expression null for the Hallmark masks: for every Hallmark set, a module of the same size built only from
# gene-gene correlation in training cells (a random seed gene and its most correlated genes). It asks whether
# curated knowledge adds anything beyond groups of genes that are expressed together.

import os
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from gene_sets import data, file_tag, names_for, tag
from scoring import cell_type, split, stim, types  # noqa: F401


def gene_set_masks(collection, alpha_name, gamma_name):
    adata = ad.read_h5ad(f"data/processed/{data}{tag}.h5ad")
    genes = adata.var_names  # Ensembl IDs; all QC genes are the gene universe

    if data == "kang":
        # today's official symbol for each Ensembl ID (HGNC); genes HGNC does not know keep their 2017 symbol
        hgnc_file = sorted(Path("data/raw").glob("hgnc_complete_set_*.txt"))[-1]
        hgnc = pd.read_csv(hgnc_file, sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
        hgnc = hgnc.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["symbol"]
        symbol = pd.Series(genes.map(hgnc), index=genes).fillna(adata.var["gene_symbols"])
        print("HGNC table:", hgnc_file.name)
        print("genes whose symbol changed since 2017:", (symbol != adata.var["gene_symbols"]).sum())
    else:  # the cytokine data name genes by symbol
        symbol = pd.Series(genes, index=genes)

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


    # readable set names for printing: "HALLMARK_IL6_JAK_STAT3_SIGNALING" -> "IL6 JAK STAT3 signaling";
    # gene symbols, words with digits and a few acronyms stay upper case, other acronyms end up lower case
    acronyms = {"DNA", "RNA", "UV", "NFKB", "TNFA", "PI3K", "TGF", "JAK", "AKT", "MTORC1", "G2M", "E2F", "UPR", "ER", "NK"}
    gene_symbols = set(symbol)


    def readable(set_id):
        words = [w if w in gene_symbols or w in acronyms or any(c.isdigit() for c in w) else w.lower()
                 for w in set_id.split("_")[1:]]
        text = " ".join(words)
        return text[0].upper() + text[1:]


    levels = [0, 0.25, 0.5, 0.75, 1]
    seeds = range(10)



    def disease_sets(names):
        """The Reactome sets (MSigDB names) to leave out as disease or pathogen sets: those under the Disease branch of Reactome's
        pathway hierarchy and, for the sets that cannot be matched to a Reactome pathway by name (pathways renamed or removed
        since MSigDB 2024.1), those with a disease or pathogen word in the name."""
        pathways = pd.read_csv(sorted(Path("data/raw").glob("ReactomePathways_*.txt"))[-1], sep="\t", header=None, names=["id", "name", "species"])
        pathways = pathways[pathways["species"] == "Homo sapiens"]
        relation = pd.read_csv(sorted(Path("data/raw").glob("ReactomePathwaysRelation_*.txt"))[-1], sep="\t", header=None, names=["parent", "child"])
        children = relation.groupby("parent")["child"].apply(list).to_dict()
        todo, disease = [pathways.loc[pathways["name"] == "Disease", "id"].item()], set()
        while todo:  # every pathway below Disease
            node = todo.pop()
            if node not in disease:
                disease.add(node)
                todo += children.get(node, [])
        ids = {re.sub("[^A-Z0-9]+", "_", n.upper()).strip("_"): i for i, n in zip(pathways["id"], pathways["name"])}
        matched = {c: ids.get(c.split("_", 1)[1]) for c in names}
        unmatched = {c for c, i in matched.items() if i is None}
        return {c for c, i in matched.items() if i in disease} | {c for c in unmatched if re.search("SARS|DEFECT|INFECTION|HIV|DISEASE|CANCER|MUTANT", c)}, unmatched


    # the Reactome set that stands in for each Hallmark set that is the known answer of a stimulus
    forced = {"HALLMARK_INTERFERON_ALPHA_RESPONSE": "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE": "REACTOME_INTERFERON_GAMMA_SIGNALING",
              "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "REACTOME_TNF_SIGNALING",
              "HALLMARK_IL6_JAK_STAT3_SIGNALING": "REACTOME_INTERLEUKIN_6_FAMILY_SIGNALING",
              "HALLMARK_IL2_STAT5_SIGNALING": "REACTOME_INTERLEUKIN_2_FAMILY_SIGNALING"}


    def select_reactome(mask, raw_size, seed=0, tolerance=0.10, max_jaccard=0.34, min_coverage=0.50, verbose=True):
        """The 50 Reactome sets of the mask, chosen from names, sizes, overlaps and the Reactome hierarchy only, never from
        expression or stimulation. Each Hallmark set gets one Reactome counterpart. The five Hallmark sets that are the known
        answer of a stimulus (interferon-alpha and -gamma, TNF-alpha, IL-6, IL-2) get the Reactome sets in `forced`, so the
        selection does not depend on which stimulus is analysed. The other 45 slots take a Reactome set of the size of the
        Hallmark set (within `tolerance`), filled from the largest Hallmark set down and drawn at random (fixed seed) among
        the sets whose Jaccard index with every set chosen so far is below `max_jaccard` (0.34 is the largest overlap between
        two Hallmark sets). Sets under the Disease branch are left out, as Hallmark has no disease or pathogen sets, and so
        are sets named after the five stimuli (including the interferon regulatory factors, IRF), as Hallmark has no second set for them. A set needs at least `min_coverage` of
        its genes among the QC genes (Hallmark's 10th percentile is 0.50), which leaves out sets of genes that PBMCs do not express."""
        hallmark = np.load("data/processed/masks_hallmark.npz")
        h_size = pd.Series(hallmark["real"].sum(axis=0), index=hallmark["sets"])
        size = mask.sum()
        disease, unmatched = disease_sets(list(mask.columns))
        coverage = size / raw_size[mask.columns]
        if verbose:
            print(f"left out as disease sets: {len(disease)} (of which {len(disease & unmatched)} unmatched to a Reactome pathway, left out by name); "
                  f"coverage below {min_coverage}: {int((coverage < min_coverage).sum())} sets; coverage of the forced sets:", coverage[list(forced.values())].round(2).to_dict())
        named_after_stimulus = "INTERFERON|IFN|TNF|INTERLEUKIN_6|IL_6|INTERLEUKIN_2(_|$)|IL2|_IRF"
        pool = sorted(c for c in mask.columns if c not in forced.values() and c not in disease and coverage[c] >= min_coverage
                      and not re.search(named_after_stimulus, c))

        def jaccard(x, y):
            return (mask[x].values & mask[y].values).sum() / (mask[x].values | mask[y].values).sum()

        chosen = list(forced.values())
        assert all(jaccard(x, y) < max_jaccard for i, x in enumerate(chosen) for y in chosen[i + 1:]), "forced sets overlap too much"
        rng = np.random.default_rng(seed)
        slots = [(h, h_size[h], forced[h], size[forced[h]], 1) for h in forced]
        for h in h_size.drop(list(forced)).sort_values(ascending=False).index:
            fits = [c for c in pool if c not in chosen and abs(size[c] - h_size[h]) <= tolerance * h_size[h]
                    and all(jaccard(c, k) < max_jaccard for k in chosen)]
            assert fits, f"no Reactome set of {h_size[h]} genes (+-{tolerance:.0%}) is left for {h}"
            chosen.append(fits[rng.integers(len(fits))])
            slots.append((h, h_size[h], chosen[-1], size[chosen[-1]], len(fits)))
        if not verbose:
            return mask[sorted(chosen)]
        print("slot: Hallmark set (genes) -> Reactome set (genes), candidates left")
        for h, hs, c, cs, n in slots:
            print(f"  {h[9:]:<40} {hs:>4} -> {c[9:]:<60} {cs:>4}  ({n})")
        print(f"reactome selection: {len(chosen)} sets, sizes {size[chosen].min()}-{size[chosen].max()} (median {size[chosen].median():.1f}); "
              f"Hallmark {h_size.min()}-{h_size.max()} (median {h_size.median():.1f}); memberships {size[chosen].sum()} (Hallmark {h_size.sum()}); "
              f"largest Jaccard {max(jaccard(x, y) for i, x in enumerate(chosen) for y in chosen[i + 1:]):.2f}")
        return mask[sorted(chosen)]


    def trainable_weights(mask, n_free=8):
        """Weights the decoder can train: the memberships, plus the free latents' links to the genes in no set."""
        real = mask.values if hasattr(mask, "values") else mask
        return int(real.sum() + n_free * (~real.any(axis=1)).sum())


    def select_reactome_matched(mask, raw_size, tolerance=0.01, tries=100):
        """The first selection seed (0, 1, 2, ...) whose Reactome mask has as many trainable decoder weights as the Hallmark mask,
        within `tolerance`: sparsity is the quantity under study, so the two collections should not differ in it. Only mask
        structure is used, never expression or stimulation."""
        hallmark = trainable_weights(np.load("data/processed/masks_hallmark.npz")["real"])
        for seed in range(tries):
            try:
                chosen = select_reactome(mask, raw_size, seed=seed, verbose=False)
            except AssertionError:  # no Reactome set of some size is left for a slot
                continue
            weights = trainable_weights(chosen)
            print(f"selection seed {seed}: {weights} trainable weights (Hallmark {hallmark}, {weights / hallmark - 1:+.1%})")
            if abs(weights / hallmark - 1) <= tolerance:
                return select_reactome(mask, raw_size, seed=seed)
        raise RuntimeError(f"no selection seed below {tries} matches the Hallmark weights within {tolerance:.0%}")


    for name, gmt in [("hallmark", "h.all.v2024.1.Hs.symbols.gmt"), ("reactome", "c2.cp.reactome.v2024.1.Hs.symbols.gmt")]:
        if name != collection:
            continue
        sets = {}
        for line in open("data/raw/" + gmt):
            set_name, _, *members = line.rstrip("\n").split("\t")
            sets[set_name] = set(members)

        mask = pd.DataFrame({s: symbol.isin(m).values for s, m in sets.items()}, index=genes)
        all_members = set().union(*sets.values())
        missing = len(all_members - set(symbol))
        print(f"\n{name}: {len(sets)} sets, {len(all_members)} genes, {missing} not among our QC genes (mostly not expressed in PBMCs)")

        if tag or data != "kang":  # other benchmark or dataset: the sets chosen on the Kang data, even if they shrank below 12 genes
            mask = mask.loc[:, list(np.load(f"data/processed/masks_{name}.npz")["sets"])]
        elif name == "reactome":
            mask = select_reactome_matched(mask, pd.Series({k: len(v) for k, v in sets.items()}))
        else:
            mask = mask.loc[:, mask.sum() >= 12]
        real = mask.values
        print(f"trainable decoder weights: {int(real.sum())} memberships + 8 free latents x {int((~real.any(axis=1)).sum())} genes in no set = "
              f"{int(real.sum() + 8 * (~real.any(axis=1)).sum())}")
        print("sets with >= 12 genes:", real.shape[1], "| genes in at least one set:", real.any(axis=1).sum(), "of", len(genes))
        print("set genes whose symbol another gene also has:", symbol[real.any(axis=1)].duplicated(keep=False).sum())

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
        for set_id, d in diff.head(5).items():
            print(f"    {readable(set_id):<55} {d:.3f}")
        ifn = mask.columns.get_loc(alpha_name)
        shuffled = [stim_minus_ctrl(real[perms[-1, s]].astype(np.float32))[ifn] for s in seeds]
        print(f"  IFN-alpha set, stim minus ctrl: real {diff[alpha_name]:.3f}, "
              f"100% shuffled {np.mean(shuffled):.3f} (range {min(shuffled):.3f} to {max(shuffled):.3f})")

        np.savez_compressed(
            f"data/processed/masks_{name}{file_tag}.npz",
            genes=np.array(genes, dtype=str), sets=np.array(mask.columns, dtype=str), real=real, perms=perms, levels=levels, bins=bins,
        )

        # Sets-only shuffle (the cleaner random control): only genes in at least one set are permuted, and only among themselves
        # inside each expression bin. Set genes stay in sets and no-set genes stay in no set, so the free latents reach the same
        # genes as in the real mask. The all-gene shuffle above moves about half of the interferon-alpha genes into the pool of
        # genes in no set. Written to its own file (perms: seed x genes) with its own random seeds, so the other masks do not change.
        in_set = real.any(axis=1)
        set_perms = np.zeros((len(seeds), len(genes)), dtype=np.int32)
        for s in seeds:
            rng = np.random.default_rng(1000 + s)
            perm = np.arange(len(genes))
            for b in range(25):
                idx = np.flatnonzero((bins == b) & in_set)
                perm[idx] = rng.permutation(idx)
            assert (real[perm].any(axis=1) == in_set).all()
            set_perms[s] = perm
        kept = [(real & real[p]).sum(axis=0) / real.sum(axis=0) for p in set_perms]
        ifn = mask.columns.get_loc(alpha_name)
        print(f"\nsets-only shuffle: {1 - np.mean(kept):.0%} of set members changed; {int((~in_set).sum())} genes in no set (unchanged);",
              "interferon-alpha genes in the free pool:", [int((real[:, ifn] & ~real[p].any(axis=1)).sum()) for p in set_perms])
        np.savez_compressed(f"data/processed/masks_{name}_setshuffle{file_tag}.npz", perms=set_perms)


def expression_matched_masks(collection, alpha_name):
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


def coexpression_groups(collection, alpha_name, gamma_name):
    masks = np.load(f"data/processed/masks_{collection}{file_tag}.npz")
    real = masks["real"]
    sizes = real.sum(axis=0)

    adata = ad.read_h5ad(f"data/processed/{data}{tag}.h5ad")
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
    ifn = real[:, sets.index(alpha_name)] | real[:, sets.index(gamma_name)]
    for s in seeds:
        m = modules[s].astype(np.float32)
        score = (x.X @ m) / m.sum(axis=0)
        diff = np.asarray(score[stim].mean(axis=0) - score[~stim].mean(axis=0)).ravel()
        best = diff.argmax()
        members = modules[s][:, best]
        print(f"  seed {s}: most raised module +{diff[best]:.3f}, {members.sum()} genes, "
              f"{ifn[members].mean():.0%} in the {collection} interferon sets")

    np.savez_compressed(
        f"data/processed/coexpression_{collection}{file_tag}.npz",
        genes=masks["genes"], sets=masks["sets"], modules=modules, seed_genes=seed_genes,
    )


for collection in ["hallmark", "reactome"]:
    alpha_name, gamma_name = names_for(collection)
    print(f"== {collection}")
    gene_set_masks(collection, alpha_name, gamma_name)
    expression_matched_masks(collection, alpha_name)
    coexpression_groups(collection, alpha_name, gamma_name)
