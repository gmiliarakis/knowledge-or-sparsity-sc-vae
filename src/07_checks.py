# Checks on the scored runs. Run after 06_evaluate.py (and 06_evaluate.py --learning-curves for the slopes). Rules and
# reasons are in the preprint (Robustness checks). Tables go to reports/tables/<data>/<collection>/ (reports/tables/kang/hallmark/ for the main experiment).
#   robustness_best_latent.csv, robustness_donors*.csv, robustness_slopes.csv
#     A. every variant scored on its best active latent (chosen on validation cells) and on its named latent;
#     B. per-donor scores of the named latents; C. per-seed slopes of the learning-curve differences (main benchmark only).
#   diagnostics_random_mask.csv   D. where the all-gene random mask puts the response, and where the interferon-alpha genes land
#                                    in it (main benchmark only).
#   diagnostics_reference.csv     E. scores without a VAE: the best principal component (validation labels) and the mean of
#                                    the interferon-alpha genes still present.
#   diagnostics_gene_removal.csv  F. best principal component after removing sets of genes (main benchmark only).
#   diagnostics_probe.csv         G. logistic regression on all latents (stimulation).
#   diagnostics_fit.csv           H. test negative log-likelihood per cell, the model's own objective on the test cells.
#   diagnostics_response.csv      I. decoded response: the stimulation shift in gene expression that the decoder reproduces,
#                                    rank-correlated (Spearman) across genes with the observed shift (no latent is chosen);
#                                    lfc_pearson is the same with the Pearson correlation.
#   diagnostics_celltype.csv      J. what the latents learned about cell type, which the model never saw: balanced accuracy of a
#                                    linear probe (all validation cells, and 10/50/200 per type) and two hard pairs.
# Hypothesis tests, intervals and multiplicity correction are in 09_tests.py; nothing here carries a p-value.
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

from gene_sets import collection, data, file_tag
from scoring import (alpha_name, auroc, donor, gamma_name, min_variance, main, results_dir, rows, split, stim, tables_dir, tag,
                     types, val_all, cell_type)

out = tables_dir
out.mkdir(parents=True, exist_ok=True)
main_benchmark = main
variants = [("plain", "vanilla"), ("real", "real"), ("shuffled", "shuffled"), ("setshuffle", "setshuffle"), ("coexpression", "coexpression"),
            ("responsematch", "responsematch")]  # only where its runs exist


def runs(folder):
    return sorted((results_dir / folder).glob("seed?"))


# A and B ---------------------------------------------------------------------------------------------------------
scored, per_donor = [], []
for folder, _ in variants:
    for run in runs(folder):
        info = json.load(open(run / "run.json"))
        model = np.load(run / "model.npz")
        z, w, ids = model["latents"], model["w"], info["latent_ids"]
        active = z[val_all].var(axis=0) > min_variance
        val_auc = np.array([auroc(z[:, k], "val") for k in range(z.shape[1])])
        best = [k for k in np.argsort(-np.abs(val_auc - 0.5)) if active[k]][0]
        best_sign = 1.0 if val_auc[best] > 0.5 else -1.0
        if folder in ("real", "shuffled", "setshuffle", "responsematch"):
            named = ids.index(alpha_name)
            named_sign = np.sign(w[w[:, named] != 0, named].mean())
        else:
            named, named_sign = best, best_sign  # no names: the chosen latent is the score of the main analysis
        scored.append({"variant": folder, "seed": info["seed"],
                       "named_latent": auroc(named_sign * z[:, named], "test") if active[named] else 0.5,
                       "best_latent": auroc(best_sign * z[:, best], "test"),
                       "best_is_named": best == named, "best_latent_name": ids[best],
                       "named_latent_active": bool(active[named])})
        for d in sorted(set(donor)):
            per_donor.append({"variant": folder, "seed": info["seed"], "donor": d,
                              "auroc": auroc(named_sign * z[:, named], "test", donor == d) if active[named] else 0.5})

scored = pd.DataFrame(scored)
scored.to_csv(out / "robustness_best_latent.csv", index=False, float_format="%.6g")
print(scored.groupby("variant")[["named_latent", "best_latent"]].agg(["mean", "std"]).round(3))
per_donor = pd.DataFrame(per_donor)
per_donor.to_csv(out / "robustness_donors.csv", index=False, float_format="%.6g")
by_donor = per_donor.groupby(["donor", "variant"])["auroc"].mean().unstack()
by_donor["real - shuffled"] = by_donor["real"] - by_donor["shuffled"]
by_donor.round(4).to_csv(out / "robustness_donors_mean.csv")

# C ---------------------------------------------------------------------------------------------------------------
if main_benchmark:
    # Learning-curve slopes (not reported in the paper): uncomment after running src/extras/learning_curves and 06_evaluate.py --learning-curves.
    # lc = pd.read_csv(out / "learning_curves_scores.csv")
    # wide = lc.pivot_table(index=["seed", "train_size"], columns="variant", values="ifn_alpha")
    # slopes = []  # per-seed slopes (descriptive)
    # for other in ["shuffled", "vanilla", "coexpression"]:
    #     per_seed = []
    #     for seed, g in wide.groupby(level=0):
    #         diff = (g["real"] - g[other]).droplevel(0)
    #         per_seed.append(np.polyfit(np.log2(diff.index.values), diff.values, 1)[0])
    #     slopes.append({"comparison": f"real - {other}", "seeds": len(per_seed), "mean_slope_per_doubling": np.mean(per_seed)})
    # pd.DataFrame(slopes).to_csv(out / "robustness_slopes.csv", index=False, float_format="%.6g")
    pass

# Expression for E, F, G ------------------------------------------------------------------------------------------
adata = ad.read_h5ad(f"data/processed/{data}{tag}.h5ad")
masks = np.load(f"data/processed/masks_{collection}{file_tag}.npz", allow_pickle=True)
assert (masks["genes"] == adata.var_names).all()
sets, real = list(masks["sets"]), masks["real"]
alpha_genes, gamma_genes = real[:, sets.index(alpha_name)], real[:, sets.index(gamma_name)]

# D ---------------------------------------------------------------------------------------------------------------
if main_benchmark:
    level = list(masks["levels"]).index(1.0)
    placement = []
    for run in runs("shuffled"):
        info = json.load(open(run / "run.json"))
        seed, ids = info["seed"], info["latent_ids"]
        model = np.load(run / "model.npz")
        z, w = model["latents"], model["w"]
        active = z[val_all].var(axis=0) > min_variance
        val_auc = np.array([auroc(z[:, k], "val") for k in range(z.shape[1])])
        best = [k for k in np.argsort(-np.abs(val_auc - 0.5)) if active[k]][0]
        mask = real[masks["perms"][level, seed]]
        in_set = mask.any(axis=1)
        reach = w[:, best] != 0  # genes the best latent can reach
        placement.append({
            "seed": seed, "best_latent": ids[best], "best_latent_kind": "free" if ids[best].startswith("free") else "named set",
            "best_latent_genes": int(reach.sum()), "best_latent_alpha_genes": int((reach & alpha_genes).sum()),
            "best_latent_gamma_genes": int((reach & gamma_genes).sum()),
            "alpha_genes_in_free_pool": int((alpha_genes & ~in_set).sum()), "alpha_genes_in_a_set": int((alpha_genes & in_set).sum()),
            "named_latent_alpha_genes": int((alpha_genes & mask[:, sets.index(alpha_name)].astype(bool)).sum()),
            "best_val_auroc": val_auc[best]})
    pd.DataFrame(placement).to_csv(out / "diagnostics_random_mask.csv", index=False, float_format="%.6g")

# E ---------------------------------------------------------------------------------------------------------------
counts = adata.layers["counts"].tocsr().astype(np.float32)
total = np.asarray(counts.sum(axis=1)).ravel()
train = np.flatnonzero(split == "train")
logx = counts.multiply(1 / total[:, None]).tocsr()
logx.data = np.log1p(logx.data * np.median(total[train]))  # the shifted log of the encoder


def best_component(matrix):
    """Best of 58 principal components of the training cells, chosen on validation labels, scored on test cells."""
    components = TruncatedSVD(n_components=58, random_state=0).fit(matrix[train]).transform(matrix)
    val_auc = np.array([auroc(components[:, k], "val") for k in range(58)])
    top = int(np.argmax(np.abs(val_auc - 0.5)))
    return auroc(np.sign(val_auc[top] - 0.5) * components[:, top], "test")


isg_score = np.asarray(logx[:, np.flatnonzero(alpha_genes)].mean(axis=1)).ravel()
reference = pd.DataFrame([
    {"reference": f"mean log expression of the {int(alpha_genes.sum())} interferon-alpha genes present (no model, no labels)",
     "test_auroc": auroc(isg_score, "test")},
    {"reference": "best of 58 principal components (chosen on validation labels)", "test_auroc": best_component(logx)},
])
reference.to_csv(out / "diagnostics_reference.csv", index=False, float_format="%.6g")
print(reference.to_string(index=False))

# F. Does the separation survive removing the genes that respond to the stimulation? ------------------------------------
# Principal components are fitted on training cells without (a) the genes of the Hallmark interferon-alpha and -gamma sets,
# and (b) the n genes with the largest stimulation effect in training cells (|t| of stimulated vs control log expression,
# megakaryocytes excluded). Each condition was sequenced in its own 10x run, so a remaining separation can be a run effect
# or a wider biological response; the data cannot tell them apart.
if main_benchmark:
    train_rows = np.flatnonzero((split == "train") & np.isin(cell_type, types))
    x_train, s_train = logx[train_rows], stim[train_rows]
    mean_stim = np.asarray(x_train[s_train].mean(axis=0)).ravel()
    mean_ctrl = np.asarray(x_train[~s_train].mean(axis=0)).ravel()
    sq = lambda m: np.asarray(m.multiply(m).mean(axis=0)).ravel()
    var_stim, var_ctrl = sq(x_train[s_train]) - mean_stim**2, sq(x_train[~s_train]) - mean_ctrl**2
    t_stat = (mean_stim - mean_ctrl) / np.sqrt(var_stim / s_train.sum() + var_ctrl / (~s_train).sum() + 1e-12)
    ranking = np.argsort(-np.abs(t_stat))
    annotated = np.flatnonzero(alpha_genes | gamma_genes)
    removals = [("none", np.array([], dtype=int)), (f"interferon-alpha and -gamma sets ({len(annotated)} genes)", annotated)]
    removals += [(f"top {n} stimulation-responsive genes", ranking[:n]) for n in (100, 500, 1000, 2000)]
    removals += [("annotated + top 1000", np.union1d(annotated, ranking[:1000]))]
    ablation = []
    for name, drop in removals:
        keep = np.setdiff1d(np.arange(logx.shape[1]), drop)
        ablation.append({"genes_removed": name, "n_removed": len(drop), "n_kept": len(keep),
                         "test_auroc": best_component(logx[:, keep])})
    pd.DataFrame(ablation).to_csv(out / "diagnostics_gene_removal.csv", index=False, float_format="%.6g")

# G. Linear probe on all latents ---------------------------------------------------------------------------------------
# A logistic regression on all 58 latent means (fitted on validation cells, scored within cell types on test cells) shows
# whether a variant holds the response anywhere, not only in the single latent the main score reads.
probe = []
for folder, variant in variants:
    for run in runs(folder):
        info = json.load(open(run / "run.json"))
        z = np.load(run / "model.npz")["latents"]
        scaler = StandardScaler().fit(z[rows["val"]])
        model = LogisticRegression(max_iter=2000).fit(scaler.transform(z[rows["val"]]), stim[rows["val"]])
        probe.append({"variant": variant, "seed": info["seed"],
                      "probe_auroc": auroc(model.decision_function(scaler.transform(z)), "test")})
probe = pd.DataFrame(probe)
probe.to_csv(out / "diagnostics_probe.csv", index=False, float_format="%.6g")
print(probe.groupby("variant")["probe_auroc"].agg(["mean", "std", "min", "max"]).round(4))

# H. Test fit ------------------------------------------------------------------------------------------------------
# The loss each run stored (negative log-likelihood per cell, the KL term included at weight 1), for validation and test cells.
fit = []
for folder, variant in variants:
    for run in runs(folder):
        info = json.load(open(run / "run.json"))
        fit.append({"variant": variant, "seed": info["seed"], "val_loss": info["val_loss"], "test_loss": info["test_loss"],
                    "best_epoch": info["best_epoch"]})
fit = pd.DataFrame(fit)
fit.to_csv(out / "diagnostics_fit.csv", index=False, float_format="%.6g")
print(fit.groupby("variant")[["val_loss", "test_loss"]].agg(["mean", "std"]).round(1))

# I. Decoded response ----------------------------------------------------------------------------------------------------
# For every test cell the decoder gives the expected share of each gene, from the cell's latent means and donor. Within each
# cell type, the mean decoded log expression of stimulated minus control cells is the shift the model reproduces; it is
# correlated across genes (detected in >= 1% of training cells) with the observed shift, and averaged over cell types. No
# latent is chosen and the gene sets play no part; the labels only define the observed shift. The second column uses the
# interferon-alpha genes alone.
detected = np.asarray((counts[train] > 0).mean(axis=0)).ravel() >= 0.01
donor_onehot = pd.get_dummies(adata.obs["donor"]).values.astype(np.float32)  # as in 04_fit.py
median_total = float(np.median(total[train]))
groups = {(t, c): rows["test"][(cell_type[rows["test"]] == t) & (stim[rows["test"]] == c)] for t in types for c in (True, False)}
observed = {key: np.asarray(logx[idx].mean(axis=0)).ravel() for key, idx in groups.items()}
observed_shift = {t: observed[(t, True)] - observed[(t, False)] for t in types}
ifn_present = alpha_genes & detected


def decoded_means(z, w, v):
    means = {}
    for key, idx in groups.items():
        total_log = np.zeros(w.shape[0])
        for i in range(0, len(idx), 1024):
            part = idx[i:i + 1024]
            logits = z[part] @ w.T + donor_onehot[part] @ v.T
            logits -= logits.max(axis=1, keepdims=True)
            share = np.exp(logits)
            share /= share.sum(axis=1, keepdims=True)
            total_log += np.log1p(share * median_total).sum(axis=0)  # the same scale as the encoder input
        means[key] = total_log / len(idx)
    return means


response = []
for folder, variant in variants:
    for run in runs(folder):
        info = json.load(open(run / "run.json"))
        model = np.load(run / "model.npz")
        means = decoded_means(model["latents"], model["w"], model["v"])
        r_all, r_pearson, r_ifn = [], [], []
        for t in types:
            shift = means[(t, True)] - means[(t, False)]
            r_all.append(spearmanr(shift[detected], observed_shift[t][detected])[0])
            r_pearson.append(np.corrcoef(shift[detected], observed_shift[t][detected])[0, 1])
            r_ifn.append(spearmanr(shift[ifn_present], observed_shift[t][ifn_present])[0])
        response.append({"variant": variant, "seed": info["seed"], "lfc_spearman": np.mean(r_all), "lfc_pearson": np.mean(r_pearson),
                         "lfc_spearman_ifn": np.mean(r_ifn)})
response = pd.DataFrame(response)
response.to_csv(out / "diagnostics_response.csv", index=False, float_format="%.6g")
print(response.groupby("variant")[["lfc_spearman", "lfc_pearson", "lfc_spearman_ifn"]].agg(["mean", "std"]).round(3))

# J. Cell type ------------------------------------------------------------------------------------------------------------
# Cell-type labels were never used in training (only for QC thresholds and for scoring within types). A linear probe on the
# latent means, fitted on validation cells and scored on test cells, shows what the latents learned beyond the stimulation.
# Few-shot: the probe sees only k validation cells per type (drawn with the run's seed). Hard pairs: AUROC between two similar types.
if data == "kang":
    pairs = [("CD4 T cells", "CD8 T cells", "auroc_cd4_vs_cd8"), ("CD14+ Monocytes", "FCGR3A+ Monocytes", "auroc_cd14_vs_fcgr3a")]
else:  # the Cytokine Dictionary data name the types differently; the same two pairs (T cells, monocytes)
    pairs = [("CD4 Memory", "CD8 Memory", "auroc_cd4_vs_cd8"), ("CD14 Mono", "CD16 Mono", "auroc_cd14_vs_fcgr3a")]
celltype = []
for folder, variant in variants:
    for run in runs(folder):
        info = json.load(open(run / "run.json"))
        z = np.load(run / "model.npz")["latents"]
        scaler = StandardScaler().fit(z[rows["val"]])
        zs = scaler.transform(z)
        record = {"variant": variant, "seed": info["seed"]}

        def probe_accuracy(train_rows):
            clf = LogisticRegression(max_iter=2000).fit(zs[train_rows], cell_type[train_rows])
            return balanced_accuracy_score(cell_type[rows["test"]], clf.predict(zs[rows["test"]]))

        record["balanced_acc_full"] = probe_accuracy(rows["val"])
        for k in (10, 50, 200):
            rng = np.random.default_rng(info["seed"])
            pools = [rows["val"][cell_type[rows["val"]] == t] for t in types]  # a type with fewer than k validation cells gives all of them
            drawn = np.concatenate([rng.choice(pool, min(k, len(pool)), replace=False) for pool in pools])
            record[f"balanced_acc_{k}"] = probe_accuracy(drawn)
        for a, b, column in pairs:
            tr = rows["val"][np.isin(cell_type[rows["val"]], [a, b])]
            te = rows["test"][np.isin(cell_type[rows["test"]], [a, b])]
            clf = LogisticRegression(max_iter=2000).fit(zs[tr], cell_type[tr])
            record[column] = roc_auc_score(cell_type[te] == clf.classes_[1], clf.decision_function(zs[te]))
        celltype.append(record)
celltype = pd.DataFrame(celltype)
celltype.to_csv(out / "diagnostics_celltype.csv", index=False, float_format="%.6g")
print(celltype.drop(columns="seed").groupby("variant").agg(["mean", "std"]).round(3).T)
