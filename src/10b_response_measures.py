# Decoded response, alternative measures and noise ceiling. Run after 10_checks.py. Reads the saved runs
# (model.npz), trains nothing, and writes to reports/tables/<data>/<collection>/:
#   diagnostics_response_alt.csv      one row per run: the decoded response of 10_checks.py (Pearson) and five other measures
#   diagnostics_response_ceiling.csv  the same measures between two random halves of the test cells, with and without the
#                                     Spearman-Brown correction to full size, mean and SD over 50 random splits
#   python src/10b_response_measures.py --learning-curves   (main benchmark) learning_curves_response.csv: the decoded response of
#                                     every run on a training subsample, and of the full-size runs
# The decoded response is the one of 10_checks.py (I): within a cell type, the mean decoded log expression of stimulated minus
# control test cells, correlated across the genes detected in at least 1% of training cells with the observed shift, and averaged
# over the scored cell types. The labels define only the observed shift. No p-values are computed here.
#   pearson        the measure of the README
#   spearman       rank correlation
#   pearson_top    Pearson on the 10% of genes with the largest absolute observed shift in each cell type
#   pearson_rest   Pearson on the other 90%
#   pearson_w      Pearson weighted by the absolute observed shift
#   pearson_250    Pearson on the 250 strongest stimulation genes of the main benchmark (the genes removed in kang_drop250)
#   pearson_no250  Pearson on all other genes (main benchmark only)
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

from gene_sets import data, main
from scoring import cell_type, results_dir, split, stim, tables_dir, tag, types

parser = argparse.ArgumentParser()
parser.add_argument("--learning-curves", action="store_true")
parser.add_argument("--splits", type=int, default=50, help="random half splits for the noise ceiling")
args = parser.parse_args()

out = tables_dir
adata = ad.read_h5ad(f"data/processed/{data}{tag}.h5ad")
counts = adata.layers["counts"].tocsr().astype(np.float32)
total = np.asarray(counts.sum(axis=1)).ravel()
train = np.flatnonzero(split == "train")
logx = counts.multiply(1 / total[:, None]).tocsr()
logx.data = np.log1p(logx.data * np.median(total[train]))  # the shifted log of the encoder
detected = np.asarray((counts[train] > 0).mean(axis=0)).ravel() >= 0.01
donor_onehot = pd.get_dummies(adata.obs["donor"]).values.astype(np.float32)  # as in 06_fit.py
median_total = float(np.median(total[train]))
test_rows = np.flatnonzero(split == "test")
groups = {(t, c): test_rows[(cell_type[test_rows] == t) & (stim[test_rows] == c)] for t in types for c in (True, False)}


def mean_log(idx):
    return np.asarray(logx[idx].mean(axis=0)).ravel()


observed_shift = {t: mean_log(groups[(t, True)]) - mean_log(groups[(t, False)]) for t in types}
# the genes of each subset are fixed from the observed shift of all test cells of the cell type
top = {}
for t in types:
    magnitude = np.where(detected, np.abs(observed_shift[t]), -np.inf)
    cutoff = np.quantile(magnitude[detected], 0.9)
    top[t] = detected & (magnitude >= cutoff)
removed_file = Path("reports/tables/kang_drop250/removed_genes_drop250.csv")
if main and removed_file.exists():
    removed = np.isin(adata.var_names, pd.read_csv(removed_file)["gene_id"])
else:
    removed = None


def weighted_pearson(a, b, w):
    mw = w.sum()
    ma, mb = (w * a).sum() / mw, (w * b).sum() / mw
    cov = (w * (a - ma) * (b - mb)).sum() / mw
    return cov / np.sqrt((w * (a - ma) ** 2).sum() / mw * (w * (b - mb) ** 2).sum() / mw)


def pearson(a, b):
    return np.corrcoef(a, b)[0, 1]


def measures(a, b, t):
    """The measures between two shifts over the genes (a decoded or half-sample shift, b the observed shift of the cell type)."""
    d = detected
    m = {"pearson": pearson(a[d], b[d]), "spearman": stats.spearmanr(a[d], b[d])[0],
         "pearson_top": pearson(a[top[t]], b[top[t]]), "pearson_rest": pearson(a[d & ~top[t]], b[d & ~top[t]]),
         "pearson_w": weighted_pearson(a[d], b[d], np.abs(observed_shift[t][d]))}
    if removed is not None:
        m["pearson_250"] = pearson(a[d & removed], b[d & removed])
        m["pearson_no250"] = pearson(a[d & ~removed], b[d & ~removed])
    return m


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


def run_measures(path):
    model = np.load(path / "model.npz")
    means = decoded_means(model["latents"], model["w"], model["v"])
    per_type = [measures(means[(t, True)] - means[(t, False)], observed_shift[t], t) for t in types]
    return {k: float(np.mean([m[k] for m in per_type])) for k in per_type[0]}


variants = [("plain", "vanilla"), ("real", "real"), ("shuffled", "shuffled"), ("setshuffle", "setshuffle"), ("coexpression", "coexpression")]

if args.learning_curves:
    rows = []
    for folder, variant in variants:
        for path in sorted((results_dir / folder).glob("seed*")):
            info = json.load(open(path / "run.json"))
            m = run_measures(path)
            rows.append({"variant": variant, "seed": info["seed"], "train_size": info.get("train_size") or int((split == "train").sum()),
                         "lfc_corr": m["pearson"], "lfc_spearman": m["spearman"]})
    table = pd.DataFrame(rows).sort_values(["train_size", "variant", "seed"], kind="stable")
    table.to_csv(out / "learning_curves_response.csv", index=False, float_format="%.6g")
    print(table.groupby(["train_size", "variant"])["lfc_corr"].agg(["mean", "std", "count"]).round(3))
    raise SystemExit

rows = []
for folder, variant in variants:
    for path in sorted((results_dir / folder).glob("seed?")):
        info = json.load(open(path / "run.json"))
        rows.append({"variant": variant, "seed": info["seed"], **run_measures(path)})
alt = pd.DataFrame(rows)
alt.to_csv(out / "diagnostics_response_alt.csv", index=False, float_format="%.6g")
print(alt.groupby("variant").mean(numeric_only=True).drop(columns="seed").round(3))

# noise ceiling: the shift of one random half of the test cells of a group against the other half
rng = np.random.default_rng(0)
ceiling = []
for rep in range(args.splits):
    per_type = []
    for t in types:
        half = {}
        for c in (True, False):
            idx = rng.permutation(groups[(t, c)])
            half[c] = (idx[: len(idx) // 2], idx[len(idx) // 2:])
        a = mean_log(half[True][0]) - mean_log(half[False][0])
        b = mean_log(half[True][1]) - mean_log(half[False][1])
        per_type.append(measures(a, b, t))
    ceiling.append({k: float(np.mean([m[k] for m in per_type])) for k in per_type[0]})
ceiling = pd.DataFrame(ceiling)
summary = []
for k in ceiling:
    r = ceiling[k]
    corrected = 2 * r / (1 + r)  # Spearman-Brown: the correlation expected between full-size shifts
    summary.append({"measure": k, "split_half": r.mean(), "split_half_sd": r.std(), "ceiling": corrected.mean(), "ceiling_sd": corrected.std()})
summary = pd.DataFrame(summary)
summary.to_csv(out / "diagnostics_response_ceiling.csv", index=False, float_format="%.6g")
print(summary.round(3).to_string(index=False))
