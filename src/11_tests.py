# Every hypothesis test of the project, with multiplicity correction. Run after 08_evaluate.py (both modes) and
# 10_checks.py. Writes reports/tables/tests.csv; the README quotes this table and no other p-value.
#
# All tests are paired by seed (10 seeds). The p-value is the exact two-sided Wilcoxon signed-rank test; zero differences
# are dropped. Its smallest possible value with 10 pairs is 2/2^10 = 0.00195, so a Holm or Bonferroni correction over more
# than 25 tests could never reject. That fixes how the tests are grouped:
#
# Family Q1 (full-size experiment, confirmatory, 10 tests): Holm step-down, which controls the family-wise error rate
#   at 0.05 under any dependence between the tests (they share seeds). Intervals are the distribution-free intervals that
#   belong to the signed-rank test (Hodges-Lehmann estimate, exact interval from the Walsh averages) at the Bonferroni level
#   1 - 0.05/m, which are simultaneous for the m tests. Holm is less conservative than Bonferroni, so an interval can touch 0
#   for a test that Holm rejects.
# Family Q1b (learning curves, exploratory, 18 tests): Benjamini-Hochberg, which controls the false discovery rate at 0.05
#   and is valid for positively dependent tests (the same seeds and nested subsamples run through all sizes). Intervals are
#   Benjamini-Yekutieli false-coverage-rate intervals, level 1 - R*q/m with R the number of rejections, of the same kind.
# Family Q1c (latent-free panel, confirmatory, 12 tests): Holm, as Q1. Four measures that choose no latent (test loss, decoded
#   response, cell-type probe, stimulation probe), each with three contrasts: Hallmark minus each of the two random masks and
#   minus the co-expression mask. The no-mask VAE is a descriptive reference, and the few-shot and hard-pair cell-type scores
#   are descriptive.
# Per-cell-type and per-donor scores are descriptive and get no p-values.
#
# Equivalence of the Hallmark and co-expression masks is tested with two one-sided Wilcoxon tests against +-margin;
# the p-value is the larger of the two (intersection-union), so it needs no further correction within the test.
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# Equivalence margin on the AUROC scale: +-0.005, half a percentage point in the probability that a random stimulated cell
# ranks above a random control cell. It was set in commit 8599a48 (2026-10-06), before the review and before the sets-only
# mask. It is one fifth of the seed-to-seed SD of the unmasked VAE (0.025) and one sixth of the smallest effect reported
# as substantial (Hallmark minus unmasked, +0.031). The best variants score 0.997-0.998, so the headroom to 1 is 0.003:
# within +-0.005 means both are at the ceiling, and the test cannot reveal a difference there.
margin = 0.005
alpha = 0.05
suffix = "_" + os.environ["SC_TAG"] if os.environ.get("SC_TAG") else ""  # SC_TAG=drop250: weaker-signal benchmark
out = Path(f"reports/tables{suffix}")  # with SC_TAG only the latent-free family is tested, on that benchmark

if not suffix:
    scores = pd.read_csv(out / "scores.csv")
    best = pd.read_csv(out / "robustness_best_latent.csv").replace({"variant": {"plain": "vanilla"}})  # that table says plain
    curves = pd.read_csv(out / "learning_curves_scores.csv")
    slopes_input = curves.pivot_table(index=["seed", "train_size"], columns="variant", values="ifn_alpha")
    full = int(scores["train_size"].max())


def wilcoxon(d, **kwargs):
    d = np.asarray(d)
    d = d[d != 0]
    return stats.wilcoxon(d, method="exact", **kwargs).pvalue if len(d) else 1.0


def signed_rank_cdf(n):
    # exact null distribution of the signed-rank statistic T (sum of the ranks of the positive differences)
    counts = np.zeros(n * (n + 1) // 2 + 1)
    counts[0] = 1
    for rank in range(1, n + 1):
        counts[rank:] = counts[rank:] + counts[:-rank].copy()
    return np.cumsum(counts) / 2**n


def interval(d, level):
    """Hodges-Lehmann estimate and exact distribution-free interval of at least the given level."""
    d = np.asarray(d, dtype=float)
    n = len(d)
    walsh = np.sort([(d[i] + d[j]) / 2 for i in range(n) for j in range(i, n)])
    cdf = signed_rank_cdf(n)
    k = int(np.searchsorted(cdf, (1 - level) / 2, side="right")) - 1  # largest k with P(T <= k) <= alpha/2
    estimate = float(np.median(walsh))
    if k < 0:
        return estimate, walsh[0], walsh[-1]
    return estimate, walsh[k], walsh[len(walsh) - 1 - k]


def holm(p):
    order = np.argsort(p)
    adjusted = np.empty(len(p))
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(p) - rank) * p[i]))
        adjusted[i] = running
    return adjusted


def benjamini_hochberg(p):
    order = np.argsort(p)[::-1]
    adjusted = np.empty(len(p))
    running = 1.0
    for rank, i in enumerate(order):
        running = min(running, p[i] * len(p) / (len(p) - rank))
        adjusted[i] = running
    return adjusted


def diff(frame, column, a, b):
    wide = frame.pivot(index="seed", columns="variant", values=column)
    return (wide[a] - wide[b]).dropna().values


tests = []  # (family, test, differences, kind)
if not suffix:
    at_full = scores[scores["train_size"] == full]
    for other in ["shuffled", "setshuffle", "vanilla", "coexpression"]:
        tests.append(("Q1", f"fixed latent: real - {other}", diff(at_full, "ifn_alpha", "real", other), "difference"))
    for a, b in [("real", "shuffled"), ("real", "setshuffle"), ("real", "vanilla"), ("shuffled", "vanilla"),
                 ("setshuffle", "vanilla")]:
        tests.append(("Q1", f"best latent: {a} - {b}", diff(best, "best_latent", a, b), "difference"))
    tests.append(("Q1", f"equivalence (+-{margin}): real - coexpression", diff(at_full, "ifn_alpha", "real", "coexpression"),
                  "equivalence"))
    for size in sorted(curves["train_size"].unique()):
        if size == full:
            continue
        for other in ["shuffled", "vanilla", "coexpression"]:
            tests.append(("Q1b", f"{size} cells: real - {other}", diff(curves[curves["train_size"] == size], "ifn_alpha", "real", other), "difference"))
    for other in ["shuffled", "vanilla", "coexpression"]:
        per_seed = []
        for seed, g in slopes_input.groupby(level=0):
            d = (g["real"] - g[other]).droplevel(0)
            per_seed.append(np.polyfit(np.log2(d.index.values), d.values, 1)[0])
        tests.append(("Q1b", f"slope per doubling of cells: real - {other}", np.array(per_seed), "difference"))

latent_free = [("test loss per cell (lower is better)", "diagnostics_fit.csv", "test_loss"),
               ("decoded response correlation", "diagnostics_response.csv", "lfc_corr"),
               ("cell-type balanced accuracy", "diagnostics_celltype.csv", "balanced_acc_full"),
               ("stimulation probe AUROC", "diagnostics_probe.csv", "probe_auroc")]
for label, file, column in latent_free:
    frame = pd.read_csv(out / file)
    for other in ["shuffled", "setshuffle", "coexpression"]:
        tests.append(("Q1c", f"{label}: real - {other}", diff(frame, column, "real", other), "difference"))

rows = []
for family, name, d, kind in tests:
    if kind == "equivalence":
        p = max(wilcoxon(d + margin, alternative="greater"), wilcoxon(d - margin, alternative="less"))
    else:
        p = wilcoxon(d)
    rows.append({"family": family, "test": name, "pairs": len(d), "mean_difference": d.mean(), "p_raw": p, "_d": d, "_kind": kind})
result = pd.DataFrame(rows)

for family, method in [("Q1", "Holm"), ("Q1b", "Benjamini-Hochberg"), ("Q1c", "Holm")]:
    sel = result["family"] == family
    if not sel.any():
        continue
    p = result.loc[sel, "p_raw"].values
    m = len(p)
    adj = holm(p) if method == "Holm" else benjamini_hochberg(p)
    reject = adj < alpha
    result.loc[sel, "method"] = method
    result.loc[sel, "p_adj"] = adj
    result.loc[sel, "reject"] = reject
    level = 1 - alpha / m if method == "Holm" else 1 - max(reject.sum(), 1) * alpha / m
    result.loc[sel, "ci_level"] = level
    for i in result.index[sel]:
        d, kind = result.at[i, "_d"], result.at[i, "_kind"]
        est, lo, hi = interval(d, level if kind == "difference" else 1 - 2 * (1 - level))
        result.at[i, "hodges_lehmann"], result.at[i, "ci_low"], result.at[i, "ci_high"] = est, lo, hi
result["m_in_family"] = result.groupby("family")["test"].transform("count")
result = result.drop(columns=["_d", "_kind"])
result.to_csv(out / "tests.csv", index=False, float_format="%.6g")
print(result[["family", "test", "mean_difference", "hodges_lehmann", "ci_low", "ci_high", "p_raw", "p_adj", "reject"]].round(4).to_string(index=False))
