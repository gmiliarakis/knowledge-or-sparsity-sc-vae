# Every hypothesis test of the project, with multiplicity correction. Run after 06_evaluate.py (both modes) and
# 07_checks.py. Writes tests.csv in reports/tables/<data>/<collection>/; the preprint quotes this table and no other p-value.
#
# All tests are paired by seed (10 seeds). The p-value is the exact two-sided Wilcoxon signed-rank test; zero differences
# are dropped. Its smallest possible value with 10 pairs is 2/2^10 = 0.00195, so a Holm or Bonferroni correction over more
# than 25 tests could never reject. That fixes how the tests are grouped:
#
# Family S1 (interferon-α latent at full size, confirmatory, 10 tests): Holm step-down, which controls the family-wise error rate
#   at 0.05 under any dependence between the tests (they share seeds). Intervals are the distribution-free intervals that
#   belong to the signed-rank test (Hodges-Lehmann estimate, exact interval from the Walsh averages) at the Bonferroni level
#   1 - 0.05/m, which are simultaneous for the m tests. Holm is less conservative than Bonferroni, so an interval can touch 0
#   for a test that Holm rejects.
# Family S3 (learning curves, exploratory, 18 tests, not reported in the preprint; the code is commented out below): Benjamini-Hochberg, which controls the false discovery rate at 0.05
#   and is valid for positively dependent tests (the same seeds and nested subsamples run through all sizes). Intervals are
#   Benjamini-Yekutieli false-coverage-rate intervals, level 1 - R*q/m with R the number of rejections, of the same kind.
# Family latent-free (Q1 to Q1c; 15 tests per benchmark, confirmatory): Holm, as S1. Five measures that choose no latent (test loss, decoded
#   response with the Spearman and with the Pearson correlation, cell-type probe, stimulation probe), each with three contrasts: Hallmark minus each of the two random masks and
#   minus the co-expression mask. The no-mask VAE is a descriptive reference, and the few-shot and hard-pair cell-type scores
#   are descriptive.
# Family responsematch (response-matched random mask, control for Q1; 5 tests per benchmark, confirmatory): Holm, as S1. The five measures
#   that choose no latent, Hallmark minus the response-matched mask. It is its own family so that the 15 tests of the latent-free
#   family keep their correction; it exists only where the response-matched runs exist (the four benchmarks, not the 250-gene repeat).
# Per-cell-type and per-donor scores are descriptive and get no p-values.
# Symmetry check: every row also carries the exact sign-test interval at the same level (no symmetry assumption), the skewness of
# the paired differences, and whether the sign-test interval and the Hodges-Lehmann interval agree on excluding zero.
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
from gene_sets import main, run  # SC_DATA, SC_SETS and SC_TAG select the benchmark
out = Path(f"reports/tables/{run}")  # outside the main experiment only the latent-free family is tested

if main:
    scores = pd.read_csv(out / "scores.csv")
    best = pd.read_csv(out / "robustness_best_latent.csv").replace({"variant": {"plain": "vanilla"}})  # that table says plain
    # Learning-curve inputs of the S3 tests (not reported in the paper): uncomment with the S3 tests below after running src/extras/learning_curves.
    # curves = pd.read_csv(out / "learning_curves_scores.csv")
    # slopes_input = curves.pivot_table(index=["seed", "train_size"], columns="variant", values="ifn_alpha")
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


def sign_interval(d, level):
    """Exact interval of at least the given level from the sign test (order statistics of the differences); needs no symmetry."""
    d = np.sort(np.asarray(d, dtype=float))
    n = len(d)
    cdf = np.cumsum([stats.binom.pmf(i, n, 0.5) for i in range(n + 1)])
    k = int(np.searchsorted(cdf, (1 - level) / 2, side="right")) - 1  # largest k with P(B <= k) <= alpha/2
    if k < 0:
        return d[0], d[-1]
    return d[k], d[n - 1 - k]


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
if main:
    at_full = scores[scores["train_size"] == full]
    for other in ["shuffled", "setshuffle", "vanilla", "coexpression"]:
        tests.append(("S1", f"fixed latent: real - {other}", diff(at_full, "ifn_alpha", "real", other), "difference"))
    for a, b in [("real", "shuffled"), ("real", "setshuffle"), ("real", "vanilla"), ("shuffled", "vanilla"),
                 ("setshuffle", "vanilla")]:
        tests.append(("S1", f"best latent: {a} - {b}", diff(best, "best_latent", a, b), "difference"))
    tests.append(("S1", f"equivalence (+-{margin}): real - coexpression", diff(at_full, "ifn_alpha", "real", "coexpression"),
                  "equivalence"))
    # Learning-curve tests (family S3, not reported in the paper): uncomment with the two lines that read `curves` above.
    # for size in sorted(curves["train_size"].unique()):
    #     if size == full:
    #         continue
    #     for other in ["shuffled", "vanilla", "coexpression"]:
    #         tests.append(("S3", f"{size} cells: real - {other}", diff(curves[curves["train_size"] == size], "ifn_alpha", "real", other), "difference"))
    # for other in ["shuffled", "vanilla", "coexpression"]:
    #     per_seed = []
    #     for seed, g in slopes_input.groupby(level=0):
    #         d = (g["real"] - g[other]).droplevel(0)
    #         per_seed.append(np.polyfit(np.log2(d.index.values), d.values, 1)[0])
    #     tests.append(("S3", f"slope per doubling of cells: real - {other}", np.array(per_seed), "difference"))

latent_free = [("test loss per cell (lower is better)", "diagnostics_fit.csv", "test_loss"),
               ("decoded response, Spearman correlation", "diagnostics_response.csv", "lfc_spearman"),
               ("decoded response, Pearson correlation", "diagnostics_response.csv", "lfc_pearson"),
               ("cell-type balanced accuracy", "diagnostics_celltype.csv", "balanced_acc_full"),
               ("stimulation probe AUROC", "diagnostics_probe.csv", "probe_auroc")]
for label, file, column in latent_free:
    frame = pd.read_csv(out / file)
    for other in ["shuffled", "setshuffle", "coexpression"]:
        tests.append(("latent-free", f"{label}: real - {other}", diff(frame, column, "real", other), "difference"))
    if "responsematch" in set(frame["variant"]):
        tests.append(("responsematch", f"{label}: real - responsematch", diff(frame, column, "real", "responsematch"), "difference"))

rows = []
for family, name, d, kind in tests:
    if kind == "equivalence":
        p = max(wilcoxon(d + margin, alternative="greater"), wilcoxon(d - margin, alternative="less"))
    else:
        p = wilcoxon(d)
    rows.append({"family": family, "test": name, "pairs": len(d), "mean_difference": d.mean(), "p_raw": p, "_d": d, "_kind": kind})
result = pd.DataFrame(rows)

for family, method in [("S1", "Holm"), ("S3", "Benjamini-Hochberg"), ("latent-free", "Holm"), ("responsematch", "Holm")]:
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
        s_lo, s_hi = sign_interval(d, level if kind == "difference" else 1 - 2 * (1 - level))
        result.at[i, "sign_ci_low"], result.at[i, "sign_ci_high"] = s_lo, s_hi
        result.at[i, "skew"] = stats.skew(d) if len(d) > 2 and np.ptp(d) > 0 else 0.0
        result.at[i, "sign_same_decision"] = bool((lo > 0 or hi < 0) == (s_lo > 0 or s_hi < 0))
result["m_in_family"] = result.groupby("family")["test"].transform("count")
result = result.drop(columns=["_d", "_kind"])
result.to_csv(out / "tests.csv", index=False, float_format="%.6g")
print(result[["family", "test", "mean_difference", "hodges_lehmann", "ci_low", "ci_high", "p_raw", "p_adj", "reject", "skew", "sign_same_decision"]].round(4).to_string(index=False))


# Share of the gap between the fully random mask and the vanilla VAE that the gene set mask recovers, for the decoded response with the Spearman
# and with the Pearson correlation: the ratio of two Hodges-Lehmann estimates of the paired differences (gene set mask minus the other mask,
# vanilla minus the other mask), the other mask being the fully random or the within-set random mask. The interval is the percentile interval of
# 2,000 resamples of the 10 seeds, drawn with the same seeds for every variant; it is conditional on the fixed split. Writes shares.csv.
rng = np.random.default_rng(0)
tri = np.triu_indices(10)


def hodges_lehmann(d):
    return float(np.median(((d[:, None] + d[None, :]) / 2)[tri]))


share_rows = []
response = pd.read_csv(out / "diagnostics_response.csv")
for measure, column in [("spearman", "lfc_spearman"), ("pearson", "lfc_pearson")]:
    w = response.pivot(index="seed", columns="variant", values=column).sort_index()
    for other, label in [("shuffled", "fully random"), ("setshuffle", "within-set random")]:
        if other not in w:
            continue
        num, den = (w["real"] - w[other]).to_numpy(), (w["vanilla"] - w[other]).to_numpy()
        boots = []
        for _ in range(2000):
            i = rng.integers(0, len(num), len(num))
            boots.append(hodges_lehmann(num[i]) / hodges_lehmann(den[i]) * 100)
        lo, hi = np.percentile(boots, [2.5, 97.5])
        share_rows.append({"quantity": f"share of the gap to the vanilla VAE recovered, against the {label} mask (%)", "measure": measure,
                           "estimate": hodges_lehmann(num) / hodges_lehmann(den) * 100, "ci_low": lo, "ci_high": hi})
pd.DataFrame(share_rows).to_csv(out / "shares.csv", index=False, float_format="%.4f")
