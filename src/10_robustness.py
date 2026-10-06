# Checks on the main result. Run after 08_evaluate.py --learning-curves.
# A. The random mask scored like the no-mask VAE: the best active latent, chosen on validation cells.
# B. Per-donor scores of the same latents (test cells), to see whether the result holds in every donor.
# C. Slope of the real - other difference on log2(training size), per seed, for the learning curves (Q1b).
# D. Equivalence of the Hallmark and co-expression masks: 90% paired t-interval against a margin of +-0.005.
# Tables go to reports/tables/robustness_*.csv.
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

alpha_name = "HALLMARK_INTERFERON_ALPHA_RESPONSE"
min_variance, margin = 0.01, 0.005
out = Path("reports/tables")

obs = ad.read_h5ad("data/processed/kang.h5ad", backed="r").obs
stim = (obs["condition"] == "stim").values
cell_type = obs["cell_type"].astype(str).values
donor = obs["donor"].astype(str).values
types = sorted(set(cell_type) - {"Megakaryocytes"})
rows = {s: np.flatnonzero((obs["split"] == s).values & (cell_type != "Megakaryocytes")) for s in ["val", "test"]}
val_all = np.flatnonzero((obs["split"] == "val").values)


def auroc(score, split, mask=None):
    r = rows[split] if mask is None else rows[split][mask[rows[split]]]
    values = []
    for t in types:
        k = r[cell_type[r] == t]
        if 0 < stim[k].sum() < len(k):
            values.append(roc_auc_score(stim[k], score[k]))
    return np.mean(values) if values else np.nan


# A and B ---------------------------------------------------------------------------------------------------------
scored, per_donor = [], []
for folder in ["plain", "real", "shuffled", "coexpression"]:
    for run in sorted(Path("results", folder).glob("seed?")):
        info = json.load(open(run / "run.json"))
        model = np.load(run / "model.npz")
        z, w, ids = model["latents"], model["w"], info["latent_ids"]
        active = z[val_all].var(axis=0) > min_variance
        val_auc = np.array([auroc(z[:, k], "val") for k in range(z.shape[1])])
        best = [k for k in np.argsort(-np.abs(val_auc - 0.5)) if active[k]][0]
        best_sign = 1.0 if val_auc[best] > 0.5 else -1.0
        if folder in ("real", "shuffled"):
            named = ids.index(alpha_name)
            named_sign = np.sign(w[w[:, named] != 0, named].mean())
        else:
            named, named_sign = best, best_sign  # no names: the chosen latent is the score of the main analysis
        scored.append({"variant": folder, "seed": info["seed"],
                       "named_latent": auroc(named_sign * z[:, named], "test") if active[named] else 0.5,
                       "best_latent": auroc(best_sign * z[:, best], "test"),
                       "best_is_named": best == named})
        for d in sorted(set(donor)):
            per_donor.append({"variant": folder, "seed": info["seed"], "donor": d,
                              "auroc": auroc(named_sign * z[:, named], "test", donor == d) if active[named] else 0.5})

scored = pd.DataFrame(scored)
scored.to_csv(out / "robustness_best_latent.csv", index=False, float_format="%.6g")
print(scored.groupby("variant")[["named_latent", "best_latent"]].agg(["mean", "std"]).round(3))
print("shuffled seeds where the best latent is the named one:", int(scored.query("variant == 'shuffled'")["best_is_named"].sum()), "of 10")
wide = scored.pivot(index="seed", columns="variant", values="best_latent")
d = wide["real"] - wide["shuffled"]
half = stats.t.ppf(0.975, 9) * d.std() / np.sqrt(10)
print(f"real - shuffled, both scored by the best active latent: {d.mean():+.4f} ({d.mean() - half:+.4f} to {d.mean() + half:+.4f}), "
      f"Wilcoxon p = {stats.wilcoxon(d).pvalue:.3g}, real higher in {(d > 0).sum()} of 10 seeds")

per_donor = pd.DataFrame(per_donor)
per_donor.to_csv(out / "robustness_donors.csv", index=False, float_format="%.6g")
by_donor = per_donor.groupby(["donor", "variant"])["auroc"].mean().unstack()
by_donor["real - shuffled"] = by_donor["real"] - by_donor["shuffled"]
print(by_donor.round(3))
by_donor.round(4).to_csv(out / "robustness_donors_mean.csv")

# C ---------------------------------------------------------------------------------------------------------------
lc = pd.read_csv(out / "learning_curves_scores.csv")
wide = lc.pivot_table(index=["seed", "train_size"], columns="variant", values="ifn_alpha")
slopes = []
for other in ["shuffled", "vanilla", "coexpression"]:
    per_seed = []
    for seed, g in wide.groupby(level=0):
        diff = (g["real"] - g[other]).droplevel(0)
        per_seed.append(np.polyfit(np.log2(diff.index.values), diff.values, 1)[0])
    s = np.array(per_seed)
    half = stats.t.ppf(0.975, len(s) - 1) * s.std(ddof=1) / np.sqrt(len(s))
    slopes.append({"comparison": f"real - {other}", "seeds": len(s), "slope_per_doubling": s.mean(),
                   "ci_low": s.mean() - half, "ci_high": s.mean() + half, "wilcoxon_p": stats.wilcoxon(s).pvalue})
slopes = pd.DataFrame(slopes)
slopes.to_csv(out / "robustness_slopes.csv", index=False, float_format="%.6g")
print(slopes.round(4).to_string(index=False))

# D ---------------------------------------------------------------------------------------------------------------
equivalence = []
for size, g in wide.groupby(level=1):
    diff = (g["real"] - g["coexpression"]).dropna()
    half = stats.t.ppf(0.95, len(diff) - 1) * diff.std() / np.sqrt(len(diff))
    lo, hi = diff.mean() - half, diff.mean() + half
    equivalence.append({"train_size": size, "mean_difference": diff.mean(), "ci90_low": lo, "ci90_high": hi,
                        "equivalent_within_margin": bool(lo > -margin and hi < margin)})
equivalence = pd.DataFrame(equivalence)
equivalence.to_csv(out / "robustness_equivalence.csv", index=False, float_format="%.6g")
print(equivalence.round(4).to_string(index=False))
