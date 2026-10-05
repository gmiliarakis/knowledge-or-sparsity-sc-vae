# Score every run: how well its interferon latents separate stimulated from control test cells. Rules are in the README
# (Evaluation). Writes one row per run, the per-cell-type AUROCs and the paired comparisons to reports/tables/. With
# --learning-curves, the runs on training subsamples are scored by the same rules, together with the main runs as the
# full-size point, into separate tables.
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

parser = argparse.ArgumentParser()
parser.add_argument("--learning-curves", action="store_true")
args = parser.parse_args()

variants = {"plain": "vanilla", "real": "real", "shuffled": "shuffled", "coexpression": "coexpression"}
alpha_name, gamma_name = "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"
min_variance = 0.01  # a latent whose posterior means vary less across validation cells is inactive

obs = ad.read_h5ad("data/processed/kang.h5ad", backed="r").obs
stim = (obs["condition"] == "stim").values
cell_type = obs["cell_type"].astype(str).values
types = sorted(set(cell_type) - {"Megakaryocytes"})  # no nucleus, so no interferon response to find
rows = {s: np.flatnonzero((obs["split"] == s).values & (cell_type != "Megakaryocytes")) for s in ["val", "test"]}
val_all = np.flatnonzero((obs["split"] == "val").values)
full = int((obs["split"] == "train").sum())


def auroc(score, split):
    # stim vs ctrl within each cell type, so a latent cannot score by separating cell types
    r = rows[split]
    return pd.Series({t: roc_auc_score(stim[r][cell_type[r] == t], score[r][cell_type[r] == t]) for t in types})


scores, per_type = [], []
for folder, variant in variants.items():
    for run in sorted(Path("results", folder).glob("seed*")):
        if "_train" in run.name and not args.learning_curves:
            continue
        info = json.load(open(run / "run.json"))
        model = np.load(run / "model.npz")
        z, w, ids = model["latents"], model["w"], info["latent_ids"]
        activity = z[val_all].var(axis=0)
        active = activity > min_variance

        if folder in ("real", "shuffled"):
            # named latents, oriented without labels: mean weight on the set's own genes made positive
            alpha, gamma = ids.index(alpha_name), ids.index(gamma_name)
            sign = np.sign(w[w[:, alpha] != 0, alpha].mean())
            pair = [alpha, gamma]
        else:
            # no meaningful names: rank active latents by validation AUROC, either direction
            val_auc = np.array([auroc(z[:, k], "val").mean() for k in range(z.shape[1])])
            order = [k for k in np.argsort(-np.abs(val_auc - 0.5)) if active[k]]
            alpha, pair = order[0], order[:2]
            sign = 1.0 if val_auc[alpha] > 0.5 else -1.0

        measured = auroc(sign * z[:, alpha], "test")
        alpha_auc = measured if active[alpha] else pd.Series(0.5, index=types)

        scaler = StandardScaler().fit(z[rows["val"]][:, pair])
        classifier = LogisticRegression().fit(scaler.transform(z[rows["val"]][:, pair]), stim[rows["val"]])
        pair_auc = auroc(classifier.decision_function(scaler.transform(z[:, pair])), "test")

        seed, size = info["seed"], info.get("train_size") or full
        scores.append({
            "variant": variant, "train_size": size, "seed": seed,
            "ifn_alpha": alpha_auc.mean(), "ifn_alpha_measured": measured.mean(), "ifn_alpha_latent": ids[alpha],
            "ifn_alpha_sign": int(sign), "ifn_alpha_activity": activity[alpha], "ifn_alpha_active": active[alpha],
            "ifn_pair": pair_auc.mean(), "pair_latents": "; ".join(ids[k] for k in pair),
            "pair_activity": "; ".join(f"{activity[k]:.3g}" for k in pair),
            "active_latents": int(active.sum()), "best_epoch": info["best_epoch"],
            "val_loss": info["val_loss"], "test_loss": info["test_loss"],
        })
        for score, values in [("ifn_alpha", alpha_auc), ("ifn_pair", pair_auc)]:
            per_type += [{"variant": variant, "train_size": size, "seed": seed, "score": score, "cell_type": t, "auroc": a}
                         for t, a in values.items()]

scores = pd.DataFrame(scores).sort_values("train_size", kind="stable")
print(scores.groupby(["train_size", "variant"], sort=False)[["ifn_alpha", "ifn_pair", "active_latents"]].agg(["mean", "count"]).round(3))

# Paired by seed. One primary test (interferon-alpha score, real vs shuffled, all training cells); everything else is
# descriptive.
comparisons = []
for size, score in [(n, s) for n in scores["train_size"].unique() for s in ["ifn_alpha", "ifn_pair"]]:
    wide = scores[scores["train_size"] == size].pivot(index="seed", columns="variant", values=score)
    for other in ["shuffled", "vanilla", "coexpression"]:
        if "real" not in wide or other not in wide:
            continue
        d = (wide["real"] - wide[other]).dropna()
        n = len(d)
        half = stats.t.ppf(0.975, n - 1) * d.std() / np.sqrt(n) if n > 1 else np.nan
        p = stats.wilcoxon(d).pvalue if (d != 0).sum() > 1 else np.nan  # zero differences are dropped
        comparisons.append({"train_size": size, "score": score, "comparison": f"real - {other}", "pairs": n,
                            "mean_difference": d.mean(), "ci_low": d.mean() - half, "ci_high": d.mean() + half, "wilcoxon_p": p,
                            "primary": score == "ifn_alpha" and other == "shuffled" and size == full})
comparisons = pd.DataFrame(comparisons)
print(comparisons.round(4).to_string(index=False))

primary = comparisons[comparisons["primary"]]
if len(primary):
    r = primary.iloc[0]
    verdict = "real beats shuffled" if r["ci_low"] > 0 else "shuffled beats real" if r["ci_high"] < 0 else "no difference shown"
    print(f"primary test ({r['pairs']} pairs): real - shuffled {r['mean_difference']:+.4f}, "
          f"95% CI [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}], Wilcoxon p = {r['wilcoxon_p']:.3g}: {verdict}")

out = Path("reports/tables")
out.mkdir(parents=True, exist_ok=True)
prefix = "learning_curves_" if args.learning_curves else ""
scores.to_csv(out / f"{prefix}scores.csv", index=False, float_format="%.6g")
pd.DataFrame(per_type).to_csv(out / f"{prefix}cell_types.csv", index=False, float_format="%.6g")
comparisons.to_csv(out / f"{prefix}comparisons.csv", index=False, float_format="%.6g")
