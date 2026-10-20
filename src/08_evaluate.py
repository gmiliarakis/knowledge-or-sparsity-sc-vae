# Score every run: how well its interferon latents separate stimulated from control test cells. Rules are in the README
# (Evaluation). Writes one row per run, the per-cell-type AUROCs and the paired comparisons to reports/tables/. With
# --learning-curves, the runs on training subsamples are scored by the same rules, together with the main runs as the
# full-size point, into separate tables.
import argparse
import json

import numpy as np
import pandas as pd

from scoring import (alpha_name, auroc_by_type as auroc, full, min_variance, results_dir, rows, tables_dir, types, val_all)

parser = argparse.ArgumentParser()
parser.add_argument("--learning-curves", action="store_true")
args = parser.parse_args()

variants = {"plain": "vanilla", "real": "real", "shuffled": "shuffled", "setshuffle": "setshuffle", "coexpression": "coexpression"}


scores, per_type = [], []
for folder, variant in variants.items():
    for run in sorted((results_dir / folder).glob("seed*")):
        if "_train" in run.name and not args.learning_curves:
            continue
        info = json.load(open(run / "run.json"))
        model = np.load(run / "model.npz")
        z, w, ids = model["latents"], model["w"], info["latent_ids"]
        activity = z[val_all].var(axis=0)
        active = activity > min_variance

        if folder in ("real", "shuffled", "setshuffle"):
            # named latent, oriented without labels: mean weight on the set's own genes made positive
            alpha = ids.index(alpha_name)
            sign = np.sign(w[w[:, alpha] != 0, alpha].mean())
        else:
            # no meaningful names: rank active latents by validation AUROC, either direction
            val_auc = np.array([auroc(z[:, k], "val").mean() for k in range(z.shape[1])])
            order = [k for k in np.argsort(-np.abs(val_auc - 0.5)) if active[k]]
            alpha = order[0]
            sign = 1.0 if val_auc[alpha] > 0.5 else -1.0

        measured = auroc(sign * z[:, alpha], "test")
        alpha_auc = measured if active[alpha] else pd.Series(0.5, index=types)

        seed, size = info["seed"], info.get("train_size") or full
        scores.append({
            "variant": variant, "train_size": size, "seed": seed,
            "ifn_alpha": alpha_auc.mean(), "ifn_alpha_measured": measured.mean(), "ifn_alpha_latent": ids[alpha],
            "ifn_alpha_sign": int(sign), "ifn_alpha_activity": activity[alpha], "ifn_alpha_active": active[alpha],
            "active_latents": int(active.sum()), "best_epoch": info["best_epoch"],
            "val_loss": info["val_loss"], "test_loss": info["test_loss"],
        })
        per_type += [{"variant": variant, "train_size": size, "seed": seed, "score": "ifn_alpha", "cell_type": t, "auroc": a}
                     for t, a in alpha_auc.items()]

scores = pd.DataFrame(scores).sort_values("train_size", kind="stable")
print(scores.groupby(["train_size", "variant"], sort=False)[["ifn_alpha", "active_latents"]].agg(["mean", "count"]).round(3))

# Paired differences by seed (descriptive). Hypothesis tests, intervals and multiplicity correction are in src/11_tests.py.
comparisons = []
for size in scores["train_size"].unique():
    wide = scores[scores["train_size"] == size].pivot(index="seed", columns="variant", values="ifn_alpha")
    for other in ["shuffled", "setshuffle", "vanilla", "coexpression"]:
        if "real" not in wide or other not in wide:
            continue
        d = (wide["real"] - wide[other]).dropna()
        comparisons.append({"train_size": size, "comparison": f"real - {other}", "pairs": len(d),
                            "mean_difference": d.mean(), "seeds_real_higher": int((d > 0).sum())})
comparisons = pd.DataFrame(comparisons)
print(comparisons.round(4).to_string(index=False))

out = tables_dir
out.mkdir(parents=True, exist_ok=True)
prefix = "learning_curves_" if args.learning_curves else ""
scores.to_csv(out / f"{prefix}scores.csv", index=False, float_format="%.6g")
pd.DataFrame(per_type).to_csv(out / f"{prefix}cell_types.csv", index=False, float_format="%.6g")
comparisons.to_csv(out / f"{prefix}comparisons.csv", index=False, float_format="%.6g")
