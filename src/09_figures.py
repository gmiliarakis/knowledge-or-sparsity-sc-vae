# Figures from the evaluation tables (src/08_evaluate.py) and the training logs, written to reports/figures/. The
# learning-curve figure is drawn once its tables exist (src/08_evaluate.py --learning-curves).
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt

surface, ink, secondary, muted, grid, baseline = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
colors = {"vanilla": "#999999", "real": "#0072B2", "shuffled": "#D55E00", "coexpression": "#009E73"}  # Okabe-Ito, colourblind-safe
names = {"vanilla": "vanilla VAE", "real": "real mask", "shuffled": "shuffled mask", "coexpression": "co-expression"}
ticks = {"vanilla": "vanilla\nVAE", "real": "real\nmask", "shuffled": "shuffled\nmask", "coexpression": "co-\nexpression"}
folders = {"vanilla": "plain", "real": "real", "shuffled": "shuffled", "coexpression": "coexpression"}
order = list(colors)
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "axes.titleweight": "medium", "axes.titlelocation": "left",
    "axes.titlepad": 10, "axes.labelsize": 8, "text.color": ink, "axes.labelcolor": secondary,
    "xtick.color": baseline, "ytick.color": baseline, "xtick.labelcolor": secondary, "ytick.labelcolor": secondary,
    "xtick.major.size": 0, "ytick.major.size": 0, "xtick.major.pad": 4, "ytick.major.pad": 4,
    "axes.edgecolor": baseline, "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": grid, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "figure.facecolor": surface, "axes.facecolor": surface, "savefig.facecolor": surface,
    "legend.frameon": False, "legend.fontsize": 8, "lines.solid_capstyle": "round",
    "savefig.dpi": 200, "savefig.bbox": "tight", "pdf.fonttype": 42,
})

tables, out = Path("reports/tables"), Path("reports/figures")
out.mkdir(parents=True, exist_ok=True)
scores = pd.read_csv(tables / "scores.csv")
cell_types = pd.read_csv(tables / "cell_types.csv")
comparisons = pd.read_csv(tables / "comparisons.csv")
offset = dict(zip(range(10), np.linspace(-0.18, 0.18, 10)))  # each seed keeps its position in every panel


def strip(ax, column, hollow=None):
    # one dot per seed, the mean as a short bar, and grey lines joining real and shuffled runs of the same seed
    x = {v: order.index(v) for v in order}
    wide = scores.pivot(index="seed", columns="variant", values=column)
    if {"real", "shuffled"} <= set(wide):
        for seed, (a, b) in wide[["real", "shuffled"]].dropna().iterrows():
            ax.plot([x["real"] + offset[seed], x["shuffled"] + offset[seed]], [a, b], color=baseline, lw=0.7, zorder=1)
    for v in order:
        runs = scores[scores["variant"] == v]
        if runs.empty:
            continue
        xs = x[v] + runs["seed"].map(offset)
        empty = runs[hollow].eq(False) if hollow else pd.Series(False, index=runs.index)
        ax.scatter(xs[~empty], runs.loc[~empty, column], s=24, color=colors[v], edgecolor=surface, lw=0.8, zorder=3)
        ax.scatter(xs[empty], runs.loc[empty, column], s=24, color=surface, edgecolor=colors[v], lw=1.2, zorder=3)
        ax.plot([x[v] - 0.3, x[v] + 0.3], [runs[column].mean()] * 2, color=ink, lw=1.6, zorder=4)
    ax.set_xticks(range(len(order)), [ticks[v] for v in order])
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.grid(axis="x", visible=False)


def chance(ax):
    # AUROC 0.5 only when the data come near it, so close values are not squashed
    low, high = ax.get_ylim()
    if low < 0.55:
        ax.axhline(0.5, color=muted, lw=0.8, zorder=0)
        ax.set_ylim(min(low, 0.48), high)


# 1. Scores per seed and the paired differences
panels = plt.figure(figsize=(7.4, 5.4)).subplot_mosaic([["a", "b"], ["c", "c"]], gridspec_kw={"hspace": 0.55})
fig, axes = panels["a"].figure, [panels["a"], panels["b"], panels["c"]]
strip(axes[0], "ifn_alpha", hollow="ifn_alpha_active")
axes[0].set(title="a  Interferon-α score", ylabel="AUROC stim vs ctrl\n(mean over 7 cell types)")
strip(axes[1], "ifn_pair")
axes[1].set(title="b  Interferon pair score")
for ax in axes[:2]:
    chance(ax)
low = min(a.get_ylim()[0] for a in axes[:2])
for ax in axes[:2]:
    ax.set_ylim(low, 1.0 + (1.0 - low) * 0.04)

ax = axes[2]
groups = {"ifn_alpha": "interferon-α", "ifn_pair": "interferon pair"}
y, rows = 0, []
for score in groups:
    ax.text(-0.04, y - 0.9, groups[score], transform=ax.get_yaxis_transform(), ha="right", color=ink,
            fontweight="medium")
    for _, r in comparisons[comparisons["score"] == score].iterrows():
        color, size = (ink, 36) if r["primary"] else (muted, 24)
        ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=color, lw=2)
        ax.scatter(r["mean_difference"], y, s=size, color=color, edgecolor=surface, lw=0.8, zorder=3)
        rows.append((y, r["comparison"].replace(" - ", " − ").replace("coexpression", "co-expression")
                     + (" (primary)" if r["primary"] else "")))
        y += 1
    y += 1
ax.axvline(0, color=baseline, lw=0.8, zorder=0)
ax.set_yticks(*zip(*rows))
ax.set_ylim(y - 1.5, -1.5)
ax.grid(axis="y", visible=False)
ax.set(title="c  Paired difference, 95% CI", xlabel="difference in AUROC")
fig.savefig(out / "scores.png")
plt.close(fig)

# 2. Per cell type
types = sorted(cell_types["cell_type"].unique())
fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.6), sharey=True, gridspec_kw={"wspace": 0.08})
for ax, (score, title) in zip(axes, [("ifn_alpha", "a  Interferon-α score"), ("ifn_pair", "b  Interferon pair score")]):
    part = cell_types[cell_types["score"] == score]
    for i, v in enumerate(order):
        runs = part[part["variant"] == v]
        if runs.empty:
            continue
        ys = runs["cell_type"].map(types.index) + (i - 1.5) * 0.17
        ax.scatter(runs["auroc"], ys, s=8, color=colors[v], alpha=0.45, lw=0, zorder=2)
        means = runs.groupby("cell_type")["auroc"].mean()
        ax.scatter(means.values, [types.index(t) + (i - 1.5) * 0.17 for t in means.index], s=34, color=colors[v],
                   edgecolor=surface, lw=0.8, zorder=3, label=names[v])
    ax.set_yticks(range(len(types)), types)
    ax.set_ylim(len(types) - 0.5, -0.5)
    ax.grid(axis="y", visible=False)
    ax.set(title=title, xlabel="AUROC stim vs ctrl, test cells")
axes[0].legend(loc="upper left", bbox_to_anchor=(0, -0.16), ncol=4, handletextpad=0.2, columnspacing=1.4)
fig.savefig(out / "cell_types.png")
plt.close(fig)

# 3. Training curves: validation loss per epoch, best epoch marked
fig, axes = plt.subplots(1, 4, figsize=(7.4, 2.4), sharey=True, gridspec_kw={"wspace": 0.08})
curves = {v: [] for v in order}
for v in order:
    for run in sorted(Path("results", folders[v]).glob("seed*")):
        if "_train" in run.name:
            continue
        curves[v].append((pd.read_csv(run / "losses.csv"), json.load(open(run / "run.json"))["best_epoch"]))
warm = next(h["epoch"][h["kl_weight"] >= 1].iloc[0] for c in curves.values() for h, _ in c)
top = max(h["val_loss"][h["epoch"] >= warm].max() for c in curves.values() for h, _ in c)
bottom = min(h["val_loss"].min() for c in curves.values() for h, _ in c)
for ax, v in zip(axes, order):
    ax.axvspan(0, warm, color=grid, alpha=0.5, lw=0, zorder=0)
    for h, best in curves[v]:
        ax.plot(h["epoch"], h["val_loss"], color=colors[v], lw=1, alpha=0.5)
        ax.scatter(best, h["val_loss"][h["epoch"] == best], s=14, color=colors[v], edgecolor=surface, lw=0.6, zorder=3)
    ax.set(title=names[v], xlabel="epoch", xlim=(0, None))
    ax.set_ylim(bottom - 0.1 * (top - bottom), top + 0.1 * (top - bottom))
    ax.grid(axis="x", visible=False)
axes[0].set_ylabel("validation loss per cell")
axes[0].text(warm / 2, 0.03, "KL\nwarm-up", transform=axes[0].get_xaxis_transform(), ha="center", color=muted,
             fontsize=7)
fig.savefig(out / "training.png")
plt.close(fig)

# 4. Active latents per run (posterior means vary by more than 0.01 across validation cells)
fig, ax = plt.subplots(figsize=(3.4, 2.6))
strip(ax, "active_latents")
ax.set(title="Active latents per run", ylabel="active latents (of 58)")
fig.savefig(out / "activity.png")
plt.close(fig)

# 5. Learning curves: mean and 95% t-interval over seeds at each training size, and the real - shuffled difference
if (tables / "learning_curves_scores.csv").exists():
    all_scores = pd.read_csv(tables / "learning_curves_scores.csv")
    all_comparisons = pd.read_csv(tables / "learning_curves_comparisons.csv")
    full = scores["train_size"].max()
    sizes = sorted(all_scores["train_size"].unique())
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7), gridspec_kw={"wspace": 0.4})
    for ax, (score, title) in zip(axes, [("ifn_alpha", "a  Interferon-α score"), ("ifn_pair", "b  Interferon pair score")]):
        for v in order:
            g = all_scores[all_scores["variant"] == v].groupby("train_size")[score]
            mean, half = g.mean(), stats.t.ppf(0.975, g.count() - 1) * g.sem()
            ax.fill_between(mean.index, mean - half, mean + half, color=colors[v], alpha=0.15, lw=0)
            ax.plot(mean.index, mean, color=colors[v], lw=2, marker="o", ms=4.5, mec=surface, mew=0.8, label=names[v])
        ax.set(title=title)
    axes[0].set_ylabel("AUROC stim vs ctrl\n(mean over 7 cell types)")
    axes[0].legend(loc="upper left", bbox_to_anchor=(0, -0.24), ncol=4, handletextpad=0.4, columnspacing=1.4)
    ax = axes[2]
    for (score, color), shift in zip([("ifn_alpha", ink), ("ifn_pair", muted)], [0.97, 1.03]):
        d = all_comparisons[(all_comparisons["score"] == score) & (all_comparisons["comparison"] == "real - shuffled")]
        x = d["train_size"] * shift  # side by side on the log axis
        ax.vlines(x, d["ci_low"], d["ci_high"], color=color, lw=2)
        ax.scatter(x, d["mean_difference"], s=24, color=color, edgecolor=surface, lw=0.8, zorder=3,
                   label={"ifn_alpha": "interferon-α", "ifn_pair": "interferon pair"}[score])
    ax.axhline(0, color=baseline, lw=0.8, zorder=0)
    ax.set(title="c  Real − shuffled, 95% CI", ylabel="difference in AUROC")
    ax.legend(loc="best", handletextpad=0.2)
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xticks(sizes, [f"{n:,}" if n < full else f"all\n{n:,}" for n in sizes], fontsize=7)
        ax.minorticks_off()
        ax.set_xlabel("training cells")
        ax.grid(axis="x", visible=False)
    fig.savefig(out / "learning_curves.png")
    plt.close(fig)

print("saved", *sorted(p.name for p in out.glob("*.png")))
