# Figures from the evaluation tables (src/08_evaluate.py) and the training logs, written to reports/figures/. The
# learning-curve figure is drawn once its tables exist (src/08_evaluate.py --learning-curves). The panel of measures that
# choose no latent needs src/10_checks.py, and the tests in src/11_tests.py supply the intervals.
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt

surface, ink, secondary, muted, grid, baseline = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
colors = {"vanilla": "#cfcfcf", "real": "#0072B2", "shuffled": "#D55E00", "setshuffle": "#E69F00", "coexpression": "#009E73"}  # Okabe-Ito, colourblind-safe
names = {"vanilla": "no mask", "real": "Hallmark", "shuffled": "fully random", "setshuffle": "within-set random", "coexpression": "co-expression"}
ticks = {"vanilla": "no\nmask", "real": "Hallmark", "shuffled": "fully\nrandom", "setshuffle": "within-set\nrandom", "coexpression": "co-\nexpression"}
folders = {"vanilla": "plain", "real": "real", "shuffled": "shuffled", "setshuffle": "setshuffle", "coexpression": "coexpression"}
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
    "savefig.dpi": 400, "savefig.bbox": "tight", "pdf.fonttype": 42,
})

tables, out = Path("reports/tables"), Path("reports/figures")
out.mkdir(parents=True, exist_ok=True)
scores = pd.read_csv(tables / "scores.csv")
tests = pd.read_csv(tables / "tests.csv").set_index("test")  # intervals and corrected p-values (src/11_tests.py)
offset = dict(zip(range(10), np.linspace(-0.18, 0.18, 10)))  # each seed keeps its position in every panel


def strip(ax, column, hollow=None, pairs=True, frame=None):
    # one dot per seed, the mean as a short bar, and grey lines joining Hallmark and random runs of the same seed
    data = scores if frame is None else frame
    x = {v: order.index(v) for v in order}
    wide = data.pivot(index="seed", columns="variant", values=column)
    if pairs and {"real", "shuffled"} <= set(wide):
        for seed, (a, b) in wide[["real", "shuffled"]].dropna().iterrows():
            ax.plot([x["real"] + offset[seed], x["shuffled"] + offset[seed]], [a, b], color=baseline, lw=0.5, zorder=1)
    for v in order:
        runs = data[data["variant"] == v]
        if runs.empty:
            continue
        xs = x[v] + runs["seed"].map(offset)
        empty = runs[hollow].eq(False) if hollow else pd.Series(False, index=runs.index)
        ax.scatter(xs[~empty], runs.loc[~empty, column], s=24, color=colors[v], edgecolor=surface, lw=0.8, zorder=3)
        ax.scatter(xs[empty], runs.loc[empty, column], s=24, color=surface, edgecolor=colors[v], lw=1.2, zorder=3)
        ax.plot([x[v] - 0.3, x[v] + 0.3], [runs[column].mean()] * 2, color=ink, lw=1.1, zorder=4)
    ax.set_xticks(range(len(order)), [ticks[v] for v in order])
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.grid(axis="x", visible=False)


def chance(ax):
    # AUROC 0.5 only when the data come near it, so close values are not squashed
    low, high = ax.get_ylim()
    if low < 0.55:
        ax.axhline(0.5, color=muted, lw=0.8, zorder=0)
        ax.set_ylim(min(low, 0.48), high)


# 0. Overview banner: data, VAE with and without a mask, the five decoder masks as trained with seed 0, and the
# decoded response of every seed
W, H = 8.2, 2.7
fig = plt.figure(figsize=(W, H))
bg = fig.add_axes([0, 0, 1, 1])
bg.set(xlim=(0, W), ylim=(0, H))
bg.axis("off")
rng = np.random.default_rng(3)


def inset(x, y, w, h):
    return fig.add_axes([x / W, y / H, w / W, h / H])


def label(x, y, text, size=6.5, color=secondary, weight="normal", ha="center", va="center"):
    bg.text(x, y, text, fontsize=size, color=color, fontweight=weight, ha=ha, va=va, linespacing=1.15)


def arrow(x0, x1, y=1.35):
    bg.annotate("", (x1, y), (x0, y), arrowprops=dict(arrowstyle="-|>", color=baseline, lw=1.1, mutation_scale=9))


def swarm(y, d_x, d_y):
    # x offsets that keep dots of diameter (d_x, d_y) from overlapping, as close to the centre as possible
    xs, placed = np.zeros(len(y)), []
    for i in np.argsort(y):
        k = 0
        while True:
            found = None
            for c in ([0.0] if k == 0 else [k * d_x * 0.2, -k * d_x * 0.2]):
                if all(((c - px) / d_x) ** 2 + ((y[i] - py) / d_y) ** 2 >= 1 for px, py in placed):
                    found = c
                    break
            if found is not None:
                break
            k += 1
        xs[i] = found
        placed.append((found, y[i]))
    return xs


for x, head in [(0.7, "1  Data"), (2.3, "2  Model"), (4.55, "3  Decoder masks"), (7.2, "4  Decoded response")]:
    label(x, 2.53, head, size=8.5, color=ink, weight="bold")  # each title centred over its panel

# 1 data: two clouds of cells
for (cy, c, text) in [(1.8, baseline, "control"), (0.9, "#CC79A7", "IFN-β, 6 h")]:
    xy = np.clip(rng.normal([0.5, cy], [0.15, 0.11], size=(80, 2)), [0.2, cy - 0.24], [0.8, cy + 0.24])
    bg.scatter(xy[:, 0], xy[:, 1], s=5.5, color=c, lw=0, zorder=2)
    label(0.88, cy, text, ha="left", color=ink)
label(0.66, 0.3, "23,919 cells\nPBMCs from 8 donors")
arrow(1.2, 1.4)

# 2 model: genes -> encoder -> latents -> decoder -> genes; a mask keeps some decoder links and removes the others
by = np.linspace(0.5, 2.2, 15)  # as tall as the encoder
for y in by:
    bg.plot([1.5, 1.65], [y, y], color=muted, lw=0.7, solid_capstyle="round")
    bg.plot([2.98, 3.13], [y, y], color=muted, lw=0.7, solid_capstyle="round")
bg.add_patch(plt.Polygon([(1.75, 0.5), (1.75, 2.2), (2.05, 1.76), (2.05, 0.94)], color=grid, lw=0))
ly = np.linspace(0.98, 1.72, 5)  # as tall as the narrow end of the encoder
kept = [(11,), (2, 8), (5,), (6, 13), (3,)]  # a mask keeps a few links per latent, to any gene


def link(y0, y1):
    return [2.4, 2.98], [y0, y1]


for i, y in enumerate(ly):
    for j, yb in enumerate(by):
        if j not in kept[i]:
            bg.plot(*link(y, yb), color="#a09e96", lw=0.4, ls=(0, (2, 2)), zorder=1)
for i, y in enumerate(ly):
    for j in kept[i]:
        bg.plot(*link(y, by[j]), color=ink, lw=0.9, solid_capstyle="round", zorder=3)
    bg.scatter(2.3, y, s=30, color=surface, edgecolor=secondary, lw=0.8, zorder=4)
label(1.9, 2.3, "encoder")
label(3.0, 2.3, "linear decoder")
label(1.57, 0.3, "12,034\ngenes")
label(2.3, 0.3, "58\nlatents")
bg.plot([1.58, 1.73], [0.04, 0.04], color=ink, lw=0.8)
label(1.78, 0.04, "link kept", ha="left", size=6)
bg.plot([2.32, 2.47], [0.04, 0.04], color="#a09e96", lw=0.5, ls=(0, (2, 2)))
label(2.52, 0.04, "link removed", ha="left", size=6)
arrow(3.29, 3.49)

# 3 the five decoder masks of seed 0, one dot per kept link: rows are the 12,034 genes sorted by Hallmark membership,
# columns the 58 latents (50 named, then 8 free)
hm = np.load("data/processed/masks_hallmark.npz")
co = np.load("data/processed/coexpression_hallmark.npz")
ss = np.load("data/processed/masks_hallmark_setshuffle.npz")
sets = hm["real"]


def full(named):
    m = np.zeros((named.shape[0], 58))
    m[:, :50] = named
    m[:, 50:] = ~named.any(axis=1, keepdims=True)
    return m


rows = np.argsort(np.where(sets.any(axis=1), sets.argmax(axis=1), 50), kind="stable")
position = np.argsort(rows)  # row of each gene in the sorted display
shuffle_row = list(hm["levels"]).index(1.0)
masks = {"real": full(sets), "coexpression": full(co["modules"][0]), "shuffled": full(sets[hm["perms"][shuffle_row, 0]]),
         "setshuffle": full(sets[ss["perms"][0]]), "vanilla": np.ones((sets.shape[0], 58))}
mask_names = {"real": "Hallmark", "coexpression": "co-expression", "shuffled": "fully random", "setshuffle": "within-set random",
              "vanilla": "no mask"}
for k, v in enumerate(["real", "coexpression", "setshuffle", "shuffled", "vanilla"]):
    if v == "vanilla":  # third column, centred on the arrows
        x, y = 4.92, 1.35 - 0.85 / 2
    else:  # the four masks in a 2 x 2 block
        x, y = 3.52 + (k % 2) * 0.65, 1.5 - (k // 2) * 1.1
    ax = inset(x, y, 0.58, 0.85)
    ax.set(xlim=(-0.5, 57.5), ylim=(len(rows) - 0.5, -0.5), xticks=[], yticks=[])
    if v == "vanilla":
        ax.add_patch(plt.Rectangle((-0.5, -0.5), 58, len(rows), color=colors[v], lw=0))  # every gene reaches every latent
    else:
        genes_, latents_ = np.nonzero(masks[v])
        ax.scatter(latents_, position[genes_], s=0.5, marker="s", color=colors[v], lw=0)
    ax.grid(False)
    for side in ax.spines.values():
        side.set_visible(False)
    label(x + 0.29, y - 0.12, mask_names[v], color=ink, size=5.2)
arrow(5.7, 5.9)

# 4 result: decoded response (correlation of the decoded and observed stimulation shift across genes) of every seed as a
# beeswarm, the mean printed under each bar
response = pd.read_csv(tables / "diagnostics_response.csv")
variants = ["real", "coexpression", "setshuffle", "shuffled", "vanilla"]  # same order as the masks in stage 3
width_in, height_in, low, high = 1.9, 1.9, 0.66, 1.0
ax = inset(6.3, 0.45, width_in, height_in)
dot_in = 0.04
for i, v in enumerate(variants):
    y_ = response[response["variant"] == v]["lfc_corr"].to_numpy()
    x_ = i + swarm(y_, dot_in / (width_in / len(variants)), dot_in / (height_in / (high - low)))
    ax.scatter(x_, y_, s=8, color=colors[v], edgecolor=surface, lw=0.3, zorder=3)
    ax.scatter(i, y_.mean(), marker="D", s=14, color=ink, edgecolor=surface, lw=0.6, zorder=4)
    ax.text(i, y_.mean() - 0.025, f"{y_.mean():.2f}", ha="center", va="top", fontsize=6.5, color=ink, zorder=5)
ax.set(xlim=(-0.5, len(variants) - 0.5), ylim=(low, high), yticks=[0.7, 0.8, 0.9, 1.0], yticklabels=["0.7", "0.8", "0.9", "1.0"])
ax.set_ylabel("correlation")
ax.yaxis.set_label_coords(-0.12, 0.5)
ax.tick_params(axis="y", labelsize=6.5)
for side in ("left", "bottom"):
    ax.spines[side].set_color(ink)
ax.set_xticks(range(len(variants)), [names[v] for v in variants], fontsize=5.5, rotation=45, ha="right", rotation_mode="anchor")
ax.grid(axis="x", visible=False)
fig.savefig(out / "overview.png")
plt.close(fig)

# 1. Interferon-α score and active latents per seed, and the paired differences
panels = plt.figure(figsize=(7.4, 5.4)).subplot_mosaic([["a", "b"], ["c", "c"]], gridspec_kw={"hspace": 0.55})
fig, axes = panels["a"].figure, [panels["a"], panels["b"], panels["c"]]
strip(axes[0], "ifn_alpha", hollow="ifn_alpha_active")
chance(axes[0])
axes[0].set(title="a  Interferon-α score", ylabel="AUROC stim vs ctrl\nmean over 7 cell types")
axes[0].set_ylim(axes[0].get_ylim()[0], 1.0 + (1.0 - axes[0].get_ylim()[0]) * 0.04)
strip(axes[1], "active_latents", pairs=False)
axes[1].set(title="b  Active latents", ylabel="active latents out of 58")

ax = axes[2]
y, rows = 0, []
wide = scores.pivot(index="seed", columns="variant", values="ifn_alpha")
for other, label_ in [("shuffled", "fully random"), ("setshuffle", "within-set random"), ("vanilla", "no mask"),
                      ("coexpression", "co-expression")]:
    r = tests.loc[f"fixed latent: real - {other}"]
    primary = other == "shuffled"  # the comparison fixed before the 10-seed runs, in black
    color, size = (ink, 22) if primary else (muted, 16)
    diff = (wide["real"] - wide[other]).dropna()
    ax.scatter(diff, y + diff.index.map(offset), s=7, color=color, alpha=0.45, lw=0, zorder=2)
    ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=color, lw=0.9, zorder=3)
    ax.scatter(r["hodges_lehmann"], y, s=size, color=color, edgecolor=surface, lw=0.6, zorder=4)
    rows.append((y, f"Hallmark−{label_}"))
    y += 1
ax.axvline(0, color=baseline, lw=0.8, zorder=0)
ax.set_yticks(*zip(*rows))
ax.set_ylim(y - 0.5, -1.0)
ax.grid(axis="y", visible=False)
ax.set(title="c  Paired difference, Hodges–Lehmann estimate and 99.5% interval", xlabel="difference in AUROC")
fig.savefig(out / "scores.png")
plt.close(fig)

# 1b. The four measures that choose no latent, one dot per seed
measures = [("a  Held-out loss per cell", "diagnostics_fit.csv", "test_loss", "negative log-likelihood + KL (lower is better)"),
            ("b  Decoded response", "diagnostics_response.csv", "lfc_corr", "correlation across genes with the observed shift"),
            ("c  Cell-type probe", "diagnostics_celltype.csv", "balanced_acc_full", "balanced accuracy"),
            ("d  Stimulation probe", "diagnostics_probe.csv", "probe_auroc", "AUROC stim vs ctrl")]
fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.4), gridspec_kw={"hspace": 0.5, "wspace": 0.3})
for ax, (title, file, column, ylabel) in zip(axes.ravel(), measures):
    strip(ax, column, pairs=False, frame=pd.read_csv(tables / file))
    ax.set(title=title, ylabel=ylabel)
    ax.tick_params(axis="x", labelsize=6.5)
fig.savefig(out / "panel.png")
plt.close(fig)

# 2. Learning curves: mean and 95% t-interval over seeds at each training size, and the real - shuffled difference
if (tables / "learning_curves_scores.csv").exists():
    all_scores = pd.read_csv(tables / "learning_curves_scores.csv")
    full = scores["train_size"].max()
    sizes = sorted(all_scores["train_size"].unique())
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0), gridspec_kw={"wspace": 0.3})
    ax = axes[0]
    for v in [v for v in order if v != "setshuffle"]:  # the within-set random mask was run at full size only
        g = all_scores[all_scores["variant"] == v].groupby("train_size")["ifn_alpha"]
        mean, half = g.mean(), stats.t.ppf(0.975, g.count() - 1) * g.sem()
        ax.fill_between(mean.index, mean - half, mean + half, color=colors[v], alpha=0.12, lw=0)
        ax.plot(mean.index, mean, color=colors[v], lw=1, marker="o", ms=3, mec=surface, mew=0.5, label=names[v])
    ax.set(title="a  Interferon-α score", ylabel="AUROC stim vs ctrl\nmean over 7 cell types")
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.2), ncol=2, handletextpad=0.4, columnspacing=1.4)
    ax = axes[1]
    for n in sizes:  # one dot per seed (Hallmark - random of the same seed), the estimate and the corrected interval
        wide = all_scores[all_scores["train_size"] == n].pivot(index="seed", columns="variant", values="ifn_alpha")
        diff = (wide["real"] - wide["shuffled"]).dropna()
        r = tests.loc["fixed latent: real - shuffled" if n == full else f"{n} cells: real - shuffled"]
        ax.scatter(n * (1 + diff.index.map(offset) * 0.5), diff, s=7, color=ink, alpha=0.4, lw=0, zorder=2)
        ax.vlines(n, r["ci_low"], r["ci_high"], color=ink, lw=0.9, zorder=3)
        ax.scatter(n, r["hodges_lehmann"], s=16, color=ink, edgecolor=surface, lw=0.5, zorder=4)
    ax.axhline(0, color=baseline, lw=0.8, zorder=0)
    ax.set(title="b  Hallmark−random, corrected interval", ylabel="difference in AUROC")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xticks(sizes, [f"{n:,}" for n in sizes], fontsize=6, rotation=45, ha="right")
        ax.minorticks_off()
        ax.set_xlabel("training cells")
        ax.grid(axis="x", visible=False)
    fig.savefig(out / "learning_curves.png")
    plt.close(fig)

# 6. Pipeline: the phases of the analysis, one box per step
phases = {"Data": ("#e3eef7", "#5b9bd5"), "Design": ("#fbe6cf", "#e69f00"), "Model": ("#ecdcf0", "#a86fb5"),
          "Evaluation": ("#d6eeea", "#2a9d8f")}
steps = [("Data", "Labelled\nsingletons"), ("Data", "Cell QC"), ("Data", "Gene\nfilter"), ("Design", "Split"),
         ("Design", "Masks"), ("Model", "VAE"), ("Evaluation", "Four\nmeasures"), ("Evaluation", "Tests")]
bw, bh, gap, gap_phase = 0.82, 0.52, 0.17, 0.3
fig_w = 0.2 + len(steps) * bw + sum(gap_phase if steps[k][0] != steps[k + 1][0] else gap for k in range(len(steps) - 1))
fig = plt.figure(figsize=(fig_w, 1.0))
bg = fig.add_axes([0, 0, 1, 1])
bg.set(xlim=(0, fig_w), ylim=(0, 1.0))
bg.axis("off")
arrow = dict(arrowstyle="-|>", color=secondary, lw=0.9, mutation_scale=8, shrinkA=0, shrinkB=0)
x, y, span = 0.1, 0.08, {}
for k, (phase, name) in enumerate(steps):
    fill, edge = phases[phase]
    bg.add_patch(matplotlib.patches.FancyBboxPatch((x, y), bw, bh, boxstyle="round,pad=0,rounding_size=0.05", fc=fill, ec=edge,
                                                   lw=1.0))
    bg.text(x + bw / 2, y + bh / 2, name, ha="center", va="center", fontsize=7, color=ink, linespacing=1.1)
    lo, hi = span.get(phase, (x, x))
    span[phase] = (min(lo, x), x + bw)
    if k < len(steps) - 1:
        step = gap_phase if steps[k + 1][0] != phase else gap
        strong = phase == "Model"  # the evaluation reads the trained VAE itself
        bg.annotate("", (x + bw + step, y + bh / 2), (x + bw, y + bh / 2),
                    arrowprops=dict(arrow, lw=1.8, color=ink, mutation_scale=10) if strong else arrow)
        x += bw + step
for phase, (lo, hi) in span.items():  # a label and a line above each phase
    bg.plot([lo, hi], [y + bh + 0.1] * 2, color=phases[phase][1], lw=1.4, solid_capstyle="round")
    bg.text((lo + hi) / 2, y + bh + 0.2, phase, ha="center", va="center", fontsize=7.5, fontweight="bold", color=ink)
fig.savefig(out / "pipeline.png")
plt.close(fig)

print("saved", *sorted(p.name for p in out.glob("*.png")))
