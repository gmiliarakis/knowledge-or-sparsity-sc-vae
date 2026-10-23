# Figures from the evaluation tables (src/08_evaluate.py) and the training logs, written to reports/figures/. The
# learning-curve figure is drawn once its tables exist (src/08_evaluate.py --learning-curves). The comparison figures need src/10_checks.py,
# and the tests in src/11_tests.py supply the intervals.
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

from gene_sets import collection, data, main, run

matplotlib.use("Agg")
import matplotlib.pyplot as plt

surface, ink, secondary, muted, grid, baseline = "#ffffff", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
colors = {"vanilla": "#cfcfcf", "real": "#0072B2", "shuffled": "#882255", "setshuffle": "#56B4E9", "coexpression": "#009E73"}  # Okabe-Ito blue, bluish green and sky blue with wine: distinct under red-green colour blindness
names = {"vanilla": "vanilla", "real": "Hallmark", "shuffled": "fully random", "setshuffle": "within-set random", "coexpression": "co-expression"}
ticks = {"vanilla": "vanilla", "real": "Hallmark", "shuffled": "fully\nrandom", "setshuffle": "within-set\nrandom", "coexpression": "co-\nexpression"}
folders = {"vanilla": "plain", "real": "real", "shuffled": "shuffled", "setshuffle": "setshuffle", "coexpression": "coexpression"}
order = list(colors)
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "axes.titleweight": "medium", "axes.titlelocation": "left",
    "axes.titlepad": 10, "axes.labelsize": 8, "axes.labelweight": "bold", "text.color": ink, "axes.labelcolor": "#1f1f1f",
    "xtick.color": baseline, "ytick.color": baseline, "xtick.labelcolor": "#1f1f1f", "ytick.labelcolor": "#1f1f1f",
    "xtick.major.size": 0, "ytick.major.size": 0, "xtick.major.pad": 4, "ytick.major.pad": 4,
    "axes.edgecolor": baseline, "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": grid, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "figure.facecolor": surface, "axes.facecolor": surface, "savefig.facecolor": surface,
    "legend.frameon": False, "legend.fontsize": 8, "lines.solid_capstyle": "round",
    "savefig.dpi": 400, "savefig.bbox": "tight", "pdf.fonttype": 42,
})

SEED_S, MEAN_S = 9, 15  # point area (pt²) of a seed and of the black mean: the same in every figure


def seed_dots(ax, x, y, color, filled=True, zorder=3):
    # a seed: one opaque dot in the colour of its variant (an open ring when `filled` is False, an inactive latent)
    if filled:
        ax.scatter(x, y, s=SEED_S, color=color, lw=0, zorder=zorder)
    else:
        ax.scatter(x, y, s=SEED_S, facecolor="none", edgecolor=color, lw=0.8, zorder=zorder)


def mean_dot(ax, x, y, zorder=5):
    ax.scatter(x, y, s=MEAN_S, color=ink, lw=0, zorder=zorder)


mask_name = "Reactome" if collection == "reactome" else "Hallmark"
if True:  # one figure per benchmark (SC_DATA, SC_SETS); the main experiment (kang/hallmark) continues below, every other benchmark stops here
    # the four measures that choose no latent, one dot per seed
    ticks["real"] = mask_name
    bench_label = {"kang": "Kang et al. IFN-β", "oesinghaus_ifn_beta": "Theis et al. IFN-β", "oesinghaus_il_2": "Theis et al. IL-2"}[data]
    tables, out = Path(f"reports/tables/{run}"), Path("reports/figures")
    out.mkdir(parents=True, exist_ok=True)
    offset = dict(zip(range(10), np.linspace(-0.18, 0.18, 10)))  # each seed keeps its position in every panel


    def strip(ax, frame, column):
        # one dot per seed and the mean as a short bar
        for i, v in enumerate(order):
            runs = frame[frame["variant"] == v]
            if runs.empty:
                continue
            xs = i + runs["seed"].map(offset)
            seed_dots(ax, xs, runs[column], colors[v])
            mean_dot(ax, i, runs[column].mean())
        ax.set_xticks(range(len(order)), [ticks[v] for v in order])
        ax.set_xlim(-0.6, len(order) - 0.4)
        ax.grid(axis="x", visible=False)
        ax.tick_params(axis="x", labelsize=6.5)


    measures = [("a  Test loss per cell", "diagnostics_fit.csv", "test_loss", "negative log-likelihood + KL\n(lower is better)"),
                ("b  Decoded response", "diagnostics_response.csv", "lfc_corr", "correlation across genes with\nthe observed shift"),
                ("c  Cell-type classifier", "diagnostics_celltype.csv", "balanced_acc_full", "balanced accuracy"),
                ("d  Stimulation classifier", "diagnostics_probe.csv", "probe_auroc", "AUROC stim vs ctrl")]
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.6), gridspec_kw={"hspace": 0.55, "wspace": 0.35})
    for ax, (title, file, column, ylabel) in zip(axes.ravel(), measures):
        strip(ax, pd.read_csv(tables / file), column)
        ax.set(title=title, ylabel=ylabel)
    fig.suptitle(f"{bench_label}, {mask_name} sets", x=0.02, ha="left", fontsize=9)
    fname = f"{data.replace('oesinghaus', 'theis')}_{collection}.png"  # figure files carry the author name used in the text
    fig.savefig(out / fname)
    plt.close(fig)
    print("saved", out / fname)
    if not main:
        raise SystemExit

tables, out = Path("reports/tables/kang/hallmark"), Path("reports/figures")
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
        seed_dots(ax, xs[~empty], runs.loc[~empty, column], colors[v])
        seed_dots(ax, xs[empty], runs.loc[empty, column], colors[v], filled=False)
        mean_dot(ax, x[v], runs[column].mean())
    ax.set_xticks(range(len(order)), [ticks[v] for v in order])
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.grid(axis="x", visible=False)


def chance(ax):
    # AUROC 0.5 only when the data come near it, so close values are not squashed
    low, high = ax.get_ylim()
    if low < 0.55:
        ax.axhline(0.5, color=muted, lw=0.8, zorder=0)
        ax.set_ylim(min(low, 0.48), high)


# 0. Overview banner: data, VAE with and without a mask, the six decoder masks as trained with seed 0, and the
# decoded response of every seed
ov_colors = {**colors, "responsematch": "#999933"}  # olive of the Tol muted palette, for the expression-matched random mask of this figure only
ov_names = {**names, "responsematch": "expression-matched random"}
W, H = 9.1, 2.7
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


for x, head in [(0.7, "1  Data"), (2.3, "2  Model"), (4.53, "3  Decoder masks"), (7.75, "4  Decoded response")]:
    label(x, 2.53, head, size=8.5, color=ink, weight="bold")  # each title centred over its panel

# 1 data: two clouds of cells
for (cy, c, text) in [(1.8, "#6f6e68", "control"), (0.9, "#CC79A7", "IFN-β, 6 h")]:
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
rm = np.load("data/processed/masks_hallmark_responsematch.npz")  # the expression-matched random masks (src/04b_build_responsematch.py)
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
         "setshuffle": full(sets[ss["perms"][0]]), "responsematch": full(sets[rm["perms"][0]]), "vanilla": np.ones((sets.shape[0], 58))}
mask_names = {"real": "Hallmark", "coexpression": "co-expression", "shuffled": "fully\nrandom", "setshuffle": "within-set\nrandom",
              "responsematch": "expression-\nmatched random", "vanilla": "vanilla"}
cells = {"real": (0, 0), "shuffled": (1, 0), "setshuffle": (2, 0), "responsematch": (0, 1), "coexpression": (1, 1), "vanilla": (2, 1)}  # column and row, in the order of the tables with vanilla last
for v in ["real", "shuffled", "setshuffle", "responsematch", "coexpression", "vanilla"]:
    x, y = 3.52 + cells[v][0] * 0.72, 1.5 - cells[v][1] * 1.1  # the six masks in a 3 x 2 block
    ax = inset(x, y, 0.58, 0.85)
    ax.set(xlim=(-0.5, 57.5), ylim=(len(rows) - 0.5, -0.5), xticks=[], yticks=[])
    if v == "vanilla":
        ax.add_patch(plt.Rectangle((-0.5, -0.5), 58, len(rows), color=ov_colors[v], lw=0))  # every gene reaches every latent
    else:
        genes_, latents_ = np.nonzero(masks[v])
        ax.scatter(latents_, position[genes_], s=0.5, marker="s", color=ov_colors[v], lw=0)
    ax.grid(False)
    for side in ax.spines.values():
        side.set_visible(False)
    label(x + 0.29, y - 0.14, mask_names[v], color=ink, size=5.2)
arrow(5.72, 5.92)

# 4 result: decoded response (correlation of the decoded and observed stimulation shift across genes) of every seed as a
# beeswarm, the mean printed under each bar
response = pd.read_csv(tables / "diagnostics_response.csv")
variants = ["real", "shuffled", "setshuffle", "responsematch", "coexpression", "vanilla"]  # same order as the masks in stage 3
width_in, height_in, low, high = 2.3, 1.9, 0.66, 1.0
ax = inset(6.6, 0.45, width_in, height_in)
dot_in = 0.04
for i, v in enumerate(variants):
    y_ = response[response["variant"] == v]["lfc_corr"].to_numpy()
    x_ = i + swarm(y_, dot_in / (width_in / len(variants)), dot_in / (height_in / (high - low)))
    seed_dots(ax, x_, y_, ov_colors[v])
    mean_dot(ax, i, y_.mean())
ax.set(xlim=(-0.5, len(variants) - 0.5), ylim=(low, high), yticks=[0.7, 0.8, 0.9, 1.0], yticklabels=["0.7", "0.8", "0.9", "1.0"])
ax.set_ylabel("correlation")
ax.yaxis.set_label_coords(-0.12, 0.5)
ax.tick_params(axis="y", labelsize=6.5)
for side in ("left", "bottom"):
    ax.spines[side].set_color(ink)
ax.set_xticks(range(len(variants)), [ov_names[v] for v in variants], fontsize=5.5, rotation=45, ha="right", rotation_mode="anchor")
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
for other, label_ in [("shuffled", "fully random"), ("setshuffle", "within-set random"), ("vanilla", "vanilla"),
                      ("coexpression", "co-expression")]:
    r = tests.loc[f"fixed latent: real - {other}"]
    primary = other == "shuffled"  # the comparison fixed before the 10-seed runs, in black
    color, size = (ink, 22) if primary else (muted, 16)
    diff = (wide["real"] - wide[other]).dropna()
    seed_dots(ax, diff, y + diff.index.map(offset), color, zorder=2)
    ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=color, lw=0.9, zorder=3)
    ax.scatter(r["hodges_lehmann"], y, s=size, color=color, edgecolor=surface, lw=0.6, zorder=4)
    rows.append((y, f"Hallmark−{label_}"))
    y += 1
ax.axvline(0, color=baseline, lw=0.8, zorder=0)
ax.set_yticks(*zip(*rows))
ax.set_ylim(y - 0.5, -1.0)
ax.grid(axis="y", visible=False)
ax.set(title="c  Paired difference, median and 99.5% interval", xlabel="difference in AUROC")
fig.savefig(out / "scores.png")
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
        seed_dots(ax, n * (1 + diff.index.map(offset) * 0.5), diff, ink, zorder=2)
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

# 6. Pipeline: the phases of the analysis, one box per step; what varies between benchmarks under the step that varies it,
# and the two extra analyses as dashed boxes under the step they change
phases = {"Data": ("#e3eef7", "#5b9bd5"), "Design": ("#f1dfe8", "#882255"), "Model": ("#ecdcf0", "#a86fb5"),
          "Evaluation": ("#d6eeea", "#2a9d8f")}
steps = [("Data", "Single-cell\ndata"), ("Data", "Cell QC"), ("Data", "Gene\nfilter"), ("Design", "Split"),
         ("Design", "Gene sets"), ("Design", "Six\nmasks"), ("Model", "VAE"), ("Evaluation", "Measures"), ("Evaluation", "Tests")]
stacks = {0: [("Kang et al. IFN-β", 0), ("Theis et al.", 0), ("IFN-β", 1), ("IL-2", 1)], 4: [("Hallmark", 0), ("Reactome", 0)],
          5: [("gene sets", 0), ("co-expression", 0), ("random", 0), ("fully random", 1), ("within-set random", 1), ("expression-\nmatched", 1), ("vanilla", 0)],
          2: [("Without the 250\nstrongest genes", 0)],
          3: [("Learning curves", 0), ("500 cells", 1), ("1,000 cells", 1), ("2,000 cells", 1), ("4,000 cells", 1),
              ("8,000 cells", 1), ("all 16,749", 1)],
          7: [("test loss", 0), ("decoded response", 0), ("classifiers", 0), ("cell-type", 1), ("stimulation", 1),
              ("labelled-latent and\nreference scores", 0)]}  # (label, level): level 1 hangs under the box above it
bw, bh, gap, gap_phase = 1.0, 0.52, 0.17, 0.3
fig_w = 0.2 + len(steps) * bw + sum(gap_phase if steps[k][0] != steps[k + 1][0] else gap for k in range(len(steps) - 1))
fig_h = 2.95
fig = plt.figure(figsize=(fig_w, fig_h))
bg = fig.add_axes([0, 0, 1, 1])
bg.set(xlim=(0, fig_w), ylim=(0, fig_h))
bg.axis("off")
arrow = dict(arrowstyle="-|>", color=secondary, lw=0.9, mutation_scale=8, shrinkA=0, shrinkB=0)
x, y, span, centre = 0.1, 2.05, {}, {}
for k, (phase, name) in enumerate(steps):
    fill, edge = phases[phase]
    bg.add_patch(matplotlib.patches.FancyBboxPatch((x, y), bw, bh, boxstyle="round,pad=0,rounding_size=0.05", fc=fill, ec=edge,
                                                   lw=1.0))
    bg.text(x + bw / 2, y + bh / 2, name, ha="center", va="center", fontsize=7, color=ink, linespacing=1.1)
    centre[k] = x + bw / 2
    lo, hi = span.get(phase, (x, x))
    span[phase] = (min(lo, x), x + bw)
    if k in stacks:  # small boxes under the step, as a tree: a spine on the left, children indented under their parent
        h, pad, left0, indent, right = 0.2, 0.06, x + 0.12, 0.17, x + bw - 0.03
        yt, boxes = y - 0.1, []
        for item, level in stacks[k]:
            box_h = h * 1.7 if "\n" in item else h
            yt -= box_h
            boxes.append((item, level, left0 + level * indent, yt, box_h))
            yt -= pad
        top_level = [b for b in boxes if b[1] == 0]
        bg.plot([x + 0.04] * 2, [y, top_level[-1][3] + top_level[-1][4] / 2], color=muted, lw=0.8, zorder=0)
        for n, (item, level, left, yb, box_h) in enumerate(boxes):
            bg.add_patch(matplotlib.patches.FancyBboxPatch((left, yb), right - left, box_h,
                                                           boxstyle="round,pad=0,rounding_size=0.04",
                                                           fc="#f3f2ee" if k in (2, 3) else surface,
                                                           ec=muted if k in (2, 3) else phases[phase][1], lw=0.8, zorder=2))
            bg.text((left + right) / 2, yb + box_h / 2, item, ha="center", va="center", fontsize=5.4 if level else 5.8,
                    color=ink, linespacing=1.1, zorder=3)
            if level == 0:  # tick from the spine
                bg.plot([x + 0.04, left], [yb + box_h / 2] * 2, color=muted, lw=0.8, zorder=0)
            else:  # elbow from the parent box (the nearest level-0 box above) to this child
                parent = next(b for b in reversed(boxes[:n]) if b[1] == 0)
                elbow = parent[2] + 0.08
                bg.plot([elbow, elbow, left], [parent[3], yb + box_h / 2, yb + box_h / 2], color=muted, lw=0.8, zorder=0)
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

# Stages of the pipeline, one figure per part of the Methods (the paper places each beside its subsection)
def stage_figure(first, last, name, width=6.3, scale=1.0):
    bw_, bh_, gap_, gp_ = 1.55, 0.6, 0.35, 0.35
    ks = list(range(first, last + 1))
    fw = 0.2 + len(ks) * bw_ + sum(gp_ if steps[k][0] != steps[k + 1][0] else gap_ for k in ks[:-1])
    tall = max(len(stacks.get(k, [])) for k in ks)
    def stack_h(k):  # the height of the small boxes under a step, as drawn below
        return sum((0.26 * 1.7 if "\n" in i else 0.26) + 0.07 for i, _ in stacks.get(k, []))
    fh = 0.95 + 0.12 + max(stack_h(k) for k in ks) + 0.1
    f = plt.figure(figsize=(fw, fh))
    ax = f.add_axes([0, 0, 1, 1]); ax.set(xlim=(0, fw), ylim=(0, fh)); ax.axis("off")
    x, y = 0.1, fh - 0.95
    first_x, spanx = x, {}
    for k in ks:
        phase, nm = steps[k]
        fill, edge = phases[phase]
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x, y), bw_, bh_, boxstyle="round,pad=0,rounding_size=0.06", fc=fill, ec=edge, lw=1.1))
        ax.text(x + bw_ / 2, y + bh_ / 2, nm, ha="center", va="center", fontsize=9, color=ink, linespacing=1.1)
        lo, hi = spanx.get(phase, (x, x)); spanx[phase] = (min(lo, x), x + bw_)
        if k in stacks:
            h, pad, left0, indent, right = 0.26, 0.07, x + 0.16, 0.22, x + bw_ - 0.04
            yt, boxes = y - 0.12, []
            for item, level in stacks[k]:
                bh2 = h * 1.7 if "\n" in item else h
                yt -= bh2; boxes.append((item, level, left0 + level * indent, yt, bh2)); yt -= pad
            top_level = [b for b in boxes if b[1] == 0]
            ax.plot([x + 0.05] * 2, [y, top_level[-1][3] + top_level[-1][4] / 2], color=muted, lw=0.9, zorder=0)
            for n, (item, level, left, yb, bh2) in enumerate(boxes):
                ax.add_patch(matplotlib.patches.FancyBboxPatch((left, yb), right - left, bh2, boxstyle="round,pad=0,rounding_size=0.05", fc="#f3f2ee" if k in (2, 3) else surface, ec=muted if k in (2, 3) else phases[phase][1], lw=0.9, zorder=2))
                ax.text((left + right) / 2, yb + bh2 / 2, item, ha="center", va="center", fontsize=7.5 if level else 8, color=ink, linespacing=1.1, zorder=3)
                if level == 0:
                    ax.plot([x + 0.05, left], [yb + bh2 / 2] * 2, color=muted, lw=0.9, zorder=0)
                else:
                    parent = next(b for b in reversed(boxes[:n]) if b[1] == 0)
                    el = parent[2] + 0.1
                    ax.plot([el, el, left], [parent[3], yb + bh2 / 2, yb + bh2 / 2], color=muted, lw=0.9, zorder=0)
        if k != ks[-1]:
            st = gp_ if steps[k + 1][0] != phase else gap_
            ax.annotate("", (x + bw_ + st, y + bh_ / 2), (x + bw_, y + bh_ / 2), arrowprops=dict(arrow, lw=1.8, color=ink, mutation_scale=10) if phase == "Model" else arrow)
            x += bw_ + st
    for ph, (lo, hi) in spanx.items():
        ax.plot([lo, hi], [y + bh_ + 0.12] * 2, color=phases[ph][1], lw=1.5, solid_capstyle="round")
        ax.text((lo + hi) / 2, y + bh_ + 0.25, ph, ha="center", va="center", fontsize=9.5, fontweight="bold", color=ink)
    f.savefig(out / name, dpi=300); plt.close(f)


stage_figure(0, 2, "stage_data.png")
stage_figure(3, 5, "stage_design.png")
stage_figure(6, 8, "stage_evaluation.png")


# 7. Overarching figures: the six benchmarks side by side (Q1c to Q1f). One line weight and one point size for every line plot.
import matplotlib.transforms as mtransforms

LW, PT = 1.0, 3  # line weight and point size
grey = "#9a9a9a"
vcolor = {**colors, "vanilla": grey}
vnames = {**names, "real": "gene sets"}
bench_rows, bench_groups, y_ = [], [], 0.0
for d_, label_ in [("kang", "Kang et al. IFN-β"), ("oesinghaus_ifn_beta", "Theis et al. IFN-β"), ("oesinghaus_il_2", "Theis et al. IL-2")]:
    first = y_
    for c_ in ("hallmark", "reactome"):
        bench_rows.append((f"{d_}/{c_}", c_.capitalize(), y_))
        y_ += 1
    bench_groups.append((label_, first, y_ - 1))
    y_ += 0.6


def group_names(ax, groups, inches=0.7):
    # the shared part of the row labels (the contrast) once on the left, instead of repeated on every row
    shift = mtransforms.offset_copy(ax.get_yaxis_transform(), fig=ax.figure, x=-inches, units="inches")
    for text, ya, yb in groups:
        ax.text(0, (ya + yb) / 2, text, ha="right", va="center", fontsize=8, fontweight="bold", color=ink, transform=shift,
                linespacing=1.1)


def advantage_figure(measures, name, size):
    # advantage of the gene sets over each other mask, with the interval of the paired test
    others = [(v, names[v], vcolor[v]) for v in ("shuffled", "setshuffle", "coexpression")]  # each line in the colour of the mask it is compared with
    rows, groups, ya = [], [], 0.0
    for d_, label_ in [("kang", "Kang et al. IFN-β"), ("oesinghaus_ifn_beta", "Theis et al. IFN-β"), ("oesinghaus_il_2", "Theis et al. IL-2")]:
        first = ya
        for c_ in ("hallmark", "reactome"):
            rows.append((f"{d_}/{c_}", c_.capitalize(), ya))
            ya += 1
        groups.append((label_, first, ya - 1))
        ya += 0.6
    with plt.rc_context({"lines.linewidth": LW, "lines.markersize": PT, "axes.grid": False, "axes.spines.left": False}):
        fig, axs = plt.subplots(1, len(measures), figsize=size, sharey=True, gridspec_kw={"wspace": 0.12})
        for ax, (key, title, sign) in zip(axs, measures):
            ax.axvline(0, color=ink, lw=0.6)
            for path, _, y in rows:
                t = pd.read_csv(f"reports/tables/{path}/tests.csv")
                for j, (o, _, c) in enumerate(others):
                    r = t[t["test"] == f"{key}: real - {o}"].iloc[0]
                    yy = y + (j - 1) * 0.26
                    lo, hi = sorted([sign * r["ci_low"], sign * r["ci_high"]])
                    ax.plot([lo, hi], [yy, yy], color=c, zorder=2)
                    if r["reject"]:  # filled: adjusted p below 0.05 and the interval excludes zero; open: the interval reaches zero
                        excludes = lo > 0 or hi < 0
                        ax.plot(sign * r["hodges_lehmann"], yy, "o", color=ink, markerfacecolor=ink if excludes else surface, zorder=3)
            ax.set_title(title, fontsize=8, loc="left")
            ax.tick_params(axis="y", length=0)
            ax.grid(axis="x", color=grid, lw=0.6)
        axs[0].set_yticks([y for _, _, y in rows], [l for _, l, _ in rows], fontsize=7.5)
        axs[0].set_ylim(rows[-1][2] + 0.7, -0.7)
        group_names(axs[0], groups)
        for _, label, c in others:
            axs[0].plot([], [], color=c, label=f"gene sets minus {label}")
        fig.legend(*axs[0].get_legend_handles_labels(), frameon=False, fontsize=7.5, ncol=3, loc="upper center", bbox_to_anchor=(0.55, -0.02))
        fig.savefig(out / name)
        plt.close(fig)


def share_scores(path):
    # decoded response of each variant, as the position between the fully random mask (0) and vanilla VAE (100)
    w = pd.read_csv(f"reports/tables/{path}/diagnostics_response.csv").pivot(index="seed", columns="variant", values="lfc_corr")
    zero = w["shuffled"].mean()
    return {v: (w[v] - zero) / (w["vanilla"].mean() - zero) * 100 for v in ("real", "coexpression", "setshuffle")}


def share_gap(path):
    # the denominator of the relative scale: decoded response of the vanilla VAE minus that of the fully random mask
    w = pd.read_csv(f"reports/tables/{path}/diagnostics_response.csv").pivot(index="seed", columns="variant", values="lfc_corr")
    return w["vanilla"].mean() - w["shuffled"].mean()


def share_panel(ax, rows, groups, title):
    # one panel: each variant placed between the fully random mask (0) and vanilla VAE (100)
    top = rows[0][2] - 0.8
    ax.axvline(0, color=colors["shuffled"], lw=0.6, zorder=1)
    ax.axvline(100, color=grey, lw=0.6, zorder=1)
    for path, _, y in rows:
        s = share_scores(path)
        for j, v in enumerate(("real", "coexpression", "setshuffle")):
            yy = y + (j - 1) * 0.28
            ax.plot([s[v].min(), s[v].max()], [yy, yy], color=vcolor[v], zorder=2)
            ax.plot(s[v].mean(), yy, "o", color=vcolor[v], zorder=3)
    ax.text(1.5, top, "fully random mask", va="center", ha="left", fontsize=7, color=colors["shuffled"])
    ax.text(101.5, top, "vanilla", va="center", ha="left", fontsize=7, color=grey)
    ax.set_yticks([y for _, _, y in rows], [l for _, l, _ in rows], fontsize=7.5)
    # the difference of each row (vanilla minus fully random decoded response, the scale of the row) in a column between the row labels and the plot
    diff_at = mtransforms.offset_copy(ax.get_yaxis_transform(), fig=ax.figure, x=-0.08, units="inches")
    ax.text(0, top, "difference", ha="right", va="center", fontsize=7.5, color=ink, transform=diff_at)
    for p, _, y in rows:
        ax.text(0, y, f"{share_gap(p):.2f}", ha="right", va="center", fontsize=7.5, color=ink, transform=diff_at)
    ax.tick_params(axis="y", length=0, pad=27)  # room for the difference column
    ax.set_xticks([-50, -25, 0, 25, 50, 75, 100])
    ax.set_xlim(-70, 125)
    ax.set_ylim(max(y for _, _, y in rows) + 0.7, top - 0.5)
    ax.grid(axis="x", color=grid, lw=0.6)
    ax.set_title(title, fontsize=8, loc="left")
    if groups:
        group_names(ax, groups, inches=0.95)


def share_figure(name):
    # a: the six benchmarks (Q1d to Q1f); b: the main experiment with and without the 250 strongest genes (Q1c)
    with plt.rc_context({"lines.linewidth": LW, "lines.markersize": PT, "axes.grid": False, "axes.spines.left": False}):
        fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(6.4, 5.0), gridspec_kw={"height_ratios": [3.4, 1.6], "hspace": 0.35})
        share_panel(ax_a, bench_rows, bench_groups, "a  Six combinations")
        share_panel(ax_b, [("kang/hallmark", "all 12,034 genes", 0.0), ("kang_drop250/hallmark", "250 strongest\nresponse genes\ndropped", 1.2)], [],
                    "b  Kang et al. data, Hallmark sets, without the 250 strongest genes")
        ax_b.set_xlabel("relative decoded response (%)")
        for v in ("real", "coexpression", "setshuffle"):
            ax_a.plot([], [], color=vcolor[v], label=vnames[v])
        fig.legend(*ax_a.get_legend_handles_labels(), frameon=False, fontsize=7.5, ncol=3, loc="upper center", bbox_to_anchor=(0.55, -0.005))
        fig.savefig(out / name)
        plt.close(fig)


def levels_figures():
    # decoded response of every seed (small point) and its mean (black point), for every variant
    order_ = ["real", "coexpression", "setshuffle", "shuffled", "vanilla"]
    offsets = np.linspace(-0.17, 0.17, 10)
    with plt.rc_context({"lines.linewidth": LW, "axes.grid": False}):
        fig, ax = plt.subplots(figsize=(7.4, 3.4))
        for gi, (path, label, _) in enumerate(bench_rows):
            d = pd.read_csv(f"reports/tables/{path}/diagnostics_response.csv")
            for k, v in enumerate(order_):
                vals = d[d["variant"] == v]["lfc_corr"]
                ax.plot([gi * 6 + k] * 2, [vals.min(), vals.max()], color=vcolor[v], lw=LW)
                ax.plot(gi * 6 + k, vals.mean(), "o", color=vcolor[v], markersize=PT)
        ax.set_xticks([gi * 6 + 2 for gi in range(6)], [l for _, l, _ in bench_rows], fontsize=7.5)
        ax.set_xlim(-1, 35)
        ax.set_ylim(0.45, 1.0)
        ax.set_ylabel("decoded response (correlation)")
        ax.grid(axis="y", color=grid, lw=0.6)
        belows = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        for ci, (_, label, _) in enumerate(zip(range(3), [g[0] for g in bench_groups], range(3))):
            ax.text(ci * 12 + 5, -0.16, label, ha="center", va="top", fontsize=8, fontweight="bold", transform=belows)
        for v in order_:
            ax.plot([], [], "o", color=vcolor[v], markersize=PT, label=vnames[v])
        ax.legend(frameon=False, fontsize=7.5, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.3))
        fig.savefig(out / "levels_benchmarks.png")
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(4.6, 3.0))
        for gi, path in enumerate(["kang/hallmark", "kang_drop250/hallmark"]):
            d = pd.read_csv(f"reports/tables/{path}/diagnostics_response.csv")
            for k, v in enumerate(order_):
                vals = d[d["variant"] == v].sort_values("seed")["lfc_corr"].to_numpy()
                seed_dots(ax, gi * 6 + k + offsets, vals, vcolor[v])
                mean_dot(ax, gi * 6 + k, vals.mean())
        ax.set_xticks([2, 8], ["all 12,034 genes", "250 strongest\nresponse genes dropped"], fontsize=7.5)
        ax.set_xlim(-0.8, 10.8)
        ax.set_ylim(0.45, 1.0)
        ax.set_ylabel("decoded response (correlation)")
        ax.grid(axis="y", color=grid, lw=0.6)
        for v in order_:
            ax.plot([], [], "o", color=vcolor[v], markersize=PT, label=vnames[v])
        ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.2), columnspacing=1, handletextpad=0.2)
        fig.savefig(out / "levels_Q1c.png")
        plt.close(fig)


def named_latent_figure():
    # score of the latent at the position of the expected-response set (real and random masks) or of the highest-scoring latent (co-expression, vanilla)
    order_ = ["real", "coexpression", "setshuffle", "shuffled", "vanilla"]
    with plt.rc_context({"lines.linewidth": LW, "lines.markersize": PT, "axes.grid": False}):
        fig, ax = plt.subplots(figsize=(7.4, 3.4))
        for gi, (path, _, _) in enumerate(bench_rows):
            d = pd.read_csv(f"reports/tables/{path}/scores.csv")
            d = d[d["train_size"] == d["train_size"].max()]
            for k, v in enumerate(order_):
                vals = d[d["variant"] == v]["ifn_alpha"]
                ax.plot([gi * 6 + k] * 2, [vals.min(), vals.max()], color=vcolor[v], lw=LW)
                ax.plot(gi * 6 + k, vals.mean(), "o", color=vcolor[v], markersize=PT)
        ax.set_xticks([gi * 6 + 2 for gi in range(6)], [l for _, l, _ in bench_rows], fontsize=7.5)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlim(-1, 35)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel("score of the latent (AUROC)")
        ax.grid(axis="y", color=grid, lw=0.6)
        belows = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        for ci, (label, _, _) in enumerate(bench_groups):
            ax.text(ci * 12 + 5, -0.16, label, ha="center", va="top", fontsize=8, fontweight="bold", transform=belows)
        for v in order_:
            ax.plot([], [], "o", color=vcolor[v], markersize=PT, label=vnames[v])
        ax.legend(frameon=False, fontsize=7.5, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.3))
        fig.savefig(out / "named_latent.png")
        plt.close(fig)


MEASURES = {"loss": ("test loss per cell (lower is better)", "test loss per cell\n(loss units)", -1),
            "response": ("decoded response correlation", "decoded response\n(correlation)", 1),
            "celltype": ("cell-type balanced accuracy", "cell-type classifier\n(balanced accuracy)", 1),
            "probe": ("stimulation probe AUROC", "stimulation classifier\n(AUROC)", 1)}
advantage_figure([MEASURES["loss"], MEASURES["response"]], "advantage.png", (6.6, 3.9))
advantage_figure([MEASURES["celltype"], MEASURES["probe"]], "advantage_probes.png", (6.6, 3.9))
share_figure("share.png")
levels_figures()
named_latent_figure()

print("saved", *sorted(p.name for p in out.glob("*.png")))
