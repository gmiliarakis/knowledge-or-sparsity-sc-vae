# Schema of the VAE (Methods, Model). Written to reports/figures/model.png and paper/figures/model.png, or to the directories
# given as arguments. Standalone: it needs no tables. Each matrix is a solid: the height suggests its rows (the 12,034-row
# plates are drawn shorter than to scale) and the depth its columns or cells (on a log scale).
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mp
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

root = Path(__file__).resolve().parent.parent
outs = [Path(a) for a in sys.argv[1:]] or [root / "reports" / "figures", root / "paper" / "figures"]

ink, secondary, flow = "#0b0b0b", "#6b6a66", "#4d4d4b"  # text, brackets, arrows between layers
C = dict(data=("#e3eef7", "#5b9bd5"), donor=("#f2f7fc", "#7aa6d0"), offsets=("#c4dbef", "#2f6fa8"), linear=("#dfe0f5", "#5b5fc7"),
         relu=("#d3eeea", "#2a9d8f"), dropout=("#e4e4e0", "#8a8983"), mean=("#dcefdc", "#4c9a4c"), logvar=("#fde7c8", "#d98e04"),
         latent=("#dde6ee", "#5f7a96"), free=("#f7f9fb", "#5f7a96"), mask=("#f1dfe8", "#882255"), weights=("#fafaf8", "#6b6a66"),
         softmax=("#ecdcf0", "#a86fb5"))  # as in the pipeline figures where the parts coincide
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
                     "text.color": ink, "savefig.dpi": 300, "mathtext.default": "it"})

W, H = 17.45, 8.95
fig = plt.figure(figsize=(W * 0.52, H * 0.52))
ax = fig.add_axes([0, 0, 1, 1]); ax.set(xlim=(0, W), ylim=(0, H)); ax.axis("off")

FACE, UNIT, PLATE = 0.12, 0.01, 2.4  # width of a front face, height per unit of a small layer, height of a 12,034-row plate
CELLS, LATENTS, DONORS = (0.15 * math.log10(n) for n in (16749, 58, 8))  # depth of a solid with that many columns
cy = 5.8  # centre line of the encoder


def mix(c, t, to=(1, 1, 1)):
    return tuple(a + (b - a) * t for a, b in zip(to_rgb(c), to))


def block(x, ymid, h, colour, depth, top=True, shade=0.3):
    """A solid with sharp edges: front face in the fill colour, lighter top, darker side."""
    fill, edge = colour
    y0, dx, dy = ymid - h / 2, depth * 1.15, depth * 0.45
    if top:
        ax.add_patch(mp.Polygon([(x, y0 + h), (x + dx, y0 + h + dy), (x + FACE + dx, y0 + h + dy), (x + FACE, y0 + h)], fc=mix(fill, 0.55), ec=edge, lw=0.8, zorder=3))
    ax.add_patch(mp.Polygon([(x + FACE, y0), (x + FACE + dx, y0 + dy), (x + FACE + dx, y0 + h + dy), (x + FACE, y0 + h)], fc=mix(fill, shade, to=to_rgb(edge)), ec=edge, lw=0.8, zorder=3))
    ax.add_patch(mp.Rectangle((x, y0), FACE, h, fc=fill, ec=edge, lw=0.8, zorder=3.1))
    return dict(x=x, r=x + FACE + dx, m=ymid, b=y0, t=y0 + h, dx=dx, dy=dy)


def text(x, y, s, ha="center", va="center", fs=8, color=ink, **kw):
    return ax.text(x, y, s, ha=ha, va=va, fontsize=fs, color=color, **kw)


def mid(b):  # the middle of the horizontal projection of a solid
    return b["x"] + (FACE + b["dx"]) / 2


def above(b, s, fs=8):
    text(mid(b), b["t"] + b["dy"] + 0.12, s, va="bottom", fs=fs)


def below(b, s):
    text(mid(b), b["b"] - 0.1, s, va="top")


def beside(b, s, gap=0.1):  # next to the long edge, at mid-height
    text(b["x"] - gap, b["m"], s, ha="right")


def on_side(b, s):  # on the side face
    text(b["x"] + FACE + b["dx"] / 2, b["m"] + b["dy"] / 2, s)


def arrow(p0, p1, colour=flow, lw=0.7):
    ax.annotate("", p1, p0, arrowprops=dict(arrowstyle="-|>", color=colour, lw=lw, mutation_scale=6, shrinkA=0, shrinkB=0), zorder=1)


def line(xs, ys, colour=flow, lw=0.7, **kw):
    ax.plot(xs, ys, color=colour, lw=lw, zorder=1, **kw)


def equation(xc, y, parts, fs=9.5, gap=0.07):
    """One equation centred on xc, each term in its own colour."""
    r, inv = fig.canvas.get_renderer(), ax.transData.inverted()
    objs = [text(0, y, s, ha="left", va="top", fs=fs, color=c) for s, c in parts]
    ws = [inv.transform((o.get_window_extent(r).x1, 0))[0] - inv.transform((o.get_window_extent(r).x0, 0))[0] for o in objs]
    x = xc - (sum(ws) + gap * (len(parts) - 1)) / 2
    for o, w in zip(objs, ws):
        o.set_x(x)
        x += w + gap


def matrix(letter):  # matrices in bold, counts and scalars upright
    return "$\\mathbf{%s}$" % letter


X = {k: v + 0.8 for k, v in dict(linear=3.05, relu=3.9, dropout=4.75, head=6.55, latent=10.15, weights=11.45).items()}
X["data"] = X["mask"] = 1.9  # the input and the mask sit further in, to leave room for their labels
h128, h58 = 128 * UNIT, 58 * UNIT
spread = 1.15  # the two heads sit this far above and below the main axis
axis = cy + CELLS * 0.45  # height of the main axis: the arrow from the input, the fork, the sampling and the decoder share it

# input: the counts and, as a separate matrix below, the donor code that is appended to them
data = block(X["data"], cy, PLATE, C["data"], CELLS)
donor = block(X["data"], data["b"] - 0.55 - CELLS * 0.45 - 0.15, 0.3, C["donor"], CELLS)
text(data["x"] - 0.4, data["m"], "input", ha="right")
beside(data, "G"); below(data, "N")
for x, y0, y1 in ((data["x"], data["b"], donor["t"]), (data["r"], data["b"] + data["dy"], donor["t"] + donor["dy"])):  # concatenated
    line([x, x], [y0, y1], C["offsets"][1], 0.8, ls=(0, (1, 2.2)), alpha=0.7)
beside(donor, "D", 0.15); on_side(donor, matrix("S"))
yv = donor["m"] + donor["dy"]  # height of the line from the donor matrix: its rear middle
cm = donor["b"] - 0.25 - PLATE / 2 - 0.07  # the mask sits just under the donor matrix
mask = block(X["mask"], cm, PLATE, C["mask"], LATENTS)
text(mask["x"] - 0.4, mask["m"], "mask " + matrix("M"), ha="right")
beside(mask, "G"); below(mask, "58")

# encoder: the log transform on the arrow, the hidden layer and the two heads
linear, relu, dropout = (block(X[k], axis, h128, C[k], CELLS) for k in ("linear", "relu", "dropout"))
for b, s in ((linear, "linear"), (relu, "ReLU"), (dropout, "dropout")):
    above(b, s)
below(linear, "128")
arrow((data["r"], axis), (linear["x"], axis))
text((data["r"] + linear["x"]) / 2, axis + 0.08, "log(1+x)", va="bottom")
lift = 0.12  # arrows meet and leave the heads this far above their middle, on one horizontal
mean = block(X["head"], axis + spread - lift, h58, C["mean"], CELLS)
logvar = block(X["head"], axis - spread - lift, h58, C["logvar"], CELLS)
above(mean, "$\\bar{\\mu}$", 9); above(logvar, "$\\log\\,\\sigma^{2}$", 9)
below(mean, "58"); below(logvar, "58")
arrow((dropout["r"], axis), (mean["x"], axis + spread)); arrow((dropout["r"], axis), (logvar["x"], axis - spread))

# latent space: the named latents and the free ones, entered and left on the main axis
named, zt = h58 - 8 * UNIT, axis - 0.12 + h58 / 2
z = block(X["latent"], zt - named / 2, named, C["latent"], CELLS)
block(X["latent"], zt - named - 4 * UNIT, 8 * UNIT, C["free"], CELLS, top=False)
above(dict(z, t=zt), matrix("Z"), 9)
text(mid(z), zt - h58 - 0.1, "58", va="top")

# sampling: z = mu + sigma * epsilon with sigma = exp(log variance / 2); the two branches mirror each other around the axis
r_node, r_exp = 0.2, 0.3
y_exp, y_mean = axis - spread, axis + spread
x_mult = logvar["r"] + 1.95
x_exp = (logvar["r"] + x_mult - r_node) / 2  # midway between the log-variance block and the multiplication


def node(x, y, s):
    ax.add_patch(mp.Circle((x, y), r_node, fc="white", ec=ink, lw=0.9, zorder=3))
    text(x, y - 0.02, s, fs=10, zorder=4)


arrow((logvar["r"], y_exp), (x_exp - r_exp, y_exp))
ax.add_patch(mp.Circle((x_exp, y_exp), r_exp, fc="white", ec=C["logvar"][1], lw=0.9, zorder=3))
text(x_exp, y_exp - 0.02, "$e^{x/2}$", fs=8.5, zorder=4)
arrow((x_exp + r_exp, y_exp), (x_mult - r_node, y_exp))
node(x_mult, y_exp, "$\\odot$")
text(x_mult, y_exp - r_node - 0.1, "$\\epsilon\\sim\\mathcal{N}(0,1)$", va="top", fs=9)
arrow((x_mult, y_exp + r_node), (x_mult, axis - r_node)); text(x_mult + 0.1, (y_exp + axis) / 2, "$\\sigma$", ha="left", fs=9)
line([mean["r"], x_mult], [y_mean, y_mean]); arrow((x_mult, y_mean), (x_mult, axis + r_node))
node(x_mult, axis, "+")
arrow((x_mult + r_node, axis), (z["x"], axis))

# decoder: the masked weights, the donor offsets and the softmax, in one line
weights = block(X["weights"], cy, PLATE, C["weights"], LATENTS, shade=0.1)
x_mask, x_offsets = weights["r"] + 0.8, weights["r"] + 1.9  # where the mask and the offsets join the line
softmax = block(x_offsets + 0.9, cy, PLATE, C["softmax"], CELLS)
above(weights, "weights " + matrix("W")); above(softmax, "softmax")
beside(weights, "G"); below(weights, "58"); below(softmax, "N")
text(softmax["r"] + 0.1, softmax["m"] + softmax["dy"], "G", ha="left")
arrow((z["r"], axis), (weights["x"], axis))
arrow((weights["r"], axis), (softmax["x"], axis), ink, 0.9)
equation(x_mask, axis - 0.16, [(matrix("W"), ink), ("$\\odot$", ink), (matrix("M"), C["mask"][1])])
equation(x_offsets, axis - 0.16, [("+", ink), (matrix("V"), C["offsets"][1]), ("·", ink), (matrix("S"), C["offsets"][1])])
y_join = axis - 0.45

# the mask restricts the weights
line([mask["r"] + 0.05, x_mask], [cm, cm], C["mask"][1], 0.9)
arrow((x_mask, cm), (x_mask, y_join), C["mask"][1], 0.9)

# the offsets V (genes x donors) take the donor matrix D straight from its rear middle; a gap where the mask arrow crosses
offsets = block(x_offsets - (FACE + DONORS * 1.15) / 2, yv, PLATE, C["offsets"], DONORS)
below(offsets, matrix("V") + " = G × D")
line([donor["r"], x_mask - 0.14], [yv, yv], C["offsets"][1], 0.9); line([x_mask + 0.14, offsets["x"]], [yv, yv], C["offsets"][1], 0.9)
text((donor["r"] + offsets["x"]) / 2, yv + 0.08, matrix("S"), va="bottom", color=C["offsets"][1])
arrow((x_offsets, offsets["t"] + offsets["dy"] / 2), (x_offsets, y_join), C["offsets"][1], 0.9)

# stages: brackets over the parts, with the same small gap between neighbours
y_br, gap = axis + 2.25, 0.12


def bracket(x0, x1, label):
    line([x0, x1], [y_br, y_br], secondary, 0.9)
    for x in (x0, x1):
        line([x, x], [y_br, y_br - 0.12], secondary, 0.9)
    text((x0 + x1) / 2, y_br + 0.1, label, va="bottom", fs=9, fontweight="bold")


bracket(data["r"] + 0.15, logvar["r"] + 0.1, "encoder")
bracket(logvar["r"] + 0.1 + gap, z["x"] - 0.1 - gap, "sampling")
bracket(z["x"] - 0.1, z["r"] + 0.1, "latent space")
bracket(z["r"] + 0.1 + gap, softmax["r"] + 0.1, "decoder")

for d in outs:
    d.mkdir(parents=True, exist_ok=True)
    fig.savefig(d / "model.png")
plt.close(fig)
