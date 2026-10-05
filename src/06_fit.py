# Fit one VAE with a linear decoder and a negative binomial likelihood (LDVAE-style). Choices are explained in the README.
import argparse
import copy
import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from torch import nn
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--mask", choices=["none", "real", "shuffled", "coexpression"], default="none")
parser.add_argument("--seed", type=int, default=0)
parser.add_argument("--width", type=int, default=128)
parser.add_argument("--lr", type=float, default=1e-3)
parser.add_argument("--max-epochs", type=int, default=400)
args = parser.parse_args()

# GPU when available (NVIDIA, then Apple), otherwise CPU, so the script runs on any machine
device = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
torch.manual_seed(args.seed)
rng = np.random.default_rng(args.seed)

n_latent = 50 + 8  # one per Hallmark set + 8 free latents
batch_size = 128
patience = 3
warmup_cells = 10_000 * 64  # Eltager et al.: 10,000 steps at batch 64, here as cells seen

adata = ad.read_h5ad("data/processed/kang.h5ad")
counts = adata.layers["counts"].tocsr().astype(np.float32)
total = np.asarray(counts.sum(axis=1)).ravel().astype(np.float32)
donor = pd.get_dummies(adata.obs["donor"]).values.astype(np.float32)
split = adata.obs["split"].values
idx = {s: np.flatnonzero(split == s) for s in ["train", "val", "test"]}
median_total = float(np.median(total[idx["train"]]))
print(f"{len(idx['train'])} train, {len(idx['val'])} val, {len(idx['test'])} test cells; device {device}")

# Mask: the first 50 latents reach only the genes of one set each: a Hallmark set (real), a Hallmark set with 100% of
# gene labels permuted (shuffled), or a co-expression module of the same size (coexpression); shuffles and modules use
# the same seed as training. The last 8 latents reach only genes in no set, so set genes can only be explained by
# named latents. The plain VAE has no mask: every latent reaches every gene.
mask = np.ones((counts.shape[1], n_latent), dtype=np.float32)
latent_ids = [f"latent {k}" for k in range(n_latent)]
if args.mask != "none":
    masks = np.load("data/processed/masks_hallmark.npz")
    assert (masks["genes"] == adata.var_names).all()
    sets = masks["real"]
    latent_ids = list(masks["sets"])
    if args.mask == "shuffled":
        sets = sets[masks["perms"][list(masks["levels"]).index(1.0), args.seed]]
    if args.mask == "coexpression":
        sets = np.load("data/processed/coexpression_hallmark.npz")["modules"][args.seed]
        latent_ids = [f"module {k + 1}" for k in range(50)]
    mask[:, :50] = sets
    mask[:, 50:] = ~sets.any(axis=1, keepdims=True)
    latent_ids += [f"free {k + 1}" for k in range(8)]


class VAE(nn.Module):
    def __init__(self, n_genes, n_donors, width, mask):
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(n_genes + n_donors, width), nn.ReLU(), nn.Dropout(0.1))
        self.mean = nn.Linear(width, mask.shape[1])
        self.log_var = nn.Linear(width, mask.shape[1])
        self.w = nn.Parameter(0.01 * torch.randn(n_genes, mask.shape[1]))  # latent -> gene weights
        self.v = nn.Parameter(torch.zeros(n_genes, n_donors))  # donor offsets, also each gene's baseline
        self.log_theta = nn.Parameter(torch.zeros(n_genes))  # NB dispersion per gene
        self.register_buffer("mask", torch.from_numpy(mask))

    def encode(self, x, d, total):
        # shifted logarithm (depth normalised to the median total, then ln(1 + x)) plus donor one-hot
        h = self.encoder(torch.cat([torch.log1p(x / total * median_total), d], dim=1))
        return self.mean(h), self.log_var(h)

    def loss(self, x, d, total, kl_weight, noise=None):
        mean, log_var = self.encode(x, d, total)
        if noise is None:
            noise = torch.randn_like(mean)
        z = mean + noise * torch.exp(0.5 * log_var)
        share = torch.softmax(z @ (self.w * self.mask).T + d @ self.v.T, dim=1)
        mu = total * share  # NB mean = observed total x predicted share
        theta = torch.exp(self.log_theta)
        nb = torch.distributions.NegativeBinomial(total_count=theta, logits=torch.log(mu + 1e-8) - torch.log(theta))
        log_nb = nb.log_prob(x).sum(dim=1)
        kl = 0.5 * (mean**2 + log_var.exp() - log_var - 1).sum(dim=1)
        return (-log_nb + kl_weight * kl).mean(), -log_nb.mean(), kl.mean()


def batch(rows):
    # sparse counts are sliced on the CPU, then moved to the device as dense tensors
    return (
        torch.from_numpy(counts[rows].toarray()).to(device),
        torch.from_numpy(donor[rows]).to(device),
        torch.from_numpy(total[rows, None]).to(device),
    )


def evaluate(rows):
    # the same random latent draws in every evaluation, so epoch-to-epoch changes reflect the model, not the draws
    generator = torch.Generator().manual_seed(args.seed)
    sums = np.zeros(3)
    with torch.no_grad():
        for i in range(0, len(rows), 1024):
            part = rows[i : i + 1024]
            noise = torch.randn(len(part), n_latent, generator=generator).to(device)
            sums += np.array([t.item() for t in model.loss(*batch(part), 1.0, noise)]) * len(part)
    return sums / len(rows)  # loss, NB, KL per cell


model = VAE(counts.shape[1], donor.shape[1], args.width, mask).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

cells_seen, best_loss, stale, history = 0, np.inf, 0, []
for epoch in range(1, args.max_epochs + 1):
    start = time.time()
    model.train()
    order = rng.permutation(idx["train"])
    train_loss = 0.0
    for i in tqdm(range(0, len(order), batch_size), desc=f"epoch {epoch}", leave=False):
        rows = order[i : i + batch_size]
        kl_weight = min(1.0, cells_seen / warmup_cells)
        loss, _, _ = model.loss(*batch(rows), kl_weight)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        cells_seen += len(rows)
        train_loss += loss.item() * len(rows)
    train_loss /= len(order)

    model.eval()
    val_loss, val_nb, val_kl = evaluate(idx["val"])
    seconds = time.time() - start
    history.append({"epoch": epoch, "train_loss": train_loss, "val_loss": val_loss, "val_nb": val_nb,
                    "val_kl": val_kl, "kl_weight": kl_weight, "seconds": seconds})
    print(f"epoch {epoch:3d}  train {train_loss:8.1f}  val {val_loss:8.1f} (NB {val_nb:8.1f}, KL {val_kl:5.1f})  "
          f"KL weight {kl_weight:.2f}  {seconds:.1f}s")

    # early stopping counts only after warm-up, when the loss has its final form
    if kl_weight < 1:
        continue
    if val_loss < best_loss:
        best_loss, best_state, best_epoch, stale = val_loss, copy.deepcopy(model.state_dict()), epoch, 0
    else:
        stale += 1
        if stale >= patience:
            break

assert best_loss < np.inf, "training ended before warm-up finished (~38 epochs); raise --max-epochs"
model.load_state_dict(best_state)
model.eval()
print(f"best epoch {best_epoch}, validation loss {best_loss:.1f}")

# latent means for every cell
latents = np.zeros((adata.n_obs, n_latent), dtype=np.float32)
with torch.no_grad():
    for i in range(0, adata.n_obs, 1024):
        rows = np.arange(i, min(i + 1024, adata.n_obs))
        latents[rows] = model.encode(*batch(rows))[0].cpu().numpy()

# test loss is stored for the evaluation and not printed, so it cannot steer choices made now
test_loss = evaluate(idx["test"])[0]

variant = {"none": "plain", "real": "real", "shuffled": "shuffled", "coexpression": "coexpression"}[args.mask]
out = Path("results") / variant / f"seed{args.seed}_width{args.width}_lr{args.lr:g}"
out.mkdir(parents=True, exist_ok=True)
np.savez(out / "model.npz", latents=latents, w=(model.w * model.mask).detach().cpu().numpy(),
         v=model.v.detach().cpu().numpy(), theta=torch.exp(model.log_theta).detach().cpu().numpy())
pd.DataFrame(history).to_csv(out / "losses.csv", index=False)
with open(out / "run.json", "w") as f:
    json.dump({**vars(args), "latent_ids": latent_ids, "best_epoch": best_epoch, "val_loss": best_loss,
               "test_loss": test_loss, "seconds_per_epoch": float(np.mean([h["seconds"] for h in history]))}, f, indent=2)
print("saved", out)
