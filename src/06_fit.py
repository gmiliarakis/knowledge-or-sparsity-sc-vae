# Fit one VAE with a linear decoder and a negative binomial likelihood (LDVAE-style).
# Choices are explained in the README.
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

from gene_sets import collection, file_tag, run, tag
from gene_sets import data as dataset  # `data` is the training tensors in main()

n_sets = 50  # one latent per Hallmark set
n_free = 8  # free latents, reaching only genes in no set
n_latent = n_sets + n_free
eval_batch = 1024
config = {"width": 128,  # encoder width, chosen on the vanilla VAE (README)
          "lr": 1e-3, "batch_size": 128, "patience": 3}  # patience as in Eltager et al. (README)


class VAE(nn.Module):
    def __init__(self, n_genes, n_donors, width, mask, median_total):
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(n_genes + n_donors, width), nn.ReLU(), nn.Dropout(0.1))
        self.mean = nn.Linear(width, mask.shape[1])
        self.log_var = nn.Linear(width, mask.shape[1])
        self.w = nn.Parameter(0.01 * torch.randn(n_genes, mask.shape[1]))  # latent -> gene weights
        self.v = nn.Parameter(torch.zeros(n_genes, n_donors))  # donor offsets, also each gene's baseline
        self.log_theta = nn.Parameter(torch.zeros(n_genes))  # NB dispersion per gene
        self.register_buffer("mask", torch.from_numpy(mask))
        self.median_total = median_total

    def encode(self, x, d, total):
        # shifted logarithm (depth normalised to the median total, then ln(1 + x)) plus donor one-hot
        h = self.encoder(torch.cat([torch.log1p(x / total * self.median_total), d], dim=1))
        return self.mean(h), self.log_var(h)

    def loss(self, x, d, total, kl_weight, noise=None):
        mean, log_var = self.encode(x, d, total)
        if noise is None:
            noise = torch.randn_like(mean)
        z = mean + noise * torch.exp(0.5 * log_var)
        share = torch.softmax(z @ (self.w * self.mask).T + d @ self.v.T, dim=1)
        mu = total * share  # NB mean = observed total x predicted share
        theta = torch.exp(self.log_theta)
        nb = torch.distributions.NegativeBinomial(
            total_count=theta, logits=torch.log(mu + 1e-8) - torch.log(theta), validate_args=False
        )
        log_nb = nb.log_prob(x).sum(dim=1)
        kl = 0.5 * (mean**2 + log_var.exp() - log_var - 1).sum(dim=1)
        return (-log_nb + kl_weight * kl).mean(), -log_nb.mean(), kl.mean()


def get_device():
    # GPU when available (NVIDIA, then Apple), otherwise CPU, so the script runs on any machine
    if torch.cuda.is_available():
        return "cuda"
    return "mps" if torch.backends.mps.is_available() else "cpu"


def build_mask(kind, var_names, seed, train_size, subsamples):
    """Return the gene x latent mask and the latent names.

    The first 50 latents reach only the genes of one set each: a Hallmark set (real), a Hallmark set with 100% of gene
    labels permuted (shuffled), or a co-expression module of the same size (coexpression); shuffles and modules use the
    same seed as training. The last 8 latents reach only genes in no set, so set genes can only be explained by named
    latents. The plain VAE has no mask: every latent reaches every gene.
    """
    if kind == "none":
        return np.ones((len(var_names), n_latent), dtype=np.float32), [f"latent {k + 1}" for k in range(n_latent)]

    masks = np.load(f"data/processed/masks_{collection}{file_tag}.npz")
    assert (masks["genes"] == var_names).all(), "gene order in the mask file does not match the data"
    sets = masks["real"]
    latent_ids = list(masks["sets"])
    # on a subsample, shuffle bins and co-expression modules are rebuilt from the subsample's cells
    if kind == "shuffled":
        if train_size:
            perm = subsamples[f"perms_{train_size}"][seed]
        else:
            perm = masks["perms"][list(masks["levels"]).index(1.0), seed]
        sets = sets[perm]
    elif kind == "setshuffle":
        # labels permuted among set genes only: the genes in no set, and so the free latents' reach, stay as in the real mask
        assert not train_size, "setshuffle masks are built for the full training set only"
        sets = sets[np.load(f"data/processed/masks_{collection}_setshuffle{file_tag}.npz")["perms"][seed]]
    elif kind == "responsematch":
        # labels permuted among set genes only, inside bins of expression level x stimulation effect (04b_build_responsematch.py): every set keeps
        # its size, overlaps, expression level and number of strongly responding genes. It uses the stimulation labels of the training cells
        assert not train_size, "responsematch masks are built for the full training set only"
        sets = sets[np.load(f"data/processed/masks_{collection}_responsematch{file_tag}.npz")["perms"][seed]]
    elif kind == "coexpression":
        if train_size:
            modules = subsamples[f"modules_{train_size}"]
        else:
            modules = np.load(f"data/processed/coexpression_{collection}{file_tag}.npz")["modules"]
        sets = modules[seed]
        latent_ids = [f"module {k + 1}" for k in range(n_sets)]

    mask = np.ones((len(var_names), n_latent), dtype=np.float32)
    mask[:, :n_sets] = sets
    mask[:, n_sets:] = ~sets.any(axis=1, keepdims=True)
    return mask, latent_ids + [f"free {k + 1}" for k in range(n_free)]


def batch(data, rows):
    # data = (counts, donor, total), all dense tensors already on the device (counts: ~1.2 GB for Kang et al.)
    rows = torch.as_tensor(rows, device=data[0].device)
    return tuple(t[rows] for t in data)


def evaluate(model, data, rows, seed, device):
    # the same random latent draws in every evaluation, so epoch-to-epoch changes reflect the model, not the draws
    model.eval()
    generator = torch.Generator().manual_seed(seed)
    sums = np.zeros(3)
    with torch.no_grad():
        for i in range(0, len(rows), eval_batch):
            part = rows[i : i + eval_batch]
            noise = torch.randn(len(part), n_latent, generator=generator).to(device)
            losses = torch.stack(model.loss(*batch(data, part), 1.0, noise)).cpu().numpy()  # one sync per chunk
            sums += losses * len(part)
    return sums / len(rows)  # loss, NB, KL per cell


def train_epoch(model, optimizer, data, order, batch_size, cells_seen, warmup_cells):
    model.train()
    total_loss = torch.zeros((), device=next(model.parameters()).device)  # summed on the device: no sync per step
    for i in range(0, len(order), batch_size):
        rows = order[i : i + batch_size]
        kl_weight = min(1.0, cells_seen / warmup_cells) if warmup_cells else 1.0
        loss, _, _ = model.loss(*batch(data, rows), kl_weight)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        cells_seen += len(rows)
        total_loss += loss.detach() * len(rows)
    return total_loss.item() / len(order), kl_weight, cells_seen


def encode_all(model, data, n_obs):
    """Latent means for every cell."""
    model.eval()  # dropout off
    latents = np.zeros((n_obs, n_latent), dtype=np.float32)
    with torch.no_grad():
        for i in range(0, n_obs, eval_batch):
            rows = np.arange(i, min(i + eval_batch, n_obs))
            latents[rows] = model.encode(*batch(data, rows))[0].cpu().numpy()
    return latents


def save_run(out, model, latents, history, info):
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / "model.npz", latents=latents, w=(model.w * model.mask).detach().cpu().numpy(),
             v=model.v.detach().cpu().numpy(), theta=torch.exp(model.log_theta).detach().cpu().numpy())
    pd.DataFrame(history).to_csv(out / "losses.csv", index=False)
    with open(out / "run.json", "w") as f:
        json.dump(info, f, indent=2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mask", choices=["none", "real", "shuffled", "setshuffle", "coexpression", "responsematch"], default="none")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-epochs", type=int, default=2000)
    parser.add_argument("--train-size", type=int, help="learning curves: train on a subsample")
    parser.add_argument("--warmup-epochs", type=float,
                        help="epochs for the KL weight to rise from 0 to 1 (default: 10,000 steps at batch 64 as in "
                             "Eltager et al., scaled to the training set: ~38 epochs at every size)")
    parser.add_argument("--min-delta", type=float, default=0.0,
                        help="relative drop in validation loss that counts as an improvement (0: any drop, as in the "
                             "reported results; 1e-3 stops about 40%% earlier)")
    parser.add_argument("--results-dir", default=f"results/{run}")
    parser.add_argument("--verbose", action="store_true", help="show a progress bar over epochs")
    args = parser.parse_args()

    device = get_device()
    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)

    adata = ad.read_h5ad(f"data/processed/{dataset}{tag}.h5ad")
    counts = adata.layers["counts"].tocsr().astype(np.float32)
    total = np.asarray(counts.sum(axis=1)).ravel().astype(np.float32)
    donor = pd.get_dummies(adata.obs["donor"]).values.astype(np.float32)
    split = adata.obs["split"].values
    idx = {s: np.flatnonzero(split == s) for s in ["train", "val", "test"]}

    # in learning-curve mode, replace the training rows with the saved subsample for this seed and size
    subsamples = None
    if args.train_size:
        subsamples = np.load("data/processed/subsamples_hallmark.npz")
        idx["train"] = subsamples[f"rows_{args.train_size}"][args.seed]

    warmup_epochs = args.warmup_epochs
    if warmup_epochs is None:
        warmup_epochs = 1e4 * 64 / (split == "train").sum()  # in cells seen, so it lasts as many epochs at every size
    warmup_cells = warmup_epochs * len(idx["train"])
    median_total = float(np.median(total[idx["train"]]))
    if args.verbose:
        print(f"{len(idx['train'])} train, {len(idx['val'])} val, {len(idx['test'])} test cells; device {device}")

    variant = "plain" if args.mask == "none" else args.mask
    run_name = f"seed{args.seed}" + (f"_train{args.train_size}" if args.train_size else "")
    out = Path(args.results_dir) / variant / run_name

    mask, latent_ids = build_mask(args.mask, adata.var_names, args.seed, args.train_size, subsamples)
    model = VAE(counts.shape[1], donor.shape[1], config["width"], mask, median_total).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=config["lr"])
    data = (torch.from_numpy(counts.toarray()).to(device),
            torch.from_numpy(donor).to(device),
            torch.from_numpy(total[:, None]).to(device))

    cells_seen, best_loss, best_epoch, stale, history = 0, np.inf, 0, 0, []
    bar = tqdm(unit="epoch", desc=f"{variant}/{run_name}", disable=not args.verbose)  # every epoch is in losses.csv
    for epoch in range(1, args.max_epochs + 1):
        start = time.time()
        train_loss, kl_weight, cells_seen = train_epoch(
            model, optimizer, data, rng.permutation(idx["train"]), config["batch_size"], cells_seen, warmup_cells
        )

        val_loss, val_nb, val_kl = evaluate(model, data, idx["val"], args.seed, device)
        seconds = time.time() - start
        history.append({"epoch": epoch, "train_loss": train_loss, "val_loss": val_loss, "val_nb": val_nb,
                        "val_kl": val_kl, "kl_weight": kl_weight, "seconds": seconds})
        bar.update()
        bar.set_postfix_str(f"train {train_loss:.0f}, val {val_loss:.0f}, KL weight {kl_weight:.2f}, "
                            f"best epoch {best_epoch}")

        # early stopping counts only after warm-up, when the loss has its final form
        if kl_weight < 1:
            continue
        if val_loss < best_loss * (1 - args.min_delta):
            best_loss, best_state, best_epoch, stale = val_loss, copy.deepcopy(model.state_dict()), epoch, 0
        else:
            stale += 1
            if stale >= config["patience"]:
                break

    bar.close()
    assert best_loss < np.inf, "training ended before warm-up finished; raise --max-epochs"
    assert stale >= config["patience"], "reached --max-epochs before early stopping; raise --max-epochs"
    model.load_state_dict(best_state)

    latents = encode_all(model, data, adata.n_obs)
    # test loss is stored for the evaluation and not printed, so it cannot steer choices made now
    test_loss = evaluate(model, data, idx["test"], args.seed, device)[0]

    info = {**vars(args), **config, "warmup_epochs": float(warmup_epochs), "latent_ids": latent_ids,
            "best_epoch": best_epoch, "val_loss": best_loss, "test_loss": test_loss,
            "seconds_per_epoch": float(np.mean([h["seconds"] for h in history]))}
    save_run(out, model, latents, history, info)
    print(f"{variant}/{run_name}: {epoch} epochs, best epoch {best_epoch}, validation loss {best_loss:.1f}, "
          f"{sum(h['seconds'] for h in history):.0f} s, device {device}")


if __name__ == "__main__":
    main()
