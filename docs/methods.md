# Methods in detail

Supporting detail for the Methods in the [README](../README.md). Numbers in square brackets refer to its reference list.

## Co-expression mask: full recipe

- **Co-expression mask** (`05_build_coexpression_null.py`). This mask tests whether pre-specified sets add anything beyond
  co-expressed genes found in the data. The groups are built from the data alone, without condition or cell-type
  labels.

  1. **Expression** Training cells only, stimulated and control together. Counts are scaled to the median total
     count and log-transformed (ln(1 + x)), the same transform the encoder reads.
  2. **Correlation** Pearson correlation between all 12,034 genes. A gene with no variance in training cells has
     correlation 0 with every other gene.
  3. **Seed genes** For each seed, 50 distinct seed genes are drawn at random from the 6,753 genes detected in ≥ 1%
     of training cells. The nearest neighbours of rarer genes are noise, and these genes are excluded.
  4. **Groups** Group j is seed gene j plus the n − 1 genes most positively correlated with it, where n is the
     size of Hallmark set j (12–197 genes). Group j takes the size of Hallmark set j and none of its biology. Its
     latent is labelled "module j" in the code. Groups are built independently and may share genes.

  The groups have three properties that affect how the results are read.

  - **Coherence** The median within-set correlation is 0.096 for the groups and 0.003 for Hallmark sets. In single
    PBMCs most members of a Hallmark set barely co-vary, whereas a group is built to co-vary. The groups are therefore a demanding
    comparison for the Hallmark sets.
  - **Overlap** The groups have the same 5,568 memberships as the Hallmark sets but cover 2,338–3,233 genes (Hallmark:
    3,187). The overlap concentrates on genes correlated with many others. Ribosomal-protein genes (RPS6, RPL3, …)
    and PTMA are in up to 28 groups per seed, whereas no Hallmark gene is in more than 10 sets. Several groups
    partly describe the same ribosomal axis. The 8 free latents reach the 8,801–9,696 genes in no group.
  - **Interferon** The groups were checked against the stimulation labels after they were built. In every seed,
    the group most raised by IFN-β is enriched for Hallmark interferon genes: 19–89% of its genes, against a base
    rate of 1.8%. In one seed it grew around UBE2L6 and contains ISG15, IFI6, IFIT1, IFIT3, MX1, IRF7 and OAS1. The
    data alone recover the core interferon response.

## Model equations

The encoder first puts the counts of a cell on a common scale. Each count $x_g$ is divided by the
total count of the cell $\ell$, multiplied by the median total $m$ of the training cells and log-transformed. Deeply and
shallowly sequenced cells then look alike.

$$\tilde{x}_g = \ln\big(1 + x_g \cdot m / \ell\big)$$

The encoder reads $\tilde{x}$ and the donor code $\mathbf{s}$ (the column of the one-hot donor matrix $\mathbf{S}$ of Figure 2 that belongs to the cell) and
proposes where the cell lies in the latent space, as a mean $\bar{\mu}$ and a log variance $\log\sigma^2$ for each of the 58 latents, the
posterior $q(z \mid x, \mathbf{s})$. Training draws a point $z$ from it; the evaluation uses its mean. The latent matrix $\mathbf{Z}$ of Figure 2 holds the $z$ of every cell.

The decoder turns $z$ back into expression. Each latent raises or lowers the genes that the
mask $\mathbf{M}$ lets it reach, through its weights $\mathbf{W}$, and the donor adds its offsets $\mathbf{V}$. A softmax over genes converts these
scores into shares $\pi_g$ that sum to 1, the fraction of the counts of the cell that each gene should receive.

$$\pi = \mathrm{softmax}\big((\mathbf{W} \odot \mathbf{M})\,z + \mathbf{V}\,\mathbf{s}\big)$$

The expected count of gene $g$ is the total count of the cell times its share. The observed count scatters around it as a
negative binomial (NB). The dispersion $\theta_g$ allows each gene to be noisier than a Poisson count.

$$x_g \sim \mathrm{NB}\big(\mu_g = \ell \cdot \pi_g,\ \theta_g\big)$$

The loss has two parts. The reconstruction term measures how unlikely the observed counts are under
the prediction. The Kullback–Leibler (KL) term pulls the posterior of each cell towards a standard normal prior. This
keeps the latents on a common scale and lets those that the data do not need become inactive.

$$\mathcal{L} = -\log p_{\mathrm{NB}}(x \mid \mu, \theta) + \beta \cdot \mathrm{KL}\big(q(z \mid x, \mathbf{s})\ \|\ \mathcal{N}(0, I)\big)$$

where, for one cell:

- $\mathcal{L}$: the loss per cell, averaged over each batch. With $\beta = 1$ it is the negative evidence lower bound
  (ELBO).
- $-\log p_{\mathrm{NB}}(x \mid \mu, \theta)$: the reconstruction term, the negative log-likelihood summed over genes.
- $x_g$: the raw count of gene $g$. $\ell = \sum_g x_g$ is the total count of the cell.
- $\mu_g = \ell\,\pi_g$: the expected count of gene $g$.
- $\pi_g$: the share of gene $g$ from a softmax over genes; the shares are positive and sum to 1.
- $\mathbf{W} \in \mathbb{R}^{G \times 58}$: the latent → gene weights (`w` in the code), with $G$ the number of genes (12,034 in the Kang et al. data).
  $\mathbf{M} \in \lbrace 0,1\rbrace^{G \times 58}$ is the mask, with 1 where latent $k$ may reach gene $g$ and 0 elsewhere. The product $\odot$ is elementwise.
- $\mathbf{V} \in \mathbb{R}^{G \times D}$: the donor offsets (`v`), which also set the baseline of each gene, with $D$ the number of donors (8 in the Kang et al. data).
- $\theta_g$: the dispersion of gene $g$ (`exp(log_theta)`), with $\mathrm{Var}(x_g) = \mu_g + \mu_g^2 / \theta_g$.
- $\beta$: the KL weight. It rises linearly from 0 to 1 over the first ≈ 38 epochs (KL warm-up). The model first learns to reconstruct, and then stays at 1.
- $\mathrm{KL}\big(q \,\|\, \mathcal{N}(0, I)\big)$: the KL divergence of the posterior from the standard normal prior,
  summed over latents.
- $z \in \mathbb{R}^{58}$: the latent values of the cell, one column of the latent matrix $\mathbf{Z}$ ($58 \times N$) in Figure 2. Training draws them from
  $q(z \mid x, \mathbf{s}) = \mathcal{N}\big(\bar{\mu},\ \mathrm{diag}\,\sigma^2\big)$ as $z = \bar{\mu} + \sigma \odot \epsilon$, with $\sigma = \exp(\tfrac{1}{2}\log\sigma^2)$
  and $\epsilon \sim \mathcal{N}(0, I)$ (`mean`, `log_var` in the code); the evaluation uses the posterior means.
- $\mathbf{s} \in \lbrace 0,1\rbrace^{D}$: the donor one-hot code of the cell, one column of the donor matrix $\mathbf{S}$ ($D \times N$ for $N$ cells) in Figure 2.
- $\tilde{x}_g$: the encoder input, the count of gene $g$ divided by $\ell$, multiplied by $m$, the median total count of the training cells, and log-transformed.

## Architecture and parameter counts

Every variant uses the same VAE, with an encoder of one hidden layer and a linear decoder (binary mask on its weight matrix in five of the six variants). The PyTorch module summary is as follows.

```
VAE(
  (encoder): Sequential(
    (0): Linear(in_features=12042, out_features=128, bias=True)
    (1): ReLU()
    (2): Dropout(p=0.1, inplace=False)
  )
  (mean): Linear(in_features=128, out_features=58, bias=True)
  (log_var): Linear(in_features=128, out_features=58, bias=True)
)
```

The encoder input is the 12,034 genes after a shifted-log transform (see the README), plus a one-hot code for the 8 donors. These are the Kang et al. data. The Theis et al. data have about 11,340 genes and 12 donors. The encoder, `w` and `v` are smaller. The
decoder has no layers of its own, only the three parameters $\mathbf{W}$, $\mathbf{V}$ and $\log\theta$ (`w`, `v` and `log_theta` in the code).

| Parameter | Shape | Role | Count |
|---|---|---|---|
| encoder | 12,042 → 128 | hidden layer | 1,541,504 |
| mean, log_var | 128 → 58 (each) | mean and log-variance of the posterior $q(z \mid x, \mathbf{s})$ | 14,964 |
| $\mathbf{W}$ (`w`) | 12,034 × 58 | latent → gene weights, multiplied elementwise by the mask $\mathbf{M}$ | 697,972 |
| $\mathbf{V}$ (`v`) | 12,034 × 8 | donor offsets (also each gene's baseline) | 96,272 |
| $\log\theta$ (`log_theta`) | 12,034 | negative binomial dispersion per gene | 12,034 |
| **total** | | | **2,362,746** |

The mask decides how many of the 697,972 weights in `w` can be trained (table in the README, Model).

## Scoring code

**Code** (`src/08_evaluate.py`, abridged; `z` contains the posterior means, `w` the masked decoder weights):

```python
def auroc(score, split):
    # stim vs ctrl per cell type, so a latent cannot score by separating cell types
    r = rows[split]
    return pd.Series({t: roc_auc_score(stim[r][cell_type[r] == t], score[r][cell_type[r] == t]) for t in types})

activity = z[val_all].var(axis=0)
active = activity > min_variance  # 0.01

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
score = alpha_auc.mean()  # unweighted mean over the 7 cell types
```

**The score measures** how well one latent separates stimulated from control cells. Any latent
whose genes respond to the stimulation can do so. The Hallmark interferon-α and
interferon-γ sets share 71 of the 95 interferon-α genes. The score therefore does not distinguish type I from type II
interferon. A latent of the fully random mask contains few interferon-α genes (3 of 95 in seed 0) and still reached up to 0.98 in a seed. The latent at the interferon-α position of the within-set random mask scored 0.93 or higher in 6 of 10 seeds.
Near 1, the AUROC is at its maximum and says little about the size of the effect, because it measures only how the cells rank. Rare cell types have few test cells. Their
AUROC values are noisy.

## Reactome sets: full recipe

The README gives the idea. This section gives the rules that choose the 50 sets and the draw that was kept.

- **Other 45 slots** Each slot takes a Reactome set of the size of its Hallmark set (within 10%), drawn at random from the sets that
  - overlap every set chosen so far by a Jaccard index below 0.34, the largest overlap between two Hallmark sets;
  - have at least 50% of their genes among the 12,034 QC genes (the 10th percentile of the Hallmark sets, which have 0.30 to 0.98);
  - are not under the Disease branch of Reactome's pathway hierarchy, as Hallmark has no disease or pathogen sets (the 39 sets that cannot be matched to a Reactome pathway by name are excluded only when their name contains SARS, DEFECT, INFECTION, HIV, DISEASE, CANCER or MUTANT);
  - contain none of INTERFERON, IFN, IRF, TNF, IL-6 or IL-2 in the name, as Hallmark has no second set for those stimuli.
- **Sparsity** Selection seeds 0, 1, 2 and so on were tried in order. The first draw with as many trainable decoder weights as the Hallmark mask (within 1%) was kept. Seed 2 gave 76,579 against 76,344 (5,003 memberships and 3,087 genes in sets, against 5,568 and 3,187; median set size 99 against 117.5; largest overlap 0.33).
- **Draw** The draw includes developmental and tissue-specific sets that PBMCs express, as Hallmark does. The Reactome hierarchy file is Reactome's release of 20 October 2026, which differs from the release behind MSigDB v2024.1.

## Human Cytokine Dictionary: full description

The README describes the data set. These are the choices that define the subset used for Q1e and Q1f.

- **Cytokines** Interferon-β, as a replication of the stimulation in the Kang et al. data, and one cytokine that is not an interferon. Its choice was fixed from the published differential-expression table alone, before any model was fitted. Among TNF-α, IL-6 and IL-2, the three that have a Reactome set among the 50, it is the one with the most differentially expressed genes summed over cell types (adjusted p < 0.05, |log fold change| ≥ 0.5, well-biased rows excluded). IL-2 had 25,835, IL-6 1,521 and TNF-α 1,419 (interferon-β 31,780).
- **Cells** The 41 GB file is read in place from the published archive (only the chosen runs of cells are transferred, about 1 GB). Cells are read until each donor × condition group has 1,800. QC as for the Kang et al. data (5 median absolute deviations per cell type and condition) removes 0.9% for IL-2, and at most 1,200 cells per group are kept, which gives 28,800 cells for each cytokine.
- **Genes** The 11,344 of the 12,031 genes of the Kang et al. data that are detected in at least 20 cells of the stimulation (11,340 for interferon-β). The masks and the weight match of the Reactome draw therefore keep their design.
- **Split and scored cell types** 70/15/15 within every donor × condition × cell type group. A cell type is scored when the test split has at least 20 cells of both conditions: 13 of 18 types (HSPC, ILC, NKT, plasmablasts and pDC are not).

## Biology in the code

The table lists the steps at which biological knowledge or an assumption about the cells enters the pipeline.

| Step | Script | Choice | Reason |
|---|---|---|---|
| Data | `03c` | Theis et al. data for Q1e and Q1f: 12 healthy donors, cells read in place from the published archive, at most 1,200 per donor × condition group, genes shared with the Kang et al. data (`SC_DATA`) | every donor appears in both conditions; the cytokine was chosen from the published differential-expression table before any model was fitted |
| Cell identity | `03` | the authors' demuxlet singlet calls and cell-type labels, used as given | a droplet containing cells of two donors is identified from their genotypes, which is stronger evidence than an expression-based doublet score |
| Cell QC | `03` | thresholds per condition × cell type | monocytes, dendritic cells and megakaryocytes differ in RNA content, and stimulation shifts it |
| Genes | `03` | all genes detected in ≥ 20 cells; no highly-variable-gene selection | selection drops interferon-induced genes |
| Split | `03` | stratified by donor × condition × cell type | every group appears in train, validation and test |
| Gene identity | `04` | Ensembl IDs mapped to current HGNC symbols (1,256 symbols changed since 2017; 427 genes HGNC no longer lists keep their 2017 symbol) | gene set symbols are current |
| Gene sets | `04` | Hallmark, sets of ≥ 12 genes (all 50); for Q1d 50 Reactome sets chosen from names, sizes, overlaps and the pathway hierarchy (`SC_SETS=reactome`) | Hallmark has one unambiguous target, the interferon-α response; Reactome is curated pathway by pathway and not built from expression |
| Fully random | `04` | labels permuted within expression bins | random genes are as highly expressed as the genes they replace |
| Within-set random | `04` | labels permuted within expression bins, among genes that are in a set | the genes in no set, and so the free latents' reach, stay as in the Hallmark mask |
| Co-expression | `05` | groups of co-expressed genes, sizes matched to Hallmark | see Gene sets and masks |
| Free latents | `06` | reach only genes in no set | when they reached all genes, they took over the stimulation signal |
| Covariates | `06` | donor in encoder and decoder; condition not used | the stimulation is the signal the latents should find |
| Likelihood | `06` | negative binomial on raw counts, observed total counts as library size | droplet counts are not zero-inflated [22] |
| Target latent | `08` | Hallmark interferon-α response or Reactome interferon α/β signaling for IFN-β; IL-2 STAT5 signaling or interleukin-2 family signaling for IL-2 (`SC_TARGET=il2`) | IFN-β is a type I interferon; IL-2 is scored against its own set |
| Orientation | `08` | sign set so the latent's mean weight on its own genes is positive | a gene set that is switched on raises its own genes |
| Scoring | `08` | AUROC per cell type on the test set, unweighted mean over the scored types (7 for the Kang et al. data, 13 for the Theis et al. data) | the question is whether every cell type responds to IFN-β, and a latent should not score by separating cell types |
| Megakaryocytes | `08` | excluded from scores, kept in training | mostly platelets, which have no nucleus |

## Data sets used and for future work

| Question | Dataset | Access |
|---|---|---|
| Q1, Q1a to Q1e | Kang et al. 2018 [17] | GEO GSE96583 (batch 2) |
| Q1d to Q1f | Human Cytokine Dictionary of Theis et al. [41], as distributed by Parse Biosciences [42] | Parse Biosciences, "10 Million Human PBMCs in a Single Experiment" (CC BY-NC 4.0) |
| Known regulators, and B cells | Human tonsil atlas [24] | Zenodo 10.5281/zenodo.8373756 |
| Known regulators | B cell activation time course and IRF4 / PRDM1 CRISPR perturbations [23] | Zenodo 10.5281/zenodo.17984776 |
| Protein as reference | Tonsil CITE-seq [24]; 10x Genomics 10k PBMCs with TotalSeq-B antibodies | Zenodo 10.5281/zenodo.8373756; 10x Genomics |
| Imperfect gene sets | Mouse gastrulation [25] | ArrayExpress E-MTAB-6967 |
