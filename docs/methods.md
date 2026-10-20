# Methods in detail

Details moved out of the README so that it stays short. Numbers in square brackets refer to the reference list in the
README.

## Co-expression mask: full recipe

- **Co-expression mask** (`05_build_coexpression_null.py`). This null tests whether curated sets add anything beyond
  co-expressed genes. The groups are built from the data alone, without condition or cell-type
  labels.

  1. **Expression.** Training cells only, stimulated and control together. Counts are scaled to the median total
     count and log-transformed (ln(1 + x)), the same transform the encoder reads.
  2. **Correlation.** Pearson correlation between all 12,034 genes. A gene with no variance in training cells has
     correlation 0 with every other gene.
  3. **Seed genes.** For each seed, 50 distinct seed genes are drawn at random from the 6,753 genes detected in ≥ 1%
     of training cells. The nearest neighbours of rarer genes are noise, and these genes are excluded.
  4. **Groups.** Group j is seed gene j plus the n − 1 genes most positively correlated with it, where n is the
     size of Hallmark set j (12–197 genes). Group j takes the size of Hallmark set j and none of its biology. Its
     latent is labelled "module j" in the code. Groups are built independently and may share genes.

  The groups have three properties that affect how the results are read.

  - **Coherence.** The median within-set correlation is 0.096 for the groups and 0.003 for Hallmark sets. In single
    PBMCs most members of a Hallmark set barely co-vary. A group, in contrast, is built to co-vary. This null is
    therefore a strong comparison.
  - **Overlap.** The groups have the same 5,568 memberships as the Hallmark sets but cover 2,338–3,233 genes (Hallmark:
    3,187). The overlap concentrates on genes correlated with many others. Ribosomal-protein genes (RPS6, RPL3, …)
    and PTMA are in up to 28 groups per seed. No Hallmark gene is in more than 10 sets. Several groups
    partly describe the same ribosomal axis. The 8 free latents reach the 8,801–9,696 genes in no group.
  - **Interferon.** The groups were checked against the stimulation labels after they were built. In every seed,
    the group most raised by IFN-β is enriched for Hallmark interferon genes: 19–89% of its genes, against a base
    rate of 1.8%. In one seed it grew around UBE2L6 and contains ISG15, IFI6, IFIT1, IFIT3, MX1, IRF7 and OAS1. The
    data alone recover the core interferon response.

## Architecture and parameter counts

Every variant uses the same VAE, with an encoder of one hidden layer and a linear decoder (binary mask on its weight matrix in four of the five variants). The PyTorch module summary is as follows.

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

The encoder input is the 12,034 genes after a shifted-log transform (below) plus a one-hot code for the 8 donors. The
decoder has no layers of its own, only the three parameters $W$, $V$ and $\log\theta$ (`w`, `v` and `log_theta` in the code).

| Parameter | Shape | Role | Count |
|---|---|---|---|
| encoder | 12,042 → 128 | hidden layer | 1,541,504 |
| mean, log_var | 128 → 58 (each) | mean and log-variance of the posterior $q(z \mid x, d)$ | 14,964 |
| $W$ (`w`) | 12,034 × 58 | latent → gene weights, multiplied elementwise by the mask $M$ | 697,972 |
| $V$ (`v`) | 12,034 × 8 | donor offsets (also each gene's baseline) | 96,272 |
| $\log\theta$ (`log_theta`) | 12,034 | negative binomial dispersion per gene | 12,034 |
| **total** | | | **2,362,746** |

The mask decides how many of the 697,972 weights in `w` can be trained (table in the README, Model).

## Scoring code

**Code** (`src/08_evaluate.py`, abridged; `z` holds the posterior means, `w` the masked decoder weights):

```python
def auroc(score, split):
    # stim vs ctrl within each cell type, so a latent cannot score by separating cell types
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

**Evidence for the score.** The named latent tracks the stimulation within every cell type. An independent
recomputation from the saved latent values with the Mann–Whitney U statistic reproduced the scores of four runs (one per
variant) to six decimals.

**Limits of the score.** The score does not show that a latent represents interferon biology specifically. Any latent
whose genes respond to the stimulation can separate stimulated from control cells. The Hallmark interferon-α and
interferon-γ sets share 71 of the 95 interferon-α genes. The score therefore cannot tell type I from type II
interferon. A shuffled latent holds few interferon-α genes (3 of 95 in seed 0) and still reached up to 0.98 in a seed; the latent at the interferon-α position of the within-set random mask scored 0.93 or higher in 6 of 10 seeds.
Near 1, the AUROC saturates and measures ranking, not the size of the effect. Rare cell types have few test cells. Their
AUROC values are noisy.

## Biology in the code

The table lists where biological knowledge or an assumption about the cells enters the pipeline.

| Step | Script | Choice | Reason |
|---|---|---|---|
| Cell identity | `03` | the authors' demuxlet singlet calls and cell-type labels, used as given | a droplet holding two donors' cells is identified from their genotypes, which is stronger evidence than an expression-based doublet score |
| Cell QC | `03` | thresholds per condition × cell type | monocytes, dendritic cells and megakaryocytes differ in RNA content, and stimulation shifts it |
| Genes | `03` | all genes detected in ≥ 20 cells; no highly-variable-gene selection | selection drops interferon-induced genes |
| Split | `03` | stratified by donor × condition × cell type | every group appears in train, validation and test |
| Gene identity | `04` | Ensembl IDs mapped to current HGNC symbols (1,256 symbols changed since 2017; 427 genes HGNC no longer lists keep their 2017 symbol) | gene set symbols are current |
| Gene sets | `04` | Hallmark, sets of ≥ 12 genes (all 50) | one unambiguous target, the interferon-α response |
| Fully random | `04` | labels permuted within expression bins | random genes are as highly expressed as the genes they replace |
| Within-set random | `04` | labels permuted within expression bins, among genes that are in a set | the genes in no set, and so the free latents' reach, stay as in the Hallmark mask |
| Co-expression | `05` | groups of co-expressed genes, sizes matched to Hallmark | see Gene sets and masks |
| Free latents | `06` | reach only genes in no set | when they reached all genes, they absorbed the stimulation signal |
| Covariates | `06` | donor in encoder and decoder; condition never | the stimulation is the signal the latents should find |
| Likelihood | `06` | negative binomial on raw counts, observed total counts as library size | droplet counts are not zero-inflated [22] |
| Target latent | `08` | HALLMARK_INTERFERON_ALPHA_RESPONSE | IFN-β is a type I interferon |
| Orientation | `08` | sign set so the latent's mean weight on its own genes is positive | a gene set that is switched on raises its own genes |
| Scoring | `08` | AUROC within cell type, unweighted mean over 7 types | the question is whether every cell type responds to IFN-β, and a latent should not score by separating cell types |
| Megakaryocytes | `08` | excluded from scores, kept in training | mostly platelets, which have no nucleus |

## Candidate datasets for future work

| Question | Dataset | Access |
|---|---|---|
| Q1, Q1b | Kang et al. 2018 [17] | GEO GSE96583 (batch 2) |
| Q1, Q2 | Human tonsil atlas [24] | Zenodo 10.5281/zenodo.8373756 |
| Q2 | B cell activation time course and IRF4 / PRDM1 CRISPR perturbations [23] | Zenodo 10.5281/zenodo.17984776 |
| Q3 | Tonsil CITE-seq [24]; 10x Genomics 10k PBMCs with TotalSeq-B antibodies | Zenodo 10.5281/zenodo.8373756; 10x Genomics |
| Q4 | Mouse gastrulation [25] | ArrayExpress E-MTAB-6967 |
