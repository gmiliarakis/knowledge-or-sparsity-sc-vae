# Knowledge or sparsity? Gene set masks in a single-cell VAE

![From stimulated cells to a masked VAE: the four decoder masks and the interferon-α score each one reached](reports/figures/overview.png)

*Stage 3 shows the decoder mask of each variant as trained with seed 0: rows are the 12,034
genes sorted by Hallmark membership, columns are the 58 latents, and each dot is a kept link between a gene and a latent.
Stage 4 shows the interferon-α score of the scored latent in all 10 seeds, with the mean printed under each
variant.*

**Contents**

- [Background](#background)
- [Methods](#methods): [study design](#study-design), [data](#data), [gene sets and masks](#gene-sets-and-masks),
  [model](#model), [evaluation](#evaluation), [learning curves](#learning-curves-q1b)
- [Results](#results): [Q1](#q1-knowledge-or-sparsity), [Q1b](#q1b-learning-curves)
- [Limitations](#limitations) · [Reproducing](#reproducing) · [Future work](#future-work) · [References](#references)

## Summary

**At a glance**

- **Question:** in a single-cell RNA-seq variational autoencoder (VAE) restricted to gene sets, does the biological knowledge in the sets or the sparsity they impose explain how well it captures a known response?
- **Data:** blood cells (PBMCs) from 8 lupus patients, with and without interferon-β stimulation (Kang et al. [17]; 23,919 cells, 12,034 genes).
- **Response:** the interferon response switched on by IFN-β (a type I interferon). Success is one latent that separates stimulated from control cells (AUROC); the interferon-α response gene set is the matching biological target.
- **Model:** VAE with a linear decoder. A binary mask on the decoder weights limits which genes each latent variable can reach.
- **Four variants:** no mask, gene set mask, random mask of the same sparsity, mask of co-expressed gene modules from the data alone.
- **Gene sets:** Hallmark collection of MSigDB v2024.1, 50 expert-curated sets (e.g. interferon response, hypoxia, inflammation).
- **Evaluation:** how well one latent distinguishes stimulated from control cells, over 10 seeds; learning curves at smaller training sizes.
- **Limitation:** one dataset with one very strong biological response.

**Key findings**

- Sparsity alone is enough to capture the response: a random mask reaches it in its best latent.
- Knowledge makes the model capture the biology: the response lands in the latent tied to the interferon-α genes, not to unrelated ones.
- Modules built from the data alone capture it about as well, but finding their latent needs the stimulation labels.
- Random and gene set masks both capture the response slightly better than an unmasked VAE.
- Knowledge did not make learning more data-efficient on this task.

## Background

Knowledge-guided single-cell models build biological knowledge into a model so that its output can be read in
biological terms. Prior knowledge can enter a model in three ways [2]:

- as input features, such as gene set or transcription-factor activity scores;
- in the architecture, as sparse or masked layers that mirror gene sets;
- as a graph of gene interactions.

This project tests the architecture route. Each latent variable (a hidden factor that the model infers for every cell)
is tied to a gene set: a masked linear decoder lets latent k reach only the genes of set k, so its value in a cell reads
as the activity of that gene set. Such models are interpretable by design, and they are expected to need fewer cells
because they need not re-discover patterns that are already well known.

Several published models are built this way, and others use gene sets in related ways.

| Model | Kind |
|---|---|
| VEGA [3], expiMap [4] | single-cell VAE with a masked decoder |
| pmVAE [5] | single-cell VAE with one sub-network per pathway |
| f-scLVM [6], PLIER [7], Spectra [8] | factor model |
| KPNN [9] | supervised neural network |
| MOFA-FLEX [10] | factor model |
| Moullet et al. [11] | self-supervised model |

In 29 pathway-informed neural networks trained for prediction, randomized models that kept the sparsity and architecture
of the originals matched or outperformed them [1]. VEGA [3] and expiMap [4], which both evaluate on
interferon-β-stimulated PBMCs, report masked decoders without a random gene set control. The closest control, in
MOFA-FLEX [10], swaps 30% of the genes of each set with random non-members and tests whether the original programs are
recovered. It is not known whether the knowledge in the gene sets helps in unsupervised single-cell models, or only the
sparsity that the mask imposes. Earlier work leaves three questions open, and each shapes the design.

- **Knowledge or sparsity?** (Q1) A gene set mask fixes 89% of the decoder weights at zero. Any gain over a model without a mask could
  come from the sparsity alone. A random mask keeps the size and overlaps of every set and changes only its
  biology. It isolates the effect of the knowledge.
- **Curated knowledge or any co-expression?** (Q1) Curated gene sets align with the main axes of expression more than
  size-matched random sets do in bulk tumour data [12]. A Hallmark mask could beat a random one only because its genes are
  co-expressed. Modules of co-expressed genes of the same sizes test whether curation adds anything beyond that.
- **Less data needed?** (Q1b) expiMap was reported to learn more sample-efficiently than the linear-decoder VAE LDVAE [13] when fewer training samples were available [4]. That comparison did not include a random gene set control. Knowledge and sparsity were confounded. Learning curves with all four
  variants test this directly.

## Methods

Step-by-step detail, the architecture printout, symbols and scoring code are in [docs/methods.md](docs/methods.md).

### Study design

The figure at the top summarises the design: the four masks and the score each one reached.

The benchmark is four variants that share one architecture and one training setup, differ only in the decoder mask
and are each fit with 10 seeds. Baselines and structure-preserving controls are included because randomized pathways and simple baselines can match more elaborate models [1, 15, 16].

| Variant | Decoder mask | Tests |
|---|---|---|
| No mask (vanilla VAE, LDVAE [13]) | none | reference without knowledge |
| Hallmark | Hallmark gene sets | knowledge + sparsity |
| Random | Hallmark sets with gene labels randomly permuted | sparsity alone |
| Co-expression | modules of co-expressed genes, same sizes | data-derived structure |

The pipeline has eight steps: data (blue), design (orange), model (purple) and evaluation (teal).

![Pipeline in eight steps, coloured by phase: data, design, model, evaluation](reports/figures/pipeline.png)

### Data

The data are from Kang et al. [17] (GEO GSE96583, batch 2). PBMCs from 8 patients with systemic lupus erythematosus
(SLE; GEO lists "subject status: SLE patient" for both samples; called donors below) were cultured for 6 h without (ctrl)
or with interferon-β (stim), and each condition was pooled and sequenced in one 10x run. Every donor appears in both
conditions. The authors' genotype-based singlet calls (demuxlet) and cell-type labels are used as given.

| Step | What is done | Result |
|---|---|---|
| Cell quality control (QC) | cells beyond 5 median absolute deviations in log total counts, log detected genes or share of counts in the top 20 genes are removed, following [18] and computed with scanpy [19], with thresholds per condition × cell type | 24,673 labelled singlets → 23,919 cells (3.1% removed; 0.8–5.0% per cell type, 20.8% of megakaryocytes) |
| Gene filter | genes detected in ≥ 20 cells are kept, and all of them are modelled | 12,034 genes |
| Split | cells are split within every donor × condition × cell type group | 70% train, 15% validation, 15% test |

Three preprocessing choices were made for this data set.

- **QC thresholds per cell type.** Thresholds per condition alone removed 5–7% of CD14+ monocytes, 9–18% of dendritic cells and 20–27% of
  megakaryocytes but under 2% of lymphocytes. The RNA profiles of these cell types differ from those of lymphocytes.
- **All genes instead of highly variable genes.** Selecting 5,000 highly variable genes drops interferon-induced genes
  such as PSMB9 and B2M and a third of the Hallmark interferon-α set.
- **No ambient-RNA correction or mitochondrial filter.** The data set has no empty droplets to estimate ambient RNA from, and mitochondrial genes have
  no counts in the count matrices.

### Gene sets and masks

The Hallmark mask uses the Hallmark collection [20] from MSigDB v2024.1. It has 50 sets of at least 12 genes (median
117.5) and no near-duplicate sets: the largest overlap, between the early and late estrogen responses, is a Jaccard
index of 0.34. Its one unambiguous target for this data is the interferon-α response, with 95 of its 97 genes present.
Gene symbols are matched through Ensembl IDs to current symbols of the HUGO Gene Nomenclature Committee (HGNC) [21].

Two null masks (negative controls) replace the Hallmark sets with sets of the same sizes. Each is built 10 times, and the run with seed k uses
version k.

- **Random mask** (`04_build_masks.py`). Genes are put into 25 bins of equal size (481–482 genes) by their mean
  count in training cells. Gene labels are then permuted within each bin: each gene takes over the set memberships
  of another gene from its bin. Every set keeps its size, its overlaps with other sets and its expression level. Only
  its biology changes. A full shuffle changes 98% of set members; the rest land on a member of their own set.
- **Co-expression mask** (`05_build_coexpression_null.py`). This null tests whether curated sets add anything beyond
  co-expressed genes. It is built from the training cells alone, without condition or cell-type labels. For each seed,
  50 seed genes are drawn at random from the 6,753 genes detected in ≥ 1% of training cells. Module j is seed gene j plus
  the n − 1 genes most positively correlated with it (Pearson correlation of log-transformed, depth-scaled counts), where
  n is the size of Hallmark set j (12–197 genes). A module keeps the size of its set and none of its biology; modules are
  built independently and may share genes. Three properties affect how the results are read.

  - **Coherence.** The median within-set correlation is 0.096 for modules and 0.003 for Hallmark sets, so this null is
    a strong comparison.
  - **Overlap.** Modules cover 2,338–3,233 genes (Hallmark: 3,187). Ribosomal-protein genes and PTMA are in up to 28
    modules per seed; no Hallmark gene is in more than 10 sets.
  - **Interferon.** In every seed, the module most raised by IFN-β is enriched for Hallmark interferon genes (19–89% of
    its genes, against a base rate of 1.8%). The data alone recover the core interferon response.

  The full recipe is in [docs/methods.md](docs/methods.md#co-expression-mask-full-recipe).

### Model

Every variant uses the same VAE. The encoder has one hidden layer (128 units, ReLU, dropout 0.1) and reads the 12,034
genes after a shifted-log transform (below) plus a one-hot code for the 8 donors. The decoder is linear and has no
layers of its own, only the parameters $W$, $V$ and $\log\theta$ (`w`, `v` and `log_theta` in the code). In three of the
four variants a binary mask multiplies $W$. The architecture printout, parameter counts and symbol list are in
[docs/methods.md](docs/methods.md#architecture-and-parameter-counts).

The mask decides how many of the 697,972 weights in `w` can be trained.

| Variant | Trainable decoder weights |
|---|---|
| No mask | 697,972 |
| Hallmark | 76,344 (5,568 set memberships + 8 free latents × 8,847 genes in no set) |
| Random | 76,344 (same sizes and overlaps as the Hallmark mask) |
| Co-expression | 75,976–83,136 (by seed) |

The encoder first puts the counts of a cell on a common scale. Each count $x_g$ is divided by the
total count of the cell $\ell$, multiplied by the median total $m$ of the training cells and log-transformed. Deeply and
shallowly sequenced cells then look alike.

$$\tilde{x}_g = \ln\big(1 + x_g \cdot m / \ell\big)$$

The encoder reads $\tilde{x}$ and the donor code $d$ and proposes where the cell sits in the latent space, as a mean and a
variance for each of the 58 latents. This proposal is the posterior $q(z \mid x, d)$. Training draws a point $z$ from it;
the evaluation uses its mean.

The decoder turns $z$ back into expression. Each latent raises or lowers the genes that the
mask $M$ lets it reach, through its weights $W$, and the donor adds its offsets $V$. A softmax over genes converts these
scores into shares $\pi_g$ that sum to 1, the fraction of the counts of the cell that each gene should receive.

$$\pi = \mathrm{softmax}\big((W \odot M)\,z + V\,d\big)$$

The expected count of gene $g$ is the total count of the cell times its share. The observed count scatters around it as a
negative binomial (NB). The dispersion $\theta_g$ allows each gene to be noisier than a Poisson count.

$$x_g \sim \mathrm{NB}\big(\mu_g = \ell \cdot \pi_g,\ \theta_g\big)$$

The loss has two parts. The reconstruction term measures how unlikely the observed counts are under
the prediction. The Kullback–Leibler (KL) term pulls the posterior of each cell towards a standard normal prior. This
keeps the latents on a common scale and lets those that the data do not need switch off.

$$\mathcal{L} = -\log p_{\mathrm{NB}}(x \mid \mu, \theta) + \beta \cdot \mathrm{KL}\big(q(z \mid x, d)\ \|\ \mathcal{N}(0, I)\big)$$

The weight $\beta$ starts at 0, so the model first learns to reconstruct. It reaches 1 after ~38 epochs.
Symbols are defined in [docs/methods.md](docs/methods.md#symbols-of-the-model).

Four design choices concern the decoder and its inputs.

- **Linear decoder in every variant.** Each latent has one weight for each gene. The mask is the only
  architectural difference between variants.
- **Named and free latents.** In the masked variants, latents 1–50 may only use the genes of their set. The 8 free
  latents may only use the 8,847 genes in no set (74%). When the free latents
  could reach all genes, they absorbed the stimulation signal and the named latents switched off. They are
  kept off the set genes.
- **Donor, not condition.** Donor is a covariate in encoder and decoder. The stimulation is the signal the latents should find, so
  condition is never a covariate.
- **Likelihood.** Negative binomial on raw counts. The library size is the observed total count, and the
  model has no zero inflation [22].

All variants are trained with the same settings.

| Setting | Value |
|---|---|
| Optimiser | Adam |
| Learning rate | 1e-3 |
| Batch size | 128 |
| KL warm-up | $\beta$ rises linearly from 0 to 1 over 10,000 steps at batch 64 [14], counted in cells seen: ~38 epochs at every training size (`--warmup-epochs`) |
| Stopping | early stopping on validation loss with a patience of 3 epochs [14], counted after the KL warm-up; any decrease counts as an improvement (`--min-delta 0`); at most 2,000 epochs (`--max-epochs`) |
| Encoder width | 128, chosen on the no-mask VAE |
| Seeds | 10 per variant; each seed sets the initial weights, the batches and (for the null masks) the mask. The train/validation/test split is fixed (seed 0) |
| Validation loss | computed with the same latent draws at every epoch. Early stopping then follows the model and not the sampling noise |
| Software | versions under Reproducing; the last digits of the losses differ between CPU and GPU |

### Evaluation

#### The interferon-α score

**Rationale.** Type I interferons induce a largely shared set of interferon-stimulated genes (ISGs) in peripheral blood
cells, with differences between cell types [17, 27]. IFN-β is a type I interferon, and the Hallmark interferon-α response
is a curated set of such genes [20]. A latent named after this set should therefore take higher values in stimulated
cells than in control cells of the same type. The stimulation labels are never used in training; they enter only here, as
the known answer. The score itself is our own choice. VEGA [3] and expiMap [4] also use this dataset, with other
measures, and we found no published standard score for a named latent.

The score is the AUROC of that latent. It uses only the ranking of the latent values, so it depends neither on the scale
of the latent, which differs between variants and seeds, nor on a threshold. Its reading is direct: 0.5 is chance, 1 is
perfect separation, and the value is the probability that a random stimulated cell has a higher latent value than a
random control cell [28], the Wilcoxon (Mann–Whitney) statistic scaled to [0, 1]. The AUROC is computed within each cell
type, because cell types differ far more from each other than stimulation does, and a latent that only separates cell
types must not score.

**Definition.**  Let $z_{ik}$ be the posterior mean of latent $k$ for test cell $i$ and $\sigma_k = \pm 1$ its
orientation, and let $s_i = \sigma_k z_{ik}$. For cell type $t$, with stimulated test cells $S_t$ and control test cells
$C_t$:

$$\mathrm{AUROC}_t = \frac{1}{|S_t|\,|C_t|} \sum_{i \in S_t} \sum_{j \in C_t} \Big( \mathbf{1}[s_i > s_j] + \tfrac{1}{2}\,\mathbf{1}[s_i = s_j] \Big)$$

The interferon-α score $A$ is the unweighted mean over the 7 cell types (megakaryocytes excluded). A latent counts as
inactive when the variance of its posterior means over the validation cells is at most 0.01, and then $A = 0.5$:

$$A = \begin{cases} \dfrac{1}{7} \sum_{t=1}^{7} \mathrm{AUROC}_t & \mathrm{Var}_{i \in \mathrm{val}}(z_{ik}) > 0.01 \\\\ 0.5 & \end{cases}$$

The test sets per cell type range from 29 control and 33 stimulated cells (dendritic cells) to 746 and 764 (CD4 T cells).

**Which latent, and which sign.** The sign of a latent is arbitrary, so $\sigma_k$ is fixed without labels where the latent
has a name.

*Hallmark and random mask.* The latent $k$ is the one named `HALLMARK_INTERFERON_ALPHA_RESPONSE` (the same position
holds random genes in the random mask). Its orientation is the sign of the mean decoder weight on its own genes, the
genes $g$ with $M_{gk} = 1$, so that a higher $z_{ik}$ raises those genes. With $G_k$ the set of these genes:

$$\sigma_k = \mathrm{sign}\Big(\frac{1}{|G_k|} \sum_{g \in G_k} W_{gk}\Big)$$

*No mask and co-expression mask.* The latents have no names. Let $A_k$ be the mean over the 7 cell types of the
AUROC of latent $k$ on validation cells. The chosen latent is the active one furthest from 0.5:

$$k^{*} = \arg\max_{k \in \mathrm{active}} \big| A_k - 0.5 \big|$$

Its orientation is $\sigma = +1$ if $A_{k^{*}} > 0.5$ and $\sigma = -1$ otherwise. It is then scored on
test cells. This rule gives these two variants the best of their active latents, where the other two variants get one
fixed latent.

**Evidence for the score.** The named latent tracks the stimulation within every cell type. An independent
recomputation from the saved latent values with the Mann–Whitney U statistic reproduced the scores of four runs (one per
variant) to six decimals. The scoring code is `src/08_evaluate.py`, abridged in
[docs/methods.md](docs/methods.md#scoring-code).

**Limits of the score.** The score does not show that a latent represents interferon biology specifically. Any latent
whose genes respond to the stimulation can separate stimulated from control cells, and the Hallmark interferon-α and
interferon-γ sets share 71 of the 95 interferon-α genes. A shuffled latent holding 3 of 95 interferon-α genes still
reached up to 0.98 in one seed. Near 1 the AUROC saturates and measures ranking, not the size of the effect, and rare
cell types have few test cells.
#### Labels, exclusions and statistics

- **Labels.** Labels never enter training. Validation labels choose latents for the variants without meaningful names
  and test labels are used only for the final scores.
- **Megakaryocytes.** They are mostly platelets, and platelets have no nucleus. Megakaryocytes are excluded from the scores and kept in
  training.
- **Statistics.** Each variant has 10 seeds, paired by seed. The one primary test compares the interferon-α score of
  the Hallmark and random masks. Hallmark beats random if the 95% paired t-interval of the difference excludes 0, and an
  exact two-sided Wilcoxon signed-rank test is reported alongside; with 10 pairs, its smallest possible p is 0.002.
  The target latent, the score and this test were fixed before the 10-seed runs. All other comparisons are secondary and
  not corrected for multiple testing.
- **Robustness checks** (`src/10_robustness.py`). (i) The random mask is also scored like the no-mask VAE, by its best
  active latent chosen on validation cells. (ii) The named latents are scored in each of the 8 donors separately.
  (iii) For the learning curves, the slope of each difference on log2(training size) is fitted per seed.

### Learning curves (Q1b)

The learning curves test whether Hallmark knowledge helps more when cells are scarce. All four variants are refit on
smaller training sets and scored as above.

- **Training sizes.** 500, 1,000, 2,000, 4,000, 8,000 and all ~16,700 training cells, as nested subsamples stratified
  by donor × condition × cell type. Each seed has its own series of subsamples, shared by all four variants. The
  comparison stays paired by seed, and the seeds also cover the choice of cells. Validation and test sets stay fixed,
  and the full-size point is the main experiment.
- **Training.** KL warm-up over ~38 epochs, patience 3 epochs and learning rate 1e-3 at every size, with 10 seeds per
  variant and size.
- **Nulls.** Co-expression modules and shuffle expression bins are rebuilt on each subsample. With few cells,
  correlations are noisier: modules overlap less and cover ~4,100 genes at 500 cells.

## Results

All 40 runs (4 variants × 10 seeds) stopped early. The interferon-α latent was active in every run. Values are
mean ± SD over 10 seeds; differences are paired by seed.

| Variant | Interferon-α score | Active latents out of 58 | Validation loss per cell |
|---|---|---|---|
| No mask | 0.966 ± 0.025 | 34–56 | 1,920 |
| Hallmark | 0.997 ± 0.002 | 36–38 | 1,963 |
| Random | 0.758 ± 0.218 | 33–37 | 1,966 |
| Co-expression | 0.998 ± 0.0004 | 28–37 | 1,951 |

![Interferon-α score and active latents per run, and paired differences](reports/figures/scores.png)

*(a) One dot per seed, bar at the mean; grey lines join Hallmark and random runs of the same seed; hollow dots would mark
inactive latents (none here). (b) Active latents per run: a latent is active if the variance of its posterior means across
validation cells exceeds 0.01. (c) Per-seed paired differences (dots), their mean and 95% t-interval; the primary
comparison, Hallmark−random, in black.*

| Comparison | Interferon-α score |
|---|---|
| Hallmark−random | +0.239 (+0.083 to +0.395), p = 0.002 |
| Hallmark−no mask | +0.031 (+0.013 to +0.049), p = 0.002 |
| Hallmark−co-expression | −0.0009 (−0.0020 to +0.0002), p = 0.04 |

Mean difference, 95% paired t-interval and Wilcoxon signed-rank p over 10 seeds.

### Q1: knowledge or sparsity?

- **The Hallmark mask places the response in its named latent.** The interferon-α latent separated stimulated from
  control test cells almost perfectly (0.997 ± 0.002). The random mask has the same set sizes, overlaps and expression
  levels. Its latent at the same position scored 0.758 ± 0.218 (0.38 to 0.98 by seed). The Hallmark mask scored
  higher in all 10 seeds, by +0.239 (p = 0.002). In the lowest random-mask seed the latent's sign pointed away from the
  response and ranked control cells above stimulated ones.
- **The random mask encodes the response in a latent that reaches unrelated genes.** Scored like the no-mask VAE, by its best
  active latent, the random mask reached 0.989 ± 0.005 (0.976 to 0.994 by seed). In none of the 10 seeds was that latent
  the interferon-α one. The Hallmark mask was still higher in all 10 seeds, by +0.008 (+0.005 to +0.012, p = 0.002). Its best
  latent was the interferon-α latent in 9 of 10 seeds. The sparsity of the mask is enough to carry the
  response; the genes in the sets decide which latent carries it and what that latent represents.
- **A co-expression mask scores very similarly to the Hallmark mask.** Its best module scored 0.998 ± 0.0004 against 0.997
  ± 0.002 for the Hallmark latent, a difference of −0.0009 (95% interval −0.0020 to +0.0002; Wilcoxon p = 0.04,
  uncorrected). Modules built from the data alone, without labels, capture the response about as well as curated sets do. Curation adds a known meaning: the Hallmark latent reaches the interferon genes and is known before
  training, while the module and the no-mask latent have to be found with validation labels.
- **The Hallmark mask beats the no-mask VAE.** Its best latent scored 0.966 ± 0.025. The Hallmark mask was higher in all 10 seeds
  (+0.031).
- **The result holds in every cell type.** The Hallmark mask reached at least 0.991 in each of the 7 cell types (NK
  cells lowest). The random mask averaged 0.73 to 0.84 in each.
- **The result holds in every donor.** The named latent scored 0.993 or higher in each of the 8 donors for the Hallmark
  mask and 0.73 to 0.78 for the random mask. The Hallmark−random difference was +0.22 to +0.27 in every
  donor. The test cells come from the same 8 donors as the training cells. This does not show generalisation to new
  donors.
- **Masks cost reconstruction.** Validation loss per cell was 1,920 for the no-mask VAE and 1,951–1,966 for the
  masked variants. The masked variants kept 34–37 of 58 latents active on average, against 47.
- **Relation to prediction networks.** There, randomized models matched or outperformed informed ones [1]. In the
  unsupervised models here, the masks agree on the signal and differ in whether the latent that carries the response has a known biological meaning.

![AUROC per cell type](reports/figures/cell_types.png)

*Small dots: single seeds; large dots: means over 10 seeds.*

![Validation loss during training](reports/figures/training.png)

*One line per seed; dots mark the best epoch; shading marks the KL warm-up.*

### Q1b: learning curves

Interferon-α score, mean over 10 seeds at each training size.

| Training cells | No mask | Hallmark | Random | Co-expression | Hallmark−random |
|---|---|---|---|---|---|
| 500 | 0.992 | 0.996 | 0.482 | 0.993 | +0.514 |
| 1,000 | 0.987 | 0.997 | 0.636 | 0.996 | +0.362 |
| 2,000 | 0.982 | 0.997 | 0.610 | 0.997 | +0.387 |
| 4,000 | 0.975 | 0.997 | 0.672 | 0.997 | +0.325 |
| 8,000 | 0.972 | 0.997 | 0.552 | 0.997 | +0.445 |
| 16,749 | 0.966 | 0.997 | 0.758 | 0.998 | +0.239 |

![Learning curves](reports/figures/learning_curves.png)

*(a) Mean over 10 seeds with a 95% t-interval band. (b) Per-seed Hallmark−random differences (dots), mean and 95%
t-interval.*

- **On this target, Hallmark knowledge does not make learning more data-efficient.** With 500 cells the Hallmark mask already scored
  0.996. The co-expression mask (0.993) and the no-mask VAE (0.992) were within 0.004 of it. The Hallmark advantage was
  already small at 500 cells.
- **A co-expression mask scores very similarly to the Hallmark mask at every size.** The mean difference (Hallmark −
  co-expression) was +0.003 at 500 cells (95% interval −0.001 to +0.007) and between +0.002 and −0.001 at larger sizes;
  from 1,000 cells up every 95% interval lies between −0.002 and +0.004. At every size the Hallmark mask was ahead in
  only 3 to 7 of 10 seeds.
- **The named latent of the random mask stays unreliable at every size.** Its mean score was 0.48 to 0.76, and single seeds ranged from 0.01 to
  0.99. The Hallmark mask beat it in all 10 seeds at every size, by +0.24 to +0.51: most at 500 cells (+0.51) and least with all cells (+0.24). The fitted slope was −0.034 per doubling of
  the training set (95% interval −0.074 to +0.007, Wilcoxon p = 0.11). The advantage over the random mask is not shown
  to depend on the training size.
- **The no-mask VAE scores lower as cells are added.** Its score fell from 0.992 to 0.966. The Hallmark advantage
  grew from +0.004 to +0.031 (slope +0.005 per doubling, 95% interval +0.003 to +0.008) and was positive in
  all 10 seeds at each size.
- **Small training sets overfit.** Validation loss exceeded training loss by 40–90 per cell (2–4%) at 500 to 2,000
  cells and by 15–20 at 4,000, and was near zero from 8,000 cells.

### Robustness checks

Checks on the main results (`src/10_robustness.py`; tables in `reports/tables/robustness_*.csv`).

| Check | Result |
|---|---|
| Random mask scored on its best latent, as the no-mask VAE is | 0.989 ± 0.005 (0.976 to 0.994); the interferon-α latent in 0 of 10 seeds |
| Hallmark scored on its best latent | 0.997 ± 0.001; the interferon-α latent in 9 of 10 seeds |
| Hallmark − random, both on the best latent | +0.008 (+0.005 to +0.012), p = 0.002; higher in 10 of 10 seeds |
| Hallmark, per donor (mean over 10 seeds) | 0.993 to 1.000 in each of the 8 donors |
| Random mask, per donor | 0.733 to 0.784 |
| Hallmark − random, per donor | +0.215 to +0.267 |
| Slope per doubling of training cells, Hallmark − random | −0.034 (−0.074 to +0.007), p = 0.11 |
| Slope per doubling of training cells, Hallmark − no mask | +0.005 (+0.003 to +0.008), p = 0.002 |
| Slope per doubling of training cells, Hallmark − co-expression | −0.0007 (−0.0013 to −0.0001), p = 0.01 |

## Limitations

- **Scope.** One dataset, one stimulus and one target: PBMCs from 8 donors (patients with SLE, whose type I interferon
  activity is often raised at baseline [29]) stimulated with IFN-β, Hallmark as the only gene set collection and the
  interferon-α response as the only target with a known answer. The signal is very strong. The Hallmark, co-expression
  and (best-latent) random masks score between 0.989 and 0.998, and the data cannot separate masks that reach this
  ceiling. Whether the other latents are correct is not tested (Q2 and Q3 in Future work).
- **Seeds.** The 10 seeds vary the training and one draw of each null mask. The split, the 8 donors and the cells are
  fixed, the test cells share their donors with the training cells and no donor was held out, so the intervals say
  nothing about new donors. The SD of the random mask (0.22) mixes the draw of the mask with the training noise.
- **Scoring and tests.** The named latent is scored for the Hallmark and random masks; the best active latent, chosen
  with validation labels, for the no-mask and co-expression variants. This favours the last two (the random mask is also
  scored that way as a robustness check). Only Hallmark against random was a primary comparison. All others are
  uncorrected for multiple testing.
- **Target of the score.** Stimulation, not interferon biology specifically (see Evaluation). The co-expression null
  also recovers the response, so it cannot show what curation adds beyond a known meaning for each latent.
- **Hard mask and linear decoder only.** Soft masks, nonlinear decoders, gene set scores as input features and graph
  models were not tested. The 8 free latents reach only genes in no set; with access to all genes, they absorbed the
  response and the named latents switched off, so this design choice shapes the result.
- **Shared settings.** Learning rate (1e-3), encoder width (128, chosen on the no-mask VAE) and patience (3 epochs) are
  the same for all variants and were not tuned for each. Other settings could favour a different variant.
- **Preprocessing.** The 2017 gene annotation (hg19), the cell-type labels and the demuxlet calls come from the authors
  of the data. Doublets of two cells from the same donor are not detected (an estimated 2% of cells), and ambient RNA is
  not corrected.

## Reproducing

The scripts are numbered in run order and are run from the repository root in a fresh virtual environment.

| Component | Version |
|---|---|
| Python | 3.13.15 |
| PyTorch | 2.14.1 (Apple GPU through MPS) |
| scanpy | 1.12.4 |
| anndata | 0.13.4 |
| NumPy | 2.5.3 |
| pandas | 3.0.6 |
| SciPy | 1.18.1 |
| scikit-learn | 1.9.1 |
| Matplotlib | 3.11.2 |
| scikit-misc | 0.5.2 |
| tqdm | 4.70.1 |

All results were produced on an Apple M1 Pro (16 GB) under macOS 15.7. `requirements.txt` sets minimum versions only.

```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python src/01_download_kang.py
python src/02_download_genesets.py
python src/03_preprocess_kang.py
python src/04_build_masks.py
python src/05_build_coexpression_null.py
bash src/07_experiment.sh
python src/08_evaluate.py
python src/09_figures.py
```

`src/06_fit.py` trains one model per call and is run by `07_experiment.sh`. Its options and defaults:

| Option | Default | Meaning |
|---|---|---|
| `--mask` | `none` | `none`, `real`, `shuffled` or `coexpression` |
| `--seed` | `0` | initial weights, batches and (for the null masks) the mask |
| `--train-size` | all training cells | number of training cells, for the learning curves |
| `--max-epochs` | `2000` | safety ceiling; a run that reaches it fails |
| `--warmup-epochs` | ~38.2 | epochs for β to rise from 0 to 1; the default is 10,000 steps at batch 64 [14], scaled to the training set |
| `--min-delta` | `0` | relative drop in validation loss that counts as an improvement; `0` is the setting of all reported results, and a value such as `1e-3` stops earlier |
| `--results-dir` | `results` | output folder |
| `--verbose` | off | progress bar over epochs; otherwise one line per run |

The defaults reproduce the reported results. Four runs (no mask and Hallmark at full size, shuffled and co-expression at
500 cells) matched the stored runs in every epoch and in the saved latent means.

The learning curves run after the main experiment and add their own tables and figure.

```
python src/learning_curves/01_subsample.py
bash src/learning_curves/02_run.sh
python src/08_evaluate.py --learning-curves
python src/09_figures.py
python src/10_robustness.py
```

The last script reads the tables of both experiments and adds the robustness checks.

| Path | Contents |
|---|---|
| `src/` | numbered scripts in run order; `06_fit.py` trains one model per call and is run by `07_experiment.sh` |
| `src/learning_curves/` | training subsamples and the runs on them (Q1b) |
| `notebooks/hvg_inspection.py` | highly-variable-gene analysis behind the choice to model all genes |
| `data/`, `results/` | created by the scripts, not tracked |
| `reports/tables/` | scores per run, AUROCs per cell type, paired comparisons and robustness checks (`robustness_*.csv`) |
| `reports/figures/` | figures |
| `docs/methods.md` | co-expression recipe, architecture and symbols, scoring code, biological choices by pipeline step, candidate datasets |

## Future work

The same comparison could be repeated on harder questions and other data. Candidate datasets are listed in
[docs/methods.md](docs/methods.md#candidate-datasets-for-future-work).

- **Q2 (known regulators).** Do the explanations recover known regulators? B cell differentiation with IRF4 and PRDM1 perturbations offers such ground truth [23].
- **Q3 (protein as reference).** Does RNA-based pathway activity agree with CITE-seq surface protein, an independent
  measurement? 
- **Q4 (imperfect gene sets).** What is the cost of incomplete or species-transferred gene sets, for example human sets
  applied to mouse gastrulation data [25]? And does a soft mask started from random sets refine to the same
  gene sets?
- **Reactome.** Reactome [26] as a second gene set collection on the Kang data, as used there by VEGA [3].

## References
1. Caranzano I, et al. Sparsity is all you need: rethinking biologically informed neural networks. *Brief Bioinform*
   (2026). [doi:10.1093/bib/bbag425](https://doi.org/10.1093/bib/bbag425)
2. Thapa K, et al. Strategies to include prior knowledge in omics analysis with deep neural networks. *Patterns*
   (2025). [doi:10.1016/j.patter.2025.101203](https://doi.org/10.1016/j.patter.2025.101203)
3. Seninge L, et al. VEGA is an interpretable generative model for inferring biological network activity in
   single-cell transcriptomics. *Nat Commun* (2021).
   [doi:10.1038/s41467-021-26017-0](https://doi.org/10.1038/s41467-021-26017-0)
4. Lotfollahi M, et al. Biologically informed deep learning to query gene programs in single-cell atlases.
   *Nat Cell Biol* (2023). [doi:10.1038/s41556-022-01072-x](https://doi.org/10.1038/s41556-022-01072-x)
5. Gut G, et al. pmVAE: learning interpretable single-cell representations with pathway modules. *bioRxiv* (2021).
   [doi:10.1101/2021.01.28.428664](https://doi.org/10.1101/2021.01.28.428664)
6. Buettner F, et al. f-scLVM: scalable and versatile factor analysis for single-cell RNA-seq. *Genome Biol* (2017).
   [doi:10.1186/s13059-017-1334-8](https://doi.org/10.1186/s13059-017-1334-8)
7. Mao W, et al. Pathway-level information extractor (PLIER) for gene expression data. *Nat Methods* (2019).
   [doi:10.1038/s41592-019-0456-1](https://doi.org/10.1038/s41592-019-0456-1)
8. Kunes RZ, et al. Supervised discovery of interpretable gene programs from single-cell data. *Nat Biotechnol*
   (2024). [doi:10.1038/s41587-023-01940-3](https://doi.org/10.1038/s41587-023-01940-3)
9. Fortelny N, Bock C. Knowledge-primed neural networks enable biologically interpretable deep learning on
   single-cell sequencing data. *Genome Biol* (2020).
   [doi:10.1186/s13059-020-02100-5](https://doi.org/10.1186/s13059-020-02100-5)
10. Qoku A, et al. MOFA-FLEX: a factor model framework for integrating omics data with prior knowledge. *bioRxiv*
    (2025). [doi:10.1101/2025.11.03.686250](https://doi.org/10.1101/2025.11.03.686250)
11. Moullet M, et al. Self-supervised learning for a gene program-centric view of cell states. *bioRxiv* (2026).
    [doi:10.64898/2026.03.24.713961](https://doi.org/10.64898/2026.03.24.713961)
12. Zhu Y, et al. A residual-ratio framework for auditing transcriptomic gene signatures against background
    expression structure. *bioRxiv* (2026). [doi:10.64898/2026.04.11.717907](https://doi.org/10.64898/2026.04.11.717907)
13. Svensson V, et al. Interpretable factor models of single-cell RNA-seq via variational autoencoders.
    *Bioinformatics* (2020). [doi:10.1093/bioinformatics/btaa169](https://doi.org/10.1093/bioinformatics/btaa169)
14. Eltager M, et al. Benchmarking variational autoencoders on cancer transcriptomics data. *PLoS One* (2023).
    [doi:10.1371/journal.pone.0292126](https://doi.org/10.1371/journal.pone.0292126)
15. Makrodimitris S, et al. An in-depth comparison of linear and non-linear joint embedding methods for bulk and
    single-cell multi-omics. *Brief Bioinform* (2024). [doi:10.1093/bib/bbad416](https://doi.org/10.1093/bib/bbad416)
16. Ahlmann-Eltze C, et al. Deep-learning-based gene perturbation effect prediction does not yet outperform simple
    linear baselines. *Nat Methods* (2025). [doi:10.1038/s41592-025-02772-6](https://doi.org/10.1038/s41592-025-02772-6)
17. Kang HM, et al. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. *Nat Biotechnol*
    (2018). [doi:10.1038/nbt.4042](https://doi.org/10.1038/nbt.4042)
18. Heumos L, et al. Best practices for single-cell analysis across modalities. *Nat Rev Genet* (2023).
    [doi:10.1038/s41576-023-00586-w](https://doi.org/10.1038/s41576-023-00586-w)
19. Wolf FA, et al. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* (2018).
    [doi:10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)
20. Liberzon A, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst* (2015).
    [doi:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)
21. Seal RL, et al. Genenames.org: the HGNC resources in 2023. *Nucleic Acids Res* (2023).
    [doi:10.1093/nar/gkac888](https://doi.org/10.1093/nar/gkac888)
22. Svensson V. Droplet scRNA-seq is not zero-inflated. *Nat Biotechnol* (2020).
    [doi:10.1038/s41587-019-0379-5](https://doi.org/10.1038/s41587-019-0379-5)
23. Demela P, et al. Competing gene regulatory networks drive naive and memory B cell differentiation. *Mol Syst Biol*
    (2026). [doi:10.1038/s44320-026-00207-8](https://doi.org/10.1038/s44320-026-00207-8)
24. Massoni-Badosa R, et al. An atlas of cells in the human tonsil. *Immunity* (2024).
    [doi:10.1016/j.immuni.2024.01.006](https://doi.org/10.1016/j.immuni.2024.01.006)
25. Pijuan-Sala B, et al. A single-cell molecular map of mouse gastrulation and early organogenesis. *Nature* (2019).
    [doi:10.1038/s41586-019-0933-9](https://doi.org/10.1038/s41586-019-0933-9)
26. Milacic M, et al. The Reactome Pathway Knowledgebase 2024. *Nucleic Acids Res* (2024).
    [doi:10.1093/nar/gkad1025](https://doi.org/10.1093/nar/gkad1025)
27. Waddell SJ, et al. Dissecting interferon-induced transcriptional programs in human peripheral blood cells.
    *PLoS One* (2010). [doi:10.1371/journal.pone.0009753](https://doi.org/10.1371/journal.pone.0009753)
28. Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) curve.
    *Radiology* (1982). [doi:10.1148/radiology.143.1.7063747](https://doi.org/10.1148/radiology.143.1.7063747)
29. Sim TM, et al. Type I interferons in systemic lupus erythematosus: a journey from bench to bedside. *Int J Mol Sci*
    (2022). [doi:10.3390/ijms23052505](https://doi.org/10.3390/ijms23052505)
