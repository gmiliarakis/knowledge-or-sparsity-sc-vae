# Supplementary results

Results that support the main report ([README](../README.md)) but are not part of its argument. Statistical tests are described
there; all tests are in `reports/tables/tests.csv`, and only adjusted p-values are quoted.

## Statistics

The tests follow the procedure of the [README](../README.md#statistics). Two further families belong to the supplement: the 10 tests on the interferon-α score at full size (Holm; family-wise error 0.05; intervals at the Bonferroni level) and the 18 learning-curve tests (Benjamini–Hochberg; false discovery rate 0.05, valid for the positive dependence between sizes that share seeds, with false-coverage-rate intervals). One family of all 40 tests is impossible, because Holm over more than 25 tests cannot reject at p ≥ 0.002. Equivalence of the Hallmark and co-expression masks on the interferon-α score is tested with two one-sided tests against ±0.005 AUROC (half a percentage point in the probability that a stimulated cell outranks a control cell); the margin predates the within-set random mask. The 12 tests of S4 form their own family (Holm).

## S1. The interferon-α latent

The score of a latent named after the Hallmark interferon-α response set, and the reference scores without a VAE.

### Definition

A latent named after the Hallmark interferon-α response set should take higher values in stimulated cells, because IFN-β is a type I interferon and these genes are induced in blood cells [17, 20, 27]. The score is the AUROC of the latent's posterior mean for stimulated against control test cells within each cell type, averaged over the 7 types (0.5 is chance). For the Hallmark and both random masks the latent at the interferon-α position is scored, oriented without labels by the sign of its mean decoder weight; it scores 0.5 when inactive (variance of its posterior means over validation cells at most 0.01). The no-mask and co-expression variants have no names, so their active latent with the largest |AUROC − 0.5| on validation cells is scored on test cells (the best latent), which uses validation labels and favours these two. The score shows where the response sits and is not specific to interferon biology: any latent whose genes respond to the stimulation can separate the conditions, the Hallmark interferon-α and -γ sets share 71 of 95 genes, and the score saturates near 1. Scoring code: `src/08_evaluate.py`, abridged in [methods.md](methods.md#scoring-code).

**Reference scores without a VAE:** the mean log expression of the interferon-α genes (no labels) and the best of 58 principal components of the training cells (chosen on validation labels), scored like the latents. The best principal component is also refitted after removing the 216 genes of the Hallmark interferon-α and -γ sets, or the 100 to 2,000 genes with the largest stimulation effect in training cells.

### Results

The interferon-α latent of the Hallmark mask scored 0.997 ± 0.001 and was its best latent in 9 of 10 seeds (the interferon-γ latent in the tenth). It scored higher than the latent at the same position of the random masks (0.758 ± 0.218 and 0.815 ± 0.281) in all 10 seeds (adjusted p = 0.018 for both). The random masks carry the response in another latent (best latent 0.989 ± 0.005 and 0.994 ± 0.002, the interferon-α latent in 0 and 3 of 10 seeds), so the score does not show that a latent is specific: a random set of about 100 genes scored 0.93 or higher in 6 of 10 within-set random seeds. Against co-expression the Hallmark latent is equivalent within ±0.005 (difference −0.0006, adjusted p = 0.010), and every variant sits at the ceiling: the best principal component scores 0.997 and the mean of the interferon-α genes 0.996 without any VAE. After removing the 216 genes of the interferon-α and -γ sets the best component still scores 0.945, and after removing the 500, 1,000 or 2,000 strongest stimulation genes 0.744, 0.649 or 0.573. The separation rests on a few hundred genes, and because each condition was sequenced in its own 10x run, run effects among them cannot be excluded.

![Interferon-α score, active latents and paired differences](../reports/figures/scores.png)

*a, interferon-α score of every seed (lines join the Hallmark and random scores of the same seed). b, latents with an active posterior. c, paired difference of Hallmark minus each other variant, with the median and the Holm-adjusted interval; the primary comparison is in black.*

## S2. Cell-type probe, few-shot and hard pairs

Balanced accuracy of the cell-type probe with 10, 50 and 200 validation cells per type, and the AUROC between CD4 and CD8 T cells and between
CD14+ and FCGR3A+ monocytes (mean ± SD over 10 seeds; the main report gives the probe on all validation cells). The probe is described under
Evaluation in the main report; types with fewer than k validation cells contribute all of them.

| Cell-type probe, descriptive | 10 per type | 50 per type | 200 per type | CD4 vs CD8 T cells | CD14+ vs FCGR3A+ monocytes |
|---|---|---|---|---|---|
| No mask | 0.815 ± 0.028 | 0.865 ± 0.008 | 0.884 ± 0.009 | 0.992 ± 0.001 | 0.962 ± 0.003 |
| Hallmark | 0.806 ± 0.023 | 0.839 ± 0.010 | 0.858 ± 0.003 | 0.984 ± 0.001 | 0.950 ± 0.004 |
| Fully random | 0.811 ± 0.026 | 0.840 ± 0.011 | 0.860 ± 0.008 | 0.986 ± 0.002 | 0.939 ± 0.005 |
| Within-set random | 0.805 ± 0.025 | 0.841 ± 0.012 | 0.857 ± 0.010 | 0.983 ± 0.003 | 0.949 ± 0.006 |
| Co-expression | 0.827 ± 0.017 | 0.854 ± 0.009 | 0.867 ± 0.009 | 0.991 ± 0.002 | 0.949 ± 0.009 |

## S3. Learning curves (Q1b)

### Design

The learning curves test whether Hallmark knowledge helps more when cells are scarce. The four variants other than the within-set random mask are refit on
smaller training sets and scored by the interferon-α latent score (S1).

- **Training sizes.** 500, 1,000, 2,000, 4,000, 8,000 and all ~16,700 training cells, as nested subsamples stratified
  by donor × condition × cell type. Each seed has its own series of subsamples, shared by these four variants. The
  comparison stays paired by seed, and the seeds also cover the choice of cells. Validation and test sets stay fixed,
  and the full-size point is the main experiment.
- **Training.** KL warm-up over ≈ 38 epochs, patience 3 epochs and learning rate 1e-3 at every size, with 10 seeds per
  variant and size.
- **Nulls.** Groups of co-expressed genes and shuffle expression bins are rebuilt on each subsample. With few cells,
  correlations are noisier: the groups overlap less and cover ~4,100 genes at 500 cells.
- **Random masks.** Only the fully random mask was run on subsamples.

### Results

Interferon-α score, mean over 10 seeds at each training size.

| Training cells | No mask | Hallmark | Random | Co-expression | Hallmark−random |
|---|---|---|---|---|---|
| 500 | 0.992 | 0.996 | 0.482 | 0.993 | +0.514 |
| 1,000 | 0.987 | 0.997 | 0.636 | 0.996 | +0.362 |
| 2,000 | 0.982 | 0.997 | 0.610 | 0.997 | +0.387 |
| 4,000 | 0.975 | 0.997 | 0.672 | 0.997 | +0.325 |
| 8,000 | 0.972 | 0.997 | 0.552 | 0.997 | +0.445 |
| 16,749 | 0.966 | 0.997 | 0.758 | 0.998 | +0.239 |

![Learning curves](../reports/figures/learning_curves.png)

*(a) Mean over 10 seeds with a 95% t-interval band. (b) Per-seed Hallmark−random differences (dots), their median and the interval at the family level.*

- No data-efficiency benefit of Hallmark knowledge could be detected, and this task cannot detect one: with 500 cells the no-mask VAE already scores 0.992 and the Hallmark mask 0.996. Hallmark and co-expression do not differ at any size (adjusted p 0.089 to 1.0).
- The Hallmark latent beats the random mask's latent at every size (all 10 seeds, adjusted p = 0.003), and this advantage is not shown to depend on training size (slope adjusted p = 0.14). The random mask's latent is unreliable at every size (single seeds 0.01 to 0.99).

## S4. Dropping the 250 strongest stimulation genes

This analysis asks whether the Hallmark advantage over the random masks rests on the strongest stimulation genes. It drops them and repeats the comparison for test loss and decoded response. Dropping them also shrinks the Hallmark sets (32 of the 95 interferon-α genes remain), so the benchmark shows as well what happens when a gene set is incomplete.

### Design

The 250 genes with the largest |t| of stimulated against control log expression in training cells (megakaryocytes excluded) are dropped, which leaves 11,784 genes (list in `reports/tables/removed_genes_drop250.csv`). They are chosen from the data and never from the gene sets. The 50 Hallmark sets, both random masks and the groups of co-expressed genes are rebuilt on the remaining genes, and 3,032 genes are in at least one set. Architecture, training, scoring and the four measures are unchanged, with 10 seeds per variant (seeds 0 and 1 were also run as a pilot before the decision to run all 10). The criterion was fixed before the 10-seed run: the sign and the Holm-adjusted p of the contrasts of Hallmark against the random masks for test loss and decoded response. Without any VAE, the best principal component scores 0.864 and the mean of the 32 remaining interferon-α genes 0.724.

### Results

| Variant | Test loss per cell | Decoded response | Cell-type balanced accuracy | Stimulation probe AUROC |
|---|---|---|---|---|
| No mask | 1718.2 ± 1.4 | 0.766 ± 0.002 | 0.881 ± 0.004 | 0.9877 ± 0.0009 |
| Hallmark | 1747.6 ± 0.7 | 0.528 ± 0.004 | 0.857 ± 0.007 | 0.9723 ± 0.0015 |
| Fully random | 1747.7 ± 2.0 | 0.522 ± 0.020 | 0.856 ± 0.006 | 0.9742 ± 0.0028 |
| Within-set random | 1749.7 ± 1.0 | 0.506 ± 0.005 | 0.850 ± 0.006 | 0.9756 ± 0.0014 |
| Co-expression | 1739.6 ± 1.6 | 0.552 ± 0.019 | 0.864 ± 0.005 | 0.9720 ± 0.0026 |

| Measure | Hallmark minus | Median difference | Interval | Adjusted p |
|---|---|---|---|---|
| Test loss per cell (lower is better) | fully random | +0.3 | −3.6 to +2.0 | 1.000 |
| | within-set random | −2.1 | −3.3 to −1.2 | 0.023 |
| | co-expression | +8.1 | +5.4 to +10.3 | 0.023 |
| Decoded response (correlation) | fully random | +0.006 | −0.020 to +0.035 | 1.000 |
| | within-set random | +0.022 | +0.016 to +0.031 | 0.023 |
| | co-expression | −0.023 | −0.048 to −0.002 | 0.023 |

Against the fully random mask the Hallmark advantage is gone: the difference in test loss is +0.3 units and in decoded response +0.006, against −3.7 and +0.087 with all genes. Against the within-set random mask the advantage persists at about one fifth of its size (−2.1 units against −10.1, and +0.022 against +0.129). The groups of co-expressed genes still fit better than the Hallmark sets (8.1 units) and reproduce the shift better (0.023 higher). The decoded response falls for every masked variant to 0.51 to 0.55, against 0.77 without a mask. The named interferon-α latent was inactive in all Hallmark and within-set random runs, and the best latent of the Hallmark and random masks scored 0.879 to 0.904 on average. Cell-type accuracy shows no difference (adjusted p of at least 0.26), and the stimulation probe is lower for Hallmark than for the within-set random mask by 0.003 (adjusted p 0.031) and does not differ otherwise.
