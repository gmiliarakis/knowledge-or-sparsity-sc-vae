# Does prior knowledge help interpretable single-cell models?

Small public-data project (Oct 2026). Knowledge-guided single-cell models wire gene sets (pathways) into their
architecture so that each latent factor reads as a biological programme. Caranzano et al. 2026
(Brief Bioinform, doi:10.1093/bib/bbag425) showed that for 29 *supervised* pathway-informed networks, structure-matched
random pathways do as well as real ones. This project asks the same question for *unsupervised* single-cell models,
where it has not been tested.

## Questions
- **Q1 (core):** does the knowledge itself help, or only the sparsity a gene-set mask imposes? A masked VAE with the
  real Hallmark mask against the same model with a 100% shuffled mask (same set sizes and overlaps, wrong genes),
  plus a plain VAE as reference and a co-expression null (groups the data finds by itself). Known answer: IFN-β
  switches on the type I interferon programme.
- **Q1b:** is learning more data-efficient with real knowledge? Black-box models are argued to be inefficient because they must re-discover well-known patterns from
  scratch. Q1 repeated at smaller training sizes.
- **Later:** Q2 (do explanations recover known B cell regulators?) and Q3 (do they agree with CITE-seq protein?),
  each on a new dataset.

## Models
One architecture: encoder (shifted-log counts + donor → one hidden layer of 128 → 58 latents), linear decoder,
negative binomial likelihood on raw counts. Variants: plain VAE (all latents reach all genes; this is LDVAE); masked
VAE with the real mask, a 100% shuffled mask, or co-expression modules (50 named latents, one per Hallmark set,
plus 8 free latents reaching only genes in no set). Donor is a covariate.

## Evaluation
Labels never enter training. Validation labels pick the latent (and sign) for the plain VAE and co-expression
models and fit the score B/C regressions; test labels are used only for the final scores. Score A (primary): AUROC
of stim vs ctrl for the named interferon-α latent, within each cell type (megakaryocytes excluded), sign set by the
latent's own weights (mean weight on its own genes positive). Score B (secondary): the interferon-α and -γ latents
together. Score C: all latents. Plain VAE and co-expression modules (no meaningful names): the best latent (B: the
best two) of all 58 and its sign chosen on validation cells, then scored on test cells. An inactive latent (variance
of its means across validation cells ≤ 0.01) scores 0.5 in A; each run's activity is reported next to its score.
10 seeds per variant, compared pair by pair. Primary test: score A, real vs shuffled; real beats shuffled if the 95%
paired t-interval of the difference excludes 0 (Wilcoxon signed-rank alongside). Other comparisons are secondary.

## Data
Kang et al. 2018 IFN-β-stimulated PBMCs (GSE96583, batch 2). Preprocessing follows Heumos et al., single-cell best
practices.

## Layout
- `src/`: one script per stage, numbered in run order:
  `01_download_kang.py`, `02_download_genesets.py`, `03_preprocess_kang.py`, `04_build_masks.py`,
  `05_build_coexpression_null.py`, `06_fit.py` (one model run per call), `07_tune.sh` (learning-rate tuning),
  `08_choose_settings.py` (picks the rate per variant), `09_experiment.sh` (10 seeds per variant)
- `notebooks/hvg_inspection.py`: evidence for using all genes; needs the outputs of 02 and 03
- `notebooks/`: exploration and figures
- `data/raw`, `data/processed`: not in git; recreated by the scripts
- `results/`: one folder per model run (not in git for now)
- `reports/figures/`: output figures

## Setup
```
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Decisions

**General**
- Knowledge can enter models as features,
  architecture or graphs (Thapa et al. 2025), and the masked decoder is used because it is the most direct form of
  "interpretable by design".
- Python (scanpy, PyTorch, scikit-learn); own PyTorch models rather than packages. Preprocessing follows the
  sc-best-practices book (Heumos et al.) unless stated below. MSigDB pinned at v2024.1.
- One plain script per stage; download separate from processing; checks are printed by the code on every run.

**Data**
- Batch 2 of GSE96583 only (ctrl vs 6 h IFN-β, same 8 donors pooled in each 10x run); batch 1 is unstimulated lupus
  samples. Starting from the authors' Cell Ranger 1.2 / hg19 count matrices (SRA has only 2 × ~23 GB BAMs, and the
  donor, cell-type and doublet labels exist only in the GEO metadata). Only the six batch-2 files are downloaded;
  `GSE96583_RAW.tar` also bundles batch 1.
- Genes indexed by Ensembl ID: 2,697 gene symbols are repeats in this annotation.
- 313 barcodes occur in both runs; the authors' metadata names the stim copies `-11`, and we follow it.

**Cell QC**
- Doublets: demuxlet calls (genotype-based) instead of a computational doublet detector. Same-donor doublets
  (~2%, estimated) are not caught.
- No ambient-RNA correction: it needs the empty droplets, and GEO has only the cell-called matrices.
- No mitochondrial filter: the 13 MT- genes have zero counts in the deposited matrices.
- Not done: re-clustering or re-annotation (authors' labels used), scaling / regress-out, cell-cycle correction.
- Outliers at 5 MADs (sc-best-practices) on log total counts, log detected genes and % counts in the top 20 genes.
  Thresholds per condition × cell type: with thresholds per condition only, normal monocytes (much RNA in a few
  genes), dendritic cells (large, many genes) and megakaryocytes were removed as damaged (up to 27% of a type).
  Per donor as well was rejected: groups of 4–20 cells give unstable MADs, and donors shared each 10x run.
- Genes kept if detected in ≥ 20 cells.
- Result: 24,673 labelled singlets → 23,919 cells (3.1% removed; 1–6% per cell type), 12,034 genes.
- Megakaryocytes still lose ~20% with their own thresholds, so they are likely a mixed group (mostly platelets).
  Platelets have no nucleus and should not mount a transcriptional IFN response, so they are excluded from the
  interferon score and kept in training.

**Split**
- 70 / 15 / 15 train / val / test, made inside every donor × condition × cell type group (seed 0), so every group
  is in all three sets; at least one val and one test cell per group (smallest group: 4 cells).
- Fitted on train, tuned and latent-matched on val, scored once on test. Learning curves subsample train only.

**Genes and likelihood**
- Gene universe: all 12,034 QC genes, no HVG selection; the same for every model and training size. The first plan,
  5,000 `seurat_v3` HVGs, missed IFN-induced genes (PSMB9, B2M, ranked > 11,000 of 12,034) and cut the Hallmark
  IFN-α set to 65 of 95 genes; the NB likelihood handles low-count genes, so selection is not needed.
  `notebooks/hvg_inspection.py` reproduces this and doubles as an HVG robustness check.
- Negative binomial likelihood on raw counts for all models (as in Makrodimitris et al. 2024, Brief Bioinform, and
  LDVAE); no zero inflation (Svensson 2020). Not Gaussian/MSE on log data: low counts are not Gaussian. No normalised
  layer is stored; it is computed where needed (figures, gene-set scores).
- Sequencing depth: NB mean = observed total counts × predicted share; removed from the encoder input, restored in
  the decoder. IFN-β raises the interferon genes' share of all counts to up to 25% (CD14+ monocytes), so after
  scaling by totals other genes look up to ~20% lower; the softmax decoder models this as one effect.

**Gene sets and mask**
- Hallmark for the main analysis: all 50 sets keep ≥ 12 genes (median 117), no near-duplicate sets, one
  unambiguous target latent (the Hallmark interferon-α response, type I IFN, 95 of 97 genes). Reactome (used by VEGA,
  expiMap and OntoVAE on this dataset; 1,010 usable sets, 1,416 near-duplicate pairs) is built but not used for now.
- Symbols matched through Ensembl IDs to current HGNC symbols (dated HGNC table): recovers 110 Hallmark genes renamed
  since 2017, including WARS1 and TENT5A in the IFN-α set.
- Sets with < 12 genes after filtering are dropped (as in expiMap).
- Linear decoder in every model (the plain VAE is LDVAE), so the mask is the only architectural difference.
  Encoder: shifted-log counts + donor → one hidden layer of 128 → 58 latents; the loss is NB on raw counts.
  Hard mask: latent k may only use the genes of set k.
- 8 free latents in every masked model, reaching only the 8,847 genes in no set (74%), so set genes can only be
  explained by named latents. First tried with free latents reaching all genes: they absorbed the stimulation
  signal and the named latents switched off (4 of 50 active with the real mask, 0 with the shuffled one). The plain
  VAE keeps 58 latents reaching all genes. Genes in no set stay in, so the knowledge does not choose the
  gene universe.

**Nulls and models**
- Shuffled masks by permuting gene labels across the gene universe: set sizes and overlaps stay identical, only the
  biology changes. 100% of labels permuted (25 / 50 / 75% also built, not used for now), 10 seeds, within 25 equal-size
  expression bins (the finest that the precision of the gene means supports; scanpy's default for control genes).
- Co-expression modules size-matched to the Hallmark sets (random seed gene + its most correlated genes, seeds
  detected in ≥ 1% of cells), built on training cells only, 10 seeds.
- Variants: plain VAE; masked VAE with the real mask, a 100% shuffled mask, or co-expression modules. Donor is a covariate in
  encoder and decoder (it appears in both conditions, so it cannot absorb the IFN-β effect); condition is the
  signal and stays out of the covariates. No batch integration (it could remove the stimulation effect).
- Loss: NB likelihood + KL (KL weight 1, warm-up over ~38 epochs, as in Eltager, ..., Makrodimitris 2023).
- Adam, batch 128, early stopping on validation loss after warm-up (patience 3 epochs, as in Eltager et al.);
  encoder width fixed at 128 (chosen on the plain VAE before the seed runs: 1,920.8 vs 1,923.9 for 256); learning
  rate {1e-3, 1e-4} tuned separately per variant by validation loss; 10 seeds per variant.
- Dropped for simplicity: β-TCVAE, the 25/50/75% dose-response, Reactome, sensitivity checks, soft mask.
- Same tuning budget and several seeds for every model; results are reported as the spread over seeds.

## Current state
Data preprocessed (23,919 cells × 12,034 genes); masks and co-expression null built; `06_fit.py` trains all model
variants; encoder width fixed at 128. Next: tune the learning rate and run 10 seeds per variant, then the evaluation
script (scores A/B/C), then the learning curves (Q1b).
