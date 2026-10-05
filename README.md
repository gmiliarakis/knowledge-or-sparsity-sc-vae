# Does prior knowledge help interpretable single-cell models?

Small public-data project (Oct 2026). Knowledge-guided single-cell models wire gene sets (pathways) into their
architecture so that each latent factor reads as a biological programme. Caranzano et al. 2026
(Brief Bioinform, doi:10.1093/bib/bbag425) showed that for 29 *supervised* pathway-informed networks, structure-matched
random pathways do as well as real ones. This project asks the same question for *unsupervised* single-cell models,
where it has not been tested.

## Questions
- **Q1 (core):** does the knowledge itself help, or only the sparsity a gene-set mask imposes? Tested as a
  dose-response: replace 0 / 25 / 50 / 75 / 100% of each gene set's genes (100% = fully shuffled).
- **Q1b:** is learning more data-efficient with real knowledge? The dose-response repeated at several training sizes.
- **Q2 / Q3 (evaluation):** do the explanations recover known biology (interferon response; B cell regulators
  PRDM1, IRF4, XBP1, PAX5, BACH2) and agree with independently measured protein (CITE-seq)?

## Models
plain VAE -> beta-TCVAE (interpretability from statistics, no biology) -> masked VAE (shuffled mask)
-> masked VAE (real mask). Hard (VEGA-style) mask first; soft/refining mask as a later step. PCA only as an optional
supplementary line. Donors handled as a covariate.

Nulls: structure-matched shuffled gene sets, and data-derived co-expression modules size-matched to the real sets.
Fair comparison: same tuning budget and several seeds for every model; results are reported as the spread over seeds.

## Data
Start: Kang et al. 2018 IFN-β-stimulated PBMCs (GSE96583), where the right answer (interferon programme up in
stimulated cells) is known. Then: B cell compartments of the tonsil atlas (Massoni-Badosa et al. 2024), Demela et al.
2026 IRF4/PRDM1 CRISPR data, tonsil CITE-seq. Preprocessing follows Heumos et al., single-cell best practices.

## Layout
- `src/`: one script per stage (`download_kang.py`, `preprocess_kang.py`, later gene sets, models, evaluation)
- `notebooks/`: exploration and figures
- `data/raw`, `data/processed`: not in git; recreated by the scripts
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
- Open: megakaryocytes still lose ~20% with their own thresholds, so they are likely a mixed group (mostly platelets).
  Platelets have no nucleus and should not mount a transcriptional IFN response: keep, drop, or use as a
  non-responding control? Decide at evaluation.

**Split**
- 70 / 15 / 15 train / val / test, made inside every donor × condition × cell type group (seed 0), so every group
  is in all three sets; at least one val and one test cell per group (smallest group: 4 cells).
- Fitted on train, tuned and latent-matched on val, scored once on test. Learning curves subsample train only.

**Genes and likelihood**
- Gene universe: all 12,034 QC genes, no HVG selection; the same for every model and training size. The first plan,
  5,000 `seurat_v3` HVGs, missed IFN-induced genes (PSMB9, B2M, ranked > 11,000 of 12,034) and cut the Hallmark
  IFN-α set to 65 of 95 genes; the NB likelihood handles low-count genes, so selection is not needed.
  `notebooks/hvg_inspection.py` reproduces this and doubles as an HVG robustness check.
- Negative binomial likelihood on raw counts for all models (as in Makrodimitris et al. 2023, Brief Bioinform, and
  LDVAE); no zero inflation (Svensson 2020). Not Gaussian/MSE on log data: low counts are not Gaussian. No normalised
  layer is stored; it is computed where needed (figures, gene-set scores).

**Gene sets and mask**
- Hallmark for the main analysis: all 50 sets keep ≥ 12 genes (median 117), no near-duplicate sets, one
  unambiguous target latent (HALLMARK_INTERFERON_ALPHA_RESPONSE, type I IFN, 95 of 97 genes). Reactome as comparison
  with VEGA, expiMap and OntoVAE, which used it on this dataset (1,010 usable sets, 1,416 near-duplicate pairs).
- Symbols matched through Ensembl IDs to current HGNC symbols (dated HGNC table): recovers 110 Hallmark genes renamed
  since 2017, including WARS1 and TENT5A in the IFN-α set.
- Sets with < 12 genes after filtering are dropped (as in expiMap).
- Hard mask on a linear decoder: latent k may only use the genes of set k. Soft mask later.
- ≤ 16 unmasked latents for genes in no set (74% of genes for Hallmark), the same number in every masked model
  (VEGA's recommendation). Genes in no set are not dropped, so the knowledge does not choose the gene universe.

**Nulls and models**
- Shuffled masks by permuting gene labels across the gene universe: set sizes and overlaps stay identical, only the
  biology changes. 0 / 25 / 50 / 75 / 100% of labels permuted (nested, dose-response), 10 seeds, within 25 equal-size
  expression bins (the finest that the precision of the gene means supports; scanpy's default for control genes).
- Co-expression modules size-matched to the real sets, built on training cells only.
- Plain VAE → beta-TCVAE → masked VAE (shuffled) → masked VAE (real); PCA optional. Donor is a covariate;
  condition is the signal and stays out of the covariates. No batch integration (it could remove the stimulation effect).
- Same tuning budget and several seeds for every model; results are reported as the spread over seeds.

## Current state
Kang data preprocessed (`data/processed/kang.h5ad`, 23,919 cells × 12,034 genes). Masks built
(`data/processed/masks_hallmark.npz`: 50 sets; `masks_reactome.npz`: 1,010 sets), with shuffled versions at
0-100% × 10 seeds. Known-answer check passes: the Hallmark IFN-α set rises by 0.262 in stimulated cells with the real
mask and by 0.012 on average when fully shuffled. Next: co-expression null, then models.
