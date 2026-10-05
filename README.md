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
Fair comparison: same tuning budget and several seeds for every model; report the spread, not the best run.

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

## Decisions (Kang 2018)

**Data**
- Batch 2 of GSE96583 only (ctrl vs 6 h IFN-β, same 8 donors pooled in each 10x run); batch 1 is unstimulated lupus
  samples. Starting from the authors' Cell Ranger 1.2 / hg19 matrices, not from reads (SRA has only 2 × ~23 GB BAMs,
  and the donor, cell-type and doublet labels exist only in the GEO metadata).
- Genes indexed by Ensembl ID: 2,697 gene symbols are repeats in this annotation.
- 313 barcodes occur in both runs; the authors' metadata names the stim copies `-11`, and we follow it.

**Cell QC**
- Doublets: demuxlet calls (genotype-based) instead of a computational doublet detector. Same-donor doublets
  (~2%, estimated) are not caught.
- No ambient-RNA correction: GEO has only the cell-called matrices, not the empty droplets it needs.
- No mitochondrial filter: the 13 MT- genes have zero counts in the deposited matrices.
- Outliers at 5 MADs (sc-best-practices) on log total counts, log detected genes and % counts in the top 20 genes.
  Thresholds per condition × cell type: with thresholds per condition only, normal monocytes (much RNA in a few
  genes), dendritic cells (large, many genes) and megakaryocytes were removed as damaged (up to 27% of a type).
  Per donor as well was rejected: groups of 4–20 cells give unstable MADs, and donors shared each 10x run.
- Result: 24,673 labelled singlets → 23,919 cells (3.1% removed; 1–6% per cell type), 12,034 genes in ≥ 20 cells.
- Open: megakaryocytes still lose ~20% with their own thresholds, so they are likely a mixed group (mostly platelets).
  Platelets have no nucleus and should not mount a transcriptional IFN response: keep, drop, or use as a
  non-responding control? Decide at evaluation.

**Split**
- 70 / 15 / 15 train / val / test, made inside every donor × condition × cell type group (seed 0), so every group
  is in all three sets; at least one val and one test cell per group (smallest group: 4 cells).
- Fitted on train, tuned and latent-matched on val, scored once on test. Learning curves subsample train only.

**Genes and likelihood**
- HVGs: `seurat_v3` on raw counts, per donor (`batch_key`), on training cells only; top 5,000 flagged, not removed.
  The number used is decided with the gene sets. All models share one gene universe, also at small training sizes.
- Negative binomial likelihood on raw counts for all models (as in Makrodimitris et al. 2023, Brief Bioinform, and
  LDVAE). No normalised layer is stored; it is computed where needed (figures, gene-set scores).

## Current state
Kang data downloaded and preprocessed (`data/processed/kang.h5ad`). Next: gene sets and masks.
