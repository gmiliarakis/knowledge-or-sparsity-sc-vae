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
- `src/`: package code (data download, preprocessing, gene-set masks, models, evaluation)
- `notebooks/`: exploration and figures
- `data/raw`, `data/processed`: not in git; recreated by the download scripts
- `reports/figures/`: output figures

## Setup
```
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Current state
Repository skeleton only. No data downloaded, no models yet.
