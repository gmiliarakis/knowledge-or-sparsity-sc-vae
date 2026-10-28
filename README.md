# Knowledge or sparsity? Gene set masks in a single-cell VAE

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23571985.svg)](https://doi.org/10.5281/zenodo.23571985)

![From stimulated cells to a masked VAE: the six decoder masks and the decoded response each one reached](reports/figures/overview.png)

**Contents**

[Summary](#summary)<br>
[Reproducing](#reproducing)<br>
[Analyses not in the preprint](#analyses-not-in-the-preprint)<br>
[Citation and licence](#citation-and-licence)

## Summary

- **Question:** in a single-cell RNA-seq variational autoencoder (VAE) restricted to gene sets, does the knowledge in the sets or the sparsity they impose explain how well it captures an expected response?
- **Data:** blood cells stimulated with interferon-β and unstimulated controls, from two datasets: 8 lupus patients stimulated for 6 h (Kang et al.) and 12 healthy donors stimulated for 24 h (Theis et al.).
- **Response:** the decoder should reproduce the stimulation shift, the difference in gene expression between control and stimulated cells.
- **Model:** a VAE whose linear decoder is restricted by a binary mask, in six variants. These are the vanilla model, gene sets, three random masks of the same sparsity (fully random, within-set random and expression-matched), and groups of co-expressed genes found in the data. The gene sets are the Hallmark collection (50 sets refined for coherent expression) or 50 Reactome sets chosen from size, overlap, pathway hierarchy and name.
- **Benchmarks:** 4 benchmarks (2 datasets × 2 gene sets) per variant, tested over 10 seeds.
- **Evaluation:** four main measures (test loss, decoded response, cell-type classifier and stimulation classifier). Decoded response is scored by the Spearman and the Pearson correlation. A noise bound, the correlation that a decoder reproducing the true shift would reach with the observed shift, puts these scores in context. The score of the latent labelled with the expected response is in the supplementary information of the preprint.
- **Result:** Every mask worsened the fit relative to the vanilla model in all tested seeds. In all four benchmarks, gene set masks reproduced the shift in expression slightly better than the random masks with Spearman correlation and better with Pearson correlation. With the Hallmark sets on the lupus data, they recovered 10% (Spearman, 95% CI [8, 11]) and 56% (Pearson, 95% CI [49, 62]) of the gap between the fully random mask and the vanilla model. They also reached a lower test loss than the random masks in seven of eight comparisons. Co-expressed groups reproduced the shift better than the gene sets in all four benchmarks and fitted better than the Hallmark sets on the lupus data.
- **Conclusion:** In these data, random masks of the same sparsity decode the response almost on par with the gene set mask. Uncurated groups of co-expressed genes found in the data exceed the gain of the gene set mask over the random masks. The curated sets offset part of the fit lost by restricting the decoder. Only the curated sets assign each latent variable a biological name chosen before training.

### Research questions

One main question and three sub-questions. Each follow-up varies one element of Q1 and retains the design. The comparison is the paired difference between the gene set mask and another mask on four main measures (the decoded response with two correlations). The rationale, methods and results are in the preprint.

- **Q1: knowledge or sparsity?** Does the gene set mask improve on a model without a mask because of the knowledge in the sets, or because of the sparsity?
- **Q1a:** Can pre-specified sets add anything over co-expressed genes found in the data?
- **Q1b:** Is the result reproduced with a second collection (Reactome)?
- **Q1c:** Is the result replicated in a second dataset?

The four benchmarks pair each stimulation with both collections. Q1b compares the Hallmark and Reactome rows within a dataset. Q1c compares the Kang et al. data with the Theis et al. data for interferon-β.

## Reproducing

**Data** All data are public. The Kang et al. data (8 donors, interferon-β) are in the Gene Expression Omnibus under GSE96583 (batch 2). The Theis et al. data are the Parse Biosciences dataset "10 Million Human PBMCs in a Single Experiment", released under CC BY-NC 4.0. The scripts read that archive in place without redistributing it. The gene sets are the Hallmark collection and the Reactome collection (C2 CP:REACTOME) of the Molecular Signatures Database v2024.1. The Reactome pathway hierarchy is the release of 20 October 2026; the SHA-256 of `ReactomePathways_2026-10-20.txt` is `f6d7a2bf89b5bcfe0250a0bc7f51bff94641447911712b8ff129f5b55e52df3a` and that of `ReactomePathwaysRelation_2026-10-20.txt` is `fd49a624d80c14eb37ae57a02e141d574d5ede3f60022bb99edbd909448a3f1e`. Gene symbols follow the HGNC.

**Results** The tables behind every figure, test and number in the preprint (scores per run, paired differences and tests for each benchmark) are in the repository. The scripts recreate the trained models.

**Code** The code is at https://github.com/gmiliarakis/knowledge-or-sparsity-sc-vae under the MIT licence. The code is archived on Zenodo under [10.5281/zenodo.23571985](https://doi.org/10.5281/zenodo.23571985). The DOI covers all versions and resolves to the latest.

The scripts are numbered in run order and run from the repository root in a fresh virtual environment. All results were produced on an Apple M1 Pro (16 GB) under macOS 15.7 with PyTorch, scikit-learn, scanpy, AnnData, NumPy, pandas, SciPy and Matplotlib; `requirements.txt` pins these exact versions:

```
python 3.13.15
torch 2.14.1 (Apple GPU through MPS)
scanpy 1.12.4
anndata 0.13.4
numpy 2.5.3
pandas 3.0.6
scipy 1.18.1
scikit-learn 1.9.1
matplotlib 3.11.2
h5py 3.16.0
requests 2.34.2
scikit-misc 0.5.2
tqdm 4.70.1
```

```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Run order

Run the scripts from the repository root in bash or zsh. `SC_DATA` (dataset) and `SC_SETS` (collection) select the benchmark a script works on; set them in the command itself and do not export them. The defaults are the main experiment, `kang` with `hallmark`. Results go to `results/<data>/<collection>/` and tables to `reports/tables/<data>/<collection>/`.

**1. Data and masks**

`01_download.py` fetches the Kang et al. files and the gene sets with their gene tables. It also reads the counts of the cells that the study uses from the Theis et al. archive. The archive (41 GB) is read in place by range requests and never downloaded; the cells go to `data/raw/` (about 1 GB is transferred for interferon-β). `02_preprocess.py` runs the cell QC and the split of each dataset, and `03_build_masks.py` builds all masks of a dataset for both collections: the gene set masks and their shuffled versions, the expression-matched mask and the co-expression groups.

```
python src/01_download.py
python src/02_preprocess.py kang
python src/02_preprocess.py IFN-beta
SC_DATA=kang python src/03_build_masks.py
SC_DATA=cytokine_ifn_beta python src/03_build_masks.py
```

**2. Fits and scoring** (the four benchmarks; 10 seeds for each of the five variants)

```
for b in "kang hallmark" "kang reactome" "cytokine_ifn_beta hallmark" "cytokine_ifn_beta reactome"; do
  d=${b% *}; c=${b#* }
  SC_DATA=$d SC_SETS=$c bash src/05_experiment.sh
  SC_DATA=$d SC_SETS=$c python src/06_evaluate.py
done
```

`05_experiment.sh` calls `04_fit.py` for every variant and seed (options below).

**3. Expression-matched random mask** (Q1; the four benchmarks, 10 seeds each)

```
bash src/extras/response_matched_control.sh fit     # 40 runs, about 3 hours
```

**4. Checks and tests** (train nothing; the checks read every variant, so run them after step 3)

```
for b in "kang hallmark" "kang reactome" "cytokine_ifn_beta hallmark" "cytokine_ifn_beta reactome"; do
  d=${b% *}; c=${b#* }
  SC_DATA=$d SC_SETS=$c python src/07_checks.py
  SC_DATA=$d SC_SETS=$c python src/08_response_measures.py
  SC_DATA=$d SC_SETS=$c python src/09_tests.py
done
```

**5. Figures** (the main benchmark last: its run draws the figures that need all four)

```
for b in "kang reactome" "cytokine_ifn_beta hallmark" "cytokine_ifn_beta reactome" "kang hallmark"; do
  d=${b% *}; c=${b#* }
  SC_DATA=$d SC_SETS=$c python src/10_figures.py
done
python src/model_schema.py
```

`07_checks.py` reads the scored runs and adds the decoded response (Spearman and Pearson correlation), the robustness checks, the reference scores, the principal-component check after removing genes and the classifier. `09_tests.py` runs every hypothesis test with its multiplicity correction and writes them to `reports/tables/<data>/<collection>/tests.csv`. Its column `family` marks the tests of the interferon-α latent (`S1`), the 15 tests of each benchmark (`latent-free`) and the expression-matched control (`responsematch`). It also writes `shares.csv`, the share of the gap to the vanilla VAE that the gene set mask recovers, with percentile intervals over the seeds. The
tables label the vanilla VAE `vanilla`. Its results folder is `plain`.

`08_response_measures.py` (run after `07_checks.py`, trains nothing) adds the alternative decoded-response measures and the noise bound (`diagnostics_response_alt.csv`, `diagnostics_response_ceiling.csv`). Like `07_checks.py`, it reads `SC_DATA`, `SC_SETS` and `SC_TAG` from the command itself.

### Options of `04_fit.py`

`src/04_fit.py` trains one model per call and is run by `05_experiment.sh`. Its options and defaults:

| Option | Default | Meaning |
|---|---|---|
| `--mask` | `none` | `none`, `real`, `shuffled` (fully random mask), `setshuffle` (within-set random mask), `coexpression` or `responsematch` (expression-matched random mask) |
| `--seed` | `0` | initial weights, batches and (for the random and co-expression masks) the mask |
| `--train-size` | all training cells | number of training cells, used by the learning-curve scripts |
| `--max-epochs` | `2000` | safety limit; a run that reaches it fails |
| `--warmup-epochs` | ≈ 38 | epochs for β to rise from 0 to 1; the default is 10,000 steps at batch 64, scaled to the training set |
| `--min-delta` | `0` | relative drop in validation loss that counts as an improvement; `0` is the setting of all results in the preprint, and a value such as `1e-3` stops earlier |
| `--results-dir` | `results` | output folder |
| `--verbose` | off | progress bar over epochs; otherwise one line per run |

The defaults reproduce the results of the preprint: four runs (vanilla and Hallmark at full size, shuffled and co-expression at 500 cells) coincided with the original runs in every epoch and in the latent means.

## Analyses not in the preprint

The repository keeps the code and the tables of three analyses that were run and are not part of the results of the preprint. Each leaves the four research questions as they are. Their tables are in `reports/tables/`; the numbers below are paired differences between the gene set mask and the other mask over 10 seeds (Hodges–Lehmann estimate and CI at the Bonferroni level of the 15 tests of the benchmark, a simultaneous coverage of 95%).

**Learning curves** (`src/extras/learning_curves/`; tables `learning_curves_*.csv` in `reports/tables/kang/hallmark/`). The vanilla, gene set, fully random and co-expression variants are refit on nested subsamples of 500, 1,000, 2,000, 4,000, 8,000 and all 16,749 training cells. The interferon-α score of the labelled latent is close to its maximum at every size (with 500 cells the vanilla VAE scores 0.992 and the Hallmark mask 0.996), so it cannot show an effect of the number of cells. On the decoded response the vanilla VAE exceeds the gene set mask and the gene set mask exceeds the fully random mask at every size (Spearman correlation at 500 cells: 0.392, 0.292 and 0.280; with all cells: 0.555, 0.426 and 0.412).

```
python src/extras/learning_curves/01_subsample.py
bash src/extras/learning_curves/02_run.sh
python src/06_evaluate.py --learning-curves
python src/08_response_measures.py --learning-curves
```

The learning-curve parts of `07_checks.py` and `09_tests.py` are commented out; uncomment them to add the slopes and the tests of this analysis.

**Without the 250 strongest genes** (`src/extras/drop250/drop_genes.py`, benchmark `SC_TAG=drop250`; tables in `reports/tables/kang_drop250/hallmark/`). The 250 genes with the largest stimulation effect (|t| of stimulated against control cells in the training data) are removed, every mask is rebuilt on the remaining 11,784 genes and the main experiment is repeated. The gene set mask is not distinguished from the fully random mask in test loss (＋0.3 units, CI [−5.2, ＋2.2]) or in the Pearson correlation (＋0.006, CI [−0.023, ＋0.054]), and the Spearman difference is ＋0.006 (CI [−0.002, ＋0.011], adjusted p 0.035). It stays ahead of the within-set random mask (test loss −2.1 units, Spearman ＋0.010, Pearson ＋0.022). The analysis does not test whether the improvement depends on these genes.

```
python src/extras/drop250/drop_genes.py 250
SC_TAG=drop250 python src/03_build_masks.py
SC_TAG=drop250 bash src/05_experiment.sh
SC_TAG=drop250 python src/06_evaluate.py
SC_TAG=drop250 python src/07_checks.py
SC_TAG=drop250 python src/08_response_measures.py
SC_TAG=drop250 python src/09_tests.py
```

**Interleukin-2** (`SC_DATA=cytokine_il_2`, with `SC_TARGET=il2` for the expected-response set of the evaluation; tables in `reports/tables/cytokine_il_2/`). The second stimulation of the Theis et al. data, run with both collections. The gene set mask exceeds the random masks in the decoded response with both correlations (over the fully random mask: Spearman ＋0.013 for Hallmark and ＋0.009 for Reactome, Pearson ＋0.033 and ＋0.010). The sets most raised by interleukin-2 are interferon sets, and the latent labelled with the interleukin-2 set does not reflect the stimulation (AUROC 0.44 for Hallmark and 0.50 for Reactome), so the benchmark is not an independent test of a stimulus other than interferon.

```
python src/01_download.py --cytokine IL-2
python src/02_preprocess.py IL-2
SC_DATA=cytokine_il_2 python src/03_build_masks.py
for c in reactome hallmark; do
  SC_DATA=cytokine_il_2 SC_SETS=$c bash src/05_experiment.sh
  SC_DATA=cytokine_il_2 SC_SETS=$c SC_TARGET=il2 python src/06_evaluate.py
done
```

Add the two interleukin-2 benchmarks to the list in `src/extras/response_matched_control.sh` to fit their expression-matched masks, then run the checks and tests as in step 4 with `SC_TARGET=il2`.

## Citation and licence

Please cite the preprint and the code archive through [`CITATION.cff`](CITATION.cff) (Zenodo: [10.5281/zenodo.23571985](https://doi.org/10.5281/zenodo.23571985); the DOI resolves to the latest version). The code is under the MIT licence ([`LICENSE`](LICENSE)). The Parse Biosciences data are under CC BY-NC 4.0.
