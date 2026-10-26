# Knowledge or sparsity? Gene set masks in a single-cell VAE

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23571985.svg)](https://doi.org/10.5281/zenodo.23571985)

![From stimulated cells to a masked VAE: the six decoder masks and the decoded response each one reached](reports/figures/overview.png)

**Contents**

[Summary](#summary)<br>
[Repository layout](#repository-layout)<br>
[Reproducing](#reproducing)<br>
[Citation and licence](#citation-and-licence)

## Summary

**Key findings**

- **Sparsity has a cost** The vanilla VAE has the lowest test loss in all six benchmarks. Its encoder width was chosen on it.
- **The gene set mask outperforms random masks in the decoded response** in all six benchmarks, by less on the Theis et al. data.
- **The gene set mask has a lower test loss than the within-set random mask in all six benchmarks and than the fully random mask in four.**
- **Co-expressed groups match or exceed the Hallmark sets** in all six benchmarks.
- **The improvement depends largely on the strongest response genes** Without the 250 strongest genes it vanishes against the fully random mask and declines to a fifth against the within-set random mask.
- **No data-efficiency benefit over the vanilla VAE** Its decoded response exceeds that of the Hallmark mask at every training size tested.

**At a glance**

- **Question:** in a single-cell RNA-seq variational autoencoder (VAE) restricted to gene sets, does the knowledge in the sets or the sparsity they impose explain how well it captures an expected response? Six follow-up questions vary the comparison mask, the number of training cells, the genes, the collection, the dataset and the stimulus ([Research questions](#research-questions)).
- **Data:** stimulated against control blood cells in three stimulations: interferon-β for 6 h in 8 lupus patients (Kang et al.), and interferon-β and interleukin-2 for 24 h in 12 healthy donors (Theis et al.).
- **Response:** the stimulation shift, the difference in gene expression between control and stimulated cells, which the decoder should reproduce.
- **Model:** a VAE whose linear decoder is restricted by a binary mask, in six variants. These are the vanilla model, gene sets, three random masks of the same sparsity (fully random, within-set random and expression-matched, which controls for the strongest response genes), and groups of co-expressed genes found in the data. The gene sets are the Hallmark collection (50 sets refined for coherent expression) or 50 Reactome sets chosen from size, overlap and pathway hierarchy alone, without reference to expression.
- **Benchmarks:** the six variants fit with 10 seeds on six benchmarks, each a pairing of a stimulation with a gene set collection (see the preprint).
- **Evaluation:** four main measures (test loss, decoded response, cell-type classifier and stimulation classifier). The score of the latent labelled with the expected response, learning curves and a repeat without the 250 strongest stimulation genes are in the supplementary information of the preprint.
- **Conclusion:** co-expressed groups found in the data match or exceed the Hallmark sets, so the improvement of the gene set mask over random masks is as well explained by genes that respond together as by curated biology.

### Research questions

One core question and six follow-ups. Each follow-up varies one element of Q1 and retains the design: the paired difference between the gene set mask and another mask on four main measures. The rationale, methods and results are in the preprint.

- **Q1: knowledge or sparsity?** Does the gene set mask improve on a model without a mask because of the knowledge in the sets, or because of the sparsity alone?
- **Q1a:** Can pre-specified sets add anything over co-expressed genes found in the data?
- **Q1b:** Can the gene set mask work with fewer training cells?
- **Q1c:** Does the improvement depend on the strongest response genes?
- **Q1d:** Is the result reproduced with a second collection (Reactome)?
- **Q1e:** Is the result replicated in a second dataset?
- **Q1f:** Is the result consistent for a second stimulus?

The six benchmarks pair each stimulation with both collections. Q1d compares the Hallmark and Reactome rows within a stimulation, Q1e the Kang et al. data with the Theis et al. data for interferon-β, and Q1f interferon-β with IL-2 in the Theis et al. data.

## Repository layout

Where the files are and which table answers which question.

| Path | Contents |
|---|---|
| `src/` | numbered scripts in run order, with the shared set-up in `scoring.py`; `06_fit.py` trains one model per call and is run by `07_experiment.sh` |
| `src/gene_sets.py` | the dataset (`SC_DATA`), gene set collection (`SC_SETS`), benchmark (`SC_TAG`) and expected response (`SC_TARGET`) a script works on |
| `src/extras/learning_curves/` | training subsamples and the runs on them (supplementary S3) |
| `src/extras/hvg_inspection.py` | highly-variable-gene analysis behind the choice to model all genes |
| `data/`, `results/` | created by the scripts and excluded from the repository; `results/<data>/<collection>/<variant>/seed<k>/`, for example `results/kang/hallmark/real/seed0/` |
| `reports/tables/<data>/<collection>/` | the main experiment is `kang/hallmark`; scores per run, AUROCs per cell type, paired differences, every hypothesis test (`tests.csv`), robustness checks (`robustness_*.csv`) and diagnostics (`diagnostics_*.csv`) |
| `reports/tables/*/*/tests.csv`, column `family` | `S1` (interferon-α latent, supplementary S1), `S3` (learning curves, Q1b), `latent-free` (Q1, Q1a and Q1d to Q1f: one table per benchmark, with the repeat without the 250 strongest genes for Q1c) and `responsematch` (the expression-matched control, Q1c) |
| Questions and tables | Q1, Q1a, Q1b: `kang/hallmark`; Q1c: `kang_drop250/hallmark`; Q1d to Q1f: the six benchmarks `kang/{hallmark,reactome}`, `cytokine_ifn_beta/{hallmark,reactome}`, `cytokine_il_2/{hallmark,reactome}` |
| `reports/figures/` | figures |

## Reproducing

**Data** All data are public. The Kang et al. data (8 donors, interferon-β) are in the Gene Expression Omnibus under GSE96583 (batch 2). The Theis et al. data are the Parse Biosciences dataset "10 Million Human PBMCs in a Single Experiment", released under CC BY-NC 4.0. The scripts read that archive in place without redistributing it. The gene sets are the Hallmark collection and the Reactome collection (C2 CP:REACTOME) of the Molecular Signatures Database v2024.1. The Reactome pathway hierarchy is the release of 20 October 2026. Gene symbols follow the HGNC.

**Results** The tables behind every figure, test and number in the preprint (scores per run, paired differences and tests for each benchmark) are in the repository. The scripts recreate the trained models.

**Code** The code is at https://github.com/gmiliarakis/knowledge-or-sparsity-sc-vae under the MIT licence. The code is archived on Zenodo under [10.5281/zenodo.23571985](https://doi.org/10.5281/zenodo.23571985), which stands for all versions and resolves to the latest.

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
python src/01_download_kang.py
python src/02_download_genesets.py
python src/03_preprocess_kang.py
python src/04_build_masks.py
python src/05_build_coexpression_null.py
bash src/07_experiment.sh
python src/08_evaluate.py
```

The commands for the other analyses follow.

### Options of `06_fit.py`

`src/06_fit.py` trains one model per call and is run by `07_experiment.sh`. Its options and defaults:

| Option | Default | Meaning |
|---|---|---|
| `--mask` | `none` | `none`, `real`, `shuffled` (fully random mask), `setshuffle` (within-set random mask), `coexpression` or `responsematch` (expression-matched random mask) |
| `--seed` | `0` | initial weights, batches and (for the random and co-expression masks) the mask |
| `--train-size` | all training cells | number of training cells, for the learning curves |
| `--max-epochs` | `2000` | safety limit; a run that reaches it fails |
| `--warmup-epochs` | ≈ 38 | epochs for β to rise from 0 to 1; the default is 10,000 steps at batch 64, scaled to the training set |
| `--min-delta` | `0` | relative drop in validation loss that counts as an improvement; `0` is the setting of all results in the preprint, and a value such as `1e-3` stops earlier |
| `--results-dir` | `results` | output folder |
| `--verbose` | off | progress bar over epochs; otherwise one line per run |

The defaults reproduce the results of the preprint: four runs (vanilla and Hallmark at full size, shuffled and co-expression at 500 cells) coincided with the original runs in every epoch and in the latent means.

### Learning curves, checks, tests and figures

The learning curves, the checks, the tests and the figures run last, in this order. The tests read the learning-curve scores and the checks, and the figures read the tests.

```
python src/extras/learning_curves/01_subsample.py
bash src/extras/learning_curves/02_run.sh
python src/08_evaluate.py --learning-curves
python src/10_checks.py
python src/10b_response_measures.py
python src/10b_response_measures.py --learning-curves
python src/11_tests.py
python src/09_figures.py
```

`10_checks.py` reads the scored runs and adds the robustness checks, the reference scores, the principal-component check after removing genes and the classifier. `11_tests.py` runs every hypothesis test with its multiplicity correction. The
tables label the vanilla VAE `vanilla`, and its results folder is `plain`.

`10b_response_measures.py` (run after `10_checks.py`, trains nothing) adds the alternative decoded-response measures and the noise bound (`diagnostics_response_alt.csv`, `diagnostics_response_ceiling.csv`). For the main experiment, `python src/10b_response_measures.py --learning-curves` also scores the decoded response of the learning-curve runs. Like `10_checks.py`, it reads `SC_DATA`, `SC_SETS` and `SC_TAG` from the command itself.

### Reactome sets

The Reactome analysis reuses the same scripts through `SC_SETS`; its results are in `results/kang/reactome/` and `reports/tables/kang/reactome/`, and only the latent-free tests are run on it. Each command reads its variables from its own line: set `SC_DATA`, `SC_SETS` and `SC_TAG` in the command itself and do not export them.

```
python src/02_download_genesets.py
SC_SETS=reactome python src/04_build_masks.py
SC_SETS=reactome python src/05_build_coexpression_null.py
SC_SETS=reactome bash src/07_experiment.sh
SC_SETS=reactome python src/08_evaluate.py
SC_SETS=reactome python src/10_checks.py
SC_SETS=reactome python src/11_tests.py
SC_SETS=reactome python src/09_figures.py
```

### Without the 250 strongest genes

The analysis that removes the 250 genes with the largest stimulation effect
reuses every script through the environment variable `SC_TAG`; its files have the suffix `_drop250` and the main files are not touched.

```
python src/03b_drop_genes.py 250
SC_TAG=drop250 python src/04_build_masks.py
SC_TAG=drop250 python src/05_build_coexpression_null.py
SC_TAG=drop250 bash src/07_experiment.sh
SC_TAG=drop250 python src/08_evaluate.py
SC_TAG=drop250 python src/10_checks.py
SC_TAG=drop250 python src/11_tests.py
```

### Theis et al. data

The Theis et al. data are read by `src/03c_preprocess_cytokine.py`, which reads the published archive in place, and `SC_DATA` selects them in every later script (`cytokine_ifn_beta`, `cytokine_il_2`). Results and tables are in a folder per dataset and collection (`results/cytokine_il_2/reactome/`, `reports/tables/cytokine_il_2/reactome/`). For IL-2, `SC_TARGET=il2` sets the expected-response set of the evaluation and the checks to the IL-2 sets; it is not used for the masks, the co-expression groups or the fits. The learning curves and the checks that depend on them are not run on these data. `python src/03c_preprocess_cytokine.py --choose DEGs.csv` applies the rule for the second cytokine to the atlas's table of differentially expressed genes.

```
python src/03c_preprocess_cytokine.py IL-2     # or IFN-beta
for c in reactome hallmark; do
  SC_DATA=cytokine_il_2 SC_SETS=$c python src/04_build_masks.py
  SC_DATA=cytokine_il_2 SC_SETS=$c python src/05_build_coexpression_null.py
  SC_DATA=cytokine_il_2 SC_SETS=$c bash src/07_experiment.sh
  SC_DATA=cytokine_il_2 SC_SETS=$c SC_TARGET=il2 python src/08_evaluate.py
  SC_DATA=cytokine_il_2 SC_SETS=$c SC_TARGET=il2 python src/10_checks.py
  SC_DATA=cytokine_il_2 SC_SETS=$c python src/11_tests.py
  SC_DATA=cytokine_il_2 SC_SETS=$c python src/09_figures.py
done
```

### Expression-matched random mask (Q1c)

Run from the project root, one stage at a time. The runs cover the six benchmarks with 10 seeds each.

```
bash src/extras/response_matched_control.sh build   # the masks, minutes
bash src/extras/response_matched_control.sh fit     # 60 runs, about 4 to 5 hours
bash src/extras/response_matched_control.sh check   # rewrites the diagnostics tables with the expression-matched runs
for c in hallmark reactome; do
  SC_DATA=kang SC_SETS=$c python src/11_tests.py
  SC_DATA=cytokine_ifn_beta SC_SETS=$c python src/11_tests.py
  SC_DATA=cytokine_il_2 SC_SETS=$c SC_TARGET=il2 python src/11_tests.py
done
```

## Citation and licence

Please cite the preprint and the code archive through [`CITATION.cff`](CITATION.cff) (Zenodo: [10.5281/zenodo.23571985](https://doi.org/10.5281/zenodo.23571985), which resolves to the latest version). The code is under the MIT licence ([`LICENSE`](LICENSE)). The Parse Biosciences data are under CC BY-NC 4.0.
