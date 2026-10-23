# Reproducing in detail

Commands beyond the main experiment, and the layout of the results and tables. The versions and the main commands are in the [README](../README.md#reproducing); numbers in square brackets refer to its reference list.

`src/06_fit.py` trains one model per call and is run by `07_experiment.sh`. Its options and defaults:

| Option | Default | Meaning |
|---|---|---|
| `--mask` | `none` | `none`, `real`, `shuffled` (fully random mask), `setshuffle` (within-set random mask) or `coexpression` |
| `--seed` | `0` | initial weights, batches and (for the random and co-expression masks) the mask |
| `--train-size` | all training cells | number of training cells, for the learning curves |
| `--max-epochs` | `2000` | safety limit; a run that reaches it fails |
| `--warmup-epochs` | ≈ 38 | epochs for β to rise from 0 to 1; the default is 10,000 steps at batch 64 [14], scaled to the training set |
| `--min-delta` | `0` | relative drop in validation loss that counts as an improvement; `0` is the setting of all reported results, and a value such as `1e-3` stops earlier |
| `--results-dir` | `results` | output folder |
| `--verbose` | off | progress bar over epochs; otherwise one line per run |

The defaults reproduce the reported results. Four runs (vanilla and Hallmark at full size, shuffled and co-expression at
500 cells) matched the stored runs in every epoch and in the saved latent means.

The learning curves, the checks, the tests and the figures run last, in this order. The tests read the learning-curve scores and the checks, and the figures read the tests.

```
python src/extras/learning_curves/01_subsample.py
bash src/extras/learning_curves/02_run.sh
python src/08_evaluate.py --learning-curves
python src/10_checks.py
python src/11_tests.py
python src/09_figures.py
```

`10_checks.py` reads the scored runs and adds the robustness checks, the reference scores, the principal-component check after removing genes and the classifier. `11_tests.py` runs every hypothesis test with its multiplicity correction. The
tables label the vanilla VAE `vanilla`, and its results folder is `plain`.

`10b_response_measures.py` (run after `10_checks.py`, trains nothing) adds the alternative decoded-response measures and the noise bound (`diagnostics_response_alt.csv`, `diagnostics_response_ceiling.csv`). For the main experiment, `python src/10b_response_measures.py --learning-curves` also scores the decoded response of the learning-curve runs. Like `10_checks.py`, it takes `SC_DATA`, `SC_SETS` and `SC_TAG` in the command itself.

The Reactome analysis reuses the same scripts through `SC_SETS`; its results go to `results/kang/reactome/` and `reports/tables/kang/reactome/`, and only the latent-free tests are run on it. Each command gets its variables from its own line: set `SC_DATA`, `SC_SETS` and `SC_TAG` in the command itself and do not export them.

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

The analysis that drops the 250 genes with the largest stimulation effect
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

The Oesinghaus data are read by `src/03c_preprocess_oesinghaus.py`, which reads the published archive in place, and `SC_DATA` selects them in every later script (`oesinghaus_ifn_beta`, `oesinghaus_il_2`). Results and tables go to a folder per data set and collection (`results/oesinghaus_il_2/reactome/`, `reports/tables/oesinghaus_il_2/reactome/`). For IL-2, `SC_TARGET=il2` sets the expected-response set of the evaluation and the checks to the IL-2 sets; it is not used for the masks, the co-expression groups or the fits. The learning curves and the checks that need them are not run on these data. `python src/03c_preprocess_oesinghaus.py --choose DEGs.csv` applies the rule for the second cytokine to the atlas's table of differentially expressed genes.

```
python src/03c_preprocess_oesinghaus.py IL-2     # or IFN-beta
for c in reactome hallmark; do
  SC_DATA=oesinghaus_il_2 SC_SETS=$c python src/04_build_masks.py
  SC_DATA=oesinghaus_il_2 SC_SETS=$c python src/05_build_coexpression_null.py
  SC_DATA=oesinghaus_il_2 SC_SETS=$c bash src/07_experiment.sh
  SC_DATA=oesinghaus_il_2 SC_SETS=$c SC_TARGET=il2 python src/08_evaluate.py
  SC_DATA=oesinghaus_il_2 SC_SETS=$c SC_TARGET=il2 python src/10_checks.py
  SC_DATA=oesinghaus_il_2 SC_SETS=$c python src/11_tests.py
  SC_DATA=oesinghaus_il_2 SC_SETS=$c python src/09_figures.py
done
```

## Expression-matched random mask (Q1c control)

Run from the project root, one stage at a time. The runs cover the six combinations with 10 seeds each.

```
bash src/extras/response_matched_control.sh build   # the masks, minutes
bash src/extras/response_matched_control.sh fit     # 60 runs, about 4 to 5 hours
bash src/extras/response_matched_control.sh check   # rewrites the diagnostics tables with the expression-matched runs
for c in hallmark reactome; do
  SC_DATA=kang SC_SETS=$c python src/11_tests.py
  SC_DATA=oesinghaus_ifn_beta SC_SETS=$c python src/11_tests.py
  SC_DATA=oesinghaus_il_2 SC_SETS=$c SC_TARGET=il2 python src/11_tests.py
done
```

## Layout of the repository, results and tables

Where the files are and which table answers which question.

| Path | Contents |
|---|---|
| `src/` | numbered scripts in run order, with the shared set-up in `scoring.py`; `06_fit.py` trains one model per call and is run by `07_experiment.sh` |
| `src/gene_sets.py` | the data set (`SC_DATA`), gene set collection (`SC_SETS`), benchmark (`SC_TAG`) and expected response (`SC_TARGET`) a script works on |
| `src/extras/learning_curves/` | training subsamples and the runs on them (supplementary S3) |
| `src/extras/hvg_inspection.py` | highly-variable-gene analysis behind the choice to model all genes |
| `data/`, `results/` | created by the scripts and left out of the repository; `results/<data>/<collection>/<variant>/seed<k>/`, for example `results/kang/hallmark/real/seed0/` |
| `reports/tables/<data>/<collection>/` | the main experiment is `kang/hallmark`; scores per run, AUROCs per cell type, paired differences, every hypothesis test (`tests.csv`), robustness checks (`robustness_*.csv`) and diagnostics (`diagnostics_*.csv`) |
| `reports/tables/*/*/tests.csv`, column `family` | `S1` (interferon-α latent, supplementary S1), `S3` (learning curves, Q1b) and `latent-free` (Q1, Q1a and Q1c to Q1f: one table per benchmark) |
| Questions and tables | Q1, Q1a, Q1b: `kang/hallmark`; Q1c: `kang_drop250/hallmark`; Q1d to Q1f: the six benchmarks `kang/{hallmark,reactome}`, `oesinghaus_ifn_beta/{hallmark,reactome}`, `oesinghaus_il_2/{hallmark,reactome}` |
| `reports/figures/` | figures |
