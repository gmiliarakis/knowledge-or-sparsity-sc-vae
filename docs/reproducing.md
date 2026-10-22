# Reproducing in detail

Commands beyond the main experiment. The versions and the main commands are in the [README](../README.md#reproducing); numbers in square brackets refer to its reference list.

`src/06_fit.py` trains one model per call and is run by `07_experiment.sh`. Its options and defaults:

| Option | Default | Meaning |
|---|---|---|
| `--mask` | `none` | `none`, `real`, `shuffled` (fully random mask), `setshuffle` (within-set random mask) or `coexpression` |
| `--seed` | `0` | initial weights, batches and (for the random and co-expression masks) the mask |
| `--train-size` | all training cells | number of training cells, for the learning curves |
| `--max-epochs` | `2000` | safety ceiling; a run that reaches it fails |
| `--warmup-epochs` | ≈ 38 | epochs for β to rise from 0 to 1; the default is 10,000 steps at batch 64 [14], scaled to the training set |
| `--min-delta` | `0` | relative drop in validation loss that counts as an improvement; `0` is the setting of all reported results, and a value such as `1e-3` stops earlier |
| `--results-dir` | `results` | output folder |
| `--verbose` | off | progress bar over epochs; otherwise one line per run |

The defaults reproduce the reported results. Four runs (no mask and Hallmark at full size, shuffled and co-expression at
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

`10_checks.py` reads the scored runs and adds the robustness checks, the reference scores, the principal-component check after removing genes and the probe. `11_tests.py` runs every hypothesis test with its multiplicity correction. The
tables label the no-mask VAE `vanilla`, and its results folder is `plain`.

`10b_response_measures.py` (run after `10_checks.py`, trains nothing) adds the alternative decoded-response measures and the noise ceiling (`diagnostics_response_alt.csv`, `diagnostics_response_ceiling.csv`). For the main benchmark, `python src/10b_response_measures.py --learning-curves` also scores the decoded response of the learning-curve runs. Like `10_checks.py`, it takes `SC_DATA`, `SC_SETS` and `SC_TAG` in the command itself.

The Reactome analysis reuses the same scripts through `SC_SETS`; its results go to `results/kang/reactome/` and `reports/tables/kang/reactome/`, and only the latent-free tests are run on it. A script never uses a variable set in an earlier command: set `SC_DATA`, `SC_SETS` and `SC_TAG` in the command itself and do not export them.

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
reuses every script through the environment variable `SC_TAG`; its files carry the suffix `_drop250` and the main files are not touched.

```
python src/03b_drop_genes.py 250
SC_TAG=drop250 python src/04_build_masks.py
SC_TAG=drop250 python src/05_build_coexpression_null.py
SC_TAG=drop250 bash src/07_experiment.sh
SC_TAG=drop250 python src/08_evaluate.py
SC_TAG=drop250 python src/10_checks.py
SC_TAG=drop250 python src/11_tests.py
```

The Oesinghaus data are read by `src/03c_preprocess_oesinghaus.py`, which reads the published archive in place, and `SC_DATA` selects them in every later script (`oesinghaus_ifn_beta`, `oesinghaus_il_2`). Results and tables go to a folder per data set and collection (`results/oesinghaus_il_2/reactome/`, `reports/tables/oesinghaus_il_2/reactome/`). For IL-2, `SC_TARGET=il2` sets the known-answer set of the evaluation and the checks to the IL-2 sets; it is not used for the masks, the co-expression groups or the fits. The learning curves and the checks that need them are not run on these data. `python src/03c_preprocess_oesinghaus.py --choose DEGs.csv` applies the rule for the second cytokine to the atlas's table of differentially expressed genes.

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

## Response-matched random mask (Q1c control)

Run from the project root, one stage at a time. The runs cover the six benchmarks with 10 seeds each.

```
bash src/extras/response_matched_control.sh build   # the masks, minutes
bash src/extras/response_matched_control.sh fit     # 60 runs, about 4 to 5 hours
bash src/extras/response_matched_control.sh check   # rewrites the diagnostics tables with the response-matched runs
for c in hallmark reactome; do
  SC_DATA=kang SC_SETS=$c python src/11_tests.py
  SC_DATA=oesinghaus_ifn_beta SC_SETS=$c python src/11_tests.py
  SC_DATA=oesinghaus_il_2 SC_SETS=$c SC_TARGET=il2 python src/11_tests.py
done
```
