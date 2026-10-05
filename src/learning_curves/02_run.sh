#!/bin/bash
# Learning curves: every variant on nested training subsamples, seeds 0-9, otherwise as src/07_experiment.sh; the full
# training set is the main experiment. Run from the project root after src/learning_curves/01_subsample.py.
set -e
for size in 500 1000 2000 4000 8000; do
  for mask in none real shuffled coexpression; do
    for seed in 0 1 2 3 4 5 6 7 8 9; do
      .venv/bin/python src/06_fit.py --mask $mask --seed $seed --train-size $size
    done
  done
done
