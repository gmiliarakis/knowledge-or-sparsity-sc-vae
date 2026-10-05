#!/bin/bash
# Learning-rate tuning: every variant at both learning rates, seed 0, width 128. Run from the project root.
set -e
for mask in none real shuffled coexpression; do
  for lr in 1e-3 1e-4; do
    .venv/bin/python src/06_fit.py --mask $mask --seed 0 --lr $lr
  done
done
