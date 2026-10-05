#!/bin/bash
# The experiment: seeds 0-9 for every variant at learning rate 1e-3, width 128. Run from the project root.
set -e
for mask in none real shuffled coexpression; do
  for seed in 0 1 2 3 4 5 6 7 8 9; do
    .venv/bin/python src/06_fit.py --mask $mask --seed $seed
  done
done
