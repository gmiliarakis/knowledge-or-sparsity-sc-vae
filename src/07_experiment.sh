#!/bin/bash
# The experiment: every variant at learning rate 1e-3, width 128. Run from the project root.
#   bash src/07_experiment.sh                       seeds 0-9 on the full data (the main experiment)
#   SC_TAG=drop250 bash src/07_experiment.sh 0 1    seeds 0 and 1 on the data without the 250 strongest stimulation genes
#   MASKS=responsematch bash src/07_experiment.sh       only the listed variants (default: none real shuffled setshuffle coexpression)
set -e
seeds=${@:-0 1 2 3 4 5 6 7 8 9}
masks=${MASKS:-none real shuffled setshuffle coexpression}
for mask in $masks; do
  for seed in $seeds; do
    .venv/bin/python src/06_fit.py --mask $mask --seed $seed
  done
done
