#!/bin/bash
# The experiment: seeds 1-9 for every variant with its chosen learning rate (seed 0 exists from tuning).
# Run from the project root after 07_tune.sh and 08_choose_settings.py.
set -eo pipefail  # pipefail: stop if results/settings.json cannot be read, instead of running nothing
.venv/bin/python -c "import json; [print(v['mask'], v['lr']) for v in json.load(open('results/settings.json')).values()]" |
while read mask lr; do
  for seed in 1 2 3 4 5 6 7 8 9; do
    .venv/bin/python src/06_fit.py --mask $mask --seed $seed --lr $lr < /dev/null
  done
done
