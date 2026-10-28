#!/bin/bash
# Response-matched random mask (the control for the advantage of the gene sets resting on the strongest response genes). Run from the project root, one stage at a time:
# The masks are built by src/03_build_masks.py, together with the gene set and co-expression masks.
#   bash src/extras/response_matched_control.sh fit      40 runs, about 3 hours (10 seeds on each of the four benchmarks)
#   bash src/extras/response_matched_control.sh check    rewrites the diagnostics tables of each benchmark with responsematch included (minutes per benchmark)
set -e
PY=.venv/bin/python
stage=${1:?stage: fit or check}
# benchmark = "<data> <collection> <target>"
# the two interleukin-2 benchmarks ("cytokine_il_2 hallmark il2" "cytokine_il_2 reactome il2") belong to an analysis the preprint does not report; add them here to run them
benchmarks=("kang hallmark -" "kang reactome -" "cytokine_ifn_beta hallmark -" "cytokine_ifn_beta reactome -")
for b in "${benchmarks[@]}"; do
  set -- $b
  d=$1; c=$2; t=$3
  export SC_DATA=$d SC_SETS=$c
  if [ "$t" = "il2" ]; then export SC_TARGET=il2; else unset SC_TARGET; fi
  case $stage in
    fit)   MASKS=responsematch bash src/05_experiment.sh ;;
    check) $PY src/07_checks.py ;;
  esac
done
