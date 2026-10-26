#!/bin/bash
# Response-matched random mask (the control for the advantage of the gene sets resting on the strongest response genes). Run from the project root, one stage at a time:
#   bash src/extras/response_matched_control.sh build    masks (minutes)
#   bash src/extras/response_matched_control.sh fit      60 runs, about 4 to 5 hours (10 seeds on each of the six benchmarks)
#   bash src/extras/response_matched_control.sh check    rewrites the diagnostics tables of each benchmark with responsematch included (minutes per benchmark)
set -e
PY=.venv/bin/python
stage=${1:?stage: build, fit or check}
# benchmark = "<data> <collection> <target>"
benchmarks=("kang hallmark -" "kang reactome -" "cytokine_ifn_beta hallmark -" "cytokine_ifn_beta reactome -" "cytokine_il_2 hallmark il2" "cytokine_il_2 reactome il2")
for b in "${benchmarks[@]}"; do
  set -- $b
  d=$1; c=$2; t=$3
  export SC_DATA=$d SC_SETS=$c
  if [ "$t" = "il2" ]; then export SC_TARGET=il2; else unset SC_TARGET; fi
  case $stage in
    build) $PY src/04b_build_responsematch.py ;;
    fit)   MASKS=responsematch bash src/07_experiment.sh ;;
    check) $PY src/10_checks.py ;;
  esac
done
