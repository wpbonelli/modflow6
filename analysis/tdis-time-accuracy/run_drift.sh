#!/bin/bash
# usage: run_drift.sh <mf6> <dir>...   runs each sim dir under int/pow/hybrid first-step modes
MF6=$1; shift
for ws in "$@"; do
  for mode in int pow hybrid; do
    ( cd "$ws" && MF6_TDIS_FIRSTSTEP=$mode MF6_TDIS_DRIFT=drift_$mode.csv "$MF6" > run_$mode.log 2>&1; echo "$ws $mode exit=$?" )
  done
done
