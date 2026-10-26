#!/usr/bin/env bash
# run_stage0b.sh -- driver for the bootstrap-efficiency Stage-0b measurement
# (GLM paths and the consistency stage).  Measurement only: no package file is
# changed.  Run from a plain terminal:
#     cd dev/glm-continuous-sims/stage0b_glm_2026-10-26 && bash run_stage0b.sh
# Each step has a hard time cap; a non-zero exit stops the driver.
set -uo pipefail
cd "$(dirname "$0")"
export STAGE0B_OUT=results
mkdir -p results logs
run() {
  local name=$1 cap=$2; shift 2
  echo "== $name (cap $cap)"
  timeout -k 15s "$cap" "$@" > "logs/$name.log" 2>&1
  local rc=$?
  tail -n 3 "logs/$name.log" | sed 's/^/   /'
  if [ $rc -ne 0 ]; then echo "STOP: $name exited $rc -- see logs/$name.log"; exit $rc; fi
}
run b0_provenance 5m  Rscript b0_provenance.R
run b1_profiles   60m Rscript b1_profiles.R
run b2_lean       60m Rscript b2_lean.R
run b3_summary    5m  Rscript b3_summary.R
echo "== all steps complete; summary at results/stage0b_summary.md"
