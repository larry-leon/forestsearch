#!/usr/bin/env bash
# run_stage0.sh -- driver for the bootstrap-efficiency Stage-0 measurement.
# Measurement only: no package file is read for writing, let alone changed.
# Run from a plain terminal (multicore is disabled under RStudio):
#     cd dev/glm-continuous-sims/stage0_bootstrap_2026-10-26 && bash run_stage0.sh
# Each step has a hard time cap; a non-zero exit stops the driver.
set -uo pipefail
cd "$(dirname "$0")"
export STAGE0_OUT=results
mkdir -p results logs
run() {
  local name=$1 cap=$2; shift 2
  echo "== $name (cap $cap)"
  timeout -k 15s "$cap" "$@" > "logs/$name.log" 2>&1
  local rc=$?
  tail -n 3 "logs/$name.log" | sed 's/^/   /'
  if [ $rc -ne 0 ]; then echo "STOP: $name exited $rc -- see logs/$name.log"; exit $rc; fi
}
run s0_provenance 5m  Rscript s0_provenance.R
run s1_profile    15m Rscript s1_profile.R
run s2_lean       15m Rscript s2_lean_identity.R
run s3_bootstrap  40m Rscript s3_bootstrap_runs.R
run s4_summary    5m  Rscript s4_summary.R
echo "== all steps complete; summary at results/stage0_summary.md"
