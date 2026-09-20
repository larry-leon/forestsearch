#!/bin/bash
# TASK_binary_stage2_v3_2026-09-19 Step 1 — sequencing driver, Mac Studio.
# It adds NOTHING to the campaign's logic: it invokes the unchanged
# scripts_or/run_or.sh once per cell, in the task's table order, under the
# 0d environment (shim PATH, W, scaled caps, ORSG_HOST, VECLIB=1).
#   * per-cell invocation is what gives the task's exact order: run_or.sh
#     iterates its own CELLS_ALL order, which puts or150_n2000 before
#     or100_n2000, and ORSG_CELLS selects a subset without reordering it;
#   * "a cell that fails writes its halt record and the run continues to the
#     next cell" (Step 1): run_or.sh's halt() exits the invocation and its
#     preflight refuses to start while HALT_or.md exists, so between cells the
#     driver files the halt record under logs_or/ and clears it exactly as
#     commit 236ef82a did.  The failed cell is NOT re-run: no retry logic.
#   * consequence, for the report: run_or.sh's cumulative CUM resets at each
#     invocation, so the ceiling is per cell.  The driver keeps the true
#     cumulative render wall itself and prints it after every cell.
set -u
SHIM=$HOME/.orsg_shim_2026-09-19
export PATH="$SHIM:$PATH"
export ORSG_WORKERS=13            # physical cores (14) - 1
export ORSG_TIMEOUT=52339         # 10800 s x 63/13, per render
export ORSG_CEILING=440000        # 90000 s x 63/13, rounded up to 10000 s
export ORSG_NSIMS=1000
export ORSG_HOST=Mac-Studio-3.local
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
REPO=/Users/larryleon/Documents/GitHub/forestsearch
QMD=$REPO/quarto/simulations/actg175/binary_020
RQ=quarto/simulations/actg175/binary_020
LOGD=$QMD/logs_or
DRV=$LOGD/driver_stage2_mac.log
CUM_ALL=0
utc () { date -u +%Y-%m-%dT%H:%M:%SZ; }
say () { echo "[$(utc)] $*" >> "$DRV"; }
# The ten cells, in the task's table order: <campaign> <cell-token>
CELLS="orgrf:or075_n500 orgrf:or100_n500 orgrf:or075_n2000 orgrf:or100_n2000 orgrf:or150_n2000 \
ordina:or075_n500 ordina:or100_n500 ordina:or075_n2000 ordina:or100_n2000 ordina:or150_n2000"
say "DRIVER START: workers $ORSG_WORKERS, timeout ${ORSG_TIMEOUT}s, ceiling ${ORSG_CEILING}s, host $ORSG_HOST"
say "DRIVER PATH: $PATH"
i=0
for spec in $CELLS; do
  i=$((i+1))
  CAMP=${spec%%:*}; CELL=${spec##*:}
  if [ -f "$QMD/HALT_or.md" ]; then
    TS=$(date -u +%Y%m%dT%H%M%SZ)
    cp "$QMD/HALT_or.md" "$LOGD/HALT_or_cleared_${TS}.md"
    HCELL=$(grep -m1 '^- Cell:' "$QMD/HALT_or.md" | sed 's/^- Cell: //')
    say "CLEARING HALT from ${HCELL:-unknown} (filed as logs_or/HALT_or_cleared_${TS}.md)"
    ( cd "$REPO" && git rm -q -- "$RQ/HALT_or.md" && \
      git commit -q -m "actg175 or: clear HALT_or.md -- ${HCELL:-a cell} halted; stage 2 continues with the next cell

TASK_binary_stage2_v3_2026-09-19 Step 1: per-cell stop-on-failure only -- the
halted cell is not re-run and no retry logic is added to the runner.  The halt
record is filed unchanged at $RQ/logs_or/HALT_or_cleared_${TS}.md.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>" -- "$RQ/HALT_or.md" ) >> "$DRV" 2>&1 \
      || say "HALT CLEAR COMMIT FAILED"
  fi
  T0=$(date +%s)
  say "CELL $i/10 START: ${CAMP}_${CELL}"
  ORSG_CAMPAIGNS="$CAMP" ORSG_CELLS="$CELL" bash "$QMD/scripts_or/run_or.sh" >> "$LOGD/runner_mac.log" 2>&1
  RC=$?
  W=$(( $(date +%s) - T0 )); CUM_ALL=$(( CUM_ALL + W ))
  SWAP=$(sysctl -n vm.swapusage)
  say "CELL $i/10 END: ${CAMP}_${CELL} rc=$RC invocation_wall_s=$W driver_cumulative_s=$CUM_ALL"
  say "CELL $i/10 SWAP: $SWAP"
done
if [ -f "$QMD/HALT_or.md" ]; then say "PASS END with HALT_or.md present (last cell); left in place for the report"; fi
say "DRIVER COMPLETE: 10 cells attempted, driver cumulative wall ${CUM_ALL} s"
