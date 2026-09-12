#!/bin/zsh
# grfmr completion driver (TASK_grfmr_completion_2026-09-12).
# usage: grfmrC.sh <cellspec-file>      cellspec lines: "<z1q-or-'-'> <n> <hr> <tag>"
#
# grfmr.sh's knob set and cell loop, UNCHANGED, with two additions the kickoff
# requires and grfmr.sh did not enforce inline:
#   1. GATE 3 PER BATCH, STOP-ON-FAILURE.  After each batch render, gate3.R is
#      run on that batch's bundle and rendered document; a failure stops the
#      whole run, not just the cell.
#   2. HARD TIMEOUT 10 h.  A watchdog kills the render and the driver at
#      36,000 s from launch.  (Ceiling 9 h is the Gate 1 planning bound; the
#      timeout is the kill bound.  Both are recorded in the report.)
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # render.sh reads $SP
cd $SP
mkdir -p logs
TIMEOUT_S=36000
KN=(FS_S7_METHOD=grf FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE
    FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE
    FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=grfmr FS_S7_WORKERS=12)
QD="${DINAMR_QMD_DIR:-$SP/..}"
MAINPID=$$
( sleep $TIMEOUT_S
  echo "=== HARD TIMEOUT ${TIMEOUT_S}s REACHED $(date): killing the render and the driver ==="
  pkill -f "quarto render sim_fs_maxeffCons_fb_mr_field_m1_template.qmd"
  kill $MAINPID ) &
WATCHDOG=$!
echo "=== GRFMR C RUN START $(date) ; hard timeout ${TIMEOUT_S}s ; watchdog pid $WATCHDOG ==="
cell () {  # $1=z1q ("" = unset) $2=n $3=hr $4=tag
  local Z=$1 N=$2 H=$3 T=$4
  local HH=$(printf "%03.0f" $(( H * 100 )))
  local ZS=""; [[ -n "$Z" ]] && ZS="_z1q60"
  local T0=$(date +%s)
  for PART in "batch 1 1000" "batch 1001 1000" "combine 1 1000"; do
    set -- ${=PART}; local MODE=$1 START=$2 NS=$3
    local OUT="grfmr_${T}_${MODE}_${START}"
    if [[ -z "$Z" ]]; then
      env -u FS_S7_Z1Q -u FS_S7_ER_JCUTS $KN \
        FS_S7_HR=$H FS_S7_N=$N FS_S7_MODE=$MODE FS_S7_START=$START FS_S7_NSIMS=$NS \
        ./render.sh $OUT || { echo "CELL FAILED: $T at $MODE/$START"; return 1; }
    else
      env -u FS_S7_ER_JCUTS $KN FS_S7_Z1Q=$Z \
        FS_S7_HR=$H FS_S7_N=$N FS_S7_MODE=$MODE FS_S7_START=$START FS_S7_NSIMS=$NS \
        ./render.sh $OUT || { echo "CELL FAILED: $T at $MODE/$START"; return 1; }
    fi
    if [[ "$MODE" == "batch" ]]; then
      local END=$(( START + NS - 1 ))
      local BUN="$QD/results/grf_effMaxSG_fb_mr_field_m1_h${HH}_knoise0_n${N}${ZS}_nb20_grfmr_res_${START}_${END}.rds"
      echo "--- GATE 3: $T batch $START ($BUN) ---"
      Rscript gate3.R "$BUN" "$QD/${OUT}.html" > "logs/gate3_${T}_batch_${START}.log" 2>&1
      local G3=$?
      grep -E '^\[(PASS|FAIL)\]|GATE 3:' "logs/gate3_${T}_batch_${START}.log"
      if [[ $G3 -ne 0 ]]; then
        echo "=== GATE 3 FAILED at $T batch $START: STOPPING THE RUN $(date) ==="
        kill $WATCHDOG 2>/dev/null; exit 3
      fi
    fi
  done
  local T1=$(date +%s)
  echo "CELL DONE: $T  wall=$((T1-T0))s"
}
while read -r Z N H T; do
  [[ -z "$T" ]] && continue
  [[ "$Z" == "-" ]] && Z=""
  echo "=== CELL START: $T (z1q='$Z' n=$N hr=$H) $(date) ==="
  cell "$Z" "$N" "$H" "$T" || echo "=== CELL ABORTED: $T, continuing ==="
done < "$1"
kill $WATCHDOG 2>/dev/null
echo "GRFMR C RUN COMPLETE $(date)"
