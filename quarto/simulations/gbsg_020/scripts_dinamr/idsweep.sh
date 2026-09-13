#!/bin/zsh
# idsweep driver (TASK_idsweep_2026-09-12): the Part B identification sweep.
# usage: idsweep.sh <cellspec-file> [--dry]
#   cellspec lines: "<z1q-or-'-'> <n> <hr> <tag>", in RUN ORDER.  --dry prints
#   every run's environment and renders nothing.
#
# Per cell, all SIXTEEN runs complete before the next cell starts:
#   consistency : effMaxSG effMinSG maxeffCons maxeff maxSG minSG
#   dina, grf   : effMaxSG effMinSG maxSG minSG maxeffCons
# On DINA and GRF `maxeffCons` stands for eff = maxeff = maxeffCons (all
# resolve to the plain effect argmax and stem-tag `eff`); it is passed as
# maxeffCons because the template's focus guard does not admit `eff`.
#
# Knobs: campaign.sh's pinned set (its METHOD / FOCUS / NBHD / CAMPAIGN are set
# per run here), plus FS_S7_MR=FALSE and FS_S7_CAMPAIGN=idsweep.  The five
# field / IJ knobs are MR-only and inert with MR off; they are pinned because
# the kickoff pins campaign.sh's set (the pBoc smoke did not set them; Gate T3's
# t3mroff rendered with them).  FS_S7_NBHD=0.20 on effMaxSG / effMinSG only,
# unset otherwise.  FS_S7_Z1Q unset (12.4%) or 0.60 (31%).  FS_S7_ER_JCUTS
# unset.  Seeds 8316951 + sim_id (template literal), sim_id 1-500, ONE batch per
# run -- 500 replicates does not need splitting, so there is no combine render.
# 12 workers; render.sh sets the three thread variables to 1.
#
# Gates: GATE A after every run (idsweep_gateA.R) and GATE I after every cell
# (idsweep_gateI.R), both stop-on-failure.  GATE 1 before every cell
# (idsweep_project.R next): the cell starts only if elapsed + its projected
# wall <= 13 h; otherwise it and every later cell are deferred.  HARD TIMEOUT
# 16 h: a watchdog kills the render and the driver.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # render.sh reads $SP
cd $SP
mkdir -p logs
TIMEOUT_S=57600
REPS=500
DRY=0; [[ "${2:-}" == "--dry" ]] && DRY=1
PIN=(FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE
     FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_WORKERS=12)
BASE=($PIN FS_S7_MR=FALSE FS_S7_CAMPAIGN=idsweep FS_S7_MODE=batch FS_S7_START=1 FS_S7_NSIMS=$REPS)
T0=$(date +%s)
WATCHDOG=""
if (( ! DRY )); then
  MAINPID=$$
  ( sleep $TIMEOUT_S
    echo "=== HARD TIMEOUT ${TIMEOUT_S}s REACHED $(date): killing the render and the driver ==="
    pkill -f "quarto render sim_fs_maxeffCons_fb_mr_field_m1_template.qmd"
    kill $MAINPID ) &
  WATCHDOG=$!
fi
stop_run () { [[ -n "$WATCHDOG" ]] && kill $WATCHDOG 2>/dev/null; exit $1; }
echo "=== IDSWEEP START $(date) ; ceiling 46800s (13 h) ; hard timeout ${TIMEOUT_S}s (16 h) ; watchdog pid ${WATCHDOG:-none} ; dry=$DRY ==="
echo "  host $(hostname) ; $(sysctl -n hw.physicalcpu) physical cores ; $(( $(sysctl -n hw.memsize) / 1073741824 )) GB ; git $(git rev-parse --short HEAD)"

run () {  # $1=engine $2=sg_focus $3=z1q ("" = unset) $4=n $5=hr $6=cell tag
  local E=$1 F=$2 Z=$3 N=$4 H=$5 T=$6
  local KN=($BASE FS_S7_HR=$H FS_S7_N=$N FS_S7_FOCUS=$F)
  local UN=(FS_S7_ER_JCUTS)
  if [[ $F == effMaxSG || $F == effMinSG ]]; then KN+=(FS_S7_NBHD=0.20); else UN+=(FS_S7_NBHD); fi
  if [[ $E == consistency ]]; then UN+=(FS_S7_METHOD); else KN+=(FS_S7_METHOD=$E); fi
  if [[ -z "$Z" ]]; then UN+=(FS_S7_Z1Q); else KN+=(FS_S7_Z1Q=$Z); fi
  local FO=$F; [[ $E != consistency && $F == maxeffCons ]] && FO=eff
  local OUT="idsweep_${T}_${E}_${FO}"
  local UARGS=(); local u; for u in $UN; do UARGS+=(-u $u); done
  echo "=== RUN START: $T $E $F -> $OUT $(date) ==="
  echo "  env set:   $KN"
  echo "  env unset: $UN"
  (( DRY )) && return 0
  env $UARGS $KN ./render.sh $OUT < /dev/null || { echo "RUN FAILED: $T $E $F"; return 1; }
  local GL="logs/idsweep_gateA_${T}_${E}_${FO}.log"
  Rscript idsweep_gateA.R $OUT $E $F "${Z:--}" $N $H "$KN" "$UN" > $GL 2>&1 < /dev/null
  local GA=$?
  grep -E '^GATE A ' $GL
  if (( GA != 0 )); then grep -E '^\[FAIL\]' $GL; echo "=== GATE A FAILED: $OUT ==="; return 2; fi
  return 0
}

CI=0; DEFER=0; DEFERRED=(); DONE=()
while read -r Z N H T; do
  [[ -z "$T" ]] && continue
  CI=$(( CI + 1 ))
  [[ "$Z" == "-" ]] && Z=""
  if (( DEFER )); then DEFERRED+=($T); echo "=== CELL DEFERRED: $T (deferred from the tail) ==="; continue; fi
  EL=$(( $(date +%s) - T0 ))
  if (( ! DRY )); then
    Rscript idsweep_project.R next $CI $EL < /dev/null | tee -a logs/idsweep_reproject.txt
    PRC=${pipestatus[1]}
    if (( PRC == 10 )); then
      DEFER=1; DEFERRED+=($T); echo "=== CELL DEFERRED: $T (Gate 1 re-projection: would exceed the 13 h ceiling) ==="; continue
    elif (( PRC != 0 )); then
      echo "=== RE-PROJECTION ERROR (rc $PRC) before $T: STOPPING $(date) ==="; stop_run 4
    fi
  fi
  echo "=== CELL START: $T (z1q='$Z' n=$N hr=$H) $(date) ; elapsed ${EL}s ==="
  (( DRY )) || { sysctl vm.swapusage; memory_pressure -Q 2>/dev/null | tail -1; }
  CT0=$(date +%s)
  for E in consistency dina grf; do
    if [[ $E == consistency ]]; then FOCI=(effMaxSG effMinSG maxeffCons maxeff maxSG minSG)
    else FOCI=(effMaxSG effMinSG maxSG minSG maxeffCons); fi
    for F in $FOCI; do
      run $E $F "$Z" $N $H $T
      RRC=$?
      if (( RRC != 0 )); then echo "=== STOPPING THE SWEEP at $T $E $F (rc $RRC) $(date) ==="; stop_run 3; fi
    done
  done
  echo "CELL DONE: $T  wall=$(( $(date +%s) - CT0 ))s"
  if (( ! DRY )); then
    Rscript idsweep_gateI.R $T "${Z:--}" $N $H > logs/idsweep_gateI_${T}.txt 2>&1 < /dev/null
    GI=$?
    grep -E '^GATE I ' logs/idsweep_gateI_${T}.txt
    if (( GI != 0 )); then grep -E '^\[FAIL\]' logs/idsweep_gateI_${T}.txt; echo "=== GATE I FAILED at $T: STOPPING THE SWEEP $(date) ==="; stop_run 5; fi
  fi
  DONE+=($T)
done < "$1"
echo "IDSWEEP COMPLETE $(date) total_wall=$(( $(date +%s) - T0 ))s ; completed ${#DONE} cells: ${DONE:-none} ; deferred ${#DEFERRED}: ${DEFERRED:-none}"
stop_run 0
