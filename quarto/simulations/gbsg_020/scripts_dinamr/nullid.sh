#!/bin/zsh
# nullid driver (TASK_null_gbsg_identification_2026-09-21): the structural-null
# counterpart of the GBSG survival design, identification and classification only.
# usage: nullid.sh <cellspec-file> [--dry]
#   cellspec lines: "<hr> <hr-tag> <n> <cell>", in RUN ORDER.  --dry prints every
#   run's environment and renders nothing.
#
# Per cell, all THREE identifier runs complete before the next cell starts:
#   consistency, dina, grf -- every one at sg_focus=effMaxSG with
#   FS_S7_NBHD=0.20, i.e. the rule max_A N(eps) at eps = 0.20.  They run on
#   IDENTICAL draws: the seed is seed_base + sim_id (template literal 8316951)
#   and the DGM build is a deterministic function of the cell, so the three
#   runs of a cell see the same 2,000 trials.  nullid_gateC.R checks that from
#   or_Hc_est, which is the Cox fit on flag_harm == 0 -- the WHOLE trial under
#   the structural null -- and so is a draw fingerprint independent of engine.
#
# Knobs: idsweep.sh's pinned set, plus FS_S7_DGM=null and FS_S7_CAMPAIGN=nullid.
# The five field / IJ knobs are MR-only and inert with MR off; pinned for
# continuity with idsweep.  FS_S7_Z1Q is UNSET and must stay unset -- the
# template rejects it under FS_S7_DGM=null, because with no planted region the
# quantile would move the z1 main effect alone.  FS_S7_ER_JCUTS unset.
# sim_id 1-2000, ONE batch per run, so there is no combine render.
# 12 workers; render.sh sets the three thread variables to 1.
#
# NO MR ANYWHERE: FS_S7_MR=FALSE and FS_S7_FB=none.  No bootstrap, no CV.
#
# Gates: GATE A after every run and GATE C after every cell, both recorded.
# A FAILING CELL DOES NOT STOP THE CAMPAIGN (task: "Gates are stop-on-failure,
# never stop-to-ask.  A failing cell writes its halt record and the run
# continues to the next cell").  The rest of that cell is abandoned, a
# HALT_nullid_<cell>.txt record is written in the study directory, and the
# driver moves to the next cell.  HARD TIMEOUT 6 h: a watchdog kills the render
# and the driver.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # render.sh reads $SP
QD="${DINAMR_QMD_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
cd $SP
mkdir -p logs
TIMEOUT_S=21600
REPS=2000
DRY=0; [[ "${2:-}" == "--dry" ]] && DRY=1
PIN=(FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE
     FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_WORKERS=12)
BASE=($PIN FS_S7_MR=FALSE FS_S7_DGM=null FS_S7_CAMPAIGN=nullid FS_S7_MODE=batch
      FS_S7_START=1 FS_S7_NSIMS=$REPS FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20)
T0=${NULLID_T0:-$(date +%s)}
FROM=${NULLID_FROM:-1}
REMAIN_S=$(( TIMEOUT_S - ($(date +%s) - T0) ))
WATCHDOG=""
if (( ! DRY )); then
  MAINPID=$$
  ( sleep $REMAIN_S
    echo "=== HARD TIMEOUT ${TIMEOUT_S}s REACHED $(date): killing the render and the driver ==="
    pkill -f "quarto render sim_fs_maxeffCons_fb_mr_field_m1_template.qmd"
    kill $MAINPID ) &
  WATCHDOG=$!
fi
stop_run () { [[ -n "$WATCHDOG" ]] && kill $WATCHDOG 2>/dev/null; exit $1; }
echo "=== NULLID START $(date) ; campaign T0 $T0 ($(date -r $T0)) ; from cell $FROM ; hard timeout ${TIMEOUT_S}s (6 h) from T0, ${REMAIN_S}s remaining ; watchdog pid ${WATCHDOG:-none} ; dry=$DRY ==="
echo "  host $(hostname) ; $(sysctl -n hw.physicalcpu) physical cores ; $(( $(sysctl -n hw.memsize) / 1073741824 )) GB ; R $(Rscript -e 'cat(paste0(R.version$major,".",R.version$minor))') ; git $(git rev-parse --short HEAD)"

halt_record () {  # $1=cell $2=stage $3=detail
  local F="$QD/HALT_nullid_$1_$(date +%Y-%m-%d).txt"
  {
    echo "HALT RECORD -- nullid, cell $1"
    echo "  task     : TASK_null_gbsg_identification_2026-09-21"
    echo "  when     : $(date)"
    echo "  host     : $(hostname) ; git $(git rev-parse --short HEAD)"
    echo "  stage    : $2"
    echo "  detail   : $3"
    echo "  action   : the rest of this cell is abandoned; the campaign continues at the next cell."
  } > "$F"
  echo "=== HALT RECORD WRITTEN: $F ==="
}

run () {  # $1=engine $2=hr $3=n $4=cell
  local E=$1 H=$2 N=$3 T=$4
  local KN=($BASE FS_S7_HR=$H FS_S7_N=$N)
  local UN=(FS_S7_ER_JCUTS FS_S7_Z1Q)
  if [[ $E == consistency ]]; then UN+=(FS_S7_METHOD); else KN+=(FS_S7_METHOD=$E); fi
  local OUT="nullid_${T}_${E}_effMaxSG"
  local UARGS=(); local u; for u in $UN; do UARGS+=(-u $u); done
  echo "=== RUN START: $T $E $(date) -> $OUT ==="
  echo "  env set:   $KN"
  echo "  env unset: $UN"
  (( DRY )) && return 0
  env $UARGS $KN ./render.sh $OUT < /dev/null || { echo "RUN FAILED (render): $T $E"; return 1; }
  local GL="logs/nullid_gateA_${T}_${E}.log"
  Rscript nullid_gateA.R $OUT $E $H $N $T "$KN" "$UN" > $GL 2>&1 < /dev/null
  local GA=$?
  grep -E '^GATE A ' $GL
  if (( GA != 0 )); then grep -E '^\[FAIL\]' $GL; echo "=== GATE A FAILED: $OUT ==="; return 2; fi
  return 0
}

CI=0; DONE=(); FAILED=()
while read -r H HT N T; do
  [[ -z "$T" ]] && continue
  CI=$(( CI + 1 ))
  if (( CI < FROM )); then DONE+=($T); echo "=== CELL $T: completed in an earlier part; skipped ==="; continue; fi
  echo "=== CELL START: $T (hr=$H n=$N) $(date) ; elapsed $(( $(date +%s) - T0 ))s ==="
  (( DRY )) || { sysctl vm.swapusage; memory_pressure -Q 2>/dev/null | tail -1; }
  CT0=$(date +%s); CELL_OK=1
  for E in consistency dina grf; do
    run $E $H $N $T
    RRC=$?
    if (( RRC != 0 )); then
      CELL_OK=0
      halt_record $T "run $E (rc $RRC)" "the $E run of cell $T failed; see logs/nullid_${T}_${E}_effMaxSG.log and logs/nullid_gateA_${T}_${E}.log"
      echo "=== CELL $T ABANDONED at $E (rc $RRC) $(date); continuing at the next cell ==="
      break
    fi
  done
  if (( ! CELL_OK )); then FAILED+=($T); continue; fi
  echo "CELL DONE: $T  wall=$(( $(date +%s) - CT0 ))s"
  if (( ! DRY )); then
    Rscript nullid_gateC.R $T $H $N > logs/nullid_gateC_${T}.txt 2>&1 < /dev/null
    GC=$?
    grep -E '^GATE C ' logs/nullid_gateC_${T}.txt
    if (( GC != 0 )); then
      grep -E '^\[FAIL\]' logs/nullid_gateC_${T}.txt
      halt_record $T "GATE C" "the per-cell gate failed; see logs/nullid_gateC_${T}.txt"
      FAILED+=($T); echo "=== GATE C FAILED at $T; continuing at the next cell ==="; continue
    fi
  fi
  DONE+=($T)
done < "$1"
echo "NULLID COMPLETE $(date) total_wall=$(( $(date +%s) - T0 ))s ; completed ${#DONE} cells: ${DONE:-none} ; failed ${#FAILED}: ${FAILED:-none}"
stop_run 0
