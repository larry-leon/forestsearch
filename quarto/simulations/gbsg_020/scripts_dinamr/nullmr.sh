#!/bin/bash
# nullmr driver (TASK_null_gbsg_mr_2026-09-21): the structural null at the
# nullid screen (c1 0.90 / c2 0.80 / p* 0.90), WITH MR post-selection
# inference.  Transplant of nullthr.sh (the nullc125 driver).
# usage: nullmr.sh <cells-file> [--dry]
#   cells file: cell lines "<hr> <hr-tag> <n> <cell>" in RUN ORDER; '#' lines
#   and blanks are skipped.  (No 'thresholds' line: the threshold knobs are
#   unset, so the template's 0.90 / 0.80 defaults apply.)
#   --dry prints every run's environment and renders nothing.
# Environment (optional; nullthr.sh's names, kept):
#   NULLTHR_CAMPAIGN  campaign tag (default nullmr)
#   NULLTHR_REPS      replicates per run (default 2000; sim_id 1..REPS, one batch)
#   NULLTHR_QUICKRUN  TRUE tags the stem _quickrun (smoke); default FALSE
#   NULLTHR_COMMIT    1 commits each completed cell (default 1); 0 for a smoke
#
# Differences from nullthr.sh (Step 2.1), and nothing else:
#   1. FS_S7_MR=TRUE, plus the per-engine MR knobs of the Step 1 record
#      (logs/nullmr_step1_record.md): common FIELD_COMPLEMENT=TRUE,
#      FIELD_SCALEC=selected, FIELD_DECOMP=TRUE, IJ_RESIDUAL=two_term, FB=none;
#      FS (campaign_p12x20.sh) adds RETURN_RESEL=TRUE and leaves FIELD_RECOV
#      unset; DINA (campaign.sh) and GRF (grfmr.sh / grfmrC.sh) add
#      FIELD_RECOV=TRUE and leave RETURN_RESEL unset.  UNIFORM / WINNER_ROWS
#      unset everywhere.
#   2. Campaign tag nullmr.
#   3. Threshold knobs unset: FS_S7_C1 / FS_S7_C2 are never passed; the gates
#      are told 0.90 / 0.80, the template defaults.
#   4. Timeouts: 6 h per render (timeout(1)), 24 h for the campaign.
#   Consequences of 1 and 3, mechanical: the bundle stem carries neither
#   _nomr nor a threshold tag (template :481, :505-506); the gates are
#   nullmr_gateA.R / nullmr_gateC.R; the halt record names this task.
# Kept from nullthr.sh: bash, the pinned focus / band / DGM / worker set
# (64 workers), every inherited FS_S7_* unset first, sim_id 1..REPS in one
# batch, halt-and-continue (a failing run abandons the rest of its cell,
# writes HALT_nullmr_<cell>_<date>.txt in the study directory, and the
# campaign continues at the next cell), and commit-per-cell by named paths.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # render.sh reads $SP
QD="$(cd "$SP/.." && pwd)"
RENDER="$SP/render.sh"
REPO=$(git -C "$QD" rev-parse --show-toplevel)
RQ=${QD#"$REPO"/}
export PATH=/usr/lib/rstudio/resources/app/bin/quarto/bin:$PATH
cd "$SP" || exit 1
mkdir -p logs
CAMP=${NULLTHR_CAMPAIGN:-nullmr}
REPS=${NULLTHR_REPS:-2000}
QR=${NULLTHR_QUICKRUN:-FALSE}
DOCOMMIT=${NULLTHR_COMMIT:-1}
WORKERS=64                 # run_p12x20.sh, the p12x20 driver that ran on pop-os
RENDER_TIMEOUT_S=21600     # 6 h per render
TIMEOUT_S=86400            # 24 h for the campaign
DRY=0; [[ "${2:-}" == "--dry" ]] && DRY=1
TRAILER=$'\n\nCo-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>'

C1=0.90; C2=0.80            # template defaults (knobs unset); for the gates only

UNSET=()
for v in $(env | grep -o '^FS_S7_[A-Z0-9_]*'); do UNSET+=(-u "$v"); done
PIN=(FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE
     FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_WORKERS=$WORKERS)
MR_FS=(FS_S7_RETURN_RESEL=TRUE)       # campaign_p12x20.sh:24-27 (FIELD_RECOV unset)
MR_DG=(FS_S7_FIELD_RECOV=TRUE)        # campaign.sh:9-11, grfmr.sh:19-21, grfmrC.sh:17-19
BASE=("${PIN[@]}" FS_S7_MR=TRUE FS_S7_DGM=null FS_S7_CAMPAIGN=$CAMP FS_S7_MODE=batch
      FS_S7_START=1 FS_S7_NSIMS=$REPS FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20)
[[ "$QR" == "TRUE" ]] && BASE+=(FS_S7_QUICKRUN=TRUE)

T0=$(date +%s)
WATCHDOG=""
if (( ! DRY )); then
  MAINPID=$$
  ( sleep $TIMEOUT_S
    echo "=== HARD TIMEOUT ${TIMEOUT_S}s REACHED $(date): killing the render and the driver ==="
    pkill -f "quarto render sim_fs_maxeffCons_fb_mr_field_m1_template.qmd"
    pkill -f "sim_fs_maxeffCons_fb_mr_field_m1_template"
    kill $MAINPID ) &
  WATCHDOG=$!
fi
stop_run () { [[ -n "$WATCHDOG" ]] && kill $WATCHDOG 2>/dev/null; exit $1; }
echo "=== NULLMR START $(date) ; campaign $CAMP ; MR on ; thresholds unset (c1 $C1 c2 $C2) ; reps $REPS ; quickrun $QR ; hard timeout ${TIMEOUT_S}s (24 h), render timeout ${RENDER_TIMEOUT_S}s (6 h) ; watchdog pid ${WATCHDOG:-none} ; dry=$DRY ==="
echo "  host $(hostname) ; $(lscpu | awk -F: '/^Core\(s\) per socket/{c=$2} /^Socket\(s\)/{s=$2} END{gsub(/ /,"",c); gsub(/ /,"",s); print c*s}') physical cores ; $(free -g | awk '/^Mem:/{print $2}') GB ; R $(Rscript -e 'cat(paste0(R.version$major,".",R.version$minor))') ; forestsearch $(Rscript -e 'cat(as.character(packageVersion("forestsearch")), packageDescription("forestsearch")$Built)') ; git $(git -C "$REPO" rev-parse --short HEAD) ; workers $WORKERS"

commit_paths () {  # $1 message, then repo-relative paths; commits ONLY those paths
  local msg=$1; shift
  (( DOCOMMIT )) || return 0
  ( cd "$REPO" && git add -- "$@" && git commit -q -m "$msg$TRAILER" -- "$@" ) \
    && echo "  committed $(git -C "$REPO" rev-parse --short HEAD): $msg" \
    || echo "  COMMIT FAILED: $msg"
}

halt_record () {  # $1=cell $2=stage $3=detail ; echoes the file name
  local F="$QD/HALT_${CAMP}_$1_$(date +%Y-%m-%d).txt"
  {
    echo "HALT RECORD -- $CAMP, cell $1"
    echo "  task     : TASK_null_gbsg_mr_2026-09-21"
    echo "  when     : $(date)"
    echo "  host     : $(hostname) ; git $(git -C "$REPO" rev-parse --short HEAD)"
    echo "  stage    : $2"
    echo "  detail   : $3"
    echo "  action   : the rest of this cell is abandoned; the campaign continues at the next cell."
  } > "$F"
  echo "=== HALT RECORD WRITTEN: $F ===" >&2
  echo "$F"
}

stem () {  # $1=engine $2=hr $3=n
  local tag=$1; [[ $1 == consistency ]] && tag=fs
  local qr=""; [[ "$QR" == "TRUE" ]] && qr="_quickrun"
  printf 'results/%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_null%03d_nb20_%s%s_res_1_%d.rds' \
    "$tag" "$(awk -v h="$2" 'BEGIN{printf "%d", h*100+0.5}')" "$3" \
    "$(awk -v h="$2" 'BEGIN{printf "%d", h*1000+0.5}')" "$CAMP" "$qr" "$REPS"
}

run () {  # $1=engine $2=hr $3=n $4=cell
  local E=$1 H=$2 N=$3 T=$4
  local KN=("${BASE[@]}" FS_S7_HR=$H FS_S7_N=$N)
  if [[ $E == consistency ]]; then KN+=("${MR_FS[@]}"); else KN+=(FS_S7_METHOD=$E "${MR_DG[@]}"); fi
  local OUT="${CAMP}_${T}_${E}_effMaxSG"
  echo "=== RUN START: $T $E $(date) -> $OUT ==="
  echo "  env set:   ${KN[*]}"
  echo "  env unset: ${UNSET[*]}"
  (( DRY )) && return 0
  env "${UNSET[@]}" "${KN[@]}" timeout $RENDER_TIMEOUT_S bash "$RENDER" "$OUT" < /dev/null
  local RC=$?
  if (( RC == 124 )); then
    pkill -f "sim_fs_maxeffCons_fb_mr_field_m1_template"
    echo "RUN FAILED (render timeout ${RENDER_TIMEOUT_S}s): $T $E"; return 3
  fi
  (( RC != 0 )) && { echo "RUN FAILED (render rc $RC): $T $E"; return 1; }
  local GL="logs/${CAMP}_gateA_${T}_${E}.log"
  Rscript nullmr_gateA.R $E $H $N $T $C1 $C2 $CAMP $REPS $QR > $GL 2>&1 < /dev/null
  local GA=$?
  grep -E '^GATE A ' $GL
  if (( GA != 0 )); then grep -E '^\[FAIL\]' $GL; echo "=== GATE A FAILED: $OUT ==="; return 2; fi
  return 0
}

CI=0; DONE=(); FAILED=()
while read -r H HT N T; do
  [[ -z "$T" || "$H" == \#* ]] && continue
  CI=$(( CI + 1 ))
  echo "=== CELL START: $T (hr=$H n=$N) $(date) ; elapsed $(( $(date +%s) - T0 ))s ==="
  (( DRY )) || free -g | awk '/^Mem:|^Swap:/'
  CT0=$(date +%s); CELL_OK=1
  for E in consistency dina grf; do
    run $E $H $N $T
    RRC=$?
    if (( RRC != 0 )); then
      CELL_OK=0
      HF=$(halt_record $T "run $E (rc $RRC)" "the $E run of cell $T failed; see logs/${CAMP}_${T}_${E}_effMaxSG.log and logs/${CAMP}_gateA_${T}_${E}.log")
      echo "=== CELL $T ABANDONED at $E (rc $RRC) $(date); continuing at the next cell ==="
      break
    fi
  done
  if (( ! CELL_OK )); then
    FAILED+=($T)
    LP=(); for E in consistency dina grf; do
      for f in "logs/${CAMP}_${T}_${E}_effMaxSG.log" "logs/${CAMP}_gateA_${T}_${E}.log"; do
        [[ -f $f ]] && LP+=("$RQ/scripts_dinamr/$f"); done; done
    commit_paths "$CAMP $T: halt record and logs (cell abandoned)" "$RQ/$(basename "$HF")" "${LP[@]}"
    continue
  fi
  echo "CELL DONE: $T  wall=$(( $(date +%s) - CT0 ))s"
  if (( ! DRY )); then
    GCL="logs/${CAMP}_gateC_${T}.txt"
    Rscript nullmr_gateC.R $T $H $N $C1 $C2 $CAMP $REPS $QR > $GCL 2>&1 < /dev/null
    GC=$?
    grep -E '^GATE C ' $GCL
    PATHS=()
    for E in consistency dina grf; do
      PATHS+=("$RQ/$(stem $E $H $N)" "$RQ/scripts_dinamr/logs/${CAMP}_${T}_${E}_effMaxSG.log"
              "$RQ/scripts_dinamr/logs/${CAMP}_gateA_${T}_${E}.log")
    done
    PATHS+=("$RQ/scripts_dinamr/$GCL")
    if (( GC != 0 )); then
      grep -E '^\[FAIL\]' $GCL
      HF=$(halt_record $T "GATE C" "the per-cell gate failed; see $GCL")
      commit_paths "$CAMP $T: GATE C FAILED -- halt record, bundles and logs" "$RQ/$(basename "$HF")" "${PATHS[@]}"
      FAILED+=($T); echo "=== GATE C FAILED at $T; continuing at the next cell ==="; continue
    fi
    commit_paths "$CAMP $T: cell complete (3 bundles, render + gate logs; Gate A x3 and Gate C PASS)" "${PATHS[@]}"
  fi
  DONE+=($T)
done < "$1"
echo "NULLMR COMPLETE $(date) total_wall=$(( $(date +%s) - T0 ))s ; completed ${#DONE[@]} cells: ${DONE[*]:-none} ; failed ${#FAILED[@]}: ${FAILED[*]:-none}"
stop_run 0
