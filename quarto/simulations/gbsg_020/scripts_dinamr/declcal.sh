#!/bin/bash
# declcal driver (TASK_declcal_CAMPAIGN_2026-09-22_v2): the declaration-
# calibration finite-sample evaluation, FS only, c1 = c2 = 1.0.
# Transplant of nullmr.sh -- the pop-os (bash, lscpu / free) descendant of the
# nullid driver nullid.sh on the same grid; nullid.sh itself reads macOS sysctl
# and cannot run here.
# usage: declcal.sh <cells-file> [<cells-file> ...] [--dry]
#   cells file: lines "<cell> <dgm> <hr> <n> <dgm-label>" in RUN ORDER; '#'
#   lines and blanks skipped.  One cells file = one block = one campaign tag
#   (the file's basename without .cells: declcal_bnull, declcal_inull,
#   declcal_power, declcal_pilot).
# Environment:
#   DECLCAL_BCAL      multiplier draws per replicate (B_cal; default 2000)
#   DECLCAL_MODE      pilot | stage2 (default stage2)
#   DECLCAL_REPS      replicates per cell (default 2000; rep 1..REPS)
#   DECLCAL_CAP_S     per-replicate hard cap in seconds (0 = none)
#   DECLCAL_HARDCAP_S campaign hard cap in seconds (default 86400)
#   DECLCAL_COMMIT    1 commits each completed block by named paths (default 1)
#
# Differences from nullmr.sh, and nothing else:
#   1. One engine (FS; FS_S7_METHOD unset), and the per-cell work is
#      declcal_run.R (the template transplant; see its header) run by Rscript,
#      not a quarto render of the template.
#   2. Knobs: FS_S7_C1=1.0 FS_S7_C2=1.0, FS_S7_DGM / FS_S7_HR / FS_S7_N per cell,
#      FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 (nullmr's pinned focus and band),
#      64 workers requested (the template caps at physical cores - 1).  No MR
#      knobs: MR post-selection correction is off; only the field is captured.
#   3. Gates: the fidelity gate (declcal_run.R exit 3) STOPS THE CAMPAIGN with
#      the disagreeing replicates recorded; the per-cell 1% abort/error gate
#      (exit 4) stops that cell, which is recorded, and the campaign continues;
#      any other failure writes a halt record and the campaign continues.
#   4. Timeouts: DECLCAL_CAP_S per replicate (inside declcal_run.R) and
#      DECLCAL_HARDCAP_S for the campaign (watchdog).  Each cell's payload is
#      checkpointed to disk per 500-replicate chunk.
#   5. Commit per BLOCK (task section 11, commit 4 "per block"), by named paths.
#   6. Refuses to start while another R / quarto campaign is running.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"
QD="$(cd "$SP/.." && pwd)"
REPO=$(git -C "$QD" rev-parse --show-toplevel)
RQ=${QD#"$REPO"/}
cd "$SP" || exit 1
mkdir -p logs
BCAL=${DECLCAL_BCAL:-2000}
MODE=${DECLCAL_MODE:-stage2}
REPS=${DECLCAL_REPS:-2000}
CAP_S=${DECLCAL_CAP_S:-0}
TIMEOUT_S=${DECLCAL_HARDCAP_S:-86400}
DOCOMMIT=${DECLCAL_COMMIT:-1}
WORKERS=64
DRY=0; FILES=()
for a in "$@"; do if [[ "$a" == "--dry" ]]; then DRY=1; else FILES+=("$a"); fi; done
TRAILER=$'\n\nCo-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>'
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # render.sh's three

LIVE=$(pgrep -af 'quarto render|Rscript|exec/R ' | grep -v -e "$$" -e pgrep || true)
if [[ -n "$LIVE" ]] && (( ! DRY )); then
  echo "=== ANOTHER R / QUARTO PROCESS IS RUNNING -- not competing for cores; stopping ==="
  echo "$LIVE"; exit 5
fi

UNSET=()
for v in $(env | grep -o '^FS_S7_[A-Z0-9_]*'); do UNSET+=(-u "$v"); done
BASE=(FS_S7_WORKERS=$WORKERS FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_C1=1.0 FS_S7_C2=1.0
      FS_S7_START=1 FS_S7_NSIMS=$REPS)

T0=$(date +%s)
WATCHDOG=""
if (( ! DRY )); then
  MAINPID=$$
  ( sleep $TIMEOUT_S
    echo "=== CAMPAIGN HARD CAP ${TIMEOUT_S}s REACHED $(date): killing the cell and the driver ==="
    pkill -f "declcal_run.R"
    kill $MAINPID ) &
  WATCHDOG=$!
fi
stop_run () { [[ -n "$WATCHDOG" ]] && kill $WATCHDOG 2>/dev/null; exit $1; }
echo "=== DECLCAL START $(date) ; mode $MODE ; B_cal $BCAL ; reps $REPS ; per-rep cap ${CAP_S}s ; campaign hard cap ${TIMEOUT_S}s ; watchdog pid ${WATCHDOG:-none} ; dry=$DRY ; blocks: ${FILES[*]} ==="
echo "  host $(hostname) ; $(lscpu | awk -F: '/^Core\(s\) per socket/{c=$2} /^Socket\(s\)/{s=$2} END{gsub(/ /,"",c); gsub(/ /,"",s); print c*s}') physical cores ; $(free -g | awk '/^Mem:/{print $2}') GB ; R $(Rscript -e 'cat(paste0(R.version$major,".",R.version$minor))') ; forestsearch $(Rscript -e 'cat(as.character(packageVersion("forestsearch")), packageDescription("forestsearch")$Built)') ; git $(git -C "$REPO" rev-parse --short HEAD) ; workers requested $WORKERS"

commit_paths () {  # $1 message, then repo-relative paths; commits ONLY those paths
  local msg=$1; shift
  (( DOCOMMIT )) || return 0
  ( cd "$REPO" && git add -- "$@" && git commit -q -m "$msg$TRAILER" -- "$@" ) \
    && echo "  committed $(git -C "$REPO" rev-parse --short HEAD): $msg" \
    || echo "  COMMIT FAILED: $msg"
}

halt_record () {  # $1=camp $2=cell $3=stage $4=detail $5=action ; echoes the file name
  local F="$QD/HALT_$1_$2_$(date +%Y-%m-%d).txt"
  {
    echo "HALT RECORD -- $1, cell $2"
    echo "  task     : TASK_declcal_CAMPAIGN_2026-09-22_v2"
    echo "  when     : $(date)"
    echo "  host     : $(hostname) ; git $(git -C "$REPO" rev-parse --short HEAD)"
    echo "  stage    : $3"
    echo "  detail   : $4"
    echo "  action   : $5"
  } > "$F"
  echo "=== HALT RECORD WRITTEN: $F ===" >&2
  echo "$F"
}

DONE=(); FAILED=()
for CF in "${FILES[@]}"; do
  CAMP=$(basename "$CF" .cells)
  BPATHS=()
  echo "=== BLOCK START: $CAMP ($CF) $(date) ; elapsed $(( $(date +%s) - T0 ))s ==="
  while read -r CELL DGM HR N LAB; do
    [[ -z "$LAB" || "$CELL" == \#* ]] && continue
    OUT="results/${CAMP}_${CELL}_res_1_${REPS}.rds"
    LOG="logs/${CAMP}_${CELL}.log"
    KN=("${BASE[@]}" FS_S7_DGM=$DGM FS_S7_HR=$HR FS_S7_N=$N FS_S7_CAMPAIGN=$CAMP
        DECLCAL_CELL=$CELL DECLCAL_DGM_LABEL=$LAB DECLCAL_BCAL=$BCAL DECLCAL_MODE=$MODE
        DECLCAL_CAP_S=$CAP_S DECLCAL_OUT=$QD/$OUT)
    echo "=== CELL START: $CAMP $CELL ($LAB, hr=$HR n=$N) $(date) ; elapsed $(( $(date +%s) - T0 ))s ==="
    echo "  env set:   ${KN[*]}"
    (( DRY )) && continue
    free -g | awk '/^Mem:|^Swap:/'
    CT0=$(date +%s)
    ( cd "$QD" && env "${UNSET[@]}" "${KN[@]}" Rscript "$SP/declcal_run.R" ) > "$LOG" 2>&1 < /dev/null
    RC=$?
    echo "CELL END: $CAMP $CELL rc=$RC wall=$(( $(date +%s) - CT0 ))s"
    grep -E '^(CELL|FIDELITY|RATES|chunk)' "$LOG"
    [[ -f "$QD/$OUT" ]] && BPATHS+=("$RQ/$OUT")
    BPATHS+=("$RQ/scripts_dinamr/$LOG")
    if (( RC == 3 )); then
      HF=$(halt_record $CAMP $CELL "FIDELITY GATE" "declared_conv != the search's own indicator; replicates on the FIDELITY line of $LOG" "CAMPAIGN STOPPED (task section 2.2): the cost design rests on this gate.")
      commit_paths "$CAMP $CELL: FIDELITY GATE FAILED -- halt record, payload and log" "$RQ/$(basename "$HF")" "${BPATHS[@]}"
      echo "=== FIDELITY GATE FAILED at $CAMP $CELL; campaign stopped ==="
      stop_run 3
    elif (( RC == 4 )); then
      HF=$(halt_record $CAMP $CELL "PER-CELL GATE" "more than 1% of replicates ended abort_time or error; see $LOG" "this cell is stopped and no rate is reported from it; the campaign continues.")
      BPATHS+=("$RQ/$(basename "$HF")"); FAILED+=("$CAMP:$CELL"); continue
    elif (( RC != 0 )); then
      HF=$(halt_record $CAMP $CELL "run (rc $RC)" "declcal_run.R failed; see $LOG" "this cell is abandoned; the campaign continues.")
      BPATHS+=("$RQ/$(basename "$HF")"); FAILED+=("$CAMP:$CELL"); continue
    fi
    DONE+=("$CAMP:$CELL")
  done < "$CF"
  (( DRY )) || commit_paths "data(sims): add declcal campaign payloads -- block $CAMP" "${BPATHS[@]}"
done
echo "DECLCAL COMPLETE $(date) total_wall=$(( $(date +%s) - T0 ))s ; completed ${#DONE[@]} cells: ${DONE[*]:-none} ; failed ${#FAILED[@]}: ${FAILED[*]:-none}"
stop_run 0
