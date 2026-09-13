#!/bin/bash
# p12x20 Part A Stage 2 runner (ADDENDUM_p12x20_stage2_unattended_2026-09-12_v2).
# Launch from this directory:
#   setsid nohup ./run_p12x20.sh > logs/p12x20_runner.log 2>&1 &
# Restart after any interruption with the same command: done cells are skipped.
#
# Transplant, not authorship.  The run and the combine are campaign_p12x20.sh
# (the committed scripts_dinamr/campaign.sh, bash port, Stage 1 §1b knobs)
# calling the committed scripts_dinamr/render.sh; the gates are gate2F.R (the
# committed scripts_dinamr/gate2.R repointed at FS, with the §5 combine
# assertions and the §6 two-direction same-draws).  This script only sequences
# the cells, discards partial output, commits, writes heartbeats and halts.

WORKERS=64          # REPORT_p12x20_gate1_2026-09-12.md, 1c (per-cell wall at 2,000 replicates)
THREADS=1           # OMP_NUM_THREADS / OPENBLAS_NUM_THREADS / VECLIB_MAXIMUM_THREADS, exported by render.sh
TAG=p12x20
DATE=2026-09-12     # the task's date: payload directory and report names carry it
BRANCH=campaign/p12x20
# Cells in task §5 order: <cell> <HR> <n>.  Prevalence 12.4% (FS_S7_Z1Q unset).
CELLS=("A1 1.50 500" "A2 1.50 1000" "A3 1.50 1500"
       "A4 1.75 500" "A5 1.75 1000" "A6 1.75 1500"
       "A7 1.00 500" "A8 1.00 1000" "A9 1.00 1500")

set -u
HERE=$(cd "$(dirname "$0")" && pwd)
QMD=$(cd "$HERE/.." && pwd)
REPO=$(git -C "$QMD" rev-parse --show-toplevel)
RQ=${QMD#"$REPO"/}
PAYDIR=${TAG}_${DATE}
PAY=$QMD/$PAYDIR
HB=$QMD/LOG_p12x20_progress.txt
OPEN=$QMD/OPEN_ITEMS_p12x20.md
HALTF=$QMD/HALT_p12x20.md
export PATH=/usr/lib/rstudio/resources/app/bin/quarto/bin:$PATH
export DINAMR_SCRATCH=$HERE P12X20_WORKERS=$WORKERS
TRAILER=$'\n\nCo-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_01WcgkteQefDFYbrj1eFqzrM'
mkdir -p "$HERE/logs" "$PAY"
cd "$REPO" || exit 1
CELL_T0=$(date +%s)

utc () { date -u +%Y-%m-%dT%H:%M:%SZ; }
log () { echo "[$(utc)] $*"; }
commit () {  # $1 message, then explicit repo-relative paths; commits ONLY those paths
  local msg=$1; shift
  git add -- "$@" && git commit -q -m "$msg$TRAILER" -- "$@"
}
heartbeat () {  # $1 cell  $2 start|done|halt  $3 elapsed minutes  [$4 commit sha]
  printf '%s\t%s\t%s\telapsed_min=%s%s\n' "$(utc)" "$1" "$2" "$3" "${4:+	commit=$4}" >> "$HB"
  commit "p12x20 $1: heartbeat $2" "$RQ/LOG_p12x20_progress.txt" || log "heartbeat commit failed ($1 $2)"
}
status_regen () {
  Rscript "$HERE/status_p12x20.R" > "$HERE/logs/p12x20_status.log" 2>&1 || log "STATUS regeneration failed (logs/p12x20_status.log)"
  [ -f "$QMD/STATUS_p12x20.md" ] && { commit "p12x20: regenerate STATUS_p12x20.md" "$RQ/STATUS_p12x20.md" || log "STATUS commit failed"; }
}
halt () {  # $1 cell  $2 assertion  [$3 file holding the observed values]
  local cell=$1 what=$2 obs=${3:-}
  {
    echo "# HALT — p12x20 Part A Stage 2"
    echo
    echo "- Halted at (UTC): $(utc)"
    echo "- Cell: $cell"
    echo "- Assertion: $what"
    echo "- HEAD at halt: $(git rev-parse --short HEAD)"
    echo "- Completed cells stay committed; nothing is rolled back."
    echo
    echo "## Observed values"
    echo
    echo '```'
    [ -n "$obs" ] && [ -f "$obs" ] && cat "$obs"
    echo '```'
  } > "$HALTF"
  commit "p12x20 $cell: HALT -- $what" "$RQ/HALT_p12x20.md" || log "HALT commit failed"
  heartbeat "$cell" halt "$(( ($(date +%s) - CELL_T0) / 60 ))"
  status_regen
  log "HALTED: $cell: $what"
  exit 1
}
open_item () {  # $1 cell  $2 text   (documentary: record and continue)
  [ -f "$OPEN" ] || printf '# OPEN ITEMS — p12x20 Part A Stage 2\n\n' > "$OPEN"
  echo "- $(utc) $1: $2" >> "$OPEN"
}
tracked_clean () { git ls-files --error-unmatch -- "$@" >/dev/null 2>&1 && git diff --quiet HEAD -- "$@"; }

# ---- preflight --------------------------------------------------------------
[ "$(git rev-parse --abbrev-ref HEAD)" = "$BRANCH" ] || { log "PREFLIGHT: not on $BRANCH"; exit 1; }
[ -f "$HALTF" ] && { log "PREFLIGHT: $RQ/HALT_p12x20.md exists; the campaign is halted"; exit 1; }
for v in OMP_NUM_THREADS OPENBLAS_NUM_THREADS VECLIB_MAXIMUM_THREADS; do
  grep -q "^export $v=$THREADS\$" "$QMD/scripts_dinamr/render.sh" \
    || { log "PREFLIGHT: render.sh does not export $v=$THREADS"; exit 1; }
done
# Worker reap (addendum §1, §3).  '[w]orkRSOCK' matches the same processes as
# 'workRSOCK' without matching a shell whose own command line carries the pattern.
pkill -f '[w]orkRSOCK'
sleep 5
NW=$(pgrep -fc '[w]orkRSOCK')
[ "$NW" = "0" ] || { log "PREFLIGHT: worker reap left $NW workRSOCK processes"; exit 1; }
log "PREFLIGHT OK: branch $BRANCH, HEAD $(git rev-parse --short HEAD), workers $WORKERS, threads $THREADS, workRSOCK 0"

# ---- cells ------------------------------------------------------------------
for spec in "${CELLS[@]}"; do
  read -r CELL HR N <<< "$spec"
  H3=$(awk -v h="$HR" 'BEGIN { printf "%03d", int(h * 100 + 0.5) }')
  CTAG=${CELL}_h${H3}_n${N}
  STEM=fs_effMaxSG_fb_mr_field_m1_h${H3}_knoise0_n${N}_nb20_${TAG}
  GREC_NAME=GATE2_${TAG}_${CELL}_${DATE}.txt
  CREP_NAME=REPORT_${TAG}_${CELL}_${DATE}.md
  GREC=$PAY/$GREC_NAME
  CREP=$PAY/$CREP_NAME

  if tracked_clean "$RQ/$PAYDIR/${STEM}_combined_1_2000.rds" "$RQ/$PAYDIR/$GREC_NAME"; then
    log "SKIP $CELL: combined payload and gate record committed"
    continue
  fi

  CELL_T0=$(date +%s)
  CELL_START_UTC=$(utc)
  heartbeat "$CELL" start 0
  # An interrupted cell restarts from its first batch: partial output is discarded.
  rm -f "$QMD/results/${STEM}_res_"*.rds "$QMD/results/${STEM}_combined_"*.rds \
        "$QMD/p12x20_${CTAG}_"*.html "$PAY/${STEM}_"* "$GREC" "$GREC.tmp" "$CREP"

  LOGC=$HERE/logs/${TAG}_${CELL}.log
  SPEC=$HERE/logs/${TAG}_${CELL}.cells
  echo "- $N $HR $CTAG" > "$SPEC"
  log "START $CELL: HR $HR n $N ($STEM), log $LOGC"
  bash "$HERE/campaign_p12x20.sh" "$CELL" "$SPEC" > "$LOGC" 2>&1
  if ! grep -q "^CELL DONE: $CTAG " "$LOGC" || [ ! -f "$QMD/results/${STEM}_combined_1_2000.rds" ]; then
    tail -60 "$LOGC" > "$HERE/logs/${TAG}_${CELL}.halt_obs"
    halt "$CELL" "run/combine did not complete (CELL DONE absent or combined payload missing)" "$HERE/logs/${TAG}_${CELL}.halt_obs"
  fi
  WALL=$(sed -n "s/^CELL DONE: $CTAG  wall=\([0-9]*\)s.*/\1/p" "$LOGC")

  # Combine assertions -> Gate 2 -> same-draws, both directions -> payload size.
  (cd "$HERE" && Rscript gate2F.R "$HR" "$N" "$CELL") > "$GREC.tmp" 2>&1
  GRC=$?
  COUNTS=$(grep '^GATE_COUNTS' "$GREC.tmp")
  [ $GRC -eq 0 ] || halt "$CELL" "gate2F.R failed: ${COUNTS:-no GATE_COUNTS line}" "$GREC.tmp"
  mv "$GREC.tmp" "$GREC"

  # Payloads and renders into the dated directory.
  FILES=()
  for f in "${STEM}_res_1_1000.rds" "${STEM}_res_1001_2000.rds" "${STEM}_combined_1_2000.rds"; do
    mv "$QMD/results/$f" "$PAY/$f" && FILES+=("$f")
  done
  for f in "p12x20_${CTAG}_batch_1.html" "p12x20_${CTAG}_batch_1001.html" "p12x20_${CTAG}_combine_1.html"; do
    if [ -f "$QMD/$f" ]; then mv "$QMD/$f" "$PAY/$f" && FILES+=("$f")
    else open_item "$CELL" "render $f not found after the run"; fi
  done
  PATHS=(); SIZES=""; OBS=$HERE/logs/${TAG}_${CELL}.size_obs; : > "$OBS"
  for f in "${FILES[@]}" "$GREC_NAME"; do
    sz=$(stat -c %s "$PAY/$f")
    echo "$f $sz B" >> "$OBS"
    [ "$sz" -gt 104857600 ] && halt "$CELL" "artifact over the 100 MB hard stop: $f ($sz B)" "$OBS"
    flag=""; [ "$sz" -gt 52428800 ] && flag="  **FLAG: over 50 MB**"
    PATHS+=("$RQ/$PAYDIR/$f")
    SIZES+="- \`$RQ/$PAYDIR/$f\`: $sz B$flag"$'\n'
  done
  REPS=$(grep -E '^  2,000 rows ' "$GREC" | grep -oE '\([0-9]+\)' | tr -d '()')

  {
    echo "# REPORT — p12x20 $CELL: HR $HR, n $N, prevalence 12.4% (FS_S7_Z1Q unset)"
    echo
    echo "- Task: \`dev/tasks/TASK_p12x20_partA_2026-09-12_v2.md\` §5-§7; addendum \`dev/tasks/ADDENDUM_p12x20_stage2_unattended_2026-09-12_v2.md\`."
    echo "- Runner: \`$RQ/scripts_p12x20/run_p12x20.sh\`; HEAD before this cell's commit: $(git rev-parse --short HEAD)."
    echo "- Stem: \`$STEM\`."
    echo "- Knobs: \`FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=p12x20 FS_S7_RETURN_RESEL=TRUE FS_S7_WORKERS=$WORKERS FS_S7_HR=$HR FS_S7_N=$N\`; unset: \`FS_S7_Z1Q FS_S7_METHOD FS_S7_MR FS_S7_FIELD_RECOV FS_S7_ER_JCUTS\`."
    echo "- Workers: $WORKERS. Threads: $THREADS (render.sh)."
    echo "- Seeds: 8316951 + sim_id; two seed-disjoint batches, sim_id 1-1000 and 1001-2000, then combine."
    echo "- Cell start (UTC): $CELL_START_UTC."
    echo "- Cell wall-clock (driver \`CELL DONE\`): ${WALL:-not found} s."
    grep -h '^WALL_SECONDS=' "$HERE/logs/p12x20_${CTAG}_"*.log 2>/dev/null | sed 's/^/- Render: /'
    echo "- Replicates combined: ${REPS:-not found}."
    echo "- Gate counts: \`$COUNTS\`."
    echo
    echo "## Artifacts (committed in this cell's commit; size in bytes)"
    echo
    printf '%s' "$SIZES"
    echo "- \`$RQ/$PAYDIR/$CREP_NAME\`: this report."
    echo
    echo "## Gate record (verbatim, \`$GREC_NAME\`)"
    echo
    echo '```'
    cat "$GREC"
    echo '```'
  } > "$CREP"
  PATHS+=("$RQ/$PAYDIR/$CREP_NAME")
  [ -f "$OPEN" ] && PATHS+=("$RQ/OPEN_ITEMS_p12x20.md")

  commit "p12x20 $CELL: HR $HR n $N" "${PATHS[@]}" || { log "CELL COMMIT FAILED: $CELL"; exit 1; }
  SHA=$(git rev-parse --short HEAD)
  heartbeat "$CELL" done "$(( ($(date +%s) - CELL_T0) / 60 ))" "$SHA"
  log "DONE $CELL: wall ${WALL:-?} s, $COUNTS, commit $SHA"
done

status_regen
log "CAMPAIGN COMPLETE"
