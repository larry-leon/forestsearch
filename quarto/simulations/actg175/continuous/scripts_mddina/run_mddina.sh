#!/bin/bash
# mddina Stage 2 runner (TASK_md_dina_campaign_2026-09-17 Stage 2): transplant of scripts_mdgrf/run_mdgrf.sh
# with the DINA identifier knob (FS_MD_METHOD=dina; FS_MD_DMIN_GRF dropped -- a GRF-only argument,
# left at the template default and recorded, inert on DINA), the dina_ stem, the tag, the date, the
# Gate 2 report name, scripts_mddina/gate2.R, and NO hard-coded model or co-author trailer on the
# commit messages (the task's carried fix).  Everything else is the GRF runner's, itself the FS runner's.
# Launch from the repo root or anywhere (paths are resolved from this file):
#   MDSG_WORKERS=<W> MDSG_TIMEOUT=<s> MDSG_CEILING=<s> setsid nohup bash <dir>/scripts_mddina/run_mddina.sh > <dir>/logs_mddina/runner.log 2>&1 &
# Restart after an interruption with the same command: done cells are skipped.
#
# Transplant of quarto/simulations/gbsg_020/scripts_p12x20/run_p12x20.sh (preflight,
# per-cell sequencing, gate, explicit-path commits, heartbeats, halt), with:
#   1. the render lines of scripts_mdf1/run_cell.sh:18-31 (batch 1-1000, batch
#      1001-2000, combine, each a `quarto render` of the MD template);
#   2. GNU `timeout` (MDSG_TIMEOUT seconds) in place of scripts_mdf1/tmo.sh;
#   3. the threads exported here (render.sh's OMP/OPENBLAS/MKL = 1);
#   4. peak memory per render from scripts_mdsgnb20/mem_sampler.sh;
#   5. the progress log <dir>/LOG_mdsgnb20_progress.txt (each render's wall and
#      peak memory) and the halt file <dir>/HALT_mdsgnb20.md;
#   6. gate2.R (the §2.3 checks) in place of gate2F.R;
#   7. commits per cell by explicit path: the cell's result directory, its combine
#      HTML (as mdf1 tracked it), its section of REPORT_md_field_rerun_gate2_2026-09-15.md,
#      and the progress log.  Raw render logs stay untracked under <dir>/logs_mdsgnb20/.
#   8. the go's ceiling (MDSG_CEILING seconds, cumulative render wall): finish the
#      current cell, then halt.
set -u
: "${MDSG_WORKERS:?MDSG_WORKERS must be set (Gate 1 go)}"
: "${MDSG_TIMEOUT:?MDSG_TIMEOUT must be set (seconds per render; Gate 1 go)}"
: "${MDSG_CEILING:?MDSG_CEILING must be set (seconds, cumulative; Gate 1 go)}"
TAG=mddina
DATE=2026-09-17
BRANCH=feature/glm-extension
# Cells in mdf1's order: <cell-label> <md> <n>
CELLS=("md40_n500 40 500" "md120_n500 120 500" "null_n500 null 500" "md40_n700 40 700")
HERE=$(cd "$(dirname "$0")" && pwd)
QMD=$(cd "$HERE/.." && pwd)
REPO=$(git -C "$QMD" rev-parse --show-toplevel)
RQ=${QMD#"$REPO"/}
LOGD=$QMD/logs_${TAG}
HB=$QMD/LOG_${TAG}_progress.txt
HALTF=$QMD/HALT_${TAG}.md
GREP=$QMD/REPORT_md_dina_gate2_${DATE}.md
export PATH=/usr/lib/rstudio/resources/app/bin/quarto/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p "$LOGD"
cd "$REPO" || exit 1
CUM=0

utc () { date -u +%Y-%m-%dT%H:%M:%SZ; }
log () { echo "[$(utc)] $*"; }
commit () {  # $1 message, then explicit repo-relative paths; commits ONLY those paths
  local msg=$1; shift
  git add -- "$@" && git commit -q -m "$msg" -- "$@"
}
heartbeat () {  # $1 cell  $2 label  $3 text
  printf '%s\t%s\t%s\t%s\n' "$(utc)" "$1" "$2" "$3" >> "$HB"
}
halt () {  # $1 cell  $2 assertion  [$3 file holding the observed values]
  local cell=$1 what=$2 obs=${3:-}
  {
    echo "# HALT — $TAG Stage 2"
    echo
    echo "- Halted at (UTC): $(utc)"
    echo "- Cell: $cell"
    echo "- Assertion: $what"
    echo "- HEAD at halt: $(git rev-parse --short HEAD)"
    echo "- Cumulative render wall: ${CUM} s (ceiling ${MDSG_CEILING} s)"
    echo "- Completed cells stay committed; nothing is rolled back."
    echo
    echo "## Observed values"
    echo
    echo '```'
    [ -n "$obs" ] && [ -f "$obs" ] && cat "$obs"
    echo '```'
  } > "$HALTF"
  heartbeat "$cell" halt "$what"
  commit "$TAG $cell: HALT -- $what" "$RQ/HALT_${TAG}.md" "$RQ/LOG_${TAG}_progress.txt" || log "HALT commit failed"
  log "HALTED: $cell: $what"
  exit 1
}
tracked_clean () { git ls-files --error-unmatch -- "$@" >/dev/null 2>&1 && git diff --quiet HEAD -- "$@"; }
render () {  # $1 cell  $2 stem  $3 label  $4 mode  $5 start  $6 nsims  $7 md  $8 n  -> sets WALL, PEAK
  local cell=$1 stem=$2 label=$3 mode=$4 start=$5 nsims=$6 md=$7 n=$8
  local out="${stem}_${label}.html" lg="$LOGD/${cell}_${label}.log" pk="$LOGD/${cell}_${label}.peak_mb"
  local t0=$(date +%s)
  ( cd "$QMD" && \
    FS_MD_METHOD=dina FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=$TAG FS_MD_FB=none \
    FS_MD_MD=$md FS_MD_N=$n FS_MD_MODE=$mode FS_MD_START=$start FS_MD_NSIMS=$nsims FS_MD_WORKERS=$MDSG_WORKERS \
    timeout "$MDSG_TIMEOUT" quarto render sim_fs_maxeffCons_mr_field_md_template.qmd --output "$out" > "$lg" 2>&1 ) &
  local rpid=$!
  bash "$HERE/mem_sampler.sh" "$rpid" "$pk" &
  local spid=$!
  wait "$rpid"; local rc=$?
  wait "$spid" 2>/dev/null
  WALL=$(( $(date +%s) - t0 )); PEAK=$(cat "$pk" 2>/dev/null || echo NA); CUM=$(( CUM + WALL ))
  echo "WALL_SECONDS=$WALL RC=$rc PEAK_MB=$PEAK OUT=$out" >> "$lg"
  heartbeat "$cell" "$label" "wall_s=$WALL peak_mb=$PEAK rc=$rc cumulative_s=$CUM"
  [ $rc -eq 0 ] || { tail -40 "$lg" > "$LOGD/${cell}_${label}.halt_obs"; halt "$cell" "render $label failed (rc=$rc; 124 = timeout ${MDSG_TIMEOUT}s)" "$LOGD/${cell}_${label}.halt_obs"; }
}

# ---- preflight --------------------------------------------------------------
[ "$(git rev-parse --abbrev-ref HEAD)" = "$BRANCH" ] || { log "PREFLIGHT: not on $BRANCH"; exit 1; }
[ -f "$HALTF" ] && { log "PREFLIGHT: $RQ/HALT_${TAG}.md exists; the campaign is halted"; exit 1; }
git diff --quiet HEAD -- "$REPO/R" "$QMD/sim_fs_maxeffCons_mr_field_md_template.qmd" || { log "PREFLIGHT: R/ or the template has uncommitted changes"; exit 1; }
pkill -f '[w]orkRSOCK'; sleep 5
NW=$(pgrep -fc '[w]orkRSOCK'); [ "$NW" = "0" ] || { log "PREFLIGHT: worker reap left $NW workRSOCK processes"; exit 1; }
BUILT=$(Rscript -e 'cat(packageDescription("forestsearch")$Built)')
log "PREFLIGHT OK: branch $BRANCH, HEAD $(git rev-parse --short HEAD), workers $MDSG_WORKERS, timeout ${MDSG_TIMEOUT}s, ceiling ${MDSG_CEILING}s, Built '$BUILT'"
heartbeat campaign start "HEAD=$(git rev-parse --short HEAD) workers=$MDSG_WORKERS timeout_s=$MDSG_TIMEOUT ceiling_s=$MDSG_CEILING built=$BUILT"

# ---- cells ------------------------------------------------------------------
for spec in "${CELLS[@]}"; do
  read -r CELL MD N <<< "$spec"
  MDTOK=$([ "$MD" = null ] && echo mdnull || echo "md$MD")
  STEM=dina_effMaxSG_mr_field_${MDTOK}_knoise0_n${N}_nb20_${TAG}
  RDIR=$QMD/mr_md_harm/${STEM}_d5000
  GREC=$LOGD/GATE2_${CELL}.txt
  if tracked_clean "$RQ/mr_md_harm/${STEM}_d5000/${STEM}_combined_1_2000.rds"; then
    log "SKIP $CELL: combined bundle committed"; continue
  fi
  if [ "$CUM" -ge "$MDSG_CEILING" ]; then
    echo "cumulative $CUM s >= ceiling $MDSG_CEILING s before $CELL" > "$LOGD/${CELL}.halt_obs"
    halt "$CELL" "ceiling reached before the cell" "$LOGD/${CELL}.halt_obs"
  fi
  CELL_T0=$(date +%s)
  heartbeat "$CELL" start "stem=$STEM"
  # An interrupted cell restarts from its first batch: partial output is discarded
  # (the save guard refuses only git-tracked paths; these are untracked).
  rm -rf "$RDIR"; rm -f "$QMD/${STEM}_batch_1_1000.html" "$QMD/${STEM}_batch_1001_2000.html" "$QMD/${STEM}_combine_1_2000.html"
  log "START $CELL: md $MD n $N ($STEM)"
  render "$CELL" "$STEM" batch_1_1000    batch   1    1000 "$MD" "$N"
  render "$CELL" "$STEM" batch_1001_2000 batch   1001 1000 "$MD" "$N"
  render "$CELL" "$STEM" combine_1_2000  combine 1    1000 "$MD" "$N"
  mv "$QMD/${STEM}_batch_1_1000.html" "$QMD/${STEM}_batch_1001_2000.html" "$LOGD/" 2>/dev/null
  # Gate 2 (§2.3)
  (cd "$HERE" && MDSG_DIR=.. MDSG_TAG=$TAG MDSG_WORKERS=$MDSG_WORKERS Rscript gate2.R "$MD" "$N" "$CELL") > "$GREC" 2>&1
  GRC=$?
  COUNTS=$(grep '^GATE_COUNTS' "$GREC")
  [ $GRC -eq 0 ] || halt "$CELL" "gate2.R failed: ${COUNTS:-no GATE_COUNTS line}" "$GREC"
  {
    [ -f "$GREP" ] || { echo "# REPORT — $TAG Gate 2 (per cell)"; echo; echo "Task: \`dev/tasks/TASK_md_dina_campaign_${DATE}.md\` Stage 2 (the GRF task's §2.3 with this task's changes). Every coverage figure of this campaign is conditional on the proposed family (DINA's family is generated from a fitted surface). Runner: \`$RQ/scripts_${TAG}/run_${TAG}.sh\`; checker: \`$RQ/scripts_${TAG}/gate2.R\`."; echo; }
    echo "## $CELL — md $MD, n $N"
    echo
    echo "- Stem: \`$STEM\`; HEAD before this cell's commit: $(git rev-parse --short HEAD); workers $MDSG_WORKERS; threads 1."
    echo "- Knobs: \`FS_MD_METHOD=dina FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=$TAG FS_MD_FB=none FS_MD_MD=$MD FS_MD_N=$N FS_MD_WORKERS=$MDSG_WORKERS\`."
    echo "- Seeds: 8316951 + sim_id; batches sim_id 1-1000 and 1001-2000, then combine."
    echo "- Cell wall: $(( $(date +%s) - CELL_T0 )) s; cumulative render wall $CUM s (ceiling $MDSG_CEILING s)."
    grep -h '^WALL_SECONDS=' "$LOGD/${CELL}_"*.log | sed 's/^/- Render: /'
    echo "- Gate counts: \`$COUNTS\`."
    echo
    echo '```'
    cat "$GREC"
    echo '```'
    echo
  } >> "$GREP"
  PATHS=("$RQ/mr_md_harm/${STEM}_d5000" "$RQ/${STEM}_combine_1_2000.html" "$RQ/REPORT_md_dina_gate2_${DATE}.md" "$RQ/LOG_${TAG}_progress.txt")
  heartbeat "$CELL" done "cell_wall_s=$(( $(date +%s) - CELL_T0 )) $COUNTS"
  commit "$TAG $CELL: md $MD n $N, 2000 replicates + combine; Gate 2 PASS ($COUNTS)" "${PATHS[@]}" || { log "CELL COMMIT FAILED: $CELL"; exit 1; }
  log "DONE $CELL: $COUNTS, commit $(git rev-parse --short HEAD)"
done
heartbeat campaign complete "cumulative_s=$CUM"
commit "$TAG: campaign complete (cumulative render wall ${CUM} s)" "$RQ/LOG_${TAG}_progress.txt" || log "final heartbeat commit failed"
log "CAMPAIGN COMPLETE"
