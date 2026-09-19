#!/bin/bash
# ACTG175 binary / OR Stage 2 runner (TASK_actg175_binary_campaign_2026-09-17 Stage 2):
# transplant of ../../continuous/scripts_mddina/run_mddina.sh with
#   1. THREE campaigns in the task's order -- orfs (consistency), then orgrf, then ordina --
#      each over the same six cells in the task's order:
#        0.75/500, 0.75/2000, 1.5/500, 1.5/2000, 1.0/500, 1.0/2000;
#   2. the FS_OR_* knobs of binary_020/sim_fs_mr_field_or_template.qmd in place of FS_MD_*,
#      with FS_OR_TARGET / FS_OR_N carrying the cell and FS_OR_METHOD the identifier;
#   3. the or_ stem, the progress log LOG_or_progress.txt, the halt file HALT_or.md, raw logs
#      untracked under logs_or/, and the Gate 2 record REPORT_actg175_or_gate2_2026-09-17.md;
#   4. scripts_or/gate2.R as the per-cell checker (it takes the campaign tag as its 4th argument
#      and, for orgrf / ordina, compares the data-level columns with orfs's bundle in the same
#      cell in both directions);
#   5. NO hard-coded model or co-author trailer on any commit message (the task's carried fix).
# Everything else -- the preflight, the per-cell sequencing, the memory sampler, the explicit-path
# commits, the heartbeats, the halt-on-failure and the cumulative ceiling -- is the MD runner's.
# Launch from anywhere (paths are resolved from this file):
#   ORSG_WORKERS=<W> ORSG_TIMEOUT=<s> ORSG_CEILING=<s> [ORSG_NSIMS=<n>] \
#   [ORSG_CAMPAIGNS="orfs"] [ORSG_CELLS="or150_n500 or075_n500"] \
#     setsid nohup bash <dir>/scripts_or/run_or.sh > <dir>/logs_or/runner.log 2>&1 &
# Restart after an interruption with the same command: done cells are skipped.
set -u
: "${ORSG_WORKERS:?ORSG_WORKERS must be set (Gate 1 go)}"
: "${ORSG_TIMEOUT:?ORSG_TIMEOUT must be set (seconds per render; Gate 1 go)}"
: "${ORSG_CEILING:?ORSG_CEILING must be set (seconds, cumulative; Gate 1 go)}"
TAG=or
DATE=2026-09-17
BRANCH=feature/glm-extension
TEMPLATE=sim_fs_mr_field_or_template.qmd
# REPLICATES PER CELL -- fixed at 1,000 by TASK_binary_launch_v2_2026-09-18 ("replicates 1,000
# per cell ... do not change").  One batch of NSIMS over GLOBAL sim_id 1..NSIMS, then the combine,
# in place of the superseded 2 x 1,000.  The seed table is the study's pre-generated one with
# MAX_SIMS = 5000, so sim_id 1..1000 is covered for every cell and the draws are identical across
# the three identifiers cell for cell, exactly as at 2,000.
NSIMS=${ORSG_NSIMS:-1000}
# The three campaigns, in the task's order: <campaign-tag> <FS_OR_METHOD> <method-stem-token>.
# ORSG_CAMPAIGNS selects a subset by tag (space-separated) so one runner drives stage 1 (orfs)
# and stage 2a (orgrf ordina) without a second copy of this logic.
CAMPAIGNS_ALL=("orfs consistency fs" "orgrf grf grf" "ordina dina dina")
# The six cells.  ORDER, per TASK_binary_launch_v2 Step 3: orfs_or150_n500 FIRST (the
# configuration that halted twice), then the other n = 500 cells, then the n = 2000 cells.
# Within each size the design points keep this directory's order (0.75, 1.5, 1.0).
# <design-token> <FS_OR_TARGET> <n>.  ORSG_CELLS selects a subset by design-token_n<N>.
CELLS_ALL=("or150 1.5 500" "or075 0.75 500" "or100 1.0 500" "or075 0.75 2000" "or150 1.5 2000" "or100 1.0 2000")
# --- subset selection -------------------------------------------------------
CAMPAIGNS=()
for cs in "${CAMPAIGNS_ALL[@]}"; do
  t=${cs%% *}
  [ -z "${ORSG_CAMPAIGNS:-}" ] || [[ " ${ORSG_CAMPAIGNS} " == *" $t "* ]] || continue
  CAMPAIGNS+=("$cs")
done
CELLS=()
for cs in "${CELLS_ALL[@]}"; do
  read -r _d _t _n <<< "$cs"
  [ -z "${ORSG_CELLS:-}" ] || [[ " ${ORSG_CELLS} " == *" ${_d}_n${_n} "* ]] || continue
  CELLS+=("$cs")
done
[ ${#CAMPAIGNS[@]} -gt 0 ] || { echo "no campaign selected (ORSG_CAMPAIGNS='${ORSG_CAMPAIGNS:-}')"; exit 1; }
[ ${#CELLS[@]} -gt 0 ] || { echo "no cell selected (ORSG_CELLS='${ORSG_CELLS:-}')"; exit 1; }
HERE=$(cd "$(dirname "$0")" && pwd)
QMD=$(cd "$HERE/.." && pwd)
REPO=$(git -C "$QMD" rev-parse --show-toplevel)
RQ=${QMD#"$REPO"/}
LOGD=$QMD/logs_${TAG}
HB=$QMD/LOG_${TAG}_progress.txt
HALTF=$QMD/HALT_${TAG}.md
GREP=$QMD/REPORT_actg175_or_gate2_${DATE}.md
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
    echo "# HALT — ACTG175 binary/OR Stage 2"
    echo
    echo "- Halted at (UTC): $(utc)"
    echo "- Cell: $cell"
    echo "- Assertion: $what"
    echo "- HEAD at halt: $(git rev-parse --short HEAD)"
    echo "- Cumulative render wall: ${CUM} s (ceiling ${ORSG_CEILING} s)"
    echo "- Completed cells stay committed; nothing is rolled back."
    echo
    echo "## Observed values"
    echo
    echo '```'
    [ -n "$obs" ] && [ -f "$obs" ] && cat "$obs"
    echo '```'
  } > "$HALTF"
  heartbeat "$cell" halt "$what"
  commit "actg175 or $cell: HALT -- $what" "$RQ/HALT_${TAG}.md" "$RQ/LOG_${TAG}_progress.txt" || log "HALT commit failed"
  log "HALTED: $cell: $what"
  exit 1
}
tracked_clean () { git ls-files --error-unmatch -- "$@" >/dev/null 2>&1 && git diff --quiet HEAD -- "$@"; }
render () {  # $1 cell  $2 stem  $3 label  $4 mode  $5 start  $6 nsims  $7 method  $8 target  $9 n  $10 camp
  local cell=$1 stem=$2 label=$3 mode=$4 start=$5 nsims=$6 meth=$7 target=$8 n=$9 camp=${10}
  local out="${stem}_${label}.html" lg="$LOGD/${cell}_${label}.log" pk="$LOGD/${cell}_${label}.peak_mb"
  local t0=$(date +%s)
  ( cd "$QMD" && \
    FS_OR_METHOD=$meth FS_OR_FOCUS=effMaxSG FS_OR_NBHD=0.20 FS_OR_RULE=neighborhood \
    FS_OR_CI=field FS_OR_FIELD_SCALEC=selected FS_OR_FIELD_COMPLEMENT=TRUE FS_OR_RESELECTION=TRUE \
    FS_OR_CAMPAIGN=$camp FS_OR_TARGET=$target FS_OR_N=$n \
    FS_OR_MODE=$mode FS_OR_START=$start FS_OR_NSIMS=$nsims FS_OR_WORKERS=$ORSG_WORKERS \
    timeout "$ORSG_TIMEOUT" quarto render "$TEMPLATE" --output "$out" > "$lg" 2>&1 ) &
  local rpid=$!
  bash "$HERE/mem_sampler.sh" "$rpid" "$pk" &
  local spid=$!
  wait "$rpid"; local rc=$?
  wait "$spid" 2>/dev/null
  WALL=$(( $(date +%s) - t0 )); PEAK=$(cat "$pk" 2>/dev/null || echo NA); CUM=$(( CUM + WALL ))
  echo "WALL_SECONDS=$WALL RC=$rc PEAK_MB=$PEAK OUT=$out" >> "$lg"
  heartbeat "$cell" "$label" "wall_s=$WALL peak_mb=$PEAK rc=$rc cumulative_s=$CUM"
  [ $rc -eq 0 ] || { tail -40 "$lg" > "$LOGD/${cell}_${label}.halt_obs"; halt "$cell" "render $label failed (rc=$rc; 124 = timeout ${ORSG_TIMEOUT}s)" "$LOGD/${cell}_${label}.halt_obs"; }
}

# ---- preflight --------------------------------------------------------------
[ "$(git rev-parse --abbrev-ref HEAD)" = "$BRANCH" ] || { log "PREFLIGHT: not on $BRANCH"; exit 1; }
[ -f "$HALTF" ] && { log "PREFLIGHT: $RQ/HALT_${TAG}.md exists; the campaign is halted"; exit 1; }
git diff --quiet HEAD -- "$REPO/R" "$QMD/$TEMPLATE" || { log "PREFLIGHT: R/ or the template has uncommitted changes"; exit 1; }
pkill -f '[w]orkRSOCK'; sleep 5
NW=$(pgrep -fc '[w]orkRSOCK'); [ "$NW" = "0" ] || { log "PREFLIGHT: worker reap left $NW workRSOCK processes"; exit 1; }
BUILT=$(Rscript -e 'cat(packageDescription("forestsearch")$Built)')
log "PREFLIGHT OK: branch $BRANCH, HEAD $(git rev-parse --short HEAD), workers $ORSG_WORKERS, replicates/cell $NSIMS, timeout ${ORSG_TIMEOUT}s, ceiling ${ORSG_CEILING}s, Built '$BUILT'"
log "PREFLIGHT SELECTION: campaigns [$(for c in "${CAMPAIGNS[@]}"; do echo -n "${c%% *} "; done)] x cells [$(for c in "${CELLS[@]}"; do read -r d t n <<< "$c"; echo -n "${d}_n${n} "; done)]"
heartbeat campaign start "HEAD=$(git rev-parse --short HEAD) workers=$ORSG_WORKERS nsims=$NSIMS timeout_s=$ORSG_TIMEOUT ceiling_s=$ORSG_CEILING built=$BUILT"

# ---- campaigns x cells ------------------------------------------------------
for cspec in "${CAMPAIGNS[@]}"; do
  read -r CAMP METH MTOK <<< "$cspec"
  heartbeat "$CAMP" campaign_start "identifier=$METH"
  log "CAMPAIGN START $CAMP (identifier $METH)"
  for spec in "${CELLS[@]}"; do
    read -r DTOK TARGET N <<< "$spec"
    CELL=${CAMP}_${DTOK}_n${N}
    STEM=${MTOK}_effMaxSG_mr_field_${DTOK}_n${N}_nb20_${CAMP}
    RDIR=$QMD/mr_or_harm/${STEM}_d5000
    GREC=$LOGD/GATE2_${CELL}.txt
    if tracked_clean "$RQ/mr_or_harm/${STEM}_d5000/${STEM}_combined_1_${NSIMS}.rds"; then
      log "SKIP $CELL: combined bundle committed"; continue
    fi
    if [ "$CUM" -ge "$ORSG_CEILING" ]; then
      echo "cumulative $CUM s >= ceiling $ORSG_CEILING s before $CELL" > "$LOGD/${CELL}.halt_obs"
      halt "$CELL" "ceiling reached before the cell" "$LOGD/${CELL}.halt_obs"
    fi
    CELL_T0=$(date +%s)
    heartbeat "$CELL" start "stem=$STEM"
    # An interrupted cell restarts from its first batch: partial output is discarded
    # (the save guard refuses only git-tracked paths; these are untracked).
    rm -rf "$RDIR"
    rm -f "$QMD/${STEM}_batch_1_${NSIMS}.html" "$QMD/${STEM}_combine_1_${NSIMS}.html"
    log "START $CELL: identifier $METH target_or_h $TARGET n $N ($STEM), ${NSIMS} replicates"
    render "$CELL" "$STEM" "batch_1_${NSIMS}"   batch   1 "$NSIMS" "$METH" "$TARGET" "$N" "$CAMP"
    render "$CELL" "$STEM" "combine_1_${NSIMS}" combine 1 "$NSIMS" "$METH" "$TARGET" "$N" "$CAMP"
    mv "$QMD/${STEM}_batch_1_${NSIMS}.html" "$LOGD/" 2>/dev/null
    # Gate 2 (Stage 2's per-cell gate)
    (cd "$HERE" && ORSG_DIR=.. ORSG_WORKERS=$ORSG_WORKERS ORSG_NSIMS=$NSIMS Rscript gate2.R "$TARGET" "$N" "$CELL" "$CAMP") > "$GREC" 2>&1
    GRC=$?
    COUNTS=$(grep '^GATE_COUNTS' "$GREC")
    [ $GRC -eq 0 ] || halt "$CELL" "gate2.R failed: ${COUNTS:-no GATE_COUNTS line}" "$GREC"
    {
      [ -f "$GREP" ] || { echo "# REPORT — ACTG175 binary/OR Gate 2 (per cell)"; echo; echo "Task: \`dev/tasks/TASK_actg175_binary_campaign_${DATE}.md\` Stage 2. GRF's and DINA's candidate families are generated from fitted surfaces, so every coverage figure of the \`orgrf\` and \`ordina\` campaigns is coverage of the estimand **conditional on the proposed family**; FS's family is the prespecified cut grid. Comparisons across the three identifiers are descriptive. Runner: \`$RQ/scripts_or/run_or.sh\`; checker: \`$RQ/scripts_or/gate2.R\`."; echo; }
      echo "## $CELL — identifier $METH, target_or_h $TARGET, n $N"
      echo
      echo "- Stem: \`$STEM\`; HEAD before this cell's commit: $(git rev-parse --short HEAD); workers $ORSG_WORKERS; threads 1."
      echo "- Knobs: \`FS_OR_METHOD=$METH FS_OR_FOCUS=effMaxSG FS_OR_NBHD=0.20 FS_OR_RULE=neighborhood FS_OR_CI=field FS_OR_FIELD_SCALEC=selected FS_OR_CAMPAIGN=$CAMP FS_OR_TARGET=$TARGET FS_OR_N=$N FS_OR_WORKERS=$ORSG_WORKERS\`."
      echo "- Seeds: the study's pre-generated table indexed by global sim_id (seed_base 8316951); one batch, sim_id 1-${NSIMS}, then combine."
      echo "- Replicates: ${NSIMS} (TASK_binary_launch_v2_2026-09-18; the superseded 2 x 1,000 layout is not used)."
      echo "- Package: $BUILT."
      echo "- Cell wall: $(( $(date +%s) - CELL_T0 )) s; cumulative render wall $CUM s (ceiling $ORSG_CEILING s)."
      grep -h '^WALL_SECONDS=' "$LOGD/${CELL}_"*.log | sed 's/^/- Render: /'
      echo "- Gate counts: \`$COUNTS\`."
      echo
      echo '```'
      cat "$GREC"
      echo '```'
      echo
    } >> "$GREP"
    PATHS=("$RQ/mr_or_harm/${STEM}_d5000" "$RQ/${STEM}_combine_1_${NSIMS}.html" "$RQ/REPORT_actg175_or_gate2_${DATE}.md" "$RQ/LOG_${TAG}_progress.txt")
    heartbeat "$CELL" done "cell_wall_s=$(( $(date +%s) - CELL_T0 )) $COUNTS"
    commit "actg175 or $CELL: $METH, target_or_h $TARGET, n $N, ${NSIMS} replicates + combine; Gate 2 PASS ($COUNTS)" "${PATHS[@]}" || { log "CELL COMMIT FAILED: $CELL"; exit 1; }
    log "DONE $CELL: $COUNTS, commit $(git rev-parse --short HEAD)"
  done
  heartbeat "$CAMP" campaign_done "cumulative_s=$CUM"
  log "CAMPAIGN DONE $CAMP (cumulative render wall ${CUM} s)"
done
heartbeat campaign complete "cumulative_s=$CUM"
CAMPLIST=$(for cs in "${CAMPAIGNS[@]}"; do echo -n "${cs%% *} "; done)
commit "actg175 or: campaigns ${CAMPLIST% } complete (${NSIMS} replicates per cell; cumulative render wall ${CUM} s)" "$RQ/LOG_${TAG}_progress.txt" || log "final heartbeat commit failed"
log "CAMPAIGN COMPLETE"
