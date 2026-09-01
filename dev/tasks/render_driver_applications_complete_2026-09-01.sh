#!/usr/bin/env bash
# Render driver — complete applications at 0.3.5 (task: cc_task_applications_complete_render_2026-09-01.md)
# Transplant of dev/tasks/render_driver_applications_0.3.5_2026-09-01.sh; changed lines only:
# document set (the two FLAGGED gbsg documents — unflagged documents are not re-rendered),
# sentinel name, log default. Run from the repo root.
# Per-document wall-clock and exit codes land in $LOG; every fresh output must be -newer $SENTINEL.
set -u

QUARTO=/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto
SENTINEL=/tmp/fs_render_sentinel_20260901_complete
LOG="${FS_RENDER_LOG:-/tmp/fs_render_complete_2026-09-01.log}"

DOCS=(
  quarto/applications/gbsg/analysis_gbsg_survival_frozen_family.qmd
  quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd
)

touch "$SENTINEL"
: > "$LOG"
echo "sentinel: $(stat -c %y "$SENTINEL")" >> "$LOG"
"$QUARTO" --version >> "$LOG" 2>&1

for qmd in "${DOCS[@]}"; do
  echo "=== RENDER $qmd  start $(date +%T)" >> "$LOG"
  /usr/bin/time -f "elapsed %E  maxRSS %MkB" "$QUARTO" render "$qmd" --to html >> "$LOG" 2>&1
  rc=$?
  echo "=== exit $rc  end $(date +%T)" >> "$LOG"
  if [ "$rc" -ne 0 ]; then
    echo "=== RENDER FAILED: $qmd (exit $rc) — stopping (per-render gate)" >> "$LOG"
    exit "$rc"
  fi
done
echo "=== ALL RENDERS DONE $(date +%T)" >> "$LOG"
