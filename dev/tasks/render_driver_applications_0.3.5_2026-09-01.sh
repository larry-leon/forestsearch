#!/usr/bin/env bash
# Render driver — applications re-render at 0.3.5 (task: cc_task_applications_render_0.3.5_2026-09-01.md)
# Transplants the 08-31 procedure verbatim: same four documents, same order (cheap first),
# same toolchain (RStudio-bundled quarto), same timing wrapper. Run from the repo root.
# Per-document wall-clock and exit codes land in $LOG; every fresh output must be -newer $SENTINEL.
set -u

QUARTO=/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto
SENTINEL=/tmp/fs_render_sentinel_20260901
LOG="${FS_RENDER_LOG:-/tmp/fs_render_0.3.5_2026-09-01.log}"

DOCS=(
  quarto/applications/actg175/template_actg175_continuous.qmd
  quarto/applications/actg175/analysis_actg175_continuous_compare_all.qmd
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
    echo "=== RENDER FAILED: $qmd (exit $rc) — stopping (G4)" >> "$LOG"
    exit "$rc"
  fi
done
echo "=== ALL RENDERS DONE $(date +%T)" >> "$LOG"
