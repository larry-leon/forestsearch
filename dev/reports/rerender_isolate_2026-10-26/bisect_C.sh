#!/usr/bin/env bash
# git bisect run script for document C (cc_task_isolate_rerender_2026-10-26 §5).
# Run from the bisect worktree. good = legacy beta-tilde 21.79 (exit 0), bad = 22.05 (exit 1), else 125.
REPO=/home/larryleon/Documents/GitHub/forestsearch
REC=$REPO/dev/reports/rerender_isolate_2026-10-26
WT=$(pwd)
sha=$(git rev-parse HEAD)
n=$(( $(ls -d ~/fs_bisect_C/step_* 2>/dev/null | wc -l) + 1 ))
step=~/fs_bisect_C/step_$n
mkdir -p "$step"
stem=analysis_actg175_continuous_intervals
log() { echo "step $n $sha $(git log -1 --format='%ad %s' --date=short | cut -c1-90) :: $*" >> "$REC/bisect_C_steps.txt"; }
Rscript -e 'devtools::install(commandArgs(TRUE)[1], dependencies = FALSE, upgrade = FALSE, quick = TRUE)' "$WT" > "$step/install.log" 2>&1 \
  || { log "install failed -> skip (125)"; exit 125; }
cp "$REPO/quarto/applications/actg175/$stem.qmd" "$step/"
export RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64
( cd "$step" && timeout 20m quarto render "$stem.qmd" > render.log 2>&1 ) \
  || { log "render failed/capped -> skip (125)"; exit 125; }
val=$(cd "$REPO/dev/reports/rerender_2026-10-26" && Rscript -e 'source("extract_lib.R"); o <- extract_outputs(commandArgs(TRUE)[1]); x <- grep("^MR \\(IJ two-term\\) \\| .*β̃\\)", o, value = TRUE)[1]; cat(sub("^MR \\(IJ two-term\\) \\| ([-0-9.]+) \\(β̃\\).*$", "\\1", x))' "$step/$stem.html" 2>/dev/null)
if [ "$val" = "21.79" ]; then log "beta-tilde=$val -> good (legacy)"; exit 0
elif [ "$val" = "22.05" ]; then log "beta-tilde=$val -> bad (new)"; exit 1
else log "beta-tilde=$val (neither 21.79 nor 22.05) -> skip (125)"; exit 125; fi
