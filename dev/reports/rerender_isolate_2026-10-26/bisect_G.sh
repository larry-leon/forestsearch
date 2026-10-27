#!/usr/bin/env bash
# git bisect run script for document G (cc_task_isolate_rerender_2026-10-26 §5).
# Run from the bisect worktree. good = legacy MR harm-confirmation HR 1.442 (exit 0), bad = 1.444 (exit 1), else 125.
REPO=/home/larryleon/Documents/GitHub/forestsearch
REC=$REPO/dev/reports/rerender_isolate_2026-10-26
WT=$(pwd)
sha=$(git rev-parse HEAD)
n=$(( $(ls -d ~/fs_bisect_G/step_* 2>/dev/null | wc -l) + 1 ))
step=~/fs_bisect_G/step_$n
mkdir -p "$step"
stem=analysis_gbsg_survival_effMaxSG
log() { echo "step $n $sha $(git log -1 --format='%ad %s' --date=short | cut -c1-90) :: $*" >> "$REC/bisect_G_steps.txt"; }
Rscript -e 'devtools::install(commandArgs(TRUE)[1], dependencies = FALSE, upgrade = FALSE, quick = TRUE)' "$WT" > "$step/install.log" 2>&1 \
  || { log "install failed -> skip (125)"; exit 125; }
cp "$REPO/quarto/applications/gbsg/$stem.qmd" "$step/"
mkdir -p "$step/_payloads"; export GH_DIR="$step/_payloads"
export RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64
( cd "$step" && timeout 20m quarto render "$stem.qmd" > render.log 2>&1 ) \
  || { log "render failed/capped -> skip (125)"; exit 125; }
val=$(cd "$REPO/dev/reports/rerender_2026-10-26" && Rscript -e 'source("extract_lib.R"); o <- extract_outputs(commandArgs(TRUE)[1]); x <- grep("^MR harm confirmation: HR = ", o, value = TRUE)[1]; cat(sub("^MR harm confirmation: HR = ([0-9.]+) .*$", "\\1", x))' "$step/$stem.html" 2>/dev/null)
if [ "$val" = "1.442" ]; then log "MR-HR=$val -> good (legacy)"; exit 0
elif [ "$val" = "1.444" ]; then log "MR-HR=$val -> bad (new)"; exit 1
else log "MR-HR=$val (neither 1.442 nor 1.444) -> skip (125)"; exit 125; fi
