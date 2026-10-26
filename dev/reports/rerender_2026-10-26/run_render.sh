#!/usr/bin/env bash
# usage: run_render.sh <id> <cap> <workdir> <cmd...>   (logs to render_<id>.log; times to timing.txt)
R=/home/larryleon/Documents/GitHub/forestsearch/dev/reports/rerender_2026-10-26
id=$1; cap=$2; wd=$3; shift 3
export RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64
cd "$wd"
t0=$(date +%s)
timeout "$cap" "$@" > "$R/render_$id.log" 2>&1
rc=$?
t1=$(date +%s)
echo "$id start=$t0 end=$t1 wall_s=$((t1-t0)) exit=$rc cap=$cap" >> "$R/timing.txt"
