#!/bin/bash
# Sample every 5 s the summed RSS (MB) of all R / Rscript / quarto / deno processes while $1 (a pid) lives;
# write the running maximum to $2.
# Transplant of scripts_mddina/mem_sampler.sh (TASK_actg175_binary_campaign_2026-09-17 §1.4), unchanged; itself scripts_mdsgnb20/mem_sampler.sh; that one was scripts_mdf1's with
# the Linux `ps -eo rss,comm` in place of macOS `ps -Ao rss,comm`.
PID=$1; OUT=$2; MAX=0
while kill -0 "$PID" 2>/dev/null; do
  S=$(ps -eo rss,comm | awk '$2 ~ /(^|\/)(R|Rscript|quarto|deno)$/ {s+=$1} END {printf "%d", s/1024}')
  [ "${S:-0}" -gt "$MAX" ] && MAX=$S
  echo "$MAX" > "$OUT"
  sleep 5
done
echo "$MAX" > "$OUT"
