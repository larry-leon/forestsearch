#!/bin/bash
# Sample every 5 s the summed RSS (MB) of all R / Rscript / quarto processes while $1 (a pid) lives;
# write the running maximum to $2.
PID=$1; OUT=$2; MAX=0
while kill -0 "$PID" 2>/dev/null; do
  S=$(ps -Ao rss,comm | awk '$2 ~ /(^|\/)(R|Rscript|quarto|deno)$/ {s+=$1} END {printf "%d", s/1024}')
  [ "${S:-0}" -gt "$MAX" ] && MAX=$S
  echo "$MAX" > "$OUT"
  sleep 5
done
echo "$MAX" > "$OUT"
