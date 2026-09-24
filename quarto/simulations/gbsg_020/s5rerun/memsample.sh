#!/bin/bash
# usage: memsample.sh <logfile> [interval_s]  -- total R RSS and system free % every interval
while true; do
  echo "$(date +%T) $(ps -axo rss,comm | awk '/exec\/R$/ {s+=$1; n++} END {printf "R_procs=%d R_rss_GB=%.2f", n, s/1048576}') free_pct=$(memory_pressure -Q 2>/dev/null | grep -o '[0-9]*%' | tr -d %)"
  sleep ${2:-10}
done >> "$1" 2>&1
