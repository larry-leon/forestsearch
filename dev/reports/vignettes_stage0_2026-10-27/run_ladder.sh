#!/usr/bin/env bash
# Timing ladder driver, vignettes stage 0 (2026-10-27). Run from this directory.
# Each rung: timeout 1200 /usr/bin/time -v Rscript ladder.R <rung> > logs/<rung>.log 2>&1
# Stops launching rungs once cumulative elapsed exceeds 3 h or a rung is killed
# by timeout. S0 is GATE A: any non-zero exit at S0 stops the ladder.
# D3 / G3 run only if S3 completed under cap (exit 0).
set -u
cd "$(dirname "$0")"
mkdir -p logs
CSV=ladder.csv
[ -f "$CSV" ] || echo "rung,exit_code,status,wall_s,time_v_elapsed,max_rss_kb,cumulative_s" > "$CSV"
CAP_TOTAL=10800
cum=0
s3_ok=0
halted=0
for rung in S0 S1 S1r S2 S3 D1 G1 D3 G3; do
  if [ "$halted" -eq 1 ]; then echo "$rung,NA,not launched (an earlier rung was killed by timeout),NA,NA,NA,$cum" >> "$CSV"; continue; fi
  if [ "$rung" = D3 ] || [ "$rung" = G3 ]; then
    if [ "$s3_ok" -ne 1 ]; then echo "$rung,NA,skipped (S3 not completed under cap),NA,NA,NA,$cum" >> "$CSV"; continue; fi
  fi
  if awk -v c="$cum" -v t="$CAP_TOTAL" 'BEGIN{exit !(c > t)}'; then
    echo "$rung,NA,not launched (cumulative cap 3 h exceeded),NA,NA,NA,$cum" >> "$CSV"; continue
  fi
  t0=$(date +%s.%N)
  timeout 1200 /usr/bin/time -v Rscript ladder.R "$rung" > "logs/$rung.log" 2>&1
  rc=$?
  t1=$(date +%s.%N)
  wall=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')
  cum=$(awk -v c="$cum" -v w="$wall" 'BEGIN{printf "%.2f", c+w}')
  tv=$(grep -m1 'Elapsed (wall clock) time' "logs/$rung.log" | sed 's/.*: //')
  rss=$(grep -m1 'Maximum resident set size' "logs/$rung.log" | awk '{print $NF}')
  if [ "$rc" -eq 124 ]; then st="exceeded cap (killed by timeout 1200 s)"
  elif [ "$rc" -eq 0 ]; then st="completed"
  elif [ "$rc" -eq 3 ]; then st="GATE A FAIL"
  else st="failed"; fi
  echo "$rung,$rc,$st,$wall,${tv:-NA},${rss:-NA},$cum" >> "$CSV"
  echo "[$(date +%T)] $rung: $st, ${wall}s, max RSS ${rss:-NA} kB, cumulative ${cum}s"
  [ "$rung" = S3 ] && [ "$rc" -eq 0 ] && s3_ok=1
  if [ "$rung" = S0 ] && [ "$rc" -ne 0 ]; then echo "STOP: S0 (GATE A) did not pass; ladder halted."; exit 3; fi
  if [ "$rc" -eq 124 ]; then
    echo "rung $rung killed by timeout; no further rungs launched."
    halted=1
  fi
done
