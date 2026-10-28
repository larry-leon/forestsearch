#!/usr/bin/env bash
# Timing ladder driver, post-selection vignette stage 0 (2026-10-28). Run from this directory.
# Each rung: timeout 1200 /usr/bin/time -v Rscript ladder.R <rung> > logs/<rung>.log 2>&1
# M-FS is GATE A: any non-zero exit at M-FS stops the ladder. A B rung runs only if
# its M rung completed. No rung is launched once cumulative elapsed exceeds 2 h.
set -u
cd "$(dirname "$0")"
mkdir -p logs
CSV=ladder.csv
[ -f "$CSV" ] || echo "rung,exit_code,status,wall_s,time_v_elapsed,max_rss_kb,cumulative_s" > "$CSV"
CAP_TOTAL=7200
cum=0
declare -A ok
for rung in M-FS M-DINA M-GRF M-FS-ij B-FS B-DINA B-GRF; do
  case "$rung" in B-*) m="M-${rung#B-}"; if [ "${ok[$m]:-0}" -ne 1 ]; then echo "$rung,NA,skipped ($m not completed),NA,NA,NA,$cum" >> "$CSV"; continue; fi;; esac
  if awk -v c="$cum" -v t="$CAP_TOTAL" 'BEGIN{exit !(c > t)}'; then
    echo "$rung,NA,not launched (cumulative cap 2 h exceeded),NA,NA,NA,$cum" >> "$CSV"; continue
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
  elif [ "$rc" -eq 0 ]; then st="completed"; ok[$rung]=1
  elif [ "$rc" -eq 3 ]; then st="GATE A FAIL"
  elif [ "$rc" -eq 4 ]; then st="bootstrap call errored"
  else st="failed"; fi
  echo "$rung,$rc,$st,$wall,${tv:-NA},${rss:-NA},$cum" >> "$CSV"
  echo "[$(date +%T)] $rung: $st, ${wall}s, max RSS ${rss:-NA} kB, cumulative ${cum}s"
  if [ "$rung" = M-FS ] && [ "$rc" -ne 0 ]; then echo "STOP: M-FS (GATE A) did not pass; ladder halted."; exit 3; fi
done
