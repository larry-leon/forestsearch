#!/bin/bash
# Section 5 full re-run driver (TASK_section5_full_rerun_2026-09-24 section 3).
# Cells from mrs5sweep/cells.txt in file order (cheapest first: n 500, 1000, 1500).
# Per cell, as scripts_p12x20/campaign_p12x20.sh: batch 1 (1000) -> batch 1001 (1000) -> combine.
# Workers: 13 for every cell (n1500 smoke s5smk1500: 13.5 GB total R RSS at 13 workers, inside 24 GB).
# Overrun guard: 3x projection (committed wall x 1.32), record and continue.
# Memory guard: system free < 20% -> kill the render, end the run.
set -u
G=/Users/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020; D=$G/s5rerun
TAG=s5rerun; W=13; FAC=1.32; MINFREE=20
LOG=$D/logs/driver.log; REC=$D/logs/cells_record.tsv
log () { echo "[$(date '+%F %T')] $*" | tee -a $LOG; }
cd $G
free_pct () { memory_pressure -Q 2>/dev/null | grep -o '[0-9]*%' | tr -d %; }
killgrp () { kill -TERM -- -$1 2>/dev/null; sleep 5; kill -KILL -- -$1 2>/dev/null; pkill -f '[w]orkRSOCK' ; pkill -f '[p]arallelly.*launch' ; true; }
bash $D/memsample.sh $D/logs/${TAG}_mem.log 30 & MS=$!
trap 'kill $MS 2>/dev/null' EXIT
log "RUN START tag=$TAG workers=$W HEAD=$(git -C $G rev-parse --short HEAD)"
while read -r ID SRC HR N Z CW CWK STEM; do
  H3=$(awk -v h=$HR 'BEGIN{printf "%03d", int(h*100+0.5)}'); ZT=""; ZA=""; [ "$Z" != "-" ] && { ZT="_z1q60"; ZA=$Z; }
  PROJ=$(awk -v w=$CW -v f=$FAC 'BEGIN{printf "%d", w*f+0.5}'); LIM=$((3*PROJ))
  log "CELL START $ID ($SRC HR $HR n $N z1q $Z): projection ${PROJ}s, guard ${LIM}s"
  T0=$(date +%s); RC=0; HALT=""
  for PART in "batch 1 1000" "batch 1001 1000" "combine 1 1000"; do
    read -r MODE START NS <<< "$PART"
    S5_WORKERS=$W perl -e 'setpgrp(0,0); exec @ARGV' bash $D/runcell.sh $TAG $ID $HR $N $MODE $START $NS "$ZA" >> $LOG 2>&1 &
    P=$!
    while kill -0 $P 2>/dev/null; do
      sleep 10
      F=$(free_pct); EL=$(( $(date +%s) - T0 ))
      if [ -n "$F" ] && [ "$F" -lt $MINFREE ]; then HALT=mem; log "MEMORY STOP $ID at $MODE/$START: free ${F}% < ${MINFREE}%, elapsed ${EL}s"; killgrp $P; break; fi
      if [ $EL -gt $LIM ]; then HALT=overrun; log "OVERRUN $ID at $MODE/$START: ${EL}s > guard ${LIM}s"; killgrp $P; break; fi
    done
    wait $P; RC=$?
    [ -n "$HALT" ] && break
    if [ $RC -ne 0 ]; then log "RENDER FAIL $ID at $MODE/$START rc=$RC"; break; fi
  done
  WALL=$(( $(date +%s) - T0 ))
  if [ "$HALT" = mem ]; then log "RUN END (memory stop) during $ID, wall ${WALL}s"; exit 4; fi
  if [ -n "$HALT" ] || [ $RC -ne 0 ]; then log "CELL NOT COMPLETE $ID (${HALT:-rc=$RC}) wall ${WALL}s; continuing"; continue; fi
  CF=results/fs_effMaxSG_fb_mr_field_m1_h${H3}_knoise0_n${N}${ZT}_nb20_${TAG}_combined_1_2000.rds
  Rscript $D/record_cell.R $ID $CF $WALL $W $REC >> $LOG 2>&1 || log "RECORD FAIL $ID"
  log "CELL DONE $ID wall=${WALL}s projection=${PROJ}s"
done < <(grep -v '^#' $G/mrs5sweep/cells.txt)
log "RUN END (all cells attempted)"
