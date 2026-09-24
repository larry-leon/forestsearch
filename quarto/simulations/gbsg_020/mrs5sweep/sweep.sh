#!/bin/bash
# Section 5 sweep driver (TASK_mr_alignment_section5_sweep_2026-09-23 section 2).
# Cheapest first; budget rule; 3x overrun timeout; stop on Gate 1 (declaration) or Gate 2 (build) failure.
set -u
G=/home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020; D=$G/mrs5sweep
START=${SWEEP_START:?}; BUDGET=9000; CAL=1.44; NS=200
LOG=$D/logs/sweep_driver.log
log () { echo "[$(date '+%F %T')] $*" | tee -a $LOG; }
cd $G
while read -r ID SRC HR N Z CW CWK STEM; do
  H3=$(awk -v h=$HR 'BEGIN{printf "%03d", int(h*100+0.5)}'); ZT=""; [ "$Z" != "-" ] && ZT="_z1q60"
  AFTER=results/fs_effMaxSG_fb_mr_field_m1_h${H3}_knoise0_n${N}${ZT}_nb20_mrs5sweep_res_1_${NS}.rds
  if [ "$ID" = A7 ]; then
    log "CELL $ID: reused from the probe (mrs5probe, 188 s)"
    Rscript $D/pair.R $ID results/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n500_nb20_mrs5probe_res_1_200.rds $STEM $D/pair/$ID.rds >> $LOG 2>&1 || { log "GATE1 FAIL $ID"; exit 2; }
    continue
  fi
  PROJ=$(awk -v w=$CW -v k=$CWK -v c=$CAL 'BEGIN{printf "%d", w*0.1*(k/64)*c + 0.5}')
  EL=$(( $(date +%s) - START ))
  if [ $(( EL + PROJ )) -gt $BUDGET ]; then log "BUDGET STOP before $ID: elapsed ${EL}s + projection ${PROJ}s > ${BUDGET}s"; echo "BUDGET_STOP $ID" >> $LOG; exit 0; fi
  log "CELL START $ID ($SRC HR $HR n $N z1q $Z): projection ${PROJ}s, elapsed ${EL}s"
  T0=$(date +%s)
  ZA=""; [ "$Z" != "-" ] && ZA=$Z
  timeout $(( 3*PROJ )) bash $D/runcell.sh mrs5sweep $ID $HR $N $NS "$ZA" >> $LOG 2>&1
  RC=$?; W=$(( $(date +%s) - T0 ))
  if [ $RC -eq 3 ]; then log "GATE2 FAIL $ID: worker build check"; exit 3; fi
  if [ $RC -eq 124 ]; then log "OVERRUN $ID: stopped at ${W}s (3x projection ${PROJ}s)"; pkill -f '[w]orkRSOCK'; continue; fi
  if [ $RC -ne 0 ]; then log "RENDER FAIL $ID rc=$RC wall ${W}s"; continue; fi
  log "CELL DONE $ID wall=${W}s projection=${PROJ}s"
  Rscript $D/pair.R $ID $AFTER $STEM $D/pair/$ID.rds >> $LOG 2>&1 || { log "GATE1 FAIL $ID: declaration differs"; exit 2; }
done < <(grep -v '^#' $D/cells.txt)
log "SWEEP END"
