#!/bin/bash
# Stage 2 driver for ONE cell: batch 1 (sim 1-1000) -> batch 2 (1001-2000) -> combine.
# Usage: run_cell.sh MD N FBJOIN THRESH_SEC SAMPLE_MEM   (MD in {40,120,null}; FBJOIN in {none,join})
MD=$1; N=$2; FB=$3; THRESH=$4; SAMPLE=$5; TAG=mdf1
SP=/private/tmp/claude-501/-Users-larryleon-Downloads/476ca7fd-7bca-4851-bd29-048be32a2d48/scratchpad
C=$SP/campaign; CEIL=14400; TMO=5400
cd ~/Documents/GitHub/forestsearch/quarto/simulations/actg175/continuous || exit 2
if [ "$MD" = "null" ]; then STEM=fs_maxeffCons_mr_field_mdnull_knoise0_n${N}_${TAG}; else STEM=fs_maxeffCons_mr_field_md${MD}_knoise0_n${N}_${TAG}; fi
CELL=md${MD}_n${N}
cum() { awk '{s+=$1} END {printf "%d", s}' $C/wall.txt 2>/dev/null; }
log() { echo "$(date '+%F %T') [$CELL] $*" | tee -a $C/status.txt; }
render() { # label mode start nsims fb output
  local label=$1 mode=$2 start=$3 nsims=$4 fb=$5 out=$6
  local cumnow=$(cum)
  if [ "${cumnow:-0}" -ge $CEIL ]; then log "CEILING reached ($cumnow s) before $label -- deferred"; echo "DEFERRED $CELL $label" >> $C/deferred.txt; return 9; fi
  log "START $label (mode=$mode start=$start nsims=$nsims fb=$fb)"
  local t0=$(date +%s)
  FS_MD_MD=$MD FS_MD_N=$N FS_MD_MODE=$mode FS_MD_START=$start FS_MD_NSIMS=$nsims FS_MD_FB=$fb \
  FS_MD_CAMPAIGN=$TAG FS_MD_WORKERS=13 FS_MD_CI=field \
    bash /private/tmp/claude-501/-Users-larryleon-Downloads/476ca7fd-7bca-4851-bd29-048be32a2d48/scratchpad/tmo.sh $TMO quarto render sim_fs_maxeffCons_mr_field_md_template.qmd --output $out > $C/${CELL}_${label}.log 2>&1
  local rc=$?; local w=$(( $(date +%s) - t0 ))
  echo "$w $CELL $label rc=$rc" >> $C/wall.txt
  log "END $label rc=$rc wall=${w}s (cumulative $(cum) s)"
  [ $rc -ne 0 ] && { log "FAILED $label (rc=$rc; 124 = hard timeout)"; return $rc; }
  if [ "$mode" = "batch" ] && [ $w -gt $THRESH ]; then log "OVER 1.5x PROJECTION: $label ${w}s > ${THRESH}s -- campaign stops after this cell's combine"; echo "OVER $CELL $label $w" >> $C/over.txt; fi
  return 0
}
if [ "$SAMPLE" = "1" ]; then bash $SP/mem_sampler.sh $$ $C/peak_mem_mb.txt & fi
render batch1 batch 1 1000 $FB ${STEM}_batch_1_1000.html || exit 1
render batch2 batch 1001 1000 none ${STEM}_batch_1001_2000.html || exit 1
render combine combine 1 1000 none ${STEM}_combine_1_2000.html || exit 1
log "CELL DONE"
