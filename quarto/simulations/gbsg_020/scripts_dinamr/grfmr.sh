#!/bin/zsh
# grfmr campaign driver (TASK_grfmr_campaign_2026-09-11, PART A).
# usage: grfmr.sh <cellspec-file>
# cellspec lines: "<z1q-or-'-'> <n> <hr> <tag>"
#
# Identical in structure to campaign.sh; the knob set swaps the DINA engine for
# the GRF one (FS_S7_METHOD=grf) and the campaign tag for `grfmr`.  Every other
# knob -- focus, band, the four field knobs, the IJ residual, FB off, 12 workers
# -- is campaign.sh's, unchanged.  FS_S7_ER_JCUTS stays UNSET: it feeds the
# consistency-only method_args and is inert on GRF, exactly as it is on DINA.
#
# dmin.grf = 0.0, grf_selection = "frontier" and grf_select_statistic =
# "effect" are TEMPLATE LITERALS (template lines 503-506), not FS_S7_* knobs,
# so this driver cannot and does not set them; gate3.R resolves all three per
# batch.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # render.sh reads $SP
cd $SP
mkdir -p logs
KN=(FS_S7_METHOD=grf FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE
    FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE
    FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=grfmr FS_S7_WORKERS=12)
cell () {  # $1=z1q ("" = unset) $2=n $3=hr $4=tag
  local Z=$1 N=$2 H=$3 T=$4
  local T0=$(date +%s)
  for PART in "batch 1 1000" "batch 1001 1000" "combine 1 1000"; do
    set -- ${=PART}; local MODE=$1 START=$2 NS=$3
    local OUT="grfmr_${T}_${MODE}_${START}"
    if [[ -z "$Z" ]]; then
      env -u FS_S7_Z1Q -u FS_S7_ER_JCUTS $KN \
        FS_S7_HR=$H FS_S7_N=$N FS_S7_MODE=$MODE FS_S7_START=$START FS_S7_NSIMS=$NS \
        ./render.sh $OUT || { echo "CELL FAILED: $T at $MODE/$START"; return 1; }
    else
      env -u FS_S7_ER_JCUTS $KN FS_S7_Z1Q=$Z \
        FS_S7_HR=$H FS_S7_N=$N FS_S7_MODE=$MODE FS_S7_START=$START FS_S7_NSIMS=$NS \
        ./render.sh $OUT || { echo "CELL FAILED: $T at $MODE/$START"; return 1; }
    fi
  done
  local T1=$(date +%s)
  echo "CELL DONE: $T  wall=$((T1-T0))s"
}
while read -r Z N H T; do
  [[ -z "$T" ]] && continue
  [[ "$Z" == "-" ]] && Z=""
  echo "=== CELL START: $T (z1q='$Z' n=$N hr=$H) $(date) ==="
  cell "$Z" "$N" "$H" "$T" || echo "=== CELL ABORTED: $T, continuing ==="
done < "$1"
echo "GRFMR RUN COMPLETE $(date)"
