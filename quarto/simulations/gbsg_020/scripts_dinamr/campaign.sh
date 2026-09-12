#!/bin/zsh
# dinamr campaign driver.  usage: campaign.sh <block-letter> <cellspec-file>
# cellspec lines: "<z1q-or-empty> <n> <hr> <tag>"
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # override with DINAMR_SCRATCH
   # EXPORTED: render.sh reads $SP for its log directory and its header
   # states the caller exports what it needs.  Unexported, a clean-environment
   # run resolved $SP to empty and died on `mkdir -p /logs`.
cd $SP
KN=(FS_S7_METHOD=dina FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE
    FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE
    FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=dinamr FS_S7_WORKERS=12)
cell () {  # $1=z1q ("" = unset) $2=n $3=hr $4=tag
  local Z=$1 N=$2 H=$3 T=$4
  local T0=$(date +%s)
  for PART in "batch 1 1000" "batch 1001 1000" "combine 1 1000"; do
    set -- ${=PART}; local MODE=$1 START=$2 NS=$3
    local OUT="dinamr_${T}_${MODE}_${START}"
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
done < "$2"
echo "BLOCK $1 COMPLETE $(date)"
