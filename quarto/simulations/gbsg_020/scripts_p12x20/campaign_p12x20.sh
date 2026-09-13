#!/bin/bash
# p12x20 campaign driver (TASK_p12x20_partA_2026-09-12_v2 §5).
# usage: campaign_p12x20.sh <block-label> <cellspec-file>
# cellspec lines: "<z1q-or-'-'> <n> <hr> <tag>"
#
# Transplant of scripts_dinamr/campaign.sh, identical in structure (per cell:
# batch 1 -> batch 1001 -> combine, via the committed scripts_dinamr/render.sh).
# Differences, and nothing else:
#   1. bash, not zsh.  zsh is not installed on pop-os, so campaign.sh cannot
#      execute on this host: ${=PART} -> read, $KN -> "${KN[@]}",
#      ./render.sh -> bash $RENDER.
#   2. The knob set is Stage 1 §1b's: the cert20 invocation with FS_S7_Z1Q
#      unset (12.4%) and campaign tag p12x20.  FS_S7_METHOD, FS_S7_MR,
#      FS_S7_FIELD_RECOV and FS_S7_ER_JCUTS stay UNSET (consistency / TRUE /
#      FALSE / 10 -- the values cert20 ran with; the first three knobs postdate
#      cert20).  Every inherited FS_S7_* is unset before the knobs are applied.
#   3. FS_S7_WORKERS comes from $P12X20_WORKERS (named in run_p12x20.sh).
#   4. Render basenames p12x20_<tag>_<mode>_<start>, as campaign.sh's dinamr_*.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # render.sh reads $SP
RENDER="$(cd "$(dirname "$0")/../scripts_dinamr" && pwd)/render.sh"
: "${P12X20_WORKERS:?P12X20_WORKERS must be set}"
cd "$SP"
mkdir -p logs
KN=(FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE
    FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE
    FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=p12x20
    FS_S7_RETURN_RESEL=TRUE FS_S7_WORKERS=$P12X20_WORKERS)
UNSET=()
for v in $(env | grep -o '^FS_S7_[A-Z0-9_]*'); do UNSET+=(-u "$v"); done
cell () {  # $1=z1q ("" = unset) $2=n $3=hr $4=tag
  local Z=$1 N=$2 H=$3 T=$4
  local T0=$(date +%s)
  local PART MODE START NS OUT
  for PART in "batch 1 1000" "batch 1001 1000" "combine 1 1000"; do
    read -r MODE START NS <<< "$PART"
    OUT="p12x20_${T}_${MODE}_${START}"
    if [[ -z "$Z" ]]; then
      env "${UNSET[@]}" "${KN[@]}" \
        FS_S7_HR=$H FS_S7_N=$N FS_S7_MODE=$MODE FS_S7_START=$START FS_S7_NSIMS=$NS \
        bash "$RENDER" "$OUT" || { echo "CELL FAILED: $T at $MODE/$START"; return 1; }
    else
      env "${UNSET[@]}" "${KN[@]}" FS_S7_Z1Q=$Z \
        FS_S7_HR=$H FS_S7_N=$N FS_S7_MODE=$MODE FS_S7_START=$START FS_S7_NSIMS=$NS \
        bash "$RENDER" "$OUT" || { echo "CELL FAILED: $T at $MODE/$START"; return 1; }
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
