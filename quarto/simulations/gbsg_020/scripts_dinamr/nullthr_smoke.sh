#!/bin/bash
# nullc125 smoke / inertness renders (TASK_null_gbsg_thresholds_2026-09-21_v2,
# Steps 1.5 and 3.1).  One cell, 20 replicates, the three identifiers, on
# nullid's pinned knob set with FS_S7_QUICKRUN=TRUE and a smoke campaign tag,
# so nothing is written to a committed stem.
# usage: nullthr_smoke.sh <campaign-tag> <hr> <n> [extra KEY=VALUE ...]
#   the extra pairs are appended to the environment (e.g. FS_S7_C1=1.25 FS_S7_C2=1.00).
# Every inherited FS_S7_* is unset first, as campaign_p12x20.sh does.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # render.sh reads $SP
RENDER="$SP/render.sh"
export PATH=/usr/lib/rstudio/resources/app/bin/quarto/bin:$PATH
TAG=$1; H=$2; N=$3; shift 3
EXTRA=("$@")
WORKERS=${NULLTHR_SMOKE_WORKERS:-20}
cd "$SP" || exit 1
mkdir -p logs
UNSET=()
for v in $(env | grep -o '^FS_S7_[A-Z0-9_]*'); do UNSET+=(-u "$v"); done
PIN=(FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE
     FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_WORKERS=$WORKERS)
BASE=("${PIN[@]}" FS_S7_MR=FALSE FS_S7_DGM=null FS_S7_CAMPAIGN=$TAG FS_S7_MODE=batch
      FS_S7_START=1 FS_S7_NSIMS=20 FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_QUICKRUN=TRUE
      FS_S7_HR=$H FS_S7_N=$N)
RC=0
for E in consistency dina grf; do
  KN=("${BASE[@]}")
  [[ $E != consistency ]] && KN+=(FS_S7_METHOD=$E)
  KN+=("${EXTRA[@]}")
  OUT="${TAG}_${E}"
  echo "=== SMOKE $TAG $E $(date) ; env: ${KN[*]}"
  env "${UNSET[@]}" "${KN[@]}" bash "$RENDER" "$OUT" < /dev/null || { echo "SMOKE RENDER FAILED: $OUT"; RC=1; }
done
exit $RC
