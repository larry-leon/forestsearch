#!/bin/zsh
# Gate T3 driver (TASK_partB_enabling_2026-09-12, Part 1), stop-on-failure.
#
# Five replicates at the STANDING IDENTITY CELL -- the same cell t2gate.sh uses:
# effMaxSG, eps 0.20, HR 1.50, n 500, FS_S7_Z1Q=0.60, seeds 8316951 + sim_id,
# sim_id 1-5.  Engine: consistency (FS_S7_METHOD unset) unless METHOD is given.
#
# usage: t3gate.sh <campaign-tag> [MR] [METHOD]
#   t3gate.sh t3pre                  BEFORE the template edit, FS_S7_MR unset
#   t3gate.sh t3post                 AFTER the edit,  FS_S7_MR unset
#   t3gate.sh t3mroff FALSE          AFTER the edit,  FS_S7_MR=FALSE
#   t3gate.sh t3dina  "" dina        supplementary: DINA, MR unset
#   t3gate.sh t3dinaoff FALSE dina   supplementary: DINA, MR off  (likewise grf)
# then compare with t3gate.R.  Pre against post is ON THIS MACHINE.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # render.sh reads $SP
cd $SP
TAG="${1:?usage: t3gate.sh <campaign-tag> [MR] [METHOD]}"
MR="${2:-}"
METHOD="${3:-}"
KN=(FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE
    FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE
    FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=$TAG FS_S7_WORKERS=12
    FS_S7_Z1Q=0.60 FS_S7_HR=1.50 FS_S7_N=500 FS_S7_MODE=batch
    FS_S7_START=1 FS_S7_NSIMS=5)
[[ -n "$MR" ]] && KN+=(FS_S7_MR=$MR)
[[ -n "$METHOD" ]] && KN+=(FS_S7_METHOD=$METHOD)
UNSET=(-u FS_S7_ER_JCUTS)
[[ -z "$MR" ]] && UNSET+=(-u FS_S7_MR)
[[ -z "$METHOD" ]] && UNSET+=(-u FS_S7_METHOD)
env $UNSET $KN ./render.sh t3gate_$TAG
echo "T3 GATE RENDER DONE: $TAG MR='${MR:-unset}' METHOD='${METHOD:-unset}' $(date)"
