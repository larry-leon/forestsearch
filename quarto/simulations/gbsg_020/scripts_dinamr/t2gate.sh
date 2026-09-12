#!/bin/zsh
# Gate T2 driver (TASK_grfmr_campaign_2026-09-11, Part T2).
#
# Five replicates at the STANDING IDENTITY CELL with FS_S7_METHOD UNSET
# (so the consistency path, the one the add-only recorder change must leave
# untouched): effMaxSG, eps 0.20, HR 1.50, n 500, FS_S7_Z1Q=0.60,
# seeds 8316951 + sim_id, sim_id 1-5.
#
# usage: t2gate.sh <campaign-tag>      e.g.  t2gate.sh t2pre   /  t2gate.sh t2post
#
# Run it once BEFORE the template edit and once AFTER, with the same argument
# style but different tags, then compare with t2gate.R.  The comparison is
# pre-change against post-change ON THIS MACHINE, which is what the task
# specifies; it is not a comparison against a committed bundle.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"   # render.sh reads $SP
cd $SP
TAG="${1:?usage: t2gate.sh <campaign-tag>}"
KN=(FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE
    FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE
    FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=$TAG FS_S7_WORKERS=12
    FS_S7_Z1Q=0.60 FS_S7_HR=1.50 FS_S7_N=500 FS_S7_MODE=batch
    FS_S7_START=1 FS_S7_NSIMS=5)
env -u FS_S7_METHOD -u FS_S7_ER_JCUTS $KN ./render.sh t2gate_$TAG
echo "T2 GATE RENDER DONE: $TAG $(date)"
