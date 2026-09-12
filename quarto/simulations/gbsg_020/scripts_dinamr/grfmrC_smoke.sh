#!/bin/zsh
# grfmr completion (TASK_grfmr_completion_2026-09-12), STAGE 1 smoke:
# 5 replicates at (12.4%, HR 1.00, n 500) on the GRF path, stop-on-failure.
# Knob set is grfmr.sh's exactly, with the campaign tag `grfmrsmk` so the
# bundle cannot collide with a campaign cell.  FS_S7_Z1Q unset = 12.4%.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"
cd $SP
KN=(FS_S7_METHOD=grf FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE
    FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE
    FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=grfmrsmk FS_S7_WORKERS=12
    FS_S7_HR=1.00 FS_S7_N=500 FS_S7_MODE=batch FS_S7_START=1 FS_S7_NSIMS=5)
env -u FS_S7_Z1Q -u FS_S7_ER_JCUTS $KN ./render.sh grfmrsmk_C124_h100_n500
echo "GRFMR C SMOKE DONE $(date)"
