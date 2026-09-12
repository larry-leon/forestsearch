#!/bin/zsh
# Post-T2 GRF smoke: 5 replicates on the GRF path, to confirm admitted_n is
# actually POPULATED (Gate T2 only proves it is NA off the GRF path) and to
# read Gate 3's three resolved alignment values off the render log before any
# 2,000-replicate cell is launched.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"
cd $SP
KN=(FS_S7_METHOD=grf FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE
    FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE
    FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=t2grf FS_S7_WORKERS=12
    FS_S7_HR=1.50 FS_S7_N=500 FS_S7_MODE=batch FS_S7_START=1 FS_S7_NSIMS=5)
env -u FS_S7_Z1Q -u FS_S7_ER_JCUTS $KN ./render.sh t2grf_smoke
echo "T2 GRF SMOKE DONE $(date)"
