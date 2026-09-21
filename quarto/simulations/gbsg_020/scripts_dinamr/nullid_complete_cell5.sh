#!/bin/zsh
# Completion run for cell null0657_n1500 (TASK_null_gbsg_identification_2026-09-21).
#
# That cell was ABANDONED by nullid.sh after a FALSE Gate A failure on its
# consistency run: the gate asserted a 95% coverage threshold on max_g T_g,
# which is not an invariant -- maxT is absent when no candidate cleared the
# effect floor, and that happens on 330 of 2000 replicates at this cell.  See
# HALT_nullid_null0657_n1500_2026-09-21.txt (withdrawn) and the corrected
# nullid_gateA.R.
#
# The consistency bundle is CORRECT and is NOT re-run.  This runs only the two
# identifiers the driver never reached, with byte-identical knobs to nullid.sh,
# then gates the cell.
export SP="${DINAMR_SCRATCH:-$(cd "$(dirname "$0")" && pwd)}"
cd $SP
mkdir -p logs
T=null0657_n1500; H=0.657; N=1500
PIN=(FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE
     FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_WORKERS=12)
BASE=($PIN FS_S7_MR=FALSE FS_S7_DGM=null FS_S7_CAMPAIGN=nullid FS_S7_MODE=batch
      FS_S7_START=1 FS_S7_NSIMS=2000 FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20
      FS_S7_HR=$H FS_S7_N=$N)
T0=$(date +%s)
echo "=== NULLID CELL-5 COMPLETION START $(date) ; host $(hostname) ; git $(git rev-parse --short HEAD) ==="
for E in dina grf; do
  OUT="nullid_${T}_${E}_effMaxSG"
  echo "=== RUN START: $T $E $(date) -> $OUT ==="
  env -u FS_S7_ER_JCUTS -u FS_S7_Z1Q $BASE FS_S7_METHOD=$E ./render.sh $OUT < /dev/null \
    || { echo "RUN FAILED: $T $E"; exit 1; }
  Rscript nullid_gateA.R $OUT $E $H $N $T > logs/nullid_gateA_${T}_${E}.log 2>&1 < /dev/null
  GA=$?
  grep -E '^GATE A ' logs/nullid_gateA_${T}_${E}.log
  if (( GA != 0 )); then grep -E '^\[FAIL\]' logs/nullid_gateA_${T}_${E}.log; echo "=== GATE A FAILED: $OUT ==="; exit 2; fi
done
echo "CELL DONE: $T  wall=$(( $(date +%s) - T0 ))s (completion run; the consistency render is from the main campaign)"
Rscript nullid_gateC.R $T $H $N > logs/nullid_gateC_${T}.txt 2>&1 < /dev/null
GC=$?
grep -E '^GATE C ' logs/nullid_gateC_${T}.txt
if (( GC != 0 )); then grep -E '^\[FAIL\]' logs/nullid_gateC_${T}.txt; exit 3; fi
echo "NULLID CELL-5 COMPLETION DONE $(date) total=$(( $(date +%s) - T0 ))s"
