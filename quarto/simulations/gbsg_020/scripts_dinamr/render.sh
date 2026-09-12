#!/bin/zsh
# dinamr campaign render driver.
# usage: render.sh <output-basename> [extra env already exported by caller]
set -e
cd "${DINAMR_QMD_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"   # the dir holding the template
OUT="$1"
export VECLIB_MAXIMUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
mkdir -p "$SP/logs"; LOG="$SP/logs/${OUT}.log"
T0=$(date +%s)
quarto render sim_fs_maxeffCons_fb_mr_field_m1_template.qmd --output "${OUT}.html" > "$LOG" 2>&1
RC=$?
T1=$(date +%s)
echo "WALL_SECONDS=$((T1-T0)) RC=$RC OUT=${OUT}" | tee -a "$LOG"
exit $RC
