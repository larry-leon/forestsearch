#!/bin/bash
# declcal_c0approx driver (TASK_declcal_c0_approx_2026-09-22): 60 field
# captures -- reps 1-20 of B1..B6 at B = 2000 with declaration_c0 -- then the
# read-only approximation (declcal_c0approx_findings.R).
# The per-cell environment is declcal.sh's (BASE + per-cell knobs from
# declcal_inull.cells) with FS_S7_NSIMS=20, DECLCAL_BCAL=2000, and
# FS_S7_WORKERS=10: the six cells run concurrently (60 workers in all).
# usage: declcal_c0approx.sh <celldir>   (per-cell payloads are written there;
#        the combined 60-row payload is results/declcal_c0approx_res.rds)
export SP="$(cd "$(dirname "$0")" && pwd)"
QD="$(cd "$SP/.." && pwd)"
CELLDIR=${1:?usage: declcal_c0approx.sh <celldir>}
mkdir -p "$CELLDIR" "$SP/logs"
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

LIVE=$(pgrep -af 'quarto render|Rscript|exec/R ' | grep -v -e "$$" -e pgrep || true)
if [[ -n "$LIVE" ]]; then
  echo "=== ANOTHER R / QUARTO PROCESS IS RUNNING -- not competing for cores; stopping ==="
  echo "$LIVE"; exit 5
fi
UNSET=()
for v in $(env | grep -o '^FS_S7_[A-Z0-9_]*'); do UNSET+=(-u "$v"); done
BASE=(FS_S7_WORKERS=10 FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_C1=1.0 FS_S7_C2=1.0
      FS_S7_START=1 FS_S7_NSIMS=20)

T0=$(date +%s)
echo "=== DECLCAL_C0APPROX START $(date) ; git $(git -C "$QD" rev-parse --short HEAD) ; forestsearch $(Rscript -e 'cat(as.character(packageVersion("forestsearch")), packageDescription("forestsearch")$Built)') ==="
PIDS=()
while read -r CELL DGM HR N LAB; do
  [[ -z "$LAB" || "$CELL" == \#* ]] && continue
  KN=("${BASE[@]}" FS_S7_DGM=$DGM FS_S7_HR=$HR FS_S7_N=$N FS_S7_CAMPAIGN=declcal_c0approx
      DECLCAL_CELL=$CELL DECLCAL_DGM_LABEL=$LAB DECLCAL_BCAL=2000 DECLCAL_MODE=stage2
      DECLCAL_CAP_S=0 DECLCAL_OUT=$CELLDIR/declcal_c0approx_${CELL}_res_1_20.rds)
  echo "  cell $CELL env: ${KN[*]}"
  ( cd "$QD" && env "${UNSET[@]}" "${KN[@]}" Rscript "$SP/declcal_c0approx_run.R" ) \
    > "$CELLDIR/declcal_c0approx_${CELL}.log" 2>&1 < /dev/null &
  PIDS+=("$CELL:$!")
done < "$SP/declcal_inull.cells"
RCALL=0
for p in "${PIDS[@]}"; do
  wait "${p#*:}"; RC=$?
  echo "CELL ${p%%:*} rc=$RC"; grep -E '^(CELL|FIDELITY|RATES)' "$CELLDIR/declcal_c0approx_${p%%:*}.log"
  (( RC != 0 )) && RCALL=1
done
echo "=== captures done: wall $(( $(date +%s) - T0 ))s ; any failure $RCALL ==="
(( RCALL )) && exit 1
( cd "$QD" && C0A_CELLDIR="$CELLDIR" Rscript "$SP/declcal_c0approx_findings.R" )
