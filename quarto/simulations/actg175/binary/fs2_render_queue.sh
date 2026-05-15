#!/usr/bin/env bash
# =============================================================================
# One-time setuop
# chmod +x quarto/simulations/actg175/fs2_render_queue.sh
# fs2_render_queue.sh
#
# Serially render all fs2 ACTG simulation qmds. Continues past failures
# so a single bad render does not kill an overnight queue.
#
# Run from the forestsearch package root (the directory containing
# 'quarto/').
#
# Usage:
#   ./quarto/simulations/actg175/binary/fs2_render_queue.sh
# =============================================================================

set -uo pipefail
shopt -s nullglob

QMDS=(quarto/simulations/actg175/binary/actg175_binary_m1_harm_*_fs2.qmd)

if [ "${#QMDS[@]}" -eq 0 ]; then
  echo "No fs2 qmds found.  Check the glob pattern and working directory."
  exit 1
fi

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="quarto/simulations/fs2_render_${TIMESTAMP}.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "fs2 render queue starting: $(date)"
echo "Files: ${#QMDS[@]}"
echo "Log:   $LOG_FILE"
echo ""

START_SEC=$(date +%s)
SUCCEEDED=()
FAILED=()

for i in "${!QMDS[@]}"; do
  qmd="${QMDS[$i]}"
  job_num=$((i + 1))

  echo "------------------------------------------------------------"
  echo "[$job_num/${#QMDS[@]}] $qmd"
  echo "  Started: $(date)"
  echo "------------------------------------------------------------"

  job_start=$(date +%s)
  if quarto render "$qmd"; then
    SUCCEEDED+=("$qmd")
    status="SUCCESS"
  else
    FAILED+=("$qmd")
    status="FAILED"
  fi
  job_min=$(( ($(date +%s) - job_start) / 60 ))

  echo ""
  echo "[$job_num/${#QMDS[@]}] $status  (${job_min} min)"
  echo ""
done

TOTAL_MIN=$(( ($(date +%s) - START_SEC) / 60 ))

echo "============================================================"
echo "Queue complete: $(date)"
echo "Total elapsed: ${TOTAL_MIN} min"
echo ""
echo "Succeeded: ${#SUCCEEDED[@]}/${#QMDS[@]}"
for q in "${SUCCEEDED[@]}"; do echo "  OK   $q"; done
echo ""
echo "Failed:    ${#FAILED[@]}/${#QMDS[@]}"
for q in "${FAILED[@]}"; do echo "  FAIL $q"; done
echo "============================================================"

[ "${#FAILED[@]}" -gt 0 ] && exit 1
exit 0
