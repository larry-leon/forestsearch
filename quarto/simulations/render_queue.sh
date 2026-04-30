#!/usr/bin/env bash
# =============================================================================
# render_queue.sh
#
# Serially render a list of quarto documents, one after another, with
# per-render timing and a final summary. Designed for overnight batch
# runs where you want each render to have full machine resources (no
# concurrent jobs competing for cores).
#
# Cross-platform: macOS and Linux.
#
# Usage:
#   ./quarto/simulations/render_queue.sh QMD_FILE [QMD_FILE ...]
#
# Examples:
#   # Render three configs sequentially:
#   ./quarto/simulations/render_queue.sh \
#     quarto/simulations/actg175/actg175_binary_m1_harm_hr_fs1.qmd \
#     quarto/simulations/actg175/actg175_binary_m1_harm_hrMaxSG_fs1.qmd \
#     quarto/simulations/actg175/actg175_binary_m1_harm_maxSG_fs1.qmd
#
#   # Render an entire fs2 spawn (after spawn_fs_bundle.sh):
#   ./quarto/simulations/render_queue.sh \
#     quarto/simulations/actg175/actg175_binary_m1_harm_*_fs2.qmd
#
# What it does:
#   1. Verifies each input file exists before starting.
#   2. Renders each in order, one at a time (waits for each to finish).
#   3. Captures stdout/stderr from each render to a per-job log file
#      under quarto/simulations/<study>/logs/.
#   4. Tracks success/failure per job; prints a summary at the end.
#   5. Continues to the next job even if one fails (so a single bad
#      render doesn't kill the entire overnight queue).
#
# What it does NOT do:
#   - Does not edit qmd setup chunks. Each qmd must already have its
#     desired sim_mode set (typically "full" for production renders).
#   - Does not run jobs in parallel. By design — that's what makes
#     this safe for sharing one machine's resources.
#
# Run from:
#   The forestsearch package root (the directory containing 'quarto/').
# =============================================================================

set -uo pipefail
# Note: NOT using -e here. We want the loop to continue past failed
# renders rather than abort on the first error.

# ---- Argument parsing --------------------------------------------------------
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 QMD_FILE [QMD_FILE ...]"
  echo ""
  echo "Examples:"
  echo "  $0 quarto/simulations/actg175/actg175_binary_m1_harm_hr_fs1.qmd \\"
  echo "     quarto/simulations/actg175/actg175_binary_m1_harm_hrMaxSG_fs1.qmd"
  echo ""
  echo "  $0 quarto/simulations/actg175/*_fs2.qmd"
  exit 1
fi

QMDS=("$@")

# ---- Pre-flight: verify all files exist before starting ---------------------
echo "Pre-flight check: verifying ${#QMDS[@]} input file(s)..."
MISSING=0
for qmd in "${QMDS[@]}"; do
  if [ ! -f "$qmd" ]; then
    echo "  MISSING: $qmd"
    MISSING=$((MISSING + 1))
  fi
done

if [ "$MISSING" -gt 0 ]; then
  echo ""
  echo "ERROR: $MISSING file(s) not found. Aborting before any renders."
  exit 1
fi

echo "All ${#QMDS[@]} file(s) present. Starting queue."
echo ""

# ---- Setup logging -----------------------------------------------------------
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_BASE_DIR="quarto/simulations"
QUEUE_LOG="${LOG_BASE_DIR}/queue_${TIMESTAMP}.log"

mkdir -p "$LOG_BASE_DIR"

# Tee everything to the queue log so a single file captures the run summary
exec > >(tee -a "$QUEUE_LOG") 2>&1

echo "============================================================"
echo "Render queue started: $(date)"
echo "Queue log: $QUEUE_LOG"
echo "============================================================"
echo ""

# ---- Track results -----------------------------------------------------------
TOTAL=${#QMDS[@]}
SUCCEEDED=()
FAILED=()
QUEUE_START_SEC=$(date +%s)

# ---- Run the queue -----------------------------------------------------------
for i in "${!QMDS[@]}"; do
  qmd="${QMDS[$i]}"
  job_num=$((i + 1))

  # Per-job log: lives next to the qmd, in the study's logs/ subdir if
  # one exists (or directly in the study dir otherwise).
  qmd_dir="$(dirname "$qmd")"
  qmd_stem="$(basename "$qmd" .qmd)"
  log_dir="${qmd_dir}/logs"
  mkdir -p "$log_dir"
  job_log="${log_dir}/render_${qmd_stem}_${TIMESTAMP}.log"

  echo "------------------------------------------------------------"
  echo "[$job_num/$TOTAL] Rendering: $qmd"
  echo "  Started: $(date)"
  echo "  Per-job log: $job_log"
  echo "------------------------------------------------------------"

  job_start_sec=$(date +%s)

  # Run the render; capture stdout+stderr to the per-job log.
  # We let it run inline as well so the queue log shows progress.
  if quarto render "$qmd" 2>&1 | tee "$job_log"; then
    job_status="SUCCESS"
    SUCCEEDED+=("$qmd")
  else
    job_status="FAILED"
    FAILED+=("$qmd")
  fi

  job_end_sec=$(date +%s)
  job_elapsed_min=$(( (job_end_sec - job_start_sec) / 60 ))

  echo ""
  echo "[$job_num/$TOTAL] $job_status: $qmd"
  echo "  Elapsed: ${job_elapsed_min} min"
  echo "  Finished: $(date)"
  echo ""
done

# ---- Final summary -----------------------------------------------------------
QUEUE_END_SEC=$(date +%s)
QUEUE_TOTAL_MIN=$(( (QUEUE_END_SEC - QUEUE_START_SEC) / 60 ))

echo "============================================================"
echo "Queue complete: $(date)"
echo "Total elapsed: ${QUEUE_TOTAL_MIN} min"
echo ""
echo "Succeeded: ${#SUCCEEDED[@]} / $TOTAL"
for q in "${SUCCEEDED[@]}"; do
  echo "  ✓ $q"
done
echo ""
echo "Failed: ${#FAILED[@]} / $TOTAL"
for q in "${FAILED[@]}"; do
  echo "  ✗ $q"
done
echo "============================================================"

# Exit with non-zero if any failures, so wrapper scripts can detect.
if [ "${#FAILED[@]}" -gt 0 ]; then
  exit 1
fi
exit 0
