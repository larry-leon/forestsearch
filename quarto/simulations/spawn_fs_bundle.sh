#!/usr/bin/env bash
# =============================================================================
# spawn_fs_bundle.sh
#
# Propagate a new FS bundle (e.g., fs2) across all existing sg_focus
# qmds within a study by cloning each *_<OLD>.qmd to *_<NEW>.qmd and
# substituting the bundle token throughout.
#
# Cross-platform: macOS (BSD sed) and Linux (GNU sed).
#
# Usage:
#   ./quarto/simulations/spawn_fs_bundle.sh OLD_BUNDLE NEW_BUNDLE [STUDY]
#
# Examples:
#   # Spawn fs2 qmds from existing fs1 qmds in the actg175 study:
#   ./quarto/simulations/spawn_fs_bundle.sh fs1 fs2
#
#   # Same, but explicit about the study (useful when multiple studies
#   # exist):
#   ./quarto/simulations/spawn_fs_bundle.sh fs1 fs2 actg175
#
# What it does:
#   1. Detects the OS (sed syntax differs between BSD and GNU).
#   2. Verifies the canonical bundle file exists for NEW_BUNDLE:
#        quarto/simulations/<study>/_fsparams/<NEW>_params.R
#      If absent, exits with a clear instruction. The user creates the
#      bundle file first (typically by copying <OLD>_params.R and
#      editing what differs).
#   3. Locates all *_<OLD>.qmd files in the study directory.
#   4. For each, creates the *_<NEW>.qmd sibling with OLD->NEW
#      substitution applied throughout the new file.
#
# What it does NOT do:
#   - Does not create or edit the <NEW>_params.R bundle file. You
#     define what the new bundle means yourself (rename function,
#     edit values, update header comment).
#   - Does not render the new qmds. Render each in 'full' mode after.
#
# Run from:
#   The forestsearch package root (the directory containing 'quarto/').
# =============================================================================

set -euo pipefail

# ---- Argument parsing --------------------------------------------------------
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
  echo "Usage: $0 OLD_BUNDLE NEW_BUNDLE [STUDY]"
  echo ""
  echo "Example: $0 fs1 fs2"
  echo "         $0 fs1 fs2 actg175"
  exit 1
fi

OLD_BUNDLE="$1"
NEW_BUNDLE="$2"
STUDY="${3:-actg175}"

# ---- Detect OS for sed flavour ----------------------------------------------
OS="$(uname -s)"
case "$OS" in
  Darwin) SED_INPLACE=(sed -i '') ;;
  Linux)  SED_INPLACE=(sed -i)    ;;
  *)      echo "ERROR: Unsupported OS '$OS'." ; exit 1 ;;
esac

# ---- Path setup --------------------------------------------------------------
STUDY_DIR="quarto/simulations/${STUDY}"
FSPARAMS_DIR="${STUDY_DIR}/_fsparams"
NEW_PARAMS_FILE="${FSPARAMS_DIR}/${NEW_BUNDLE}_params.R"

# ---- Sanity checks -----------------------------------------------------------
if [ ! -d "$STUDY_DIR" ]; then
  echo "ERROR: study directory not found: $STUDY_DIR"
  echo "  Run this script from the forestsearch package root."
  echo "  Current directory: $(pwd)"
  exit 1
fi

if [ ! -d "$FSPARAMS_DIR" ]; then
  echo "ERROR: _fsparams/ directory not found: $FSPARAMS_DIR"
  exit 1
fi

if [ ! -f "$NEW_PARAMS_FILE" ]; then
  echo "ERROR: bundle definition file not found: $NEW_PARAMS_FILE"
  echo ""
  echo "  Define ${NEW_BUNDLE} first by creating $NEW_PARAMS_FILE."
  echo "  Typical approach:"
  echo "    cp ${FSPARAMS_DIR}/${OLD_BUNDLE}_params.R \\"
  echo "       $NEW_PARAMS_FILE"
  echo "    # Then edit $NEW_PARAMS_FILE:"
  echo "    #   - Rename get_${OLD_BUNDLE}_params() to get_${NEW_BUNDLE}_params()"
  echo "    #   - Update header comment describing what ${NEW_BUNDLE} is for"
  echo "    #   - Edit parameter values that differ from ${OLD_BUNDLE}"
  exit 1
fi

# Verify the function name matches what the qmds will call
EXPECTED_FN="get_${NEW_BUNDLE}_params"
if ! grep -q "^${EXPECTED_FN} <- function" "$NEW_PARAMS_FILE"; then
  echo "WARNING: $NEW_PARAMS_FILE does not appear to define ${EXPECTED_FN}()."
  echo "  The cloned qmds will fail to find the function unless you rename"
  echo "  the function in $NEW_PARAMS_FILE."
  echo ""
  echo "Continue anyway? [y/N] "
  read -r REPLY
  if [ "$REPLY" != "y" ] && [ "$REPLY" != "Y" ]; then
    exit 1
  fi
fi

# ---- Locate source qmds -----------------------------------------------------
mapfile -t SOURCE_QMDS < <(
  find "$STUDY_DIR" -maxdepth 1 -type f -name "*_${OLD_BUNDLE}.qmd" \
    2>/dev/null | sort
)

if [ "${#SOURCE_QMDS[@]}" -eq 0 ]; then
  echo "ERROR: No qmds matching '*_${OLD_BUNDLE}.qmd' found in $STUDY_DIR."
  exit 1
fi

# ---- Plan: list what we'll do, check for conflicts --------------------------
echo "Spawn plan: ${OLD_BUNDLE} -> ${NEW_BUNDLE} in $STUDY_DIR"
echo ""
echo "Bundle file (must exist): $NEW_PARAMS_FILE  ✓"
echo ""
echo "Cloning ${#SOURCE_QMDS[@]} qmd(s):"

CONFLICTS=0
for src in "${SOURCE_QMDS[@]}"; do
  tgt="${src//_${OLD_BUNDLE}.qmd/_${NEW_BUNDLE}.qmd}"
  if [ -e "$tgt" ]; then
    echo "  CONFLICT: $tgt already exists"
    CONFLICTS=$((CONFLICTS + 1))
  else
    echo "  $src"
    echo "  -> $tgt"
  fi
done

if [ "$CONFLICTS" -gt 0 ]; then
  echo ""
  echo "ERROR: $CONFLICTS target file(s) already exist."
  echo "  Delete them first if you really want to re-spawn."
  exit 1
fi

# ---- Execute ----------------------------------------------------------------
echo ""
echo "OS detected: $OS"
echo "Proceeding..."
echo ""

for src in "${SOURCE_QMDS[@]}"; do
  tgt="${src//_${OLD_BUNDLE}.qmd/_${NEW_BUNDLE}.qmd}"
  cp "$src" "$tgt"
  "${SED_INPLACE[@]}" "s/${OLD_BUNDLE}/${NEW_BUNDLE}/g" "$tgt"
  echo "  wrote $tgt"
done

echo ""
echo "Done. Spawned ${#SOURCE_QMDS[@]} new qmd(s)."
echo ""
echo "Next steps:"
echo "  1. Spot-check one of the new qmds to verify substitutions."
echo "  2. Render each in 'full' mode (one render per qmd):"
echo ""
for src in "${SOURCE_QMDS[@]}"; do
  tgt="${src//_${OLD_BUNDLE}.qmd/_${NEW_BUNDLE}.qmd}"
  echo "       quarto render $tgt"
done
echo ""
echo "  3. Add ${NEW_BUNDLE} to the cross-cut summary's fs_levels vector:"
echo "       ${STUDY_DIR}/_summaries/*_focusAll_fsAll_summary.qmd"
echo ""
echo "  4. Re-render the cross-cut summary and the fsparams_summary:"
echo "       quarto render ${STUDY_DIR}/_summaries/<...>.qmd"
echo "       quarto render ${FSPARAMS_DIR}/fsparams_summary.qmd"
