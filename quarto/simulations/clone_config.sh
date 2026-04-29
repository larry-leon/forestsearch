#!/usr/bin/env bash
# =============================================================================
# clone_config.sh
#
# Clone an existing per-config qmd to a new one by renaming the file
# and replacing token references inside it. Works on macOS and Linux
# (System76 / Pop!_OS).
#
# Usage:
#   ./quarto/simulations/clone_config.sh OLD_TOKEN NEW_TOKEN
#
# Examples:
#   # Clone hrMaxSG_fs1 to hr_fs1 (focus swap):
#   ./quarto/simulations/clone_config.sh hrMaxSG hr
#
#   # Clone hr_fs1 to hr_fs2 (single-config bundle swap):
#   ./quarto/simulations/clone_config.sh fs1 fs2
#
# What it does:
#   1. Detects the OS (sed syntax differs between BSD and GNU).
#   2. Finds the qmd whose filename contains _<OLD>_ under
#      quarto/simulations/ (excluding _summaries, _fsparams, _archive).
#      If multiple qmds match, lists them and exits.
#   3. Copies it to the analogous _<NEW>_ filename.
#   4. Replaces every occurrence of OLD_TOKEN with NEW_TOKEN inside
#      the new qmd.
#
# What it does NOT do:
#   - Does not create the .rds bundle (render the new qmd in 'full'
#     mode for that).
#   - Does not edit any FS bundle parameters. If swapping fs1->fs2
#     also changes the bundle definition, edit _fsparams/fsN_params.R
#     yourself first (or use spawn_fs_bundle.sh, which handles the
#     entire fs spawn across all sg_focus configs at once).
#
# Run from:
#   The forestsearch package root (the directory containing 'quarto/').
# =============================================================================

set -euo pipefail

# ---- Argument parsing --------------------------------------------------------
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 OLD_TOKEN NEW_TOKEN"
  echo ""
  echo "Examples:"
  echo "  $0 hrMaxSG hr      # clone hrMaxSG config to hr config"
  echo "  $0 fs1 fs2         # clone one fs1 config to one fs2 config"
  exit 1
fi

OLD="$1"
NEW="$2"

# ---- Detect OS for sed flavour ----------------------------------------------
OS="$(uname -s)"
case "$OS" in
  Darwin) SED_INPLACE=(sed -i '') ;;
  Linux)  SED_INPLACE=(sed -i)    ;;
  *)      echo "ERROR: Unsupported OS '$OS'." ; exit 1 ;;
esac

# ---- Locate qmd to clone -----------------------------------------------------
SIM_ROOT="quarto/simulations"

if [ ! -d "$SIM_ROOT" ]; then
  echo "ERROR: '$SIM_ROOT' not found."
  echo "  Run this script from the forestsearch package root."
  echo "  Current directory: $(pwd)"
  exit 1
fi

# Find candidate qmds across all study directories whose filename
# contains _<OLD>_ as a stem token. Exclude any path under _summaries,
# _fsparams, or _archive.
mapfile -t CANDIDATES < <(
  find "$SIM_ROOT" -maxdepth 2 -type f -name "*_${OLD}_*.qmd" \
    -not -path "*_summaries*" -not -path "*_fsparams*" \
    -not -path "*_archive*" 2>/dev/null | sort
)

if [ "${#CANDIDATES[@]}" -eq 0 ]; then
  echo "ERROR: No qmds matching '*_${OLD}_*.qmd' found under $SIM_ROOT."
  exit 1
fi

if [ "${#CANDIDATES[@]}" -gt 1 ]; then
  echo "Multiple matching qmds found. Please specify which to clone"
  echo "(use a more specific OLD_TOKEN, or clone manually):"
  for i in "${!CANDIDATES[@]}"; do
    echo "  [$((i+1))] ${CANDIDATES[$i]}"
  done
  exit 1
fi

SOURCE_QMD="${CANDIDATES[0]}"
TARGET_QMD="${SOURCE_QMD//_${OLD}_/_${NEW}_}"

# ---- Sanity checks -----------------------------------------------------------
if [ "$SOURCE_QMD" = "$TARGET_QMD" ]; then
  echo "ERROR: source and target paths are identical."
  exit 1
fi

if [ -e "$TARGET_QMD" ]; then
  echo "ERROR: target file already exists: $TARGET_QMD"
  echo "  Delete it first if you really want to re-clone."
  exit 1
fi

# ---- Clone + replace ---------------------------------------------------------
echo "OS detected: $OS"
echo "Cloning:"
echo "  source: $SOURCE_QMD"
echo "  target: $TARGET_QMD"
echo

cp "$SOURCE_QMD" "$TARGET_QMD"
"${SED_INPLACE[@]}" "s/${OLD}/${NEW}/g" "$TARGET_QMD"

echo "Done. New file:"
ls -la "$TARGET_QMD"
echo
echo "Next steps:"
echo "  1. Verify substitutions by opening $TARGET_QMD."
echo "  2. (If this is an fs* swap) ensure _fsparams/${NEW}_params.R"
echo "     exists with a get_${NEW}_params() function."
echo "  3. Render the new qmd in 'full' mode:"
echo "       quarto render $TARGET_QMD"
