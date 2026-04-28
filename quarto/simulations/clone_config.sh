#!/usr/bin/env bash
# =============================================================================
# clone_config.sh
#
# Clone an existing config to a new one by renaming files and replacing
# token references inside them.  Works on macOS and Linux (System76 / Pop!_OS).
#
# Usage:
#   ./clone_config.sh OLD_TOKEN NEW_TOKEN
#
# Example:
#   ./clone_config.sh hrMaxSG hr
#
# What it does:
#   1. Detects whether you're on macOS or Linux (sed syntax differs).
#   2. Copies the four files of the OLD_TOKEN config to NEW_TOKEN names:
#        actg175_binary_m1_harm_OLD_fs1.qmd  -> ..._NEW_fs1.qmd
#        actg175_binary_m1_harm_OLD_fs1.R    -> ..._NEW_fs1.R
#        sim_actg175_binary_m1_harm_OLD_fs1.R-> sim_..._NEW_fs1.R
#      (does NOT copy the .rds — that's an artifact, will be created
#       by running the new runner)
#   3. Replaces every occurrence of OLD_TOKEN with NEW_TOKEN inside each
#      of the three new files.
#
# What it does NOT do:
#   - Does not create the .rds bundle (run the new runner script for that).
#   - Does not adjust FS parameter values.  If the new token implies
#     different parameter values (e.g., fs1 -> fs2 means a different
#     parameter bundle), you must edit those by hand AFTER this script.
#
# Run from:
#   The forestsearch package root (the directory containing 'quarto/').
# =============================================================================

set -euo pipefail

# ---- Argument parsing --------------------------------------------------------
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 OLD_TOKEN NEW_TOKEN"
  echo "Example: $0 hrMaxSG hr"
  exit 1
fi

OLD="$1"
NEW="$2"
DIR="quarto/simulations/actg175"
STEM_OLD="actg175_binary_m1_harm_${OLD}_fs1"
STEM_NEW="actg175_binary_m1_harm_${NEW}_fs1"

# ---- Sanity checks -----------------------------------------------------------
if [ ! -d "$DIR" ]; then
  echo "ERROR: '$DIR' not found."
  echo "  Run this script from the forestsearch package root."
  echo "  Current directory: $(pwd)"
  exit 1
fi

for ext in qmd R; do
  if [ ! -f "$DIR/${STEM_OLD}.${ext}" ]; then
    echo "ERROR: source file not found: $DIR/${STEM_OLD}.${ext}"
    exit 1
  fi
done
if [ ! -f "$DIR/sim_${STEM_OLD}.R" ]; then
  echo "ERROR: source file not found: $DIR/sim_${STEM_OLD}.R"
  exit 1
fi

# ---- Detect OS for sed flavour ----------------------------------------------
# macOS uses BSD sed (requires '' after -i); Linux uses GNU sed (no '' needed).
OS="$(uname -s)"
case "$OS" in
  Darwin)
    SED_INPLACE=(sed -i '')
    ;;
  Linux)
    SED_INPLACE=(sed -i)
    ;;
  *)
    echo "ERROR: Unsupported OS '$OS'.  This script supports macOS and Linux."
    exit 1
    ;;
esac

echo "OS detected: $OS"
echo "Cloning: $STEM_OLD  ->  $STEM_NEW"
echo

# ---- Copy + replace ---------------------------------------------------------
cp -n "$DIR/${STEM_OLD}.qmd"   "$DIR/${STEM_NEW}.qmd"
cp -n "$DIR/${STEM_OLD}.R"     "$DIR/${STEM_NEW}.R"
cp -n "$DIR/sim_${STEM_OLD}.R" "$DIR/sim_${STEM_NEW}.R"

echo "Files copied. Now substituting '$OLD' -> '$NEW' inside them..."

"${SED_INPLACE[@]}" "s/${OLD}/${NEW}/g" "$DIR/${STEM_NEW}.qmd"
"${SED_INPLACE[@]}" "s/${OLD}/${NEW}/g" "$DIR/${STEM_NEW}.R"
"${SED_INPLACE[@]}" "s/${OLD}/${NEW}/g" "$DIR/sim_${STEM_NEW}.R"

echo
echo "Done.  New files:"
ls -la "$DIR/${STEM_NEW}.qmd" "$DIR/${STEM_NEW}.R" "$DIR/sim_${STEM_NEW}.R"
echo
echo "Next steps:"
echo "  1. Verify the substitutions by opening any of the new files."
echo "  2. (Optional) Edit FS parameters if '$NEW' implies different values."
echo "  3. Run the new runner:"
echo "       Rscript $DIR/${STEM_NEW}.R"
