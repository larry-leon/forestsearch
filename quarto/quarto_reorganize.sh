#!/bin/bash
# =============================================================================
# quarto_reorganize.sh
# Reorganize quarto/ directory into topic subdirectories (Option A)
#
# Run from the quarto/ directory:
#   cd /path/to/forestsearch/quarto
#   bash quarto_reorganize.sh
#
# Prerequisites:
#   - Clean git working tree (commit or stash first)
#   - Review the _quarto.yml change (removes output-dir)
#
# What this does:
#   1. Creates topic subdirectories
#   2. Moves files with git mv (preserves history)
#   3. Updates _data paths in moved .qmd files
#   4. Removes output-dir from _quarto.yml (renders HTML alongside .qmd)
#   5. Moves existing _output/*.html to their new topic directories
#
# After running:
#   - git diff --stat  to review all changes
#   - Render one doc per directory to verify paths resolve
#   - git add -A && git commit -m "Reorganize quarto/ into topic directories"
# =============================================================================

set -euo pipefail

if [ ! -f "_quarto.yml" ]; then
  echo "ERROR: Run this script from the quarto/ directory."
  exit 1
fi

echo "=== Step 1: Create directories ==="
mkdir -p simulations calibration grf applications guides qc archive

# =============================================================================
# Step 2: Move files
# =============================================================================

echo "=== Step 2a: Production simulations ==="
git mv actg175_binary_benefit_simulations.qmd           simulations/
git mv actg175_continuous_simulations.qmd                simulations/
git mv actg175_survival_benefit_simulations.qmd          simulations/
git mv reproduce_table1_m1.qmd                           simulations/

echo "=== Step 2b: Calibration & selection ==="
git mv calibration_binary_leff.qmd                       calibration/
git mv calibration_survival_fourway.qmd                  calibration/
git mv calibration_survival_leff_grid.qmd                calibration/
git mv calibration_pipeline_guide.qmd                    calibration/
git mv selection_binary_frontier.qmd                     calibration/
git mv selection_survival_frontier.qmd                   calibration/

echo "=== Step 2c: GRF standalone ==="
git mv grf_dmin_calibration_actg175_binary.qmd           grf/
git mv grf_dmin_calibration_binary.qmd                   grf/
git mv grf_review_and_recommendations.qmd                grf/
git mv grf_scoring_binary_simulation.qmd                 grf/
git mv grf_subgroup_identification.qmd                   grf/

echo "=== Step 2d: Applications & validation ==="
git mv biomarker_comparison.qmd                          applications/
git mv count_data_demo.qmd                               applications/
git mv gbsg_analysis.qmd                                 applications/
git mv gbsg_analysis_cox_poisson.qmd                     applications/
git mv gbsg_poisson_simulations.qmd                      applications/
git mv forestsearch_scenario_validation.qmd              applications/
git mv forestsearch_scenario_validation_100.qmd          applications/
git mv validation_glm_simulation_study.qmd               applications/
git mv validation_hte_tests_crump.qmd                    applications/

echo "=== Step 2e: Theory, guides & literature ==="
git mv theory_glm_detection_probability.qmd              guides/
git mv theory_glm_leff_correction.qmd                    guides/
git mv cox_causal_review.qmd                             guides/
git mv extreme_subgroups.qmd                             guides/
git mv actg175_literature_review.qmd                     guides/
git mv forestsearch_glm_statistical_review.md            guides/
git mv glm_hte_literature_summary.qmd                    guides/
git mv count_data_hte_summary.qmd                        guides/
git mv glm_extension_document_map.md                     guides/
git mv glm_extension_render_guide.qmd                    guides/
git mv glm_simulation_suite_design.qmd                   guides/
git mv forestsearch_outcome_orientation_guide.qmd        guides/
git mv selection_glm_guide.qmd                           guides/
git mv suggest_thresholds_guide.qmd                      guides/
git mv discussion_template.qmd                           guides/
git mv leff_workflow_guide.md                            guides/
git mv claude_code_guide_forestsearch.qmd                guides/
git mv vscode_claudecode_forestsearch.qmd                guides/

echo "=== Step 2f: QC & test scripts ==="
git mv qc_benefit_search_switching.R                     qc/
git mv qc_binary_pipeline.R                              qc/
git mv qc_exhaustive_glm_pipeline.R                      qc/
git mv qc_glm_dgm_extension.R                           qc/
git mv qc_glm_sign_diagnosis.R                          qc/
git mv qc_grf_actg175_binary_benefit.R                  qc/
git mv qc_grf_binary_benefit.R                          qc/
git mv qc_grf_dmin_calibration_binary.R                 qc/
git mv qc_grf_scenario_validation.R                     qc/
git mv qc_null_classification.R                         qc/
git mv qc_poisson_offset.R                              qc/
git mv qc_search_alignment_diagnostic.R                 qc/
git mv test_adverse_outcome.R                           qc/
git mv test_glm_simulation_suite.R                      qc/

echo "=== Step 2g: Archive (working copies & superseded) ==="
git mv actg175_binary_benefit_simulations_working.qmd    archive/
git mv actg175_binary_benefit_simulations_working2.qmd   archive/
git mv actg175_survival_benefit_simulations_working.qmd  archive/
git mv reproduce_table1_m1_working.qmd                   archive/
git mv actg175_simulations_binary_working.qmd            archive/
git mv actg175_simulations_1.qmd                         archive/
git mv actg175_simulations_benefit.qmd                   archive/
git mv actg175_simulations_binary.qmd                    archive/
git mv actg175_simulations_ddi-vs-zdv.qmd                archive/
git mv actg175_simulations_review.md                     archive/
git mv COMMIT_SUMMARY.md                                 archive/
git mv COMMIT_SUMMARY_benefit_inversion.md               archive/
git mv selection_binary_frontier.rmarkdown               archive/
# Non-tracked files (regular mv)
[ -f Rplots.pdf ] && mv Rplots.pdf archive/
[ -f .Rhistory ]  && mv .Rhistory  archive/

# Also move the orphaned HTML helper directory
[ -d "cox_vs_glm_approximation_null-search_3_files" ] && \
  mv "cox_vs_glm_approximation_null-search_3_files" archive/

# =============================================================================
# Step 3: Update _data paths in moved .qmd files
# =============================================================================
# Files that move into subdirectories need _data -> ../_data
# because Quarto sets the working directory to the .qmd's location.

echo "=== Step 3: Update _data paths ==="

# Pattern 1: output_dir <- "_data"
# Pattern 2: readRDS("_data/...") or file.path("_data", ...)
# Pattern 3: literal "_data/" in path strings

for dir in simulations calibration grf applications; do
  for f in "$dir"/*.qmd; do
    [ -f "$f" ] || continue

    # Count replacements for reporting
    n=$(grep -c '"_data"' "$f" 2>/dev/null || echo 0)
    n2=$(grep -c '"_data/' "$f" 2>/dev/null || echo 0)

    if [ "$((n + n2))" -gt 0 ]; then
      # output_dir <- "_data"  ->  output_dir <- "../_data"
      sed -i 's|"_data"|"../_data"|g' "$f"
      # "_data/filename"  ->  "../_data/filename"
      sed -i 's|"_data/|"../_data/|g' "$f"
      echo "  Updated $f ($((n + n2)) replacements)"
    fi
  done
done

# =============================================================================
# Step 4: Update _quarto.yml — remove output-dir
# =============================================================================
# With output-dir removed, Quarto renders HTML alongside the .qmd source.
# Each topic directory will contain both .qmd and .html files together.

echo "=== Step 4: Update _quarto.yml ==="

# Comment out or remove the output-dir line
sed -i 's|^output-dir:.*|# output-dir: _output  # removed: HTML now renders alongside .qmd|' _quarto.yml

echo "  _quarto.yml updated (output-dir commented out)"
cat _quarto.yml

# =============================================================================
# Step 5: Relocate existing _output/*.html to topic directories
# =============================================================================

echo "=== Step 5: Move existing HTML output ==="

if [ -d "_output" ]; then
  for html in _output/*.html; do
    [ -f "$html" ] || continue
    base=$(basename "$html" .html)

    # Try to find the matching .qmd in a subdirectory
    match=$(find simulations calibration grf applications guides \
            -name "${base}.qmd" 2>/dev/null | head -1)

    if [ -n "$match" ]; then
      dest_dir=$(dirname "$match")
      mv "$html" "$dest_dir/"
      echo "  $html -> $dest_dir/"
    else
      echo "  $html -> no matching .qmd found (left in _output/)"
    fi
  done

  # Remove _output if empty
  rmdir _output 2>/dev/null && echo "  _output/ removed (empty)" || \
    echo "  _output/ still has files (check manually)"
fi

# =============================================================================
# Step 6: Summary
# =============================================================================

echo ""
echo "=== Reorganization complete ==="
echo ""
echo "Root contents:"
ls -F | head -20
echo ""
echo "Subdirectories:"
for d in simulations calibration grf applications guides qc archive; do
  n=$(find "$d" -maxdepth 1 -type f | wc -l)
  printf "  %-15s %2d files\n" "$d/" "$n"
done
echo ""
echo "Next steps:"
echo "  1. git diff --stat                    # review changes"
echo "  2. quarto render simulations/reproduce_table1_m1.qmd   # test one"
echo "  3. git add -A && git commit -m 'Reorganize quarto/ into topic directories'"
echo ""
echo "If paths break, check for hardcoded '_data' references that"
echo "were not caught by the sed patterns above.  Search with:"
echo "  grep -rn '\"_data' simulations/ calibration/ grf/ applications/"
