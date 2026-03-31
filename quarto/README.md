# quarto/README.md
#
# forestsearch --- Extended Quarto Documentation
# =============================================
#
# This directory holds standalone Quarto documents that supplement the
# pkgdown site with material requiring Quarto-specific features.
#
# DIRECTORY LAYOUT
# ----------------
#
#   quarto/
#   +-- _quarto.yml                  project config (output-dir: _output)
#   +-- README.md                    this file
#   +-- _output/                     rendered HTML (Quarto-managed, disposable)
#   +-- _data/                       .rds calibration bundles (persistent, precious)
#   |
#   |   # Theory
#   +-- theory_glm_detection_probability.qmd          (T1)
#   +-- theory_glm_leff_correction.qmd                (T2)
#   |
#   |   # Calibration (expensive, render once)
#   +-- calibration_binary_leff.qmd                   (C1)
#   +-- calibration_survival_fourway.qmd              (C2)
#   +-- calibration_survival_leff_grid.qmd            (C3)
#   |
#   |   # Selection (fast, loads .rds from _data/)
#   +-- selection_binary_frontier.qmd                 (S1)
#   +-- selection_survival_frontier.qmd               (S2)
#   +-- selection_glm_guide.qmd                       (S3)
#   |
#   |   # Validation
#   +-- validation_glm_simulation_study.qmd           (V1)
#   +-- validation_hte_tests_crump.qmd                (V2)
#   |
#   |   # Earlier standalone documents
#   +-- cox_causal_review.qmd
#   +-- biomarker_comparison.qmd
#   +-- index.qmd                    landing page
#   |
#   +-- _output/                     rendered HTMLs
#   |   +-- *.html
#   |
#   +-- _data/                       persistent .rds results
#       +-- calibration_binary_leff.rds
#       +-- calibration_survival_leff.rds
#       +-- calibration_survival_fourway_YYYYMMDD.rds
#
#
# TWO OUTPUT DIRECTORIES
# ----------------------
#
#   _output/   Quarto HTML renders. Disposable --- re-rendered any time.
#              Managed by Quarto via output-dir in _quarto.yml.
#
#   _data/     .rds calibration and simulation bundles. Persistent ---
#              hours of computation. Produced by calibration documents
#              (C1, C2, C3) via save_leff_calibration(). Consumed by
#              selection and validation documents via load_leff_calibration().
#
#
# RENDERING
# ---------
#
#   From the repo root --- all documents:
#     quarto render quarto/
#
#   Single document:
#     quarto render quarto/calibration_survival_fourway.qmd
#
#   From RStudio console:
#     quarto::quarto_render("quarto/selection_binary_frontier.qmd")
#
#   Recommended render order:
#     Step 1 (once):  calibration_binary_leff.qmd, calibration_survival_fourway.qmd
#     Step 2 (fast):  selection_binary_frontier.qmd, selection_survival_frontier.qmd
#     Step 3 (mins):  validation_glm_simulation_study.qmd, validation_hte_tests_crump.qmd
#
#   See glm_extension_render_guide.qmd for full details.
#
#
# .Rbuildignore  (add these lines so CRAN ignores this directory entirely)
# ------------------------------------------------------------------------
#
#   ^quarto$
#   ^quarto/.*
#   ^render_quarto_docs\.R$
#
#
# .gitignore options
# ------------------
#
#   Option A --- commit rendered HTMLs, ignore .rds (recommended):
#     quarto/_data/
#
#   Option B --- ignore both rendered output and data:
#     quarto/_output/
#     quarto/_data/
#
#   The _data/ directory contains large .rds files (potentially 100+ MB)
#   that should not be in version control.  Regenerate by rendering the
#   calibration documents.
