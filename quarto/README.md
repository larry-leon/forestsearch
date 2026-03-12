# quarto/README.md
#
# forestsearch — Extended Quarto Documentation
# =============================================
#
# This directory holds standalone Quarto documents that supplement the
# pkgdown site with material requiring Quarto-specific features.
#
# DIRECTORY LAYOUT
# ----------------
#
#   quarto/
#   ├── _quarto.yml                  project config (output-dir: _output)
#   ├── README.md                    this file
#   ├── index.qmd                    landing page linking both documents
#   ├── cox_causal_review.qmd        causal estimands reference
#   ├── biomarker_comparison.qmd     STEPP / MFPI / spline Cox comparison
#   └── _output/                     rendered HTMLs (see .gitignore note below)
#       ├── index.html
#       ├── cox_causal_review.html
#       └── biomarker_comparison.html
#
# RENDERING
# ---------
#
#   From the repo root — all documents:
#     quarto render quarto/
#
#   Single document:
#     quarto render quarto/cox_causal_review.qmd
#
#   From RStudio console:
#     source("render_quarto_docs.R")
#
#   Output lands in quarto/_output/.
#   Each HTML is fully self-contained (embed-resources: true).
#
# .Rbuildignore  (add these lines so CRAN ignores this directory entirely)
# ------------------------------------------------------------------------
#
#   ^quarto$
#   ^quarto/.*
#   ^render_quarto_docs\.R$
#
# .gitignore options
# ------------------
#
#   Option A — commit rendered HTMLs (simplest; no CI needed):
#     Do NOT add quarto/_output/ to .gitignore.
#     Render locally, commit the HTMLs, share the repo link.
#
#   Option B — ignore rendered output (keep git lean):
#     Add to .gitignore:
#       quarto/_output/
#     HTMLs are regenerated on demand with `quarto render quarto/`.
#
# NOTE ON biomarker_comparison.qmd
# ---------------------------------
#
#   The vignette: front-matter block has been removed from this file.
#   It is NOT a formal CRAN vignette; it lives here as a standalone
#   reference document only.
