
# Please fetch forestsearch codebase https://github.com/larry-leon/forestsearch
#
# In oc_analyses_gbsg.R and sim_aft_gbsg.R there are "gbsg" specific functions, create_gbsg_dgm() and simulate_from_gbsg_dgm().
# Can these be replaced with the more general functions generate_aft_dgm_flex() and simulate_from_dgm()?
# If so, please create a new "oc_analyses.R", including any necessary helpers, which comprehensively replaces oc_analyses_gbsg.R and sim_aft_gbsg.R.
# I would like to phase out any specific reference to "gbsg", except when using as an example dataset.   Please create a multiple step strategy so
# that we can verify any revisions at each step are working properly.
#
#
# Let's start over, this is going nowhere and in circles.   The original R codebase has now been restored.
# Let's start with create_gbsg_dgm(), and proceed in separate steps so that each change can be verified before
# proceeding to the next step.
# To check the success of migrating create_gbsg_dgm() the following legacy results
# should be created:
# dgm_alt <- create_gbsg_dgm(
#   model = "alt",
#   k_treat = 1.0,
#   k_inter = 2.0,      # Interaction effect multiplier
#   k_z3 = 1.0,
#   z1_quantile = 0.25, # ER threshold at 25th percentile
#   n_super = 5000,
#   cens_type = "weibull",
#   seed = 8316951,
#   verbose = TRUE
# )
# # Enrich DGM with CDE values (computed from super-population theta_0/theta_1)
# dgm_alt <- compute_dgm_cde(dgm_alt)
# # Examine the DGM (print method now shows both HR and AHR metrics)
# print(dgm_alt)
# GBSG-Based AFT Data Generating Mechanism (Aligned)
# ===================================================
#   Model type: alt
# Super-population size: 5000
# Effect Modifiers:
#   k_treat: 1
# k_inter: 2
# k_z3: 1
# Hazard Ratios (Cox-based, stacked PO):
#   Overall (causal): 0.7331
# Harm subgroup (H): 2.9651
# Complement (Hc): 0.6612
# Ratio HR(H)/HR(Hc): 4.4846
# Average Hazard Ratios (from loghr_po):
#   AHR (overall): 0.7431
# AHR_harm (H): 3.8687
# AHR_no_harm (Hc): 0.5848
# Ratio AHR(H)/AHR(Hc): 6.6157
# Subgroup definition: z1 == 1 & z3 == 1 (low ER & premenopausal)
# ER threshold: 8 (quantile = 0.25)
# Subgroup size: 634 (12.7%)
# Analysis variables: v1, v2, v3, v4, v5, v6, v7

# NOTE: on a file lock
# make sure no other Git processes are active
#ps aux | grep git
#rm -f .git/index.lock

# dev-workflow.R
# Development workflow helper script for ForestSearch package
# Source this file to load helper functions, or run sections interactively

# git rm -r --cached dev/private
# echo "dev/private/" >> .gitignore
# git commit -m "Stop tracking dev/private"
# git push

#
# Prior devtools-helpers file
library(usethis)
library(devtools)

# Step 1: Create the package structure
# This will create the template directory
# with Rproject setup
#usethis::create_package("../../GitHub/forestsearch")

# Copy functions into R directory
# Step 3: Add README and MIT license
usethis::use_readme_rmd(open = FALSE)
usethis::use_mit_license("Larry Leon")

# Step 4: Add dependencies to DESCRIPTION
# desc::desc_set_dep("survival", file = "DESCRIPTION")
# desc::desc_set_dep("ggplot2", file = "DESCRIPTION")
# desc::desc_set_dep("grf", file = "DESCRIPTION")
# desc::desc_set_dep("policytree", file = "DESCRIPTION")

pkgs <- c("data.table", "foreach", "future", "doFuture", "future.apply", "glmnet", "gt", "randomForest", "stringr", "survival", "grf",
          "ggplot2", "policytree", "future.callr")
for (pkg in pkgs) {
  desc::desc_set_dep(pkg, file = "DESCRIPTION")
}

#desc::desc_set_dep("weightedSurv", type = "Suggests", file = "DESCRIPTION")

#desc::desc_set_dep(package = "weightedsurv", type = "Imports")

desc::desc_set_remotes("larry-leon/weightedsurv")

usethis::use_package("patchwork", type = "Suggests")

# DiagrammeR for grf graph
usethis::use_package("DiagrammeR", type = "Suggests")

usethis::use_package("tidyr", type = "Suggests")

usethis::use_package("rlang")

usethis::use_package("dplyr")

# Required per mrct_simulation()

usethis::use_package("progressr")

usethis::use_package("cubature", type = "Suggests")


# Step 5: Generate documentation
# Also, run this if revising R files such as @importFrom

rm(list=ls())

.rs.restartR()

# Restart R (Cmd+Shift+F10)
remove.packages("forestsearch")


#source("apply_cran_fixes.R")
# removed from root after applied

# clean up old documentation
unlink("man/*.Rd")

devtools::document()

devtools::load_all(".")

devtools::check()

devtools::check(args = c("--as-cran", "--no-build-vignettes"))

devtools::clean_dll()
unlink("man/*.Rd")
unlink("Meta", recursive = TRUE)

# Restart R session (Ctrl+Shift+F10 in RStudio), then:

devtools::document()

devtools::install(dependencies = FALSE)



unlink("man/glm_cs_fit.Rd")
devtools::document()
devtools::check(args = c("--as-cran", "--no-build-vignettes"))

# CRAN-level check (stricter; use this before submission)
devtools::check(cran = TRUE)

devtools::clean_dll()

devtools::install(quick = TRUE)

devtools::install(dependencies = FALSE)

devtools::install()

devtools::build()

library(forestsearch)
args(forestsearch)
help(forestsearch)

# Testing
# Be careful: Claude tends to modify tests until they pass without anything actionable
devtools::test()


# From your forestsearch project directory:
devtools::clean_dll()
devtools::document()
devtools::install(dependencies = FALSE, upgrade = "never")

# Or
remove.packages("forestsearch")
devtools::install()


# In RStudio, from the forestsearch project directory:
devtools::clean_dll()
devtools::document()
devtools::install(upgrade = "never")

# 1. Remove the corrupted installation
remove.packages("forestsearch")
# 2. Restart R (Session > Restart R in RStudio, or quit/reopen)
# 3. Reinstall from source
library(pak)
pak::pak("larry-leon/forestsearch")


devtools::install(args = "--no-byte-compile")
library(forestsearch)
packageVersion("forestsearch")


usethis::use_testthat()
# Run just these tests
devtools::test(filter = "glm-pipeline")


# Nuclear option: manually delete the library folder, then reinstall
#unlink(find.package("forestsearch"), recursive = TRUE)
# Restart R again, then:
#pak::pak("larry-leon/forestsearch")


devtools::build()                          # creates the tarball

devtools::check_built(
  "/Users/larryleon/Documents/GitHub/forestsearch/dev/code-working/CRAN_submission_steps/version_010/forestsearch_0.1.0.tar.gz",
  args = "--as-cran"
)



# Build the tarball
#R CMD build forestsearch/
# Run the full CRAN check
#R CMD check --as-cran forestsearch_0.1.0.tar.gz


# corrupt .rdb cached.  Force remove it
# from terminal
rm -rf /home/larryleon/R/x86_64-pc-linux-gnu-library/4.5/forestsearch
# restart R (Session -> Restart R)

# Retain dev/working locally
# Remove from Git's index only (keeps local files intact)
#git rm -r --cached dev/working/
# Then commit
#git add .gitignore

# echo "dev/working/" >> .gitignore
# git rm -r --cached dev/working/
# git add .gitignore
# git commit -m "Stop tracking dev/working/, keep locally"

# End prior flow

dir.create("vignettes/articles", recursive = TRUE, showWarnings = FALSE)
usethis::use_build_ignore("vignettes/articles")

pkgdown::check_pkgdown()

# One-time setup: creates a gh-pages branch and GitHub Actions workflow
#usethis::use_pkgdown_github_pages()
#This does three things: sets the URL to https://larry-leon.github.io/forestsearch, creates .github/workflows/pkgdown.yaml for
#automatic rebuilds on push, and configures the gh-pages branch.

#After pushing, the site is live at https://larry-leon.github.io/forestsearch.

# usethis::create_github_token()
# # store the PAT
# gitcreds::gitcreds_set()
# usethis::gh_token_help()    # Should show your token status
# usethis::git_sitrep()       # Full diagnostics

# USING THIS
# Build locally, then push docs/ to gh-pages branch
#pkgdown::build_site()
# Deploy docs/ to gh-pages branch
#pkgdown::deploy_to_branch(branch = "gh-pages", commit_message = "Update pkgdown site")


# note: deploy_to_branch calls build_site(), so only run here if deploying
# otherwise will duplicate the compiling of documents

# 5-March-2026
# DO THIS!!!
usethis::use_github_action("pkgdown")

git add .github/workflows/pkgdown.yaml
git push


#git add vignettes/articles/extreme_subgroups.Rmd




# Don't do this --> use github actions
#pkgdown::deploy_to_branch()


git status          # docs/ files should show as modified
git add docs/
git commit -m "Rebuild pkgdown site"
git push

# Trying to find disconnect

#The docs/ directory isn't showing up in git status at all, which means it's either gitignored or was never tracked. Check:

cat .gitignore | grep docs



# Build only articles that haven't been built yet (skips already-built ones)
# DNW (does not work, rebuilds everything)

pkgdown::check_pkgdown()

pkgdown::build_articles(lazy = TRUE)


# Rebuild the articles index page and navbar only
pkgdown::build_articles_index()
pkgdown::build_home()


# The typical workflow when adding a single new vignette is:
# 1. Build just the new article
pkgdown::build_article("articles/biomarker_effects")
# 2. Refresh the index so it appears in the articles listing
pkgdown::build_articles_index()


#. Try this --
# Render only the new vignette(s)
pkgdown::build_article("articles/biomarker_effects")
pkgdown::build_article("articles/causal_effects_brief_review")

# Refresh the index page so they appear in the listing
pkgdown::build_articles_index()

# Rebuild the navbar (home page picks up the updated menu)
pkgdown::build_home()


# https://larry-leon.github.io/forestsearch

# Preview without building everything
pkgdown::build_reference()


pkgdown::build_site()

pkgdown::preview_site()

# Build single article for quick preview
#pkgdown::build_article("getting-started")

pkgdown::build_article("forestsearch")

pkgdown::build_article("articles/methodology")

pkgdown::build_article("articles/treatment_effect_definitions")

pkgdown::build_article("articles/paper_simulations")

pkgdown::build_article("articles/extreme_subgroups")

pkgdown::build_article("articles/biomarker_effects")

pkgdown::build_article(
  "articles/paper_simulations",
  quiet = FALSE
)


pkgdown::build_article(
  "articles/biomarker_effects",
  quiet = FALSE
)

browseURL("docs/articles/biomarker_effects.html")


pkgdown::build_article("causal_effects_brief_review")


pkgdown::build_article(
  "articles/causal_effects_brief_review",
  quiet = FALSE
)



# If you ever want a standalone PDF version later, you can render it separately
# with quarto render biomarker_effects.qmd --to pdf using a copy that has the pdf:
# block restored — but it should live outside the pkgdown build pipeline.


# Open it directly
#browseURL("docs/articles/articles/paper_simulations.html")

pkgdown::build_article(
  "articles/treatment_effect_definitions",
  quiet = FALSE
)

# Setup for storing diagrams for incorporation in pkgdown
#Step 1: Create the pkgdown/assets/figures/ Directory
#The pkgdown/ directory does not exist by default — simply create it.
#pkgdown automatically detects it and copies the contents of assets/
#to the site root on every build. Nothing needs to be declared in _pkgdown.yml.
#From R:
#dir.create("pkgdown/assets/figures", recursive = TRUE)

#Or from the terminal:
#mkdir -p pkgdown/assets/figures


# =============================================================================
# pkgdown Deployment Notes for forestsearch
# =============================================================================
#
# Step 1: Build and deploy the pkgdown site to the gh-pages branch
# -----------------------------------------------------------------------------
#
# From the package root directory in R:
#
#   pkgdown::deploy_to_branch()
#
# This calls build_site() internally, then commits and pushes the
# rendered HTML to the gh-pages branch. There is no need to call
# build_site() separately before deploying.
#
# Step 2: Enable GitHub Pages (one-time setup)
# -----------------------------------------------------------------------------
#
# 1. Navigate to: https://github.com/larry-leon/forestsearch
# 2. Go to: Settings > Pages (left sidebar)
# 3. Under "Source", select:
#      - Branch: gh-pages
#      - Folder: / (root)
# 4. Click "Save"
#
# The site will be live within 1-2 minutes at:
#
#   https://larry-leon.github.io/forestsearch
#
# Step 3: Verify deployment
# -----------------------------------------------------------------------------
#
# - Check for a green checkmark on the latest gh-pages branch commit
# - Or check Settings > Environments for an active "github-pages" deployment
# - Browse the live site to confirm pages render correctly
#
# Step 4: Remove GitHub Actions pkgdown workflow (one-time cleanup)
# -----------------------------------------------------------------------------
#
# If you receive email notifications from GitHub stating
# "pkgdown.yaml: All jobs have failed", this is a separate CI workflow
# at .github/workflows/pkgdown.yaml that tries to auto-build the site
# on every push. Since we deploy manually via deploy_to_branch(), this
# workflow is unnecessary and can be removed.
#
# From the terminal in the package root:
#
#   git rm .github/workflows/pkgdown.yaml
#   git commit -m "Remove pkgdown CI workflow (deploying manually)"
#   git push
#
# This stops the failing emails. The existing live site is unaffected.
# There is NO need to re-run deploy_to_branch() after this step.
#
# Ongoing workflow
# -----------------------------------------------------------------------------
#
# After initial setup, the only command needed is:
#
#   pkgdown::deploy_to_branch()
#
# Run this whenever you want to update the live site after making
# changes to documentation, vignettes, or package code. There is no
# need to run it after routine git pushes that don't affect the site.
#
# Notes
# -----------------------------------------------------------------------------
#
# - The GitHub Pages setting persists once configured. Each subsequent
#   call to pkgdown::deploy_to_branch() will automatically trigger a
#   rebuild without needing to revisit Settings.
#
# - To rebuild the site locally without deploying (preview only):
#
#     pkgdown::build_site()
#     pkgdown::preview_site()
#
# - To validate the _pkgdown.yml configuration:
#
#     pkgdown::check_pkgdown()
#
# =============================================================================


# re-run
#usethis::use_pkgdown_github_pages()

# IMPORTANT
# Restore the last committed version
git checkout HEAD -- _pkgdown.yml


pkgdown::check_pkgdown()    # Validate _pkgdown.yml against NAMESPACE
pkgdown::build_reference()  # Build just the function reference
pkgdown::build_articles()   # Build just the vignettes


pkgdown::as_pkgdown()$vignettes


# Alternative: Skip usethis entirely
# If you'd rather not configure a PAT right now, you can set up GitHub Pages manually:
# # 1. Build the site locally
pkgdown::build_site()

# 2. Deploy to gh-pages branch (uses git, not GitHub API)
pkgdown::deploy_to_branch(branch = "gh-pages")

# Then go to your repo on GitHub: Settings → Pages → Source → Deploy from a branch → gh-pages / / (root) → Save.
# Your site will be live at https://larry-leon.github.io/forestsearch within a minute or two.



# Confirm no scattered declarations remain
grep -rn "globalVariables" R/ | grep -v "globals.R"

# Find all globalVariables() calls outside globals.R
matches <- grep(
  "globalVariables",
  list.files("R", full.names = TRUE, pattern = "\\.R$") |>
    setdiff("R/globals.R") |>
    sapply(readLines, simplify = FALSE) |>
    unlist()
)


# Starting with fresh yaml

# Remove your current file first
file.rename("_pkgdown.yml", "_pkgdown_old.yml")

## Auto-generate from actual namespace
#pkgdown::template_reference()

# Revised yaml 9-Feb-2026
usethis::use_pkgdown()
pkgdown::build_site()

# See all exported functions
# feed output into claude
#sort(getNamespaceExports("forestsearch"))

pkgdown::check_pkgdown()

pkgdown::build_site()
# Opens in browser automatically

# Preview without building everything
pkgdown::build_reference()
pkgdown::preview_site()

# Trying to clean up yaml WTF?
# Get all exported functions
exports <- getNamespaceExports("forestsearch")

# Check each explicitly listed function in your _pkgdown.yml
yml_functions <- c(
  "forestsearch", "print.forestsearch", "summary.forestsearch",
  "forestsearch_bootstrap_dofuture", "bootstrap_results", "get_dfRes",
  "summarize_bootstrap_results", "summarize_bootstrap_subgroups",
  "format_bootstrap_table", "format_subgroup_summary_tables",
  "create_factor_summary_tables",
  "forestsearch_Kfold", "forestsearch_tenfold", "forestsearch_KfoldOut",
  "CV_sgs", "print.fs_kfold", "print.fs_tenfold",
  "subgroup.consistency",
  "grf.subg.harm.survival", "grf.estimates.out",
  "cox_summary", "cox_ahr_cde_analysis", "cox_cs_fit", "plot_subgroup_effects",
  "get_FSdata", "dummy_encode", "get_dfpred", "subgroup.search",
  "plot_subgroup_results_forestplot", "print.fs_forestplot", "plot.fs_forestplot",
  "create_subgroup_summary_df", "compute_sg_hr", "sens_text",
  "sg_tables", "sg_estimates",
  "generate_aft_dgm_flex", "create_gbsg_dgm", "create_dgm_for_mrct",
  "simulate_from_gbsg_dgm",
  "run_simulation_analysis", "default_fs_params", "default_grf_params",
  "summarize_simulation_results",
  "mrct_region_sims", "summaryout_mrct",
  "cv_summary_tables", "cv_metrics_tables",
  "get_param"
)

missing <- setdiff(yml_functions, exports)
cat("Missing from namespace:\n")
print(missing)



# =============================================================================
# SETUP - Run once when starting development
# =============================================================================

# Install development dependencies
install.packages(c("devtools", "usethis", "roxygen2", "testthat", "pkgdown", "spelling"))

# Initialize package infrastructure (run once)
if (FALSE) {

  usethis::use_testthat()
  usethis::use_vignette("getting-started")
  usethis::use_readme_rmd()
  usethis::use_news_md()
  usethis::use_github_action_check_standard()
  usethis::use_spell_check()
  usethis::use_code_of_conduct()
}

# =============================================================================
# DAILY WORKFLOW - Common development tasks
# =============================================================================

# Load all package functions into environment (like library() but for development)
devtools::load_all()

# Generate/update documentation from roxygen2 comments
devtools::document()

# Run R CMD check (do this frequently!)
devtools::check()

# Quick check - skips some slower tests
devtools::check(args = "--no-examples")

# Run tests only
devtools::test()

# Run tests for a specific file
devtools::test_active_file()  # In RStudio, tests file related to current file
testthat::test_file("tests/testthat/test-bootstrap.R")

# =============================================================================
# NAMESPACE & IMPORTS
# =============================================================================
# Fix common namespace issues

# Check for missing imports
devtools::check() |>
  # Look for: "no visible global function definition for"
  # Look for: "no visible binding for global variable"

  # Add imports to NAMESPACE via roxygen2 tags:
  # @importFrom data.table data.table .N setnames
  # @importFrom stats median quantile sd

  # For pipe operator
  usethis::use_pipe()

# For data.table .SD, .N, etc. (add to R/package-imports.R)
if (FALSE) {
  cat('
#\' @importFrom data.table .N .SD .I .GRP .BY .EACHI
#\' @importFrom rlang .data
NULL
')
}

# =============================================================================
# CRAN CHECKS
# =============================================================================

# Full CRAN-like check
devtools::check(cran = TRUE)

# Check on multiple R versions (via R-hub)
devtools::check_rhub()

# Check on Windows (if you're on Mac/Linux)
devtools::check_win_devel()

# Spell check
spelling::spell_check_package()

# Check URLs in documentation
urlchecker::url_check()

# =============================================================================
# BUILDING & INSTALLING
# =============================================================================
# Build source package
devtools::build()

# Build binary package
devtools::build(binary = TRUE)
# Install package locally
devtools::install()

# Install with vignettes
devtools::install(build_vignettes = TRUE)

# =============================================================================
# VIGNETTES
# =============================================================================

# Build vignettes
devtools::build_vignettes()

# Preview a specific vignette
rmarkdown::render("vignettes/getting-started.Rmd")

# Move working vignette to official vignettes folder
# file.copy("dev/vignettes-working/my-analysis.Rmd", "vignettes/", overwrite = TRUE)

# =============================================================================
# DOCUMENTATION WEBSITE (pkgdown)
# =============================================================================

# Initialize pkgdown (run once)
usethis::use_pkgdown()

# Build documentation website locally
pkgdown::build_site()

# Preview site
pkgdown::preview_site()

# Build single article for quick preview
pkgdown::build_article("getting-started")

# =============================================================================
# VERSION & RELEASES
# =============================================================================

# Increment version number
usethis::use_version("patch")  # 0.1.0 -> 0.1.1
usethis::use_version("minor")  # 0.1.0 -> 0.2.0
usethis::use_version("major")  # 0.1.0 -> 1.0.0

# Development version
usethis::use_dev_version()     # 0.1.0 -> 0.1.0.9000

# Update NEWS.md with changes before release

# Create GitHub release
usethis::use_github_release()

# =============================================================================
# DATA
# =============================================================================

# Save internal data (not exported, available to package functions)
# usethis::use_data(internal_object, internal = TRUE)

# Save exported data (available to users)
# usethis::use_data(exported_dataset, overwrite = TRUE)

# Document datasets in R/data.R

# =============================================================================
# HELPER FUNCTIONS FOR THIS PROJECT
# =============================================================================

#' Quick development cycle
#' @description Document, load, and optionally test
quick_dev <- function(test = FALSE) {
  devtools::document()
  devtools::load_all()
  if (test) devtools::test()
  invisible(TRUE)
}

#' Check for common issues before committing
pre_commit_check <- function() {
  message("=== Documenting ===")
  devtools::document()


  message("\n=== Running tests ===")
  test_results <- devtools::test()

  message("\n=== Quick check ===
")
  devtools::check(args = c("--no-manual", "--no-vignettes"))

  message("\n=== Spell check ===")
  spelling::spell_check_package()

  invisible(TRUE)
}

#' Full pre-release check
pre_release_check <- function() {
  message("=== Full CRAN check ===")
  devtools::check(cran = TRUE)

  message("\n=== URL check ===")
  urlchecker::url_check()

  message("\n=== Reverse dependency check ===")
  # revdepcheck::revdep_check()  # If you have reverse dependencies

  invisible(TRUE)
}

#' Clean up generated files
clean_package <- function() {
  unlink("man", recursive = TRUE)
  unlink("Meta", recursive = TRUE)
  unlink("doc", recursive = TRUE)
  unlink("inst/doc", recursive = TRUE)
  unlink(list.files(pattern = "\\.tar\\.gz$"))
  message("Cleaned generated files. Run devtools::document() to regenerate.")
}

# =============================================================================
# PROJECT-SPECIFIC WORKFLOWS
# =============================================================================

#' Regenerate synthetic data for examples
regenerate_example_data <- function() {
  source("data-raw/generate-example-data.R")
  message("Example data regenerated")
}

#' Run bootstrap simulation for testing
run_test_simulation <- function(n_boots = 10) {
  devtools::load_all()
  # Your simulation code here
  message(sprintf("Completed %d bootstrap iterations", n_boots))
}

# =============================================================================
# TROUBLESHOOTING NOTES
# =============================================================================

# Common errors and fixes:
#
# 1. "no visible global function definition for 'xyz'"
#    -> Add @importFrom pkg xyz to the function's roxygen block
#
# 2. "no visible binding for global variable 'xyz'"
#    -> For data.table: add utils::globalVariables(c("xyz")) in R/globals.R
#    -> Or use .data$xyz with @importFrom rlang .data
#
# 3. "library() call in package code"
#    -> Replace library(pkg) with pkg::function() or @importFrom
#
# 4. Unicode characters causing issues
#    -> Replace with ASCII: >= instead of ≥, use \u2265 in strings
#
# 5. Examples taking too long
#    -> Wrap slow examples in \donttest{} or \dontrun{}

# =============================================================================
# USEFUL REFERENCES
# =============================================================================

# R Packages book: https://r-pkgs.org/
# CRAN policies: https://cran.r-project.org/web/packages/policies.html
# Writing R Extensions: https://cran.r-project.org/doc/manuals/R-exts.html

message("dev-workflow.R loaded. Key functions: quick_dev(), pre_commit_check()")


library(usethis)
library(devtools)

# Step 1: Create the package structure
# This will create the template directory
# with Rproject setup
#usethis::create_package("../../GitHub/forestsearch")

# Copy functions into R directory
# Step 3: Add README and MIT license
usethis::use_readme_rmd(open = FALSE)
usethis::use_mit_license("Larry Leon")

# Step 4: Add dependencies to DESCRIPTION
# desc::desc_set_dep("survival", file = "DESCRIPTION")
# desc::desc_set_dep("ggplot2", file = "DESCRIPTION")
# desc::desc_set_dep("grf", file = "DESCRIPTION")
# desc::desc_set_dep("policytree", file = "DESCRIPTION")

pkgs <- c("data.table", "foreach", "future", "doFuture", "future.apply", "glmnet", "gt", "randomForest", "stringr", "survival", "grf",
          "ggplot2", "policytree", "future.callr")
for (pkg in pkgs) {
  desc::desc_set_dep(pkg, file = "DESCRIPTION")
}

#desc::desc_set_dep("weightedSurv", type = "Suggests", file = "DESCRIPTION")

#desc::desc_set_dep(package = "weightedsurv", type = "Imports")

desc::desc_set_remotes("larry-leon/weightedsurv")

usethis::use_package("patchwork", type = "Suggests")

# DiagrammeR for grf graph
usethis::use_package("DiagrammeR", type = "Suggests")

usethis::use_package("tidyr", type = "Suggests")

usethis::use_package("rlang")

usethis::use_package("dplyr")

# Step 5: Generate documentation
# Also, run this if revising R files such as @importFrom

rm(list=ls())

.rs.restartR()

# clean up old documentation
unlink("man/*.Rd")

devtools::document()

devtools::load_all()

devtools::check()

devtools::clean_dll()


# Run these before CRAN submission
devtools::check(cran = TRUE)     # Full CRAN validation
devtools::check_rhub()            # Check on multiple platforms
devtools::check_win_devel()       # Windows checks


tools::showNonASCIIfile("~/Documents/GitHub/forestsearch/R/summary_utility_functions.R")



# gitignore additions
#.Rproj.user
#vignettes/results/
#  dev/private/
#  *_files/
# vignettes/working/



# Stop tracking files in GitHub

#git rm -r --cached *_files/

# Issue with trying to remove weightedSurv and replace with weightedsurv
# Remove from search path and unload
library(forestsearch)  # load it first if not loaded
detach("package:forestsearch", unload = TRUE, force = TRUE)
# Clear workspace and restart
rm(list = ls())
.rs.restartR()

# Notes
# Incorporate in AI prompts when documenting
# Every \item in a \describe{} block must have both a label and a description.
# Do not leave a lone \item{...} without {...} after it.
# Do not next \item inside another \item
# Use plain text, dashes, or a single paragraph for subpoints.
#
#
# León LF, Jemielita T, Guo Z, Marceau West R, Anderson KM.
# Exploratory subgroup identification in the heterogeneous Cox model: A relatively simple procedure.
# Statistics in Medicine. 2024; 43(20): 3921-3942. doi: 10.1002/sim.10163


#roxygen2::roxygenise()


# Remove and reinstall the package
remove.packages("forestsearch")
# Clear package cache
.libPaths()  # Find your library path
# Then reinstall from your source
# Either from GitHub or local source:
#devtools::install_github("larry-leon/forestsearch")
# OR
devtools::load_all()  # if developing locally

#devtools::install()

# run this in terminal (next to console [go to tools terminal tab])
#git pull --no-rebase



gitcreds::gitcreds_set()

usethis::use_git()

usethis::use_github()

# If functions are in namespace but not directly loaded
# devtools::load_all()

# Or  access hidden files:  mypackage::my_function




#library(codetools)
#codetools::findGlobals(forestsearch_bootstrap_dofuture, merge = FALSE)$functions

codetools::findGlobals(forestsearch, merge = FALSE)$functions

codetools::findGlobals(subgroup.search, merge = FALSE)$functions

codetools::findGlobals(get_FSdata, merge = FALSE)$functions


check <- c("calc_cov",  "calculate_counts", "analyze_subgroups", "calculate_potential_hr","ci.est","count.id","CV_sgs",
           "cox_summary","df_counting","double_robust_scores", "extract_subgroup","format_results", "get_targetEst","getci_Cox",
           "getCIs","grf.estimates.out","hrCI_format","km_summary","n_pcnt","plot_subgroup","plot_weighted_km",
           "prepare_subgroup_data","quiet","rmst_calculation","sg_tables","sort_subgroups","SummaryStat","var_summary",
           "get_FSdata", "dummy","run_bootstrap",
           "forestsearch", "forestsearch_bootstrap_dofuture","get_combinations_info",
           "get_dfpred",
           "grf.subg.harm.survival",
           "subgroup.search",
           "subgroup.consistency",
           "lasso_selection",
           "get_Cox_sg",
           "get_conf_force",
           "filter_by_lassokeep",
           "is.continuous",
           "process_conf_force_expr",
           "is_flag_continuous",
           "is_flag_drop", "acm.disjctif",  "acm.util.df2", "acm.util.df", "dummy2","ztrail","one.zero",
           "get_dfRes", "get_subgroup_membership",
           "SG_tab_estimates",
           "prepare_data",
           "run_grf",
           "evaluate_subgroups",
           "summarize_results",
           "clean_data", "qlow", "qhigh","FS_labels","thiscut","get_cut_name",
           "bootstrap_results", "remove_redundant_subgroups", "sg_consistency_out","get_split_hr","cut_var",
           "bootstrap_ystar", "ensure_packages", "fit_cox_models", "build_cox_formula",
           "format_CI","setup_parallel_SGcons", "get_covs_in", "extract_idx_flagredundancy"
)

duplicated_elements <- check[duplicated(check)]
duplicated_elements


# Per Claude

# In your local forestsearch repo
# git checkout -b fix/must-fix-implementations
#
# # Copy the new files to R/ directory (as shown above)
#
# # Stage changes
# git add R/bootstrap_helpers.R
# git add R/summary_utility_functions.R
# git add R/improved_grf_functions.r
# git add R/input_validation_utils.R

# Commit
# git commit -m "Implement must-fix recommendations
#
# - Add comprehensive input validation to all functions
# - Fix division-by-zero issues throughout
# - Standardize variable naming conventions
# - Add input_validation_utils.R with validation helpers
# - Maintain 100% backward compatibility"
#
# # Push to GitHub
# git push origin fix/must-fix-implementations

# Then create a Pull Request on GitHub


# Generate documentation
devtools::document()

# Load your updated package
devtools::load_all()


# Run tests
devtools::test()

# Check package
devtools::check()

# Search all R files in your package
# grep -rn "enumerate" ~/path/to/ForestSearch/R/ --include="*.R"
#
# # Or from the package root
# cd ~/path/to/ForestSearch
# grep -rn "enumerate" R/
#
#
#   How to Update Project Files
#
# Go to your Claude Project (in the left sidebar, click on the project name)
# Open Project Settings (click the gear icon or "Project settings")
# Remove outdated files:
#
#   Find the files in the "Project knowledge" or "Files" section
# Delete the old versions
#
#
# Upload current files:
#
#   Click "Add content" or "Upload files"
# Select the current versions from your R directory
#
#
# Return to this chat - the new files will be available in /mnt/project
#

# Finding Your Project Knowledge (Files)
# You'll find the project knowledge base on the right side of your project's main page. Click on the "+" button to add content to the project.
# Step-by-Step:
#
#   Go to your Project page:
#
#   Click on "Projects" in the left sidebar (or go to claude.ai/projects)
# Click on your ForestSearch project
#
#
# Look to the RIGHT side of the screen:
#
#   You should see a "Project Knowledge" or "Knowledge" panel
# Your 6 uploaded files should be listed there
#
#
# To remove old files:
#
#   Hover over a file name
# Click the trash/delete icon or "×" that appears
#
#
# To add updated files:
#
#   Click the "+" button in that panel
# Select your current R files from your computer
