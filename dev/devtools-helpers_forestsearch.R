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

pkgs <- c("data.table", "foreach", "future", "doFuture", "future.apply", "glmnet", "gt", "randomForest", "stringr", "survival", "grf", "ggplot2", "policytree", "DiagrammeR",
          "future.callr")
for (pkg in pkgs) {
  desc::desc_set_dep(pkg, file = "DESCRIPTION")
}

#desc::desc_set_dep("weightedSurv", type = "Suggests", file = "DESCRIPTION")

desc::desc_set_dep(package = "weightedsurv", type = "Imports")

desc::desc_set_remotes("larry-leon/weightedsurv")

usethis::use_package("patchwork", type = "Suggests")


#usethis::use_package("ggplot2")

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

# Load your updated package
devtools::load_all()

# Generate documentation
devtools::document()

# Run tests
devtools::test()

# Check package
devtools::check()

