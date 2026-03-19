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


# Notes from John B.
# I agree with many of its decisions and fixes
# * I am concerned about the descriptive stubs. These are `@examples` sections
# which are comments-only. I suspect these will be flagged by the CRAN reviewer
# as insufficient. Do these need to be exported functions? For example,
# `extract_selected_tree_cuts()` seems like an internal-only function to me
# * You committed the package tarball `forestsearch_0.1.0.tar.gz`. Please remove
# this file and add the pattern `*.tar.gz` to your `.gitignore`
# * It notes that there are 4 example sections that still use `\dontrun{}`. I see
# many more when I search for it. This cannot be delayed for future releases as
# Claude suggests. You need to address all of these before resubmitting

# If the example is short and runs quickly, there is no harm in keeping it.
# But it doesn't make much sense to include complex, long-running examples that
# have to be wrapped in \dontrun{} if an end user is unlikely to ever run the example code anways
#
# In the branch cran-review-1, {forestsearch} currently has 118 exported functions in NAMESPACE and
# 255 manual pages in man/. However, only 42 out of the 118 are actually used in the vignettes/articles
# under vignettes/. This suggests that only those 42 functions are required to be used directly
# by an end user performing an analysis. I recommend focusing your attention on writing thorough
# documentation for only these 42 important functions.
# All the other functions can be internal to the package and only receive a bare minimum of documentation

# Four examples in R/oc_analyses_gbsg.R still use \dontrun{}. These call create_gbsg_dgm() without suppressWarnings().
# If CRAN asks to convert these to \donttest{}, they will need suppressWarnings() wrappers


# I am preparing an R package called forestsearch for CRAN submission. The CRAN reviewer raised 8 issues.
# I have made corrections and would like you to independently verify each fix. I will share the corrected files grouped by fix category.
# For each group, please confirm the fix is correctly applied and flag anything that looks wrong or incomplete.


I am preparing an R package called forestsearch for CRAN submission.
Please fetch the most recent forestsearch codebase https://github.com/larry-leon/forestsearch.
An extensive review has been conducted by Claude, however an independent reviewer has raised the following concerns:
(A) I am concerned about the descriptive stubs. These are `@examples` sections which are comments-only. I suspect these will
be flagged by the CRAN reviewer  as insufficient. Do these need to be exported functions? For example,
`extract_selected_tree_cuts()` seems like an internal-only function

