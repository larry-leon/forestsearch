#!/usr/bin/env Rscript
# =============================================================================
# §5A evidence probe -- stage 3
#   Q7b GLM outcomes under the DEFAULT consistency_method = "split"
#   Q9  gate family vs the identifier's admissible set
#   Q10 does the gate run inside bootstrap replicates / CV folds?
# READ-ONLY with respect to the package.  Writes only under dev/terminology-work/.
# =============================================================================

suppressMessages(library(forestsearch))
source("tests/testthat/helper-synthetic-dgm.R")
dir.create("dev/terminology-work/out", showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a
hr <- function(s) cat("\n===== ", s, " =====\n", sep = "")
ns <- asNamespace("forestsearch")

DGA  <- list(draws = 200L, seed = 11L)
CONF <- c("age", "stage", "sex", "noise")
d    <- .make_survival_data(N = 400L, HR_harm = 3.0)
mk <- function(...) .fs_args_for("survival", confounders = CONF,
  extra = utils::modifyList(
    list(use_grf = FALSE, use_lasso = FALSE, quiet = TRUE), list(...)))

# ---------------------------------------------------------------------------
hr("Q7b: GLM outcomes -- gate under consistency_method 'split' (DEFAULT) vs 'resample'")
cat("consistency_method default =",
    deparse(formals(forestsearch)$consistency_method)[1], "-> \"split\"\n\n")
specs <- list(
  binary     = list(df = .make_binary_data(N = 400L, OR_harm = 3.0),
                    conf = c("age", "biomarker_hi")),
  continuous = list(df = .make_continuous_data(N = 400L, MD_harm = 2.0),
                    conf = c("age", "biomarker_hi")),
  count      = list(df = .make_count_data(N = 400L, IRR_harm = 3.0),
                    conf = c("age", "biomarker_hi"))
)
rows <- list()
for (ot in names(specs)) {
  for (cm in c("split", "resample")) {
    sp <- specs[[ot]]
    aa <- .fs_args_for(ot, confounders = sp$conf,
      extra = list(use_grf = FALSE, use_lasso = FALSE, quiet = TRUE,
                   consistency_method = cm,
                   debias_gate = TRUE, debias_gate_args = DGA))
    rr <- tryCatch(.run_fs_capture(sp$df, aa),
      error = function(e) list(value = NULL,
        warnings = paste("ERROR:", conditionMessage(e))))
    fb <- rr$value
    rows[[length(rows) + 1L]] <- data.frame(
      outcome = ot, consistency_method = cm,
      found    = !is.null(fb$sg.harm) && length(fb$sg.harm) > 0,
      gate_ran = !is.null(fb$debias_gate),
      hfd      = format(fb$harm_flag_debiased %||% NA),
      warned_about_gate = any(grepl("debias|gate", rr$warnings, ignore.case = TRUE)),
      stringsAsFactors = FALSE)
  }
}
glm_tab <- do.call(rbind, rows); print(glm_tab, right = FALSE)
saveRDS(glm_tab, "dev/terminology-work/out/glm_split_vs_resample.rds")

# ---------------------------------------------------------------------------
hr("Q9: gate re-selection family vs the identifier's admissible set")
a9 <- utils::modifyList(mk(sg_focus = "maxSG"),
                        list(debias_gate = TRUE, debias_gate_args = DGA))
r9  <- .run_fs_capture(d, a9); fs9 <- r9$value
cat("gate n_family                       :", fs9$debias_gate$n_family, "\n")
cat("search filter_counts:\n"); print(unlist(fs9$find.grps$filter_counts))
cat("\nmax_subgroups_search (formal default):",
    deparse(formals(forestsearch)$max_subgroups_search), "| used here: 5\n")
cat("gate family applies (main.R:2839-2843): maxk and n.min ONLY",
    "-- no d0.min/d1.min, no max_subgroups_search\n")
cat("search applies (subgroup_search.R:591-601, 639): n.min AND d0.min/d1.min\n")

# ---------------------------------------------------------------------------
hr("Q10: fs_debias_gate() invocation count inside a 2-replicate bootstrap")
fsb <- fs9
cat("original analysis gate ran :", !is.null(fsb$debias_gate), "\n")
cat("args_call_all$debias_gate  :", format(fsb$args_call_all$debias_gate), "\n")
cnt <- new.env(parent = emptyenv()); cnt$n <- 0L
ok <- tryCatch({
  suppressMessages(trace("fs_debias_gate", where = ns, print = FALSE,
                         tracer = function() cnt$n <- cnt$n + 1L)); TRUE
}, error = function(e) { cat("trace failed:", conditionMessage(e), "\n"); FALSE })
if (ok) {
  cnt$n <- 0L
  bb <- tryCatch(suppressWarnings(suppressMessages(
    forestsearch_bootstrap_dofuture(fsb, nb_boots = 2L,
      parallel_args = list(plan = "sequential", workers = 1L,
                           show_message = FALSE)))),
    error = function(e) { cat("bootstrap ERROR:", conditionMessage(e), "\n"); NULL })
  n_boot <- cnt$n
  cat("calls during bootstrap(nb_boots = 2) :", n_boot, "\n")

  cnt$n <- 0L
  cv <- tryCatch(suppressWarnings(suppressMessages(
    forestsearch_Kfold(fsb, n_folds = 2L, n_repeats = 1L,
      parallel_args = list(plan = "sequential", workers = 1L,
                           show_message = FALSE)))),
    error = function(e) { cat("CV note:", conditionMessage(e), "\n"); NULL })
  cat("calls during Kfold(n_folds = 2)      :", cnt$n, "\n")
  suppressMessages(try(untrace("fs_debias_gate", where = ns), silent = TRUE))
  cat("\n0 = confined to the original analysis; >0 = also recomputed per replicate/fold\n")
}
cat("\n[stage 3 done]\n")
