#!/usr/bin/env Rscript
# =============================================================================
# §5A evidence probe -- stage 4: CV-fold recomputation (correct Kfold signature)
# and the §5A.3.1 no-downstream-branching check.
# READ-ONLY with respect to the package.
# =============================================================================

suppressMessages(library(forestsearch))
source("tests/testthat/helper-synthetic-dgm.R")

hr <- function(s) cat("\n===== ", s, " =====\n", sep = "")
ns <- asNamespace("forestsearch")

CONF <- c("age", "stage", "sex", "noise")
d <- .make_survival_data(N = 400L, HR_harm = 3.0)
a <- .fs_args_for("survival", confounders = CONF,
  extra = list(use_grf = FALSE, use_lasso = FALSE, quiet = TRUE,
               sg_focus = "maxSG", debias_gate = TRUE,
               debias_gate_args = list(draws = 200L, seed = 11L)))
fs <- .run_fs_capture(d, a)$value
cat("original gate ran:", !is.null(fs$debias_gate), "\n")

hr("Q10b: fs_debias_gate() calls inside forestsearch_Kfold(Kfolds = 2)")
cnt <- new.env(parent = emptyenv()); cnt$n <- 0L
suppressMessages(trace("fs_debias_gate", where = ns, print = FALSE,
                       tracer = function() cnt$n <- cnt$n + 1L))
cnt$n <- 0L
cv <- tryCatch(suppressWarnings(suppressMessages(
  forestsearch_Kfold(fs, Kfolds = 2L,
    parallel_args = list(plan = "sequential", workers = 1L,
                         show_message = FALSE)))),
  error = function(e) { cat("Kfold ERROR:", conditionMessage(e), "\n"); NULL })
cat("calls during Kfold(Kfolds = 2):", cnt$n, "\n")
suppressMessages(try(untrace("fs_debias_gate", where = ns), silent = TRUE))

hr("Q11 (§5A.3.1): is the reported subgroup identical with and without the gate?")
a_off <- utils::modifyList(a, list(debias_gate = FALSE))
fs_off <- .run_fs_capture(d, a_off)$value
cat("sg.harm  gate ON :", paste(fs$sg.harm,     collapse = " & "), "\n")
cat("sg.harm  gate OFF:", paste(fs_off$sg.harm, collapse = " & "), "\n")
cat("identical sg.harm      :", identical(fs$sg.harm, fs_off$sg.harm), "\n")
cat("identical sg.harm.id   :",
    identical(fs$grp.consistency$sg.harm.id, fs_off$grp.consistency$sg.harm.id), "\n")
cat("identical max_sg_est   :", identical(fs$max_sg_est, fs_off$max_sg_est), "\n")
cat("identical df.est dims  :",
    identical(dim(fs$df.est), dim(fs_off$df.est)), "\n")

hr("Q12: fields present ONLY when the gate ran (return-shape delta)")
n_on  <- names(fs); n_off <- names(fs_off)
cat("names identical:", identical(n_on, n_off), "\n")
cat("gate-only fields:", paste(setdiff(n_on, n_off), collapse = ", "), "\n")
cat("debias_gate (OFF) is NULL   :", is.null(fs_off$debias_gate), "\n")
cat("harm_flag_debiased (OFF)    :", format(fs_off$harm_flag_debiased), "\n")

hr("Q13: include_complement -- formal default vs the value callers actually pass")
cat("fs_debias_gate() formal default   :",
    deparse(formals(get("fs_debias_gate", envir = ns))$include_complement), "\n")
cat(".fs_apply_debias_gate passes      : TRUE",
    "(fs_debias_gate_methods.R:187)\n")
cat("consistency branch passes         : TRUE",
    "(forestsearch_main.R:2884)\n")
cat("observed complement present       :", !is.null(fs$debias_gate$complement), "\n")
cat("\n[stage 4 done]\n")
