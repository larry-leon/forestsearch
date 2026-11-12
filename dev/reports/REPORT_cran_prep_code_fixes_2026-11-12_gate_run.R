# gate_run.R -- TASK_cran_prep_code_fixes_2026-11-12, "After the fixes" step 1 (default-path gate).
# Usage: Rscript --vanilla gate_run.R <tree> <out_rds>
# The six equivalence cases of REPORT_cran_prep_seed_legacy_2026-11-11.md section 3.3 (script
# /tmp/seed_check_2026-11-11/equiv_run.R, step-5 block removed): load <tree> with
# pkgload::load_all(), set.seed(123) before each call, sequential plan, the package test
# helpers own synthetic data.
a <- commandArgs(TRUE); tree <- a[1]; out <- a[2]
t_all <- proc.time()[3]
suppressMessages(pkgload::load_all(tree, quiet = TRUE))
source(file.path(tree, "tests/testthat/helper-synthetic-dgm.R"))
future::plan("sequential")
grDevices::pdf(NULL)

seq_args <- list(plan = "sequential", workers = 1L, show_message = FALSE)
state <- function() if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(suppressMessages(expr)))
  out
}

# Data (built before any set.seed(123); the helpers seed themselves)
df_all  <- .make_survival_data(N = 600L, HR_harm = 2.5, seed = 42L)
df_test <- .make_survival_data(N = 300L, HR_harm = 2.5, seed = 43L)
df_bin  <- .make_binary_data(N = 600L, OR_harm = 3, seed = 44L)
surv_args <- list(df.analysis = df_all, confounders.name = c("age", "stage", "sex", "noise"),
                  outcome.name = "time", event.name = "event", treat.name = "treat",
                  id.name = "id", parallel_args = seq_args)
bin_args  <- list(df.analysis = df_bin,
                  confounders.name = c("age", "wtkg", "biomarker", "biomarker_hi"),
                  outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
                  outcome_type = "binary", parallel_args = seq_args)

calls <- list(
  c1_forestsearch_survival_defaults = function(env) do.call(forestsearch, surv_args),
  c2_forestsearch_mr_inference      = function(env) do.call(forestsearch, c(surv_args, list(mr_inference = TRUE))),
  c3_forestsearch_df_test           = function(env) do.call(forestsearch, c(surv_args, list(df.test = df_test))),
  c4_bootstrap_on_c1                = function(env) forestsearch_bootstrap_dofuture(env$c1_forestsearch_survival_defaults, nb_boots = 10L, parallel_args = seq_args),
  c5_kfold_cv_on_c1                 = function(env) forestsearch_Kfold(env$c1_forestsearch_survival_defaults, Kfolds = 5L, parallel_args = seq_args),
  c6_forestsearch_binary            = function(env) do.call(forestsearch, bin_args)
)
res <- list(); st <- list(); wall <- c()
for (nm in names(calls)) {
  set.seed(123)
  before <- state()
  t0 <- proc.time()[3]
  res[[nm]] <- quietly(calls[[nm]](res))
  wall[nm] <- proc.time()[3] - t0
  st[[nm]] <- list(before = before, after = state())
  cat(sprintf("%-36s %6.1f s  state restored: %s\n", nm, wall[nm], identical(before, state())))
}

saveRDS(list(tree = tree, res = res, state = st, wall = wall,
             versions = list(R = R.version.string, fs = as.character(packageVersion("forestsearch")))),
        out)
cat(sprintf("total %.1f s\n", proc.time()[3] - t_all))
