# vi_common.R -- shared setup for the vi.grf.min / max_n_confounders
# deprecation verification (task dev/tasks/cc_task_vi_deprecation_2026-10-26.md).
# Stage A (stageA.R) runs against the installed CRAN 0.4.0 build; Stage C
# (stageC.R) against the installed branch build.  Both source this file.
#
# The three fixtures are transplanted verbatim from the committed Stage-0b
# bundle dev/glm-continuous-sims/stage0b_glm_2026-10-26/: fixture_survival()
# from stage0b_common.R, fixture_continuous.R and fixture_binary.R whole.

suppressPackageStartupMessages({
  library(forestsearch)
  library(survival)
  library(future)
  library(doFuture)
  library(data.table)
})

OUT <- Sys.getenv("VI_OUT", ".")

# ---- fixtures (verbatim transplant) -----------------------------------------
fixture_survival <- function() {
  df.analysis <- survival::gbsg
  df.analysis <- within(df.analysis, {
    id          <- seq_len(nrow(df.analysis))
    time_months <- rfstime / 30.4375
    grade3      <- ifelse(grade == "3", 1, 0)
  })
  confounders.name <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
  args <- list(
    df.analysis = df.analysis,
    outcome.name = "time_months", event.name = "status",
    treat.name = "hormon", id.name = "id",
    confounders.name = confounders.name,
    is.RCT = TRUE, seedit = 8316951, est.scale = "hr",
    use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
    conf.cont_jcuts = list(er = 10, pgr = 10), collapse_cuts = TRUE,
    conf_force = c("er <= 0", "pgr <= 0"),
    subgroup_method = "consistency", sg_focus = "effMaxSG",
    selection_rule = "neighborhood", effect_neighborhood = 0.20,
    mr_inference = FALSE,
    n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
    hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.80,
    fs.splits = 200, max.minutes = 3,
    details = FALSE, quiet = TRUE, plot.sg = FALSE,
    parallel_args = list(plan = "sequential")
  )
  list(label = "survival: vignette gbsg config (Stage-0 fixture V)", df = df.analysis, args = args,
       source = "inst/doc/survival-post-selection.qmd (CRAN 0.4.0), mr_inference = FALSE")
}

# ---- fixture_continuous.R (verbatim) ----
# fixture_continuous.R -- the continuous (MD) fixture for Stage 0b.
#
# TRANSPLANTED (CC, 2026-10-26) from the committed ACTG175 continuous
# application: quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd,
# chunks `data-prep` (data preparation, verbatim) and `anchor` (the
# forestsearch() parameter block, verbatim), with exactly these overrides for
# the per-replicate shape:
#   mr_inference = FALSE (mr_inference_args dropped), details = FALSE,
#   quiet = TRUE, plot.sg = FALSE, plot.grf = FALSE,
#   parallel_args = list(plan = "sequential").
# seedit = params$seed is resolved to the document's YAML params value 8316951.
# use_lasso / use_grf / use_twostage / consistency_method and every threshold
# are kept as the document has them.

fixture_continuous <- function() {
  if (!requireNamespace("speff2trial", quietly = TRUE))
    stop("speff2trial is required for the ACTG175 fixtures")
  params <- list(seed = 8316951)

  # ---- data-prep (verbatim) ----
  cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
  bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
  confounders.name <- c(cont_vars, bin_vars)

  actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat <- ifelse(actg_df$arms == 1L, 1L, 0L)
  actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$y_decline <- actg_df$cd40 - actg_df$cd420
  for (v in bin_vars) actg_df[[v]] <- as.factor(actg_df[[v]])

  N <- nrow(actg_df)
  stopifnot(N == 1083L)

  # ---- anchor (verbatim parameter block, overrides marked) ----
  args <- list(
    df.analysis      = actg_df,
    confounders.name = confounders.name,
    outcome.name     = "y_decline",
    treat.name       = "treat",
    id.name          = "id",
    outcome_type     = "continuous",
    effect_measure   = "MD",
    adverse_outcome  = TRUE,
    seedit           = params$seed,
    sg_focus         = "effMaxSG",
    effect_neighborhood = 0.20,
    selection_rule   = "neighborhood",
    consistency_method = "resample",
    effect.threshold       = 10,
    consistency.threshold  = 5,
    pconsistency.threshold = 0.90,
    use_twostage     = TRUE,
    conf.cont_jcuts  = list(age = 10, preanti = 10, wtkg = 10,
                            karnof = 10, cd40 = 10, cd80 = 10),
    cut_type         = "default",
    maxk             = 2L,
    n.min            = 60L,
    d0.min           = 10L,
    d1.min           = 10L,
    fs.splits        = 500L,
    use_lasso        = FALSE,
    use_grf          = FALSE,
    use_dina         = FALSE,
    is.RCT           = TRUE,
    parallel_args    = list(plan = "sequential"),   # override (document: same)
    details          = FALSE,                       # override (document: same)
    quiet            = TRUE,                        # override (document: same)
    mr_inference     = FALSE,                       # override (document: TRUE; mr_inference_args dropped)
    plot.sg          = FALSE,                       # override (document: absent)
    plot.grf         = FALSE                        # override (document: absent)
  )
  list(label = "continuous: ACTG175 MD, committed application (effMaxSG, eps = 0.20)",
       df = actg_df, args = args,
       source = "quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd, chunks data-prep + anchor")
}

# ---- fixture_binary.R (verbatim) ----
# fixture_binary.R -- the binary (OR) fixture for Stage 0b.
#
# TRANSPLANTED (CC, 2026-10-26) from the committed ACTG175 binary application:
# quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd,
# chunks `parallel-setup` (the fs_* settings the call references),
# `analysis-roles`, `data-treat`, `data-outcomes`, `data-covariates` (data
# preparation, verbatim) and `forestsearch-main` (the forestsearch() parameter
# block, verbatim), with exactly these overrides for the per-replicate shape:
#   mr_inference = FALSE (mr_inference_args dropped), details = FALSE,
#   quiet = TRUE, plot.sg = FALSE, plot.grf = FALSE,
#   parallel_args = list(plan = "sequential").
# use_lasso / use_grf / vi.grf.min / use_twostage / consistency_method and
# every threshold are kept as the document has them.

fixture_binary <- function() {
  if (!requireNamespace("speff2trial", quietly = TRUE))
    stop("speff2trial is required for the ACTG175 fixtures")

  # ---- parallel-setup: settings referenced by the call (verbatim) ----
  fs_sg_focus              <- "effMaxSG"   # GLM-natural alias of "hrMaxSG"
  fs_selection_rule        <- "neighborhood"       # neighborhood ∩ Pareto (effMaxSG/effMinSG)
  fs_or_threshold          <- 1.0
  fs_or_consistency        <- 1.0
  fs_pconsistency          <- 0.90
  fs_consistency_method    <- "resample"
  fs_subgroup_method <- "consistency"  # "consistency" | "dina" | "grf"
  fs_grf_selection   <- "frontier"     # only if subgroup_method == "grf": "tree" | "frontier"
  fs_grf_select_statistic <- "effect"  # only if subgroup_method=="grf" & grf_selection=="frontier": "dr" | "effect".
  fs_grf_depth       <- 2L             # GRF depth (screening + GRF selection)
  fs_dmin_grf        <- 0              # GRF DR-score screening floor (RD scale for binary; 0 = permissive)
  fs_dina_args       <- list()         # extra args forwarded to DINA selection
  fs_dina_select_statistic <- "effect" # only if subgroup_method=="dina": "dina" | "effect".

  # ---- analysis-roles (verbatim) ----
  analysis_seed <- 8316951L
  set.seed(analysis_seed)
  cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
  bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
  confounders.name <- c(cont_vars, bin_vars)
  id.name   <- "id"
  treat.name   <- "treat"
  outcome.name <- "y_binary_neg"   # FS adverse-scale outcome (1 = no improvement)
  fs_outcome_type   <- "binary"
  fs_effect_measure <- "OR"

  # ---- data-treat / data-outcomes / data-covariates (verbatim) ----
  actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat <- ifelse(actg_df$arms == 1L, 1L, 0L)
  actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$y_binary_pos <- as.integer(actg_df$cd420 > actg_df$cd40)
  actg_df$y_binary_neg <- 1L - actg_df$y_binary_pos
  for (v in bin_vars) actg_df[[v]] <- as.numeric(actg_df[[v]])

  # ---- forestsearch-main (verbatim parameter block, overrides marked) ----
  args <- list(
    df.analysis               = actg_df,
    confounders.name          = confounders.name,
    outcome.name              = outcome.name,
    event.name                = outcome.name,    # binary: outcome IS the event
    treat.name                = treat.name,
    id.name                   = id.name,
    is.RCT                    = TRUE,
    seedit                    = analysis_seed,
    outcome_type              = fs_outcome_type,
    effect_measure            = fs_effect_measure,
    adverse_outcome           = TRUE,
    use_lasso                 = FALSE,
    use_grf                   = FALSE,
    grf_depth                 = fs_grf_depth,
    dmin.grf                  = fs_dmin_grf,
    frac.tau                  = 0.80,
    vi.grf.min                = -0.2,
    return_selected_cuts_only = TRUE,
    tune_grf                  = FALSE,
    use_dina                  = FALSE,
    dina_args                 = fs_dina_args,
    cut_type                  = "default",
    cont.cutoff               = 4L,
    conf.cont_jcuts           = list(cd40 = 10, wtkg = 10),
    subgroup_method           = fs_subgroup_method,        # consistency | dina | grf
    grf_selection             = fs_grf_selection,          # used when subgroup_method = "grf"
    grf_select_statistic      = fs_grf_select_statistic,   # frontier only
    dina_select_statistic     = fs_dina_select_statistic,  # used when subgroup_method = "dina"
    sg_focus                  = fs_sg_focus,
    selection_rule            = fs_selection_rule,
    effect_neighborhood       = 0.10,
    mr_inference              = FALSE,                     # override (document: run_mr = TRUE; mr_inference_args dropped)
    n.min                     = 60L,
    d0.min                    = 10L,
    d1.min                    = 10L,
    maxk                      = 2L,
    effect.threshold          = fs_or_threshold,
    consistency.threshold     = fs_or_consistency,
    pconsistency.threshold    = fs_pconsistency,
    consistency_method        = fs_consistency_method,   # "resample" (gate requires it)
    fs.splits                 = 500L,
    use_twostage              = TRUE,
    show_candidate_summary    = TRUE,
    details                   = FALSE,                     # override (document: TRUE)
    quiet                     = TRUE,                      # override (document: FALSE)
    plot.sg                   = FALSE,   # KM not applicable; use plot_sg_glm_outcomes()
    plot.grf                  = FALSE,                     # override (document: absent)
    parallel_args             = list(plan = "sequential")  # override (document: multisession, n_cores)
  )
  list(label = "binary: ACTG175 OR, committed application (fixed family, effMaxSG, eps = 0.10)",
       df = actg_df, args = args,
       source = "quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd, chunks parallel-setup + analysis-roles + data-* + forestsearch-main")
}

FIXTURES <- list(survival = fixture_survival, continuous = fixture_continuous, binary = fixture_binary)

# ---- forest-fit counter and runner (verbatim from the task document, §2) ----
# Counts grf forest fits made anywhere inside expr, via base::trace() on the
# grf namespace (forestsearch calls grf::causal_forest / grf::causal_survival_forest
# through `::`, which resolves in that namespace).  Session-only; untraced on exit.
count_grf_fits <- function(expr) {
  cnt <- new.env(); cnt$n <- 0L
  # The tracer is evaluated in the traced function's frame, so it must carry
  # the counter environment inside the call (a closure named here is not found).
  tracer <- bquote(assign("n", get("n", envir = .(cnt)) + 1L, envir = .(cnt)))
  fns <- c("causal_forest", "causal_survival_forest")
  suppressMessages(for (f in fns) trace(f, tracer = tracer, where = asNamespace("grf"), print = FALSE))
  on.exit(suppressMessages(for (f in fns) untrace(f, where = asNamespace("grf"))), add = TRUE)
  value <- force(expr)
  list(value = value, n_fits = cnt$n)
}
# Runs one forestsearch() call, collecting warnings without suppressing the result.
run_collect <- function(args) {
  w <- character(0)
  r <- count_grf_fits(withCallingHandlers(
    suppressMessages(do.call(forestsearch::forestsearch, args)),
    warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") }))
  list(fs = r$value, n_fits = r$n_fits, warnings = w)
}

# ---- run variants -------------------------------------------------------------
# (i): args as committed.  (ii): vi.grf.min = NULL, max_n_confounders absent.
variant_args <- function(args, run) {
  if (run == "ii") {
    args$max_n_confounders <- NULL
    args["vi.grf.min"] <- list(NULL)
  }
  args
}

# ---- save / compare helpers ---------------------------------------------------
# Identity components (args_call_all is stored but excluded from identity).
IDENTITY_COMPONENTS <- c("sg.harm", "membership", "out_sg_result", "hr.subgroups",
                         "filter_counts", "n_candidates")

extract_run <- function(rc, id.name) {
  fs <- rc$fs
  memb <- NULL
  if (!is.null(fs$df.est)) {
    o <- order(fs$df.est[[id.name]])
    memb <- data.frame(id = fs$df.est[[id.name]][o],
                       treat.recommend = fs$df.est$treat.recommend[o])
  }
  gc <- fs$grp.consistency
  nc <- if (!is.null(gc)) gc[grep("^n_candidates_", names(gc))] else NULL
  list(
    sg.harm       = fs$sg.harm,
    membership    = memb,
    out_sg_result = gc$out_sg$result,
    hr.subgroups  = fs$find.grps$out.found$hr.subgroups,
    filter_counts = fs$find.grps$filter_counts,
    n_candidates  = nc,
    n_fits        = rc$n_fits,
    n_warnings    = length(rc$warnings),
    warnings      = rc$warnings,
    args_call_all = fs$args_call_all
  )
}

BOOT_COMPONENTS <- c("results", "SG_CIs", "H_estimates", "Hc_estimates", "Ystar_mat", "FSsg_tab")

extract_boot <- function(bt) {
  res <- bt$results
  if (!is.null(res)) {
    drop <- intersect(c("tmins_search", "tmins_iteration"), names(res))
    if (length(drop)) {
      if (data.table::is.data.table(res)) res <- res[, setdiff(names(res), drop), with = FALSE]
      else res <- res[, setdiff(names(res), drop), drop = FALSE]
    }
    attr(res, "timing") <- NULL
  }
  found <- intersect(BOOT_COMPONENTS, names(bt))
  out <- lapply(stats::setNames(found, found), function(k) if (k == "results") res else bt[[k]])
  out$names_found <- found
  out
}

run_boot <- function(fs) {
  bt <- suppressMessages(forestsearch_bootstrap_dofuture(
    fs.est = fs, nb_boots = 20L, seed = 8316951L,
    parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE)))
  future::plan("sequential")
  bt
}

out_file <- function(stage, fixture, run) file.path(OUT, sprintf("stage%s_%s_%s.rds", stage, fixture, run))

summary_row <- function(fixture, run, x) {
  data.frame(fixture = fixture, run = run, n_fits = x$n_fits, n_warnings = x$n_warnings,
             sg.harm = paste(x$sg.harm, collapse = " & "), stringsAsFactors = FALSE)
}
