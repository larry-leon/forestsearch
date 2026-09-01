# =============================================================================
# fastpath_verify_2026-09-02.R -- Stage C gates, routing proof, realized
# recovery (task: dev/tasks/cc_task_fastpath_unadjusted_2026-09-02.md sections
# 5-7).  Run with 0.3.3 installed, AFTER the postchange battery has written
# fastpath_postchange_2026-09-02.rds.
# Writes: fastpath_verify_2026-09-02.rds, two .Rprof files, and the log.
# =============================================================================
suppressPackageStartupMessages(library(forestsearch))
suppressPackageStartupMessages(library(survival))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
D <- "dev/glm-continuous-sims"
ns <- asNamespace("forestsearch")
cat(sprintf("installed forestsearch %s (expect 0.3.3) | git HEAD %s\n",
            as.character(utils::packageVersion("forestsearch")),
            system("git rev-parse --short HEAD", intern = TRUE)))
stopifnot(utils::packageVersion("forestsearch") == "0.3.3")

# ---- section 5: the equality gates ------------------------------------------
A <- readRDS(file.path(D, "fastpath_baseline_2026-09-02.rds"))
C <- readRDS(file.path(D, "fastpath_postchange_2026-09-02.rds"))
stopifnot(A$pkg_version == "0.3.2", C$pkg_version == "0.3.3")

first_diff <- function(a, b, path = "") {
  if (identical(a, b)) return(NULL)
  if (is.list(a) && is.list(b) && !is.data.frame(a) && !is.data.frame(b) &&
      identical(names(a), names(b)) && length(a) == length(b)) {
    for (i in seq_along(a)) {
      d <- first_diff(a[[i]], b[[i]], paste0(path, "$", if (!is.null(names(a))) names(a)[i] else i))
      if (!is.null(d)) return(d)
    }
    return(list(path = paste0(path, " [attributes differ]"), a = attributes(a), b = attributes(b)))
  }
  if (is.data.frame(a) && is.data.frame(b) && identical(names(a), names(b))) {
    for (cn in names(a)) if (!identical(a[[cn]], b[[cn]])) {
      j <- which(a[[cn]] != b[[cn]] | xor(is.na(a[[cn]]), is.na(b[[cn]])))[1]
      return(list(path = paste0(path, "$", cn, if (!is.na(j)) paste0("[", j, "]")),
                  a = if (!is.na(j)) a[[cn]][j] else a[[cn]], b = if (!is.na(j)) b[[cn]][j] else b[[cn]]))
    }
    return(list(path = paste0(path, " [frame attributes differ]"), a = attributes(a), b = attributes(b)))
  }
  list(path = path, a = a, b = b)
}

fail_any <- FALSE
eq_rows <- list()
compare_case <- function(nm) {
  a <- A$results[[nm]]; c_ <- C$results[[nm]]
  if (is.null(a$pruned) || is.null(c_$pruned)) {
    err_ok <- identical(a$error, c_$error)
    cat(sprintf("%-22s errored in both stages, messages identical: %s\n", nm, err_ok))
    eq_rows[[nm]] <<- data.frame(case = nm, component = "<error message>", identical = err_ok)
    if (!err_ok) fail_any <<- TRUE
    return(invisible())
  }
  comps <- union(names(a$pruned), names(c_$pruned))
  ok <- vapply(comps, function(k) identical(a$pruned[[k]], c_$pruned[[k]]), logical(1))
  extra <- c(warnings = identical(a$warnings, c_$warnings),
             counters = identical(a$counters, c_$counters),
             sg.harm = identical(a$sg.harm, c_$sg.harm))
  eq_rows[[nm]] <<- rbind(data.frame(case = nm, component = comps, identical = unname(ok)),
                          data.frame(case = nm, component = names(extra), identical = unname(extra)))
  cat(sprintf("%-22s components %d/%d identical; warnings %s; counters %s\n",
              nm, sum(ok), length(ok), extra[["warnings"]], extra[["counters"]]))
  bad <- c(comps[!ok], names(extra)[!extra])
  if (length(bad)) {
    fail_any <<- TRUE
    for (k in comps[!ok]) {
      d <- first_diff(a$pruned[[k]], c_$pruned[[k]], paste0(nm, "$", k))
      cat("  FIRST DIFF at ", d$path, "\n   Stage A: ", sep = "")
      print(d$a); cat("   Stage C: "); print(d$b)
    }
    if (!extra[["warnings"]]) { cat("  warnings A:\n"); print(a$warnings); cat("  warnings C:\n"); print(c_$warnings) }
    if (!extra[["counters"]]) { print(a$counters); print(c_$counters) }
  }
}
cat("== GATE E-1 (fast-path fixtures) ==\n"); for (nm in c("F1", "F2", "F3")) compare_case(nm)
cat("== GATE E-2 (routing fixtures) ==\n");  for (nm in c("F4", "F5", "F6")) compare_case(nm)
cat("== GATE E-3 (edge fixtures) ==\n");     for (nm in c("E1", "E2", "E3", "E4", "E5", "E6")) compare_case(nm)
cat("== GATE E-4 (bootstrap payload) ==\n")
boot_ok <- identical(A$boot$pruned, C$boot$pruned)
cat(sprintf("bootstrap pruned payload identical: %s\n", boot_ok))
if (!boot_ok) { fail_any <- TRUE; d <- first_diff(A$boot$pruned, C$boot$pruned, "boot")
                cat("  FIRST DIFF at ", d$path, "\n"); print(d$a); print(d$b) }
eq_matrix <- do.call(rbind, eq_rows)
if (fail_any) { saveRDS(list(eq_matrix = eq_matrix), file.path(D, "fastpath_verify_2026-09-02.rds"))
                stop("EQUALITY GATE FAILED -- see log; stopping before routing/profile.") }
cat("\nALL EQUALITY GATES PASS: every retained component identical() across stages.\n\n")

# ---- fixtures for sections 6-7 (F1/F2/F4/F5 args, verbatim blocks) ----------
actg_frame <- function() {
  actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat <- 1L - ifelse(actg_df$arms == 1L, 1L, 0L)
  actg_df$cd4_change <- actg_df$cd420 - actg_df$cd40
  actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$z1  <- as.factor(ifelse(actg_df$age > 34, 1L, 0L))
  actg_df$z2  <- as.factor(ifelse(actg_df$preanti <= 744.5, 1L, 0L))
  actg_df$z3  <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L))
  actg_df$z4  <- as.factor(ifelse(actg_df$karnof <= median(actg_df$karnof), 1L, 0L))
  actg_df$z5  <- as.factor(ifelse(actg_df$cd40 <= median(actg_df$cd40), 1L, 0L))
  actg_df$z6  <- as.factor(ifelse(actg_df$cd80 <= median(actg_df$cd80), 1L, 0L))
  actg_df$z7  <- as.factor(actg_df$hemo);  actg_df$z8  <- as.factor(actg_df$homo)
  actg_df$z9  <- as.factor(actg_df$drugs); actg_df$z10 <- as.factor(actg_df$race)
  actg_df$z11 <- as.factor(actg_df$gender);actg_df$z12 <- as.factor(actg_df$symptom)
  for (v in c("hemo", "homo", "drugs", "race", "gender", "symptom", "str2")) actg_df[[v]] <- as.factor(actg_df[[v]])
  actg_df
}
PAY <- "quarto/simulations/actg175/continuous/mr_md_harm/fs_maxeffCons_mr_md40_knoise0_n500_s1000_d5000/fs_maxeffCons_mr_md40_knoise0_n500_res_1_1000.rds"
truth <- readRDS(PAY)$truth
dgm <- generate_glm_dgm(
  data = actg_frame(), factor_vars = paste0("z", 1:12), outcome_var = "cd4_change",
  treatment_var = "treat", outcome_type = "continuous", effect_measure = "MD",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_treat = 1, k_inter = truth$beta_inter, n_super = 5000L,
  seed = 8316951L, verbose = FALSE)
SEED_BASE <- 8316951L; SIM_ID <- 1L
df_cont <- simulate_from_glm_dgm(dgm, n = 500L, seed = SEED_BASE + SIM_ID)
confs_cont <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80",
                "hemo", "homo", "drugs", "race", "gender", "symptom", "str2")
args_cont <- function(extra = list()) {
  base <- list(
    df.analysis = df_cont, confounders.name = confs_cont,
    outcome.name = "y_sim", treat.name = "treat_sim", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = 30, consistency.threshold = 10, pconsistency.threshold = 0.90,
    fs.splits = 400L, n.min = 60L, d0.min = 12L, d1.min = 12L, maxk = 2L,
    vi.grf.min = NULL, sg_focus = "maxeffCons", selection_rule = "neighborhood",
    effect_neighborhood = 0.10, stop_threshold = NULL, consistency_method = "resample",
    conf.cont_jcuts = list(age = 10, preanti = 10),
    use_lasso = FALSE, use_dina = FALSE, use_grf = FALSE, use_twostage = TRUE,
    is.RCT = TRUE, adverse_outcome = FALSE, details = FALSE, quiet = TRUE,
    seedit = SEED_BASE + SIM_ID, parallel_args = list(plan = "sequential", workers = 1L),
    mr_inference = FALSE)
  utils::modifyList(base, extra)
}
df_surv <- survival::gbsg
df_surv <- within(df_surv, { id <- seq_len(nrow(df_surv)); time_months <- rfstime / 30.4375; grade3 <- ifelse(grade == "3", 1, 0) })
args_surv <- function(extra = list()) {
  base <- list(
    df.analysis = df_surv, outcome.name = "time_months", event.name = "status",
    treat.name = "hormon", id.name = "id",
    confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
    is.RCT = TRUE, seedit = 8316951, est.scale = "hr",
    use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
    conf.cont_jcuts = list(er = 10, pgr = 10), collapse_cuts = TRUE,
    conf_force = c("er <= 0", "pgr <= 0"),
    subgroup_method = "consistency", sg_focus = "effMaxSG",
    selection_rule = "neighborhood", effect_neighborhood = 0.10,
    mr_inference = FALSE,
    n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
    hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.90,
    minp = 0.025, use_twostage = TRUE, stop_threshold = NULL,
    consistency_method = "resample", fs.splits = 1000, max.minutes = 3,
    vi.grf.min = NULL,
    show_candidate_summary = FALSE, details = FALSE, quiet = TRUE, plot.sg = FALSE,
    parallel_args = list(plan = "sequential", workers = 1L))
  utils::modifyList(base, extra)
}

# ---- section 6: routing proof -----------------------------------------------
.e <- new.env()
reset6 <- function() { .e$lm_n <- 0L; .e$lm_callers <- character(0)
                       .e$sc_n <- 0L; .e$sc_callers <- character(0)
                       .e$fast_n <- 0L; .e$slow_n <- 0L; .e$attr_seen <- NA }
caller_of <- function(k = -1) { cl <- sys.call(k); if (is.null(cl)) "<top>" else deparse(cl[[1]])[1] }
set6 <- function() {
  suppressMessages(trace("lm", where = asNamespace("stats"), print = FALSE,
    tracer = quote({ .e <- get(".e", envir = globalenv())
                     .e$lm_n <- .e$lm_n + 1L
                     cl <- sys.call(-1)
                     .e$lm_callers <- c(.e$lm_callers, if (is.null(cl)) "<top>" else deparse(cl[[1]])[1]) })))
  suppressMessages(trace("summary.coxph", where = asNamespace("survival"), print = FALSE,
    tracer = quote({ .e <- get(".e", envir = globalenv())
                     .e$sc_n <- .e$sc_n + 1L
                     cl <- sys.call(-1)
                     .e$sc_callers <- c(.e$sc_callers, if (is.null(cl)) "<top>" else deparse(cl[[1]])[1]) })))
  suppressMessages(trace("fit_glm_for_subgroup_fast", where = ns, print = FALSE,
    tracer = quote({ .e <- get(".e", envir = globalenv()); .e$fast_n <- .e$fast_n + 1L })))
  suppressMessages(trace("fit_glm_for_subgroup", where = ns, print = FALSE,
    tracer = quote({ .e <- get(".e", envir = globalenv()); .e$slow_n <- .e$slow_n + 1L })))
  suppressMessages(trace("subgroup.search", where = ns, print = FALSE,
    tracer = quote({ .e <- get(".e", envir = globalenv())
                     .e$attr_seen <- isTRUE(attr(estimator_fn, "fs_fast_unadjusted")) })))
}
clear6 <- function() for (f in c("lm", "summary.coxph")) suppressMessages(try(untrace(f,
  where = asNamespace(if (f == "lm") "stats" else "survival")), silent = TRUE))
clear6b <- function() for (f in c("fit_glm_for_subgroup_fast", "fit_glm_for_subgroup", "subgroup.search"))
  suppressMessages(try(untrace(f, where = ns), silent = TRUE))

assign(".e", .e, envir = globalenv())
routing <- list()
route_case <- function(nm, args) {
  reset6(); set6()
  fs <- suppressWarnings(do.call(forestsearch, args))
  clear6(); clear6b()
  lc <- sort(table(.e$lm_callers), decreasing = TRUE)
  sc <- sort(table(.e$sc_callers), decreasing = TRUE)
  cat(sprintf("%-8s attr_on_search %-5s fast_fit_calls %5d slow_fit_calls %5d stats::lm %5d summary.coxph %4d\n",
              nm, .e$attr_seen, .e$fast_n, .e$slow_n, .e$lm_n, .e$sc_n))
  if (length(lc)) { cat("  lm callers:            "); print(lc) }
  if (length(sc)) { cat("  summary.coxph callers: "); print(sc) }
  routing[[nm]] <<- list(attr_seen = .e$attr_seen, fast_n = .e$fast_n, slow_n = .e$slow_n,
                         lm_n = .e$lm_n, sc_n = .e$sc_n,
                         lm_callers = lc, sc_callers = sc, sg = fs$sg.harm)
}
route_case("F1", args_cont())
route_case("F2", args_surv())
route_case("F4", args_cont(list(adjust_covariates = c("age", "wtkg"))))
route_case("F5", args_surv(list(adjust_covariates = "age")))
route_case("F6", args_cont(list(ps_method = "logistic", ps_adjust_method = "iptw")))
# F3 anchor: predicate value only (attribute observation), no full re-run cost concern -- run it
cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
actg_oc <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
actg_oc$id <- seq_len(nrow(actg_oc)); actg_oc$treat <- ifelse(actg_oc$arms == 1L, 1L, 0L)
actg_oc <- actg_oc[!is.na(actg_oc$cd420), ]; actg_oc$y_decline <- actg_oc$cd40 - actg_oc$cd420
for (v in bin_vars) actg_oc[[v]] <- as.factor(actg_oc[[v]])
args_anchor <- list(df.analysis = actg_oc, confounders.name = c(cont_vars, bin_vars),
  outcome.name = "y_decline", treat.name = "treat", id.name = "id",
  outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE, seedit = 8316951,
  sg_focus = "maxeffCons", selection_rule = "neighborhood", consistency_method = "resample",
  effect.threshold = 10, consistency.threshold = 5, pconsistency.threshold = 0.90,
  use_twostage = TRUE, conf.cont_jcuts = list(age = 10, preanti = 10, wtkg = 10, karnof = 10, cd40 = 10, cd80 = 10),
  cut_type = "default", maxk = 2L, n.min = 60L, d0.min = 10L, d1.min = 10L, fs.splits = 500L,
  use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE, is.RCT = TRUE,
  parallel_args = list(plan = "sequential"), details = FALSE, quiet = TRUE)
route_case("F3", args_anchor)

# ---- section 7: realized recovery -------------------------------------------
BUCKETS <- list(
  medians      = c("survfit", "survfit.formula", "survfitKM", "summary.survfit", "survfit.Surv"),
  grf          = c("causal_forest", "causal_survival_forest", "fit_causal_forest_glm", "variable_importance",
                   "regression_forest", "fit_grf_variable_importance"),
  fit          = c("fit_glm_for_subgroup", "fit_cox_for_subgroup", "fit_glm_for_subgroup_fast"),
  cuts         = c("get_FSdata"),
  consistency  = c("evaluate_subgroup_consistency", "evaluate_consistency_twostage", "consistency_resample",
                   ".consistency_via_splits"),
  future_overhead = c("future_lapply", "future_xapply", "getGlobalsAndPackages", "findGlobals", "globalsOf",
                      "makeFuture", "run.Future", "resolve", "value", "result", "Future", "future", "sequential",
                      "SequentialFuture", "getGlobalsAndPackagesXApply", "packagesOf", "cleanup"),
  dedup_select = c("remove_near_duplicate_subgroups", ".maxeff_membership_dedup", "sort_subgroups", "sort_subgroups_preview"),
  enumeration  = c("subgroup.search", "generate_combination_indices", "evaluate_combination_with_status",
                   "meets_prevalence_threshold", "extract_idx_flagredundancy", "get_covs_in", "has_positive_variance",
                   "format_search_results", "get_subgroup_membership", "calculate_max_combinations"),
  consistency_other = c("subgroup.consistency"))
parse_rprof <- function(file) {
  ln <- readLines(file)
  hdr <- grepl("sample.interval=", ln, fixed = TRUE) | grepl("^#File ", ln)
  interval <- as.numeric(sub(".*sample.interval=", "", ln[grepl("sample.interval=", ln, fixed = TRUE)][1])) / 1e6
  samples <- ln[!hdr]
  stacks <- lapply(strsplit(samples, " ", fixed = TRUE), function(tok) {
    tok <- tok[nzchar(tok)]; tok <- tok[!grepl("^[0-9]+#[0-9]+$", tok)]; gsub('"', "", tok) })
  bucket_of <- function(st) { for (b in names(BUCKETS)) if (any(st %in% BUCKETS[[b]])) return(b); "everything_else" }
  bk <- vapply(stacks, bucket_of, character(1))
  tab <- table(factor(bk, levels = c(names(BUCKETS), "everything_else")))
  data.frame(bucket = names(tab), samples = as.integer(tab), secs = as.integer(tab) * interval,
             share = as.numeric(tab) / length(bk), stringsAsFactors = FALSE) -> out
  list(table = out, n_samples = length(bk), interval = interval, total_secs = length(bk) * interval)
}
prof_case <- function(nm, args) {
  invisible(suppressWarnings(do.call(forestsearch, args)))   # warm-up
  pf <- file.path(D, sprintf("fastpath_profile_%s_2026-09-02.Rprof", nm))
  gc(); Rprof(pf, interval = 0.010, line.profiling = TRUE, memory.profiling = FALSE)
  t0 <- proc.time()[[3]]
  fs <- suppressWarnings(do.call(forestsearch, args))
  wall <- proc.time()[[3]] - t0
  Rprof(NULL)
  pb <- parse_rprof(pf)
  cat(sprintf("\n%s profiled: %.2f s wall, %d samples\n", nm, wall, pb$n_samples))
  print(pb$table, digits = 4, row.names = FALSE)
  list(wall = wall, buckets = pb$table, rprof = pf)
}
prof <- list(F1 = prof_case("F1", args_cont()), F2 = prof_case("F2", args_surv()))
old <- readRDS(file.path(D, "bootstrap_reprofile_2026-09-01.rds"))
for (pair in list(c("F1", "cont"), c("F2", "surv"))) {
  nb <- prof[[pair[1]]]$buckets; ob <- old$results[[pair[2]]]$buckets
  m <- merge(nb, ob, by = "bucket", suffixes = c("_033", "_032"), sort = FALSE)
  m <- m[match(nb$bucket, m$bucket), c("bucket", "secs_032", "secs_033", "share_032", "share_033")]
  cat(sprintf("\n---- %s buckets: 0.3.2 (%.2f s) -> 0.3.3 (%.2f s)\n",
              pair[1], old$results[[pair[2]]]$wall_profiled, prof[[pair[1]]]$wall))
  print(m, digits = 4, row.names = FALSE)
}
cat(sprintf("\nStage-C F1 bootstrap: wall %.1f s, per replicate %.3f s (0.3.2: 5.02 s); B = 1000 -> %.1f min (0.3.2: 83.7 min)\n",
            C$boot$wall, C$boot$per_rep, 1000 * C$boot$per_rep / 60))

saveRDS(list(eq_matrix = eq_matrix, boot_ok = boot_ok, routing = routing,
             prof = prof, stageC_boot = list(wall = C$boot$wall, per_rep = C$boot$per_rep),
             pkg_version = as.character(utils::packageVersion("forestsearch")),
             built_at = Sys.time()),
        file.path(D, "fastpath_verify_2026-09-02.rds"))
cat("\nwritten:", file.path(D, "fastpath_verify_2026-09-02.rds"), "\n")
