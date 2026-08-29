# =============================================================================
# vi_grf_smoke_2026-08-30.R -- vi.grf.min = -0.2 (order only) vs NULL (Section 5
# skipped): equivalence of results, and whether the ordering buys efficiency
# -----------------------------------------------------------------------------
# Task: dev/tasks/cc_task_vi_grf_smoke_2026-08-30.md.  Read-only with respect to
# the package: no R/ change, no default change.  Synthetic data only (the
# testthat helpers' generators plus fixture-specific variants below).
#
# For each fixture F1..F8 and each of SEEDS seeds, forestsearch() is run twice on
# the same data with the same seedit -- vi.grf.min = -0.2 and vi.grf.min = NULL
# -- and the pair is compared on SUBJECT MEMBERSHIP (never the definition
# string), on the candidate family as a set of rows keyed by the selected cut
# columns (the q-labels are stable identities across the two runs; only their
# column order moves), on the consistency-stage counters, and on wall-clock.
# Section 5's own cost is measured two ways: the paired wall-clock difference,
# and a direct timing of the same grf fit on the same cut matrix.
#
# Writes: dev/glm-continuous-sims/vi_grf_smoke_2026-08-30.rds (+ .log via Rscript)
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-synthetic-dgm.R")
OUT   <- "dev/glm-continuous-sims/vi_grf_smoke_2026-08-30.rds"
SEEDS <- 1:20
N     <- 500L

cat("forestsearch", as.character(utils::packageVersion("forestsearch")),
    "| formal vi.grf.min default:", deparse(formals(forestsearch)$vi.grf.min), "\n")

# ---- fixture generators -------------------------------------------------------
# F1: plain continuous DGM, 8 covariates, one real subgroup (age <= 50 & bm_hi == 1)
gen_F1 <- function(seed, n = N) {
  set.seed(seed)
  df <- data.frame(id = seq_len(n), treat = rbinom(n, 1L, 0.5),
                   age = round(rnorm(n, 55, 12)), wtkg = round(rnorm(n, 75, 12), 1),
                   biomarker = rnorm(n), cd4 = round(rnorm(n, 350, 100)),
                   x5 = rnorm(n), x6 = round(runif(n, 0, 10)),
                   bm_hi = factor(rbinom(n, 1L, 0.5), levels = 0:1),
                   sex = factor(rbinom(n, 1L, 0.5), levels = 0:1))
  inh <- df$age <= 50 & df$bm_hi == "1"
  df$y <- rnorm(n, mean = df$treat * ifelse(inh, 1.5, 0), sd = 1)
  list(df = df, conf = c("age", "wtkg", "biomarker", "cd4", "x5", "x6", "bm_hi", "sex"),
       type = "continuous")
}
# F2: F1 plus a binary whose {rare == 1} cut has prevalence 496/500 = 0.992,
#     just below (n - rmin)/n = 0.99 ... i.e. adding it FIRST removes only 4 <= rmin
#     subjects, so extract_idx_flagredundancy() flags it only when it leads.
gen_F2 <- function(seed, n = N) {
  g <- gen_F1(seed, n); set.seed(seed + 1000L)
  rare <- rep(1L, n); rare[sample.int(n, 4L)] <- 0L
  g$df$rare <- factor(rare, levels = 0:1)
  # give `rare` a real interaction so it can be selected / matter for ordering
  g$df$y <- g$df$y + g$df$treat * ifelse(g$df$rare == "1", 0.4, 0)
  g$conf <- c(g$conf, "rare"); g
}
# F3: two covariates with identical membership (x_dup == bm_hi) and one near-
#     identical (x_near differs from bm_hi in one subject)
gen_F3 <- function(seed, n = N) {
  g <- gen_F1(seed, n); set.seed(seed + 2000L)
  g$df$x_dup  <- g$df$bm_hi
  nr <- as.integer(as.character(g$df$bm_hi)); i <- sample.int(n, 1L); nr[i] <- 1L - nr[i]
  g$df$x_near <- factor(nr, levels = 0:1)
  g$conf <- c(g$conf, "x_dup", "x_near"); g
}
# F4: 10 covariates with max_n_confounders = 5 (the positive control)
gen_F4 <- function(seed, n = N) {
  g <- gen_F1(seed, n); set.seed(seed + 3000L)
  g$df$x7 <- rnorm(n); g$df$x8 <- factor(rbinom(n, 1L, 0.4), levels = 0:1)
  g$conf <- c(g$conf, "x7", "x8"); g$extra <- list(max_n_confounders = 5L); g
}
# F5: all covariates pure noise (no interaction anywhere)
gen_F5 <- function(seed, n = N) {
  g <- gen_F1(seed, n); set.seed(seed + 4000L)
  g$df$y <- rnorm(n, mean = 0.2 * g$df$treat, sd = 1); g
}
# F6: survival with zero events in the treatment arm (the survival VI guard)
gen_F6 <- function(seed, n = N) {
  df <- .make_survival_data(N = n, HR_harm = 2.0, seed = seed)
  df$event[df$treat == 1L] <- 0L
  list(df = df, conf = c("age", "stage", "sex", "noise"), type = "survival")
}
# F7: survival, normal
gen_F7 <- function(seed, n = N) {
  df <- .make_survival_data(N = n, HR_harm = 2.0, seed = seed)
  list(df = df, conf = c("age", "stage", "sex", "noise"), type = "survival")
}
# F8: binary outcome
gen_F8 <- function(seed, n = N) {
  df <- .make_binary_data(N = n, OR_harm = 2.5, seed = seed)
  list(df = df, conf = c("age", "wtkg", "biomarker", "biomarker_hi", "sex"), type = "binary")
}
GEN <- list(F1 = gen_F1, F2 = gen_F2, F3 = gen_F3, F4 = gen_F4, F5 = gen_F5,
            F6 = gen_F6, F7 = gen_F7, F8 = gen_F8)

# ---- the forestsearch() call, per outcome type --------------------------------
fs_call <- function(g, vi, seed, sg_focus = "maxeffCons", extra = list(), pcons = 0.5) {
  base <- list(df.analysis = g$df, confounders.name = g$conf, treat.name = "treat", id.name = "id",
               vi.grf.min = vi, sg_focus = sg_focus, n.min = 30L, maxk = 2L, fs.splits = 30L,
               pconsistency.threshold = pcons, use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
               is.RCT = TRUE, details = FALSE, quiet = TRUE, seedit = seed,
               parallel_args = list(plan = "sequential"))
  ot <- switch(g$type,
    continuous = list(outcome.name = "y", event.name = "y", outcome_type = "continuous",
                      effect_measure = "MD", adverse_outcome = TRUE,
                      effect.threshold = 0.5, consistency.threshold = 0.0),
    survival   = list(outcome.name = "time", event.name = "event", outcome_type = "survival",
                      hr.threshold = 1.25, hr.consistency = 1.0),
    binary     = list(outcome.name = "y", event.name = "y", outcome_type = "binary",
                      effect_measure = "OR", adverse_outcome = TRUE,
                      hr.threshold = 1.10, hr.consistency = 1.0))
  args <- c(base, ot, g$extra %||% list(), extra)
  t0 <- proc.time()[["elapsed"]]
  fit <- suppressWarnings(suppressMessages(do.call(forestsearch, args)))
  list(fit = fit, secs = proc.time()[["elapsed"]] - t0)
}

# ---- extraction and comparison ------------------------------------------------
extract <- function(r) {
  fit <- r$fit; gc_ <- fit$grp.consistency
  hs <- fit$find.grps$out.found$hr.subgroups
  cand_keys <- if (is.null(hs) || !nrow(hs)) character(0) else {
    qc <- grep("^q[0-9]+\\.[01]$", names(hs), value = TRUE)
    apply(as.matrix(hs[, qc, with = FALSE]), 1, function(z) paste(sort(qc[z == 1]), collapse = "&"))
  }
  memb <- if (!is.null(gc_$sg.harm.id)) sort(which(gc_$sg.harm.id == 1)) else integer(0)
  eff  <- if (!is.null(gc_$out_sg$result) && nrow(gc_$out_sg$result)) gc_$out_sg$result$hr[1] else NA_real_
  pc   <- if (!is.null(gc_$out_sg$result) && nrow(gc_$out_sg$result)) gc_$out_sg$result$Pcons[1] else NA_real_
  fcnt <- fit$find.grps$filter_counts
  list(declared = !is.null(fit$sg.harm), sg_string = paste(fit$sg.harm, collapse = " & "),
       memb = memb, effect = eff, pcons = pc, cand = sort(cand_keys), n_cand = length(cand_keys),
       n_total = gc_$n_candidates_total %||% NA, n_eval = gc_$n_candidates_evaluated %||% NA,
       n_passed = gc_$n_passed %||% NA, es_trig = isTRUE(gc_$early_stop_triggered),
       es_cand = gc_$early_stop_candidate %||% NA,
       n_fits = fcnt$n_passed_sample_size %||% NA,      # combinations reaching a model fit
       n_passed_hr = fcnt$n_passed_hr %||% NA,
       conf_eval = fit$confounders.evaluated, secs = r$secs)
}
classify <- function(a, b) {
  if (a$declared != b$declared) return("differs substantively")
  if (!identical(a$memb, b$memb)) return("differs substantively")
  if (identical(a$sg_string, b$sg_string)) "identical" else "clause-order only"
}

# ---- direct timing of Section 5's forest on the same cut matrix ---------------
# Rebuild X exactly as Section 5 does: the evaluated cut columns as 0/1 numerics
# (FSconfounders.name = FSdata$confs_names; df[, ...] are the {0,1} factors).
cut_matrix <- function(df, labels) {
  X <- sapply(labels, function(lb) {
    if (grepl("<=", lb, fixed = TRUE)) {
      v <- trimws(strsplit(lb, "<=", fixed = TRUE)[[1]]); as.numeric(df[[v[1]]] <= as.numeric(v[2]))
    } else as.numeric(as.character(df[[lb]]))
  })
  X <- as.matrix(X); storage.mode(X) <- "numeric"; X
}
time_forest <- function(g, labels, seed) {
  X <- cut_matrix(g$df, labels); df <- g$df
  t0 <- proc.time()[["elapsed"]]
  if (g$type == "survival") {
    Y <- df$time; Event <- df$event; Treat <- df$treat
    y_evt_t <- Y[Treat == 1 & Event == 1]; y_evt_c <- Y[Treat == 0 & Event == 1]
    if (length(y_evt_t) == 0L || length(y_evt_c) == 0L) return(c(secs = 0, skipped = 1))
    tau.rmst <- min(max(y_evt_t), max(y_evt_c))
    invisible(try(suppressWarnings(grf::causal_survival_forest(X, Y, Treat, Event, W.hat = 0.5,
      horizon = 0.7 * tau.rmst, seed = seed, tune.parameters = "none")), TRUE))
  } else {
    Y_vi <- if (g$type == "binary") 1L - df$y else -df$y
    invisible(try(suppressWarnings(fit_causal_forest_glm(X, Y_vi, df$treat, TRUE, seedit = seed, tune_grf = FALSE)), TRUE))
  }
  c(secs = proc.time()[["elapsed"]] - t0, skipped = 0)
}

# ---- run --------------------------------------------------------------------
rows <- list(); pairs <- list()
run_block <- function(fx, sg_focus, pcons = 0.5) {
  for (s in SEEDS) {
    g <- GEN[[fx]](s)
    ra <- fs_call(g, -0.2, s, sg_focus, pcons = pcons); rb <- fs_call(g, NULL, s, sg_focus, pcons = pcons)
    a <- extract(ra); b <- extract(rb)
    cl <- classify(a, b)
    tf <- time_forest(g, b$conf_eval, s)
    rows[[length(rows) + 1L]] <<- data.frame(
      fixture = fx, focus = sg_focus, pcons = pcons, seed = s, class = cl,
      declared_a = a$declared, declared_b = b$declared,
      memb_same = identical(a$memb, b$memb), string_same = identical(a$sg_string, b$sg_string),
      n_harm_a = length(a$memb), n_harm_b = length(b$memb),
      effect_a = a$effect, effect_b = b$effect, pcons_a = a$pcons, pcons_b = b$pcons,
      cand_same = identical(a$cand, b$cand), n_cand_a = a$n_cand, n_cand_b = b$n_cand,
      n_conf_eval_a = length(a$conf_eval), n_conf_eval_b = length(b$conf_eval),
      conf_eval_same_set = setequal(a$conf_eval, b$conf_eval),
      n_total_a = a$n_total, n_total_b = b$n_total, n_eval_a = a$n_eval, n_eval_b = b$n_eval,
      n_passed_a = a$n_passed, n_passed_b = b$n_passed,
      es_trig_a = a$es_trig, es_trig_b = b$es_trig, es_cand_a = a$es_cand, es_cand_b = b$es_cand,
      n_fits_a = a$n_fits, n_fits_b = b$n_fits, n_passed_hr_a = a$n_passed_hr, n_passed_hr_b = b$n_passed_hr,
      secs_a = a$secs, secs_b = b$secs, secs_diff = a$secs - b$secs,
      forest_secs = tf[["secs"]], forest_skipped = tf[["skipped"]] == 1,
      stringsAsFactors = FALSE)
    pairs[[paste(fx, sg_focus, pcons, s, sep = "_")]] <<- list(
      sg_a = a$sg_string, sg_b = b$sg_string, conf_a = a$conf_eval, conf_b = b$conf_eval,
      cand_a = a$cand, cand_b = b$cand)
    cat(sprintf("%s %-10s pcons %.1f seed %2d: %-22s  n_eval %2s/%2s  es_cand %2s/%2s  secs %.2f/%.2f  forest %.2f\n",
                fx, sg_focus, pcons, s, cl, a$n_eval, b$n_eval, a$es_cand, b$es_cand, a$secs, b$secs, tf[["secs"]]))
  }
}
t_all <- proc.time()[["elapsed"]]
for (fx in names(GEN)) run_block(fx, "maxeffCons")
run_block("F1", "eff")                 # a non-maxeffCons focus: early stopping cannot fire
# strict early-stop arm: at pconsistency = 0.9 the first candidate rarely clears the
# gate, so the scan runs deeper and n_candidates_evaluated / early_stop_candidate
# can actually differ between the two orderings if the ordering reached the loop
for (fx in c("F1", "F2", "F4", "F7")) run_block(fx, "maxeffCons", pcons = 0.9)
res <- do.call(rbind, rows); rownames(res) <- NULL
cat(sprintf("\nall runs: %d pairs in %.0f s\n", nrow(res), proc.time()[["elapsed"]] - t_all))

# ---- summaries ------------------------------------------------------------------
cat("\n== classification by fixture ==\n")
print(with(res, table(paste(fixture, focus, pcons), class)))
cat("\n== candidate family identical as a set of rows; consistency counters identical ==\n")
print(aggregate(cbind(cand_same, conf_eval_same_set,
                      counters_same = (n_total_a == n_total_b) & (n_eval_a == n_eval_b) & (n_passed_a == n_passed_b)) ~ fixture + focus + pcons,
                data = res, FUN = mean))
cat("\n== efficiency: paired differences (-0.2 minus NULL) ==\n")
eff <- aggregate(cbind(secs_a, secs_b, secs_diff, forest_secs,
                       d_n_eval = n_eval_a - n_eval_b, d_es_cand = es_cand_a - es_cand_b,
                       d_n_fits = n_fits_a - n_fits_b, es_rate_a = es_trig_a, es_rate_b = es_trig_b) ~ fixture + focus + pcons,
                 data = res, FUN = function(v) mean(v, na.rm = TRUE))
print(eff, digits = 3)
se <- function(v) { v <- v[is.finite(v)]; if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA }
cat("\n== paired-difference SEs ==\n")
print(aggregate(cbind(secs_diff, d_n_eval = n_eval_a - n_eval_b, d_es_cand = es_cand_a - es_cand_b) ~ fixture + focus + pcons,
                data = res, FUN = se), digits = 3)
saveRDS(list(res = res, pairs = pairs, seeds = SEEDS, n = N,
             built_at = Sys.time(), pkg_version = as.character(utils::packageVersion("forestsearch"))), OUT)
cat("written:", OUT, "\n")
