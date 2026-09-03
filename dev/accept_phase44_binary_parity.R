# accept_phase44_binary_parity.R ------------------------------------------
# [delivery sentinel: p44r1-6e2e3f60]
# Phase 4.4 acceptance gate: binary gate-lift for subgroup_glm().
#
# Proves, against the INSTALLED forestsearch build:
#   G1  OR parity: subgroup_glm(binary, OR) vs an independent hand chain
#       (glm logistic, model SE, qnorm UB, exp back-transform), both
#       adverse_outcome settings, across deterministic gbsg subsets.
#   G2  RD parity: subgroup_glm(binary, RD) vs a hand transcription of
#       the three-tier RD chain (identity-binomial -> logistic
#       G-computation -> raw means with converged = FALSE).
#   G3  Degenerate fixtures: single-arm, perfect separation, zero-event
#       arm, tiny n -- end-to-end hand-chain equality including the
#       wrapper's converged/NA/finite guards.
#   G4  robust_se no-op on the validated measures (OR and RD).
#   G5  Support-gate messages verbatim: count, RR, IRR/IRD, bad binary
#       measure, bad continuous measure; continuous NULL -> MD intact.
#   G6  attr(, "effect") stamps: binary OR c(0.5,1)/c(2,3); binary RD
#       and continuous MD identity stamps unchanged (D-A4).
#   G7  DGM-aware default fit (D-A5): missing(fit) == explicit
#       DGM-matched subgroup_glm() for a binary DGM; continuous default
#       path byte-identical to the pre-4.4 fixed default; schedule
#       invariance (workers 1 vs 2).
#   G8  Serialization ceiling: binary fitter ships under 512 KB.
#   G9  Display policy (D-A4/D-A3): binary summary headers
#       Pr(OR<0.5)/Pr(OR>1.0)/Pr(UB(OR)>=2)/Pr(UB(OR)>=3); plot axis
#       delegation for metadata-carrying results; HR constants retained
#       for legacy and metadata-less (override) survival summaries.
#   G10 Calibration smoke: calibrate_glm_interaction() to OR(Q) = 2.0.
#
# Stop-first-failure; each gate prints its evidence before the verdict.
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(forestsearch)
  library(survival)
  library(future)
})

z975 <- stats::qnorm(0.975)

fail <- function(...) stop(sprintf(...), call. = FALSE)
gate <- function(id, ok, detail = "") {
  cat(sprintf("[%s] %s %s\n", id, if (ok) "PASS" else "FAIL", detail))
  if (!ok) fail("Gate %s failed. %s", id, detail)
}

cat("== G0: environment ==========================================\n")
cat(sprintf("R %s | forestsearch %s\n",
            paste(R.version$major, R.version$minor, sep = "."),
            as.character(utils::packageVersion("forestsearch"))))

# -- Fixture: gbsg binary --------------------------------------------------
gb <- survival::gbsg
gb2 <- data.frame(
  y      = as.integer(gb$status),
  w      = as.integer(gb$hormon),
  meno   = as.integer(gb$meno),
  grade3 = as.integer(gb$grade == 3),
  age    = gb$age,
  size   = gb$size
)
set.seed(7)
idx40 <- sample.int(nrow(gb2), 40L)
idx18 <- sample.int(nrow(gb2), 18L)
subsets <- list(
  full        = gb2,
  meno1       = gb2[gb2$meno == 1, ],
  meno0       = gb2[gb2$meno == 0, ],
  grade3      = gb2[gb2$grade3 == 1, ],
  grade12     = gb2[gb2$grade3 == 0, ],
  young       = gb2[gb2$age <= 45, ],
  large       = gb2[gb2$size > 30, ],
  meno1grade3 = gb2[gb2$meno == 1 & gb2$grade3 == 1, ],
  rand40      = gb2[idx40, ],
  rand18      = gb2[idx18, ]
)

# -- Hand chains (independent transcriptions of the documented paths) ------
hand_or <- function(df, adverse = TRUE, level = 0.95) {
  if (!adverse) df$y <- 1L - df$y
  z <- stats::qnorm(1 - (1 - level) / 2)
  r <- tryCatch({
    fit <- stats::glm(y ~ w, data = df,
                      family = stats::binomial(link = "logit"))
    list(est = stats::coef(fit)[["w"]],
         se  = sqrt(diag(stats::vcov(fit)))[["w"]],
         conv = fit$converged)
  }, error = function(e) NULL)
  if (is.null(r) || !isTRUE(r$conv) || is.na(r$est) || is.na(r$se)) {
    return(c(NA_real_, NA_real_))
  }
  out <- c(exp(r$est), exp(r$est + z * r$se))
  if (!all(is.finite(out))) return(c(NA_real_, NA_real_))
  out
}

hand_rd <- function(df, adverse = TRUE, level = 0.95) {
  if (!adverse) df$y <- 1L - df$y
  z  <- stats::qnorm(1 - (1 - level) / 2)
  y0 <- df$y[df$w == 0]; y1 <- df$y[df$w == 1]
  p0 <- mean(y0, na.rm = TRUE); p1 <- mean(y1, na.rm = TRUE)
  n0 <- sum(df$w == 0L, na.rm = TRUE); n1 <- sum(df$w == 1L, na.rm = TRUE)
  tier1 <- tryCatch({
    fit <- stats::glm(y ~ w, data = df,
                      family = stats::binomial(link = "identity"),
                      start = c(p0, p1 - p0))
    if (!fit$converged) stop("nc")
    list(est = stats::coef(fit)[["w"]],
         se  = sqrt(diag(stats::vcov(fit)))[["w"]], conv = TRUE)
  },
  error = function(e) NULL,
  warning = function(wn) {
    if (grepl("converge", conditionMessage(wn), ignore.case = TRUE)) {
      return(NULL)
    }
    invokeRestart("muffleWarning")
  })
  r <- tier1
  if (is.null(r)) {
    r <- tryCatch({
      fitl <- stats::glm(y ~ w, data = df,
                         family = stats::binomial(link = "logit"))
      if (!fitl$converged) stop("nc")
      nd0 <- nd1 <- df; nd0$w <- 0L; nd1$w <- 1L
      pr0 <- stats::predict(fitl, newdata = nd0, type = "response")
      pr1 <- stats::predict(fitl, newdata = nd1, type = "response")
      se_t <- sqrt(diag(stats::vcov(fitl)))[["w"]]
      list(est = mean(pr1) - mean(pr0),
           se  = abs(mean(pr1 * (1 - pr1))) * se_t, conv = TRUE)
    }, error = function(e) NULL)
  }
  if (is.null(r)) {
    r <- list(est = p1 - p0,
              se  = sqrt(p1 * (1 - p1) / max(n1, 1L) +
                         p0 * (1 - p0) / max(n0, 1L)),
              conv = FALSE)
  }
  if (!isTRUE(r$conv) || is.na(r$est) || is.na(r$se)) {
    return(c(NA_real_, NA_real_))
  }
  out <- c(r$est, r$est + z * r$se)
  if (!all(is.finite(out))) return(c(NA_real_, NA_real_))
  out
}

cat("== G1: OR parity vs hand chain ==============================\n")
fit_or  <- subgroup_glm(outcome.name = "y", treat.name = "w",
                        outcome_type = "binary")
fit_orB <- subgroup_glm(outcome.name = "y", treat.name = "w",
                        outcome_type = "binary", adverse_outcome = FALSE)
for (nm in names(subsets)) {
  d <- subsets[[nm]]
  a <- fit_or(d);  b <- hand_or(d, adverse = TRUE)
  gate(paste0("G1/", nm, "/advT"), identical(a, b),
       sprintf("wrapper=(%.6g, %.6g) hand=(%.6g, %.6g)",
               a[1], a[2], b[1], b[2]))
  a2 <- fit_orB(d); b2 <- hand_or(d, adverse = FALSE)
  gate(paste0("G1/", nm, "/advF"), identical(a2, b2), "")
}

cat("== G2: RD parity vs transcribed tier chain ==================\n")
fit_rd <- subgroup_glm(outcome.name = "y", treat.name = "w",
                       outcome_type = "binary", effect_measure = "RD")
est_rd_raw <- make_effect_estimator("binary", treat.name = "w",
                                    outcome.name = "y",
                                    effect_measure = "RD")
for (nm in names(subsets)) {
  d <- subsets[[nm]]
  a <- fit_rd(d); b <- hand_rd(d, adverse = TRUE)
  mu <- est_rd_raw(d)$method_used
  gate(paste0("G2/", nm), identical(a, b),
       sprintf("method=%s wrapper=(%.6g, %.6g)", mu, a[1], a[2]))
}

cat("== G3: degenerate fixtures ==================================\n")
deg <- list(
  single_arm = within(gb2[1:30, ], w <- 1L),
  separation = within(gb2[1:24, ], y <- w),
  zero_event = within(gb2[1:24, ], y <- 0L),
  tiny_n3    = gb2[1:3, ]
)
for (nm in names(deg)) {
  d <- deg[[nm]]
  a  <- fit_or(d); b  <- hand_or(d, adverse = TRUE)
  ar <- fit_rd(d); br <- hand_rd(d, adverse = TRUE)
  gate(paste0("G3/", nm, "/OR"), identical(a, b),
       sprintf("wrapper=(%s, %s)", format(a[1]), format(a[2])))
  gate(paste0("G3/", nm, "/RD"), identical(ar, br),
       sprintf("wrapper=(%s, %s)", format(ar[1]), format(ar[2])))
}
gate("G3/single_arm/NA", all(is.na(fit_or(deg$single_arm))) &&
       all(is.na(fit_rd(deg$single_arm))),
     "aliased treatment -> c(NA, NA) on both measures")

cat("== G4: robust_se no-op on validated measures ================\n")
for (m in c("OR", "RD")) {
  eT <- make_effect_estimator("binary", treat.name = "w",
                              outcome.name = "y", effect_measure = m,
                              robust_se = TRUE)
  eF <- make_effect_estimator("binary", treat.name = "w",
                              outcome.name = "y", effect_measure = m,
                              robust_se = FALSE)
  same <- vapply(subsets, function(d) {
    a <- eT(d); b <- eF(d)
    identical(a$estimate, b$estimate) && identical(a$se, b$se)
  }, logical(1L))
  gate(paste0("G4/", m), all(same),
       sprintf("robust_se TRUE == FALSE on %d subsets", length(same)))
}

cat("== G5: support-gate messages ================================\n")
msg_of <- function(expr) tryCatch({ expr; NA_character_ },
                                  error = function(e) conditionMessage(e))
m1 <- msg_of(subgroup_glm(outcome_type = "count"))
gate("G5/count", grepl('outcome_type = "count" is not yet wrapper-validated',
                       m1, fixed = TRUE) &&
       grepl("Phase 4.5 count arc", m1, fixed = TRUE), m1)
m2 <- msg_of(subgroup_glm(outcome_type = "binary", effect_measure = "RR"))
gate("G5/RR", grepl("log-binomial -> modified-Poisson chain", m2,
                    fixed = TRUE), m2)
m3 <- msg_of(subgroup_glm(outcome_type = "binary", effect_measure = "IRR"))
gate("G5/IRR", grepl("rate measure (requires offset.name)", m3,
                     fixed = TRUE), m3)
m4 <- msg_of(subgroup_glm(outcome_type = "binary", effect_measure = "IRD"))
gate("G5/IRD", grepl("rate measure (requires offset.name)", m4,
                     fixed = TRUE), m4)
m5 <- msg_of(subgroup_glm(outcome_type = "binary", effect_measure = "MD"))
gate("G5/binary-MD",
     grepl('not available for outcome_type = "binary"', m5, fixed = TRUE),
     m5)
m6 <- msg_of(subgroup_glm(outcome_type = "continuous",
                          effect_measure = "OR"))
gate("G5/continuous-OR",
     grepl('not available for outcome_type = "continuous"', m6,
           fixed = TRUE), m6)

cat("== G6: effect-metadata stamps (D-A4) ========================\n")
effA <- attr(fit_or, "effect")
gate("G6/OR",
     identical(effA$measure, "OR") && isTRUE(effA$log_scale) &&
       identical(effA$null_value, 1) &&
       identical(effA$est_thresholds, c(0.5, 1)) &&
       identical(effA$ub_thresholds, c(2, 3)) &&
       identical(effA$est_label, "OR") &&
       identical(effA$ub_label, "UB(OR)"),
     sprintf("est_thr=%s ub_thr=%s",
             paste(effA$est_thresholds, collapse = ","),
             paste(effA$ub_thresholds, collapse = ",")))
effR <- attr(fit_rd, "effect")
gate("G6/RD",
     identical(effR$measure, "RD") && !isTRUE(effR$log_scale) &&
       identical(effR$null_value, 0) &&
       identical(effR$est_thresholds, c(NA_real_, 0)) &&
       identical(effR$ub_thresholds, c(NA_real_, NA_real_)), "")
effM <- attr(subgroup_glm(), "effect")
gate("G6/MD-unchanged",
     identical(effM$measure, "MD") &&
       identical(effM$est_thresholds, c(NA_real_, 0)) &&
       identical(effM$ub_thresholds, c(NA_real_, NA_real_)),
     "continuous stamps byte-identical to Phase 4.1")

cat("== G7: DGM-aware default fit (D-A5) =========================\n")
dgm_bin <- generate_glm_dgm(
  data = gb2, factor_vars = c("meno", "grade3"),
  outcome_var = "y", treatment_var = "w", outcome_type = "binary",
  subgroup_vars = c("meno", "grade3"),
  subgroup_cuts = list(meno = 1L, grade3 = 1L),
  model = "null", n_super = 2500L, seed = 101L)
sg <- list(
  list(id = "flag_itt == 1",             name = "All Patients", grp = "ITT"),
  list(id = "meno == 1",                 name = "Post-meno",    grp = "Clinical"),
  list(id = "meno == 0",                 name = "Pre-meno",     grp = "Clinical"),
  list(id = "grade3 == 1",               name = "Grade 3",      grp = "Clinical"),
  list(id = "meno == 1 & grade3 == 1",   name = "Post/G3",      grp = "Interaction")
)
plan(sequential)
res_def <- run_subgroup_sims(dgm_bin, sg, n_sims = 8, n = 300,
                             benchmarks = NULL, workers = 1,
                             seed_base = 0L,
                             hr_true = dgm_bin$hazard_ratios$overall)
res_exp <- run_subgroup_sims(dgm_bin, sg, n_sims = 8, n = 300,
                             fit = subgroup_glm(
                               outcome_type    = dgm_bin$outcome_type,
                               effect_measure  = dgm_bin$effect_measure,
                               adverse_outcome = isTRUE(dgm_bin$adverse_outcome)),
                             benchmarks = NULL, workers = 1,
                             seed_base = 0L,
                             hr_true = dgm_bin$hazard_ratios$overall)
gate("G7/binary-default==explicit",
     identical(res_def$sim_hrs, res_exp$sim_hrs) &&
       identical(res_def$sim_ubs, res_exp$sim_ubs) &&
       identical(res_def$sim_ns,  res_exp$sim_ns) &&
       identical(res_def$effect,  res_exp$effect),
     sprintf("effect measure carried: %s", res_def$effect$measure))
gate("G7/binary-measure-is-OR",
     identical(res_def$effect$measure, "OR"),
     "default fit analyzes on the DGM's own scale")

res_w2 <- run_subgroup_sims(dgm_bin, sg, n_sims = 8, n = 300,
                            benchmarks = NULL, workers = 2,
                            seed_base = 0L,
                            hr_true = dgm_bin$hazard_ratios$overall)
plan(sequential)
gate("G7/schedule-invariance",
     identical(res_def$sim_hrs, res_w2$sim_hrs) &&
       identical(res_def$sim_ubs, res_w2$sim_ubs),
     "workers = 1 vs 2 give identical matrices")

dgm_cont <- generate_glm_dgm(
  data = gb2, factor_vars = c("meno", "grade3"),
  outcome_var = "age", treatment_var = "w", outcome_type = "continuous",
  subgroup_vars = c("meno", "grade3"),
  subgroup_cuts = list(meno = 1L, grade3 = 1L),
  model = "null", n_super = 2500L, seed = 102L)
res_c_def <- run_subgroup_sims(dgm_cont, sg, n_sims = 8, n = 300,
                               benchmarks = NULL, workers = 1,
                               seed_base = 0L,
                               hr_true = dgm_cont$hazard_ratios$overall)
res_c_old <- run_subgroup_sims(dgm_cont, sg, n_sims = 8, n = 300,
                               fit = subgroup_glm(),
                               benchmarks = NULL, workers = 1,
                               seed_base = 0L,
                               hr_true = dgm_cont$hazard_ratios$overall)
gate("G7/continuous-byte-identity",
     identical(res_c_def$sim_hrs, res_c_old$sim_hrs) &&
       identical(res_c_def$sim_ubs, res_c_old$sim_ubs) &&
       identical(res_c_def$sim_ns,  res_c_old$sim_ns) &&
       identical(res_c_def$effect,  res_c_old$effect),
     "DGM-aware default == pre-4.4 fixed default on continuous")

cat("== G8: serialization ceiling ================================\n")
sz <- length(serialize(fit_or, NULL))
gate("G8/512KB", sz < 512 * 1024,
     sprintf("binary OR fitter serializes at %d bytes", sz))

cat("== G9: display policy (D-A4 headers, D-A3 axes) =============\n")
S_bin <- summary(res_def)
hdr <- colnames(S_bin$results_tbl)
need <- c("Pr(OR<0.5)", "Pr(OR>1.0)", "mOR",
          "Pr(UB(OR)>=2)", "Pr(UB(OR)>=3)", "mUB(OR)")
gate("G9/headers", all(need %in% hdr),
     paste(hdr, collapse = " | "))
gate("G9/highrisk-populated",
     is.list(S_bin$highrisk) && any(S_bin$highrisk$include),
     "UB tail c(2,3) drives a working high-risk panel")

S_cont <- summary(res_c_def)

gb$time_months <- gb$rfstime / 30.4375
dgm_surv <- generate_aft_dgm_flex(
  data = gb,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "time_months", event_var = "status",
  treatment_var = "hormon", subgroup_vars = NULL,
  model = "null", n_super = 1500, seed = 99, verbose = FALSE)
sg_surv <- list(
  list(id = "flag_itt == 1", name = "All Patients", grp = "ITT"),
  list(id = "meno == 1",     name = "Post-meno",    grp = "Clinical"),
  list(id = "meno == 0",     name = "Pre-meno",     grp = "Clinical"),
  list(id = "grade == 3",    name = "Grade 3",      grp = "Clinical"),
  list(id = "meno == 1 & grade == 3", name = "Post/G3", grp = "Interaction")
)
res_surv <- run_subgroup_sims(dgm_surv, sg_surv, n_sims = 4, n = 200,
                              analysis_time = 84, max_entry = 24,
                              cens_adjust = 0.0412,
                              benchmarks = NULL, workers = 1,
                              seed_base = 0L, hr_true = 0.70)
S_leg <- summary(res_surv)
S_ov  <- summary(res_surv, est_thresholds = c(0.4, 1.2))

ns <- asNamespace("forestsearch")
orig_gg <- get("gg_forest", envir = ns)
cap <- new.env()
stub <- function(...) { cap$args <- list(...); invisible(NULL) }
unlockBinding("gg_forest", ns)
assign("gg_forest", stub, envir = ns)
lockBinding("gg_forest", ns)
restore_gg <- function() {
  unlockBinding("gg_forest", ns)
  assign("gg_forest", orig_gg, envir = ns)
  lockBinding("gg_forest", ns)
}
axis_probe <- function(S, metric, panel) {
  cap$args <- NULL
  plot(S, metric = metric, panel = panel)
  list(xlim = cap$args$xlim, ticks = cap$args$ticks_at,
       xlog = cap$args$xlog)
}
ok9 <- tryCatch({
  p1 <- axis_probe(S_bin, "hr", "single")
  gate("G9/binary-hr-delegated",
       is.null(p1$xlim) && is.null(p1$ticks) && isTRUE(p1$xlog),
       "metadata + ratio scale -> data-driven axes, xlog TRUE")
  p2 <- axis_probe(S_bin, "ub", "combo")
  gate("G9/binary-ub-delegated",
       is.null(p2$xlim) && is.null(p2$ticks) && isTRUE(p2$xlog), "")
  p3 <- axis_probe(S_cont, "hr", "single")
  gate("G9/continuous-delegated-unchanged",
       is.null(p3$xlim) && identical(p3$xlog, FALSE),
       "identity delegation exactly as Phase 4.3")
  p4 <- axis_probe(S_leg, "hr", "single")
  gate("G9/survival-legacy-constants",
       identical(p4$xlim, c(0.15, 3.5)) && isTRUE(p4$xlog),
       "legacy branch untouched")
  p5 <- axis_probe(S_ov, "hr", "single")
  gate("G9/survival-override-constants",
       identical(p5$xlim, c(0.15, 3.5)) &&
         identical(p5$ticks,
                   c(0.20, 0.35, 0.50, 0.70, 1.00, 1.50, 2.50)),
       "metadata-less generic branch keeps HR constants (D-A3 guard)")
  p6 <- axis_probe(S_ov, "ub", "combo")
  gate("G9/survival-override-ub-constants",
       identical(p6$xlim, c(0.30, 9.0)), "")
  TRUE
}, finally = restore_gg())
stopifnot(isTRUE(ok9))

cat("== G10: calibration smoke, OR(Q) = 2.0 ======================\n")
dgm_cal <- calibrate_glm_interaction(
  data = gb2, factor_vars = c("meno", "grade3"),
  outcome_var = "y", treatment_var = "w", target_effect = 2.0,
  outcome_type = "binary",
  subgroup_vars = c("meno", "grade3"),
  subgroup_cuts = list(meno = 1L, grade3 = 1L),
  k_inter_range = c(0, 3), grid_step = 0.05,
  n_super = 2500L, seed = 101L, verbose = FALSE)
achieved <- dgm_cal$hazard_ratios$harm_subgroup
gate("G10/class", inherits(dgm_cal, "glm_dgm"), "")
gate("G10/achieved",
     is.finite(achieved) && abs(log(achieved / 2.0)) <= log(1.025) + 1e-8,
     sprintf("achieved OR(Q) = %.4f (target 2.0, tol_rel 2.5%%), k_inter = %.3f",
             achieved, dgm_cal$model_params$k_inter))

cat("\nALL GATES PASS -- Phase 4.4 binary parity accepted.\n")
