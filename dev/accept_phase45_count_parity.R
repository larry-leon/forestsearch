# accept_phase45_count_parity.R -------------------------------------------
# [delivery sentinel: p45r1-0f32cf97]
# Phase 4.5 acceptance gate: count gate-lift for subgroup_glm().
#
# Proves, against the INSTALLED forestsearch build:
#   G1  IRR parity vs an independent hand chain (Poisson + offset,
#       mirrored robust/model SE conditional incl. the sandwich
#       fallback and the eps exposure clamp), robust_se TRUE and FALSE.
#   G2  IRD parity vs the transcribed delta chain (mean log-exposure,
#       grad c(lam1-lam0, lam1), V[1:2,1:2]).
#   G3  Unit-vs-varying exposure contrast: t == 1 reproduces the
#       no-offset Poisson model bitwise; varying t differs.
#   G4  Degenerates characterized then asserted: zero-count arm,
#       single-arm (aliased), tiny n, zero-exposure rows (clamp path).
#   G5  Support-gate messages: count without offset.name (wrapper stop
#       pre-empting the estimator's terser construction error), bad
#       measure on count; binary rate-measure message unchanged.
#   G6  attr(, "effect") stamps: IRR ratio stamps; IRD identity stamps.
#   G7  D-B4 default fit: count DGM (offset_var set) default == explicit;
#       count DGM without offset_var stops at construction; continuous
#       and binary default paths byte-identical; schedule invariance.
#   G8  Serialization ceiling on the count fitter.
#   G9  Display: Pr(IRR<0.5)-family headers; delegated axes for IRR
#       (ratio + metadata) and IRD (identity).
#   G10 Calibration smoke: IRR(Q) = 2.0 with the offset fixture.
#
# Stop-first-failure; each gate prints its evidence before the verdict.
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(forestsearch)
  library(survival)
  library(future)
})

fail <- function(...) stop(sprintf(...), call. = FALSE)
gate <- function(id, ok, detail = "") {
  cat(sprintf("[%s] %s %s\n", id, if (ok) "PASS" else "FAIL", detail))
  if (!ok) fail("Gate %s failed. %s", id, detail)
}
msg_of <- function(expr) tryCatch({ expr; NA_character_ },
                                  error = function(e) conditionMessage(e))

cat("== G0: environment ==========================================\n")
cat(sprintf("R %s | forestsearch %s | sandwich %s\n",
            paste(R.version$major, R.version$minor, sep = "."),
            as.character(utils::packageVersion("forestsearch")),
            if (requireNamespace("sandwich", quietly = TRUE))
              "present (robust branch active)" else
              "ABSENT (model-SE fallback branch active)"))

# -- Fixture: gbsg-derived Poisson counts with varying exposure ------------
gb <- survival::gbsg
set.seed(11)
lin <- -0.6 + 0.45 * (gb$grade == 3) + 0.30 * as.integer(gb$meno) -
  0.25 * as.integer(gb$hormon) + 0.004 * (gb$age - 50)
t_exp <- gb$rfstime / 365.25
gc2 <- data.frame(
  y      = stats::rpois(nrow(gb), lambda = exp(lin) * t_exp),
  w      = as.integer(gb$hormon),
  t      = t_exp,
  meno   = as.integer(gb$meno),
  grade3 = as.integer(gb$grade == 3),
  age    = gb$age,
  size   = gb$size
)
set.seed(7)
idx40 <- sample.int(nrow(gc2), 40L)
subsets <- list(
  full        = gc2,
  meno1       = gc2[gc2$meno == 1, ],
  meno0       = gc2[gc2$meno == 0, ],
  grade3      = gc2[gc2$grade3 == 1, ],
  young       = gc2[gc2$age <= 45, ],
  large       = gc2[gc2$size > 30, ],
  meno1grade3 = gc2[gc2$meno == 1 & gc2$grade3 == 1, ],
  rand40      = gc2[idx40, ]
)

# -- Hand chains (independent transcriptions) ------------------------------
.hand_se_irr <- function(fit, robust) {
  if (isTRUE(robust)) {
    tryCatch({
      if (!requireNamespace("sandwich", quietly = TRUE)) stop("ns")
      sqrt(sandwich::sandwich(fit)["w", "w"])
    }, error = function(e) sqrt(diag(stats::vcov(fit)))[["w"]])
  } else sqrt(diag(stats::vcov(fit)))[["w"]]
}
hand_irr <- function(df, robust = TRUE, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  tv <- df$t
  if (any(tv <= 0, na.rm = TRUE)) tv <- pmax(tv, .Machine$double.eps)
  df$.log_offset <- log(tv)
  r <- tryCatch({
    fit <- stats::glm(y ~ w + offset(.log_offset), data = df,
                      family = stats::poisson(link = "log"))
    list(est = stats::coef(fit)[["w"]], se = .hand_se_irr(fit, robust),
         conv = fit$converged)
  }, error = function(e) NULL)
  if (is.null(r) || !isTRUE(r$conv) || is.na(r$est) || is.na(r$se))
    return(c(NA_real_, NA_real_))
  out <- c(exp(r$est), exp(r$est + z * r$se))
  if (!all(is.finite(out))) return(c(NA_real_, NA_real_))
  out
}
hand_ird <- function(df, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  tv <- df$t
  if (any(tv <= 0, na.rm = TRUE)) tv <- pmax(tv, .Machine$double.eps)
  df$.log_offset <- log(tv)
  r <- tryCatch({
    fit <- stats::glm(y ~ w + offset(.log_offset), data = df,
                      family = stats::poisson(link = "log"))
    beta <- stats::coef(fit); V <- stats::vcov(fit)
    mlt  <- mean(log(tv), na.rm = TRUE)
    lam0 <- exp(beta[[1L]] + mlt)
    lam1 <- exp(beta[[1L]] + beta[["w"]] + mlt)
    grad <- c(lam1 - lam0, lam1)
    list(est = lam1 - lam0,
         se  = sqrt(as.numeric(t(grad) %*% V[1:2, 1:2] %*% grad)),
         conv = fit$converged)
  }, error = function(e) NULL)
  if (is.null(r) || !isTRUE(r$conv) || is.na(r$est) || is.na(r$se))
    return(c(NA_real_, NA_real_))
  out <- c(r$est, r$est + z * r$se)
  if (!all(is.finite(out))) return(c(NA_real_, NA_real_))
  out
}

cat("== G1: IRR parity vs hand chain =============================\n")
fit_irr  <- subgroup_glm(outcome.name = "y", treat.name = "w",
                         outcome_type = "count", offset.name = "t")
fit_irrM <- subgroup_glm(outcome.name = "y", treat.name = "w",
                         outcome_type = "count", offset.name = "t",
                         robust_se = FALSE)
for (nm in names(subsets)) {
  d <- subsets[[nm]]
  a <- fit_irr(d);  b <- hand_irr(d, robust = TRUE)
  gate(paste0("G1/", nm, "/robust"), identical(a, b),
       sprintf("wrapper=(%.6g, %.6g)", a[1], a[2]))
  a2 <- fit_irrM(d); b2 <- hand_irr(d, robust = FALSE)
  gate(paste0("G1/", nm, "/model"), identical(a2, b2), "")
}

cat("== G2: IRD parity vs transcribed delta chain ================\n")
fit_ird <- subgroup_glm(outcome.name = "y", treat.name = "w",
                        outcome_type = "count", offset.name = "t",
                        effect_measure = "IRD")
for (nm in names(subsets)) {
  d <- subsets[[nm]]
  a <- fit_ird(d); b <- hand_ird(d)
  gate(paste0("G2/", nm), identical(a, b),
       sprintf("wrapper=(%.6g, %.6g)", a[1], a[2]))
}

cat("== G3: unit-vs-varying exposure contrast ====================\n")
d1 <- subsets$full; d1$t <- 1
a_unit <- fit_irr(d1)
fit0 <- stats::glm(y ~ w, data = d1,
                   family = stats::poisson(link = "log"))
gate("G3/unit==no-offset",
     identical(a_unit[1], exp(stats::coef(fit0)[["w"]])),
     "log(1)=0 offset reproduces the no-offset model bitwise")
gate("G3/varying-differs",
     !identical(fit_irr(subsets$full)[1], exp(stats::coef(fit0)[["w"]])),
     "")

cat("== G4: degenerate fixtures ==================================\n")
deg <- list(
  zero_count_arm = within(gc2[1:40, ], y <- ifelse(w == 1L, 0L, y)),
  single_arm     = within(gc2[1:30, ], w <- 1L),
  tiny_n3        = gc2[1:3, ],
  zero_exposure  = { d <- gc2[1:40, ]; d$t[1:5] <- 0; d }
)
for (nm in names(deg)) {
  d <- deg[[nm]]
  a <- fit_irr(d); b <- hand_irr(d, robust = TRUE)
  gate(paste0("G4/", nm, "/IRR"), identical(a, b),
       sprintf("wrapper=(%s, %s)", format(a[1]), format(a[2])))
  ar <- fit_ird(d); br <- hand_ird(d)
  gate(paste0("G4/", nm, "/IRD"), identical(ar, br), "")
}
gate("G4/single_arm/NA", all(is.na(fit_irr(deg$single_arm))) &&
       all(is.na(fit_ird(deg$single_arm))),
     "aliased treatment -> c(NA, NA) on both measures")

cat("== G5: support-gate messages ================================\n")
m1 <- msg_of(subgroup_glm(outcome_type = "count"))
gate("G5/count-no-offset",
     grepl("requires", m1, fixed = TRUE) &&
       grepl("offset.name", m1, fixed = TRUE) &&
       grepl("unit column of 1s", m1, fixed = TRUE), m1)
m1e <- msg_of(make_effect_estimator("count", treat.name = "w",
                                    outcome.name = "y",
                                    effect_measure = "IRR",
                                    offset.name = NULL))
gate("G5/wrapper-pre-empts-estimator",
     !identical(m1, m1e) && grepl("offset.name", m1e, fixed = TRUE),
     sprintf("estimator msg: %s", m1e))
m2 <- msg_of(subgroup_glm(outcome_type = "count", offset.name = "t",
                          effect_measure = "OR"))
gate("G5/count-OR",
     grepl('not available for outcome_type = "count"', m2, fixed = TRUE),
     m2)
m3 <- msg_of(subgroup_glm(outcome_type = "binary", effect_measure = "IRR"))
gate("G5/binary-IRR-unchanged",
     grepl("rate measure (requires offset.name)", m3, fixed = TRUE), m3)
fit_null <- subgroup_glm(outcome_type = "count", offset.name = "t")
gate("G5/count-NULL-is-IRR",
     identical(attr(fit_null, "effect")$measure, "IRR"), "")

cat("== G6: effect-metadata stamps ===============================\n")
effI <- attr(fit_irr, "effect")
gate("G6/IRR",
     identical(effI$measure, "IRR") && isTRUE(effI$log_scale) &&
       identical(effI$null_value, 1) &&
       identical(effI$est_thresholds, c(0.5, 1)) &&
       identical(effI$ub_thresholds, c(2, 3)) &&
       identical(effI$est_label, "IRR") &&
       identical(effI$ub_label, "UB(IRR)"), "")
effD <- attr(fit_ird, "effect")
gate("G6/IRD",
     identical(effD$measure, "IRD") && !isTRUE(effD$log_scale) &&
       identical(effD$null_value, 0) &&
       identical(effD$est_thresholds, c(NA_real_, 0)) &&
       identical(effD$ub_thresholds, c(NA_real_, NA_real_)), "")

cat("== G7: DGM-aware default fit with offset (D-B4) =============\n")
sg <- list(
  list(id = "flag_itt == 1",           name = "All Patients", grp = "ITT"),
  list(id = "meno == 1",               name = "Post-meno",    grp = "Clinical"),
  list(id = "meno == 0",               name = "Pre-meno",     grp = "Clinical"),
  list(id = "grade3 == 1",             name = "Grade 3",      grp = "Clinical"),
  list(id = "meno == 1 & grade3 == 1", name = "Post/G3",      grp = "Interaction")
)
dgm_cnt <- generate_glm_dgm(
  data = gc2, factor_vars = c("meno", "grade3"),
  outcome_var = "y", treatment_var = "w", outcome_type = "count",
  offset_var = "t",
  subgroup_vars = c("meno", "grade3"),
  subgroup_cuts = list(meno = 1L, grade3 = 1L),
  model = "null", n_super = 2500L, seed = 201L)
gate("G7/offset-nesting",
     is.null(dgm_cnt$offset_var) &&
       identical(dgm_cnt$model_params$offset_var, "t"),
     "offset_var lives in model_params (top-level is absent) -- D-B4 premise")
plan(sequential)
res_def <- run_subgroup_sims(dgm_cnt, sg, n_sims = 8, n = 300,
                             benchmarks = NULL, workers = 1,
                             seed_base = 0L,
                             hr_true = dgm_cnt$hazard_ratios$overall)
res_exp <- run_subgroup_sims(dgm_cnt, sg, n_sims = 8, n = 300,
                             fit = subgroup_glm(
                               outcome_type    = dgm_cnt$outcome_type,
                               effect_measure  = dgm_cnt$effect_measure,
                               offset.name     = dgm_cnt$model_params$offset_var,
                               adverse_outcome = TRUE),
                             benchmarks = NULL, workers = 1,
                             seed_base = 0L,
                             hr_true = dgm_cnt$hazard_ratios$overall)
gate("G7/count-default==explicit",
     identical(res_def$sim_hrs, res_exp$sim_hrs) &&
       identical(res_def$sim_ubs, res_exp$sim_ubs) &&
       identical(res_def$sim_ns,  res_exp$sim_ns) &&
       identical(res_def$effect,  res_exp$effect),
     sprintf("measure carried: %s", res_def$effect$measure))
gate("G7/count-measure-is-IRR",
     identical(res_def$effect$measure, "IRR"), "")

res_w2 <- run_subgroup_sims(dgm_cnt, sg, n_sims = 8, n = 300,
                            benchmarks = NULL, workers = 2,
                            seed_base = 0L,
                            hr_true = dgm_cnt$hazard_ratios$overall)
plan(sequential)
gate("G7/schedule-invariance",
     identical(res_def$sim_hrs, res_w2$sim_hrs) &&
       identical(res_def$sim_ubs, res_w2$sim_ubs), "")

dgm_cnt0 <- generate_glm_dgm(
  data = gc2, factor_vars = c("meno", "grade3"),
  outcome_var = "y", treatment_var = "w", outcome_type = "count",
  offset_var = NULL,
  subgroup_vars = c("meno", "grade3"),
  subgroup_cuts = list(meno = 1L, grade3 = 1L),
  model = "null", n_super = 1000L, seed = 202L)
m7 <- msg_of(run_subgroup_sims(dgm_cnt0, sg, n_sims = 2, n = 100,
                               benchmarks = NULL, workers = 1,
                               seed_base = 0L))
gate("G7/no-offset-DGM-stops",
     grepl("requires", m7, fixed = TRUE) &&
       grepl("offset.name", m7, fixed = TRUE) &&
       grepl("unit column of 1s", m7, fixed = TRUE),
     "count DGM without offset_var stops at fitter construction (D-B5)")

dgm_cont <- generate_glm_dgm(
  data = gc2, factor_vars = c("meno", "grade3"),
  outcome_var = "age", treatment_var = "w", outcome_type = "continuous",
  subgroup_vars = c("meno", "grade3"),
  subgroup_cuts = list(meno = 1L, grade3 = 1L),
  model = "null", n_super = 2500L, seed = 102L)
res_c_def <- run_subgroup_sims(dgm_cont, sg, n_sims = 6, n = 300,
                               benchmarks = NULL, workers = 1,
                               seed_base = 0L,
                               hr_true = dgm_cont$hazard_ratios$overall)
res_c_old <- run_subgroup_sims(dgm_cont, sg, n_sims = 6, n = 300,
                               fit = subgroup_glm(),
                               benchmarks = NULL, workers = 1,
                               seed_base = 0L,
                               hr_true = dgm_cont$hazard_ratios$overall)
gate("G7/continuous-byte-identity",
     identical(res_c_def$sim_hrs, res_c_old$sim_hrs) &&
       identical(res_c_def$sim_ubs, res_c_old$sim_ubs) &&
       identical(res_c_def$effect,  res_c_old$effect),
     "D-B4's added offset argument is inert for continuous DGMs")

dgm_bin <- generate_glm_dgm(
  data = transform(gc2, yb = as.integer(y > 0)),
  factor_vars = c("meno", "grade3"),
  outcome_var = "yb", treatment_var = "w", outcome_type = "binary",
  subgroup_vars = c("meno", "grade3"),
  subgroup_cuts = list(meno = 1L, grade3 = 1L),
  model = "null", n_super = 2500L, seed = 103L)
res_b_def <- run_subgroup_sims(dgm_bin, sg, n_sims = 6, n = 300,
                               benchmarks = NULL, workers = 1,
                               seed_base = 0L,
                               hr_true = dgm_bin$hazard_ratios$overall)
res_b_exp <- run_subgroup_sims(dgm_bin, sg, n_sims = 6, n = 300,
                               fit = subgroup_glm(
                                 outcome_type    = "binary",
                                 effect_measure  = dgm_bin$effect_measure,
                                 adverse_outcome = isTRUE(dgm_bin$adverse_outcome)),
                               benchmarks = NULL, workers = 1,
                               seed_base = 0L,
                               hr_true = dgm_bin$hazard_ratios$overall)
gate("G7/binary-byte-identity",
     identical(res_b_def$sim_hrs, res_b_exp$sim_hrs) &&
       identical(res_b_def$sim_ubs, res_b_exp$sim_ubs) &&
       identical(res_b_def$effect,  res_b_exp$effect),
     "binary default path unchanged by D-B4")

cat("== G8: serialization ceiling ================================\n")
sz <- length(serialize(fit_irr, NULL))
gate("G8/512KB", sz < 512 * 1024,
     sprintf("count IRR fitter serializes at %d bytes", sz))

cat("== G9: display policy =======================================\n")
S_irr <- summary(res_def)
hdr <- colnames(S_irr$results_tbl)
need <- c("Pr(IRR<0.5)", "Pr(IRR>1.0)", "mIRR",
          "Pr(UB(IRR)>=2)", "Pr(UB(IRR)>=3)", "mUB(IRR)")
gate("G9/IRR-headers", all(need %in% hdr), paste(hdr, collapse = " | "))
res_ird <- run_subgroup_sims(dgm_cnt, sg, n_sims = 6, n = 300,
                             fit = subgroup_glm(
                               outcome_type = "count",
                               effect_measure = "IRD",
                               offset.name = "t"),
                             benchmarks = NULL, workers = 1,
                             seed_base = 0L,
                             hr_true = 0)
S_ird <- summary(res_ird)
gate("G9/IRD-headers", "mIRD" %in% colnames(S_ird$results_tbl) &&
       "mUB(IRD)" %in% colnames(S_ird$results_tbl),
     paste(colnames(S_ird$results_tbl), collapse = " | "))

ns <- asNamespace("forestsearch")
orig_gg <- get("gg_forest", envir = ns)
cap <- new.env()
stub <- function(...) { cap$args <- list(...); invisible(NULL) }
unlockBinding("gg_forest", ns); assign("gg_forest", stub, envir = ns)
lockBinding("gg_forest", ns)
restore_gg <- function() {
  unlockBinding("gg_forest", ns); assign("gg_forest", orig_gg, envir = ns)
  lockBinding("gg_forest", ns)
}
ok9 <- tryCatch({
  cap$args <- NULL; plot(S_irr, metric = "ub", panel = "combo")
  gate("G9/IRR-delegated-axes",
       is.null(cap$args$xlim) && is.null(cap$args$ticks_at) &&
         isTRUE(cap$args$xlog),
       "ratio + metadata -> data-driven axes")
  cap$args <- NULL; plot(S_ird, metric = "hr", panel = "single")
  gate("G9/IRD-delegated-axes",
       is.null(cap$args$xlim) && identical(cap$args$xlog, FALSE),
       "identity delegation as Phase 4.3/4.4")
  TRUE
}, finally = restore_gg())
stopifnot(isTRUE(ok9))

cat("== G10: calibration smoke, IRR(Q) = 2.0 =====================\n")
dgm_cal <- calibrate_glm_interaction(
  data = gc2, factor_vars = c("meno", "grade3"),
  outcome_var = "y", treatment_var = "w", target_effect = 2.0,
  outcome_type = "count", offset_var = "t",
  subgroup_vars = c("meno", "grade3"),
  subgroup_cuts = list(meno = 1L, grade3 = 1L),
  k_inter_range = c(0, 3), grid_step = 0.05,
  n_super = 2500L, seed = 201L, verbose = FALSE)
achieved <- dgm_cal$hazard_ratios$harm_subgroup
gate("G10/class", inherits(dgm_cal, "glm_dgm"), "")
gate("G10/achieved",
     is.finite(achieved) && abs(log(achieved / 2.0)) <= log(1.025) + 1e-8,
     sprintf("achieved IRR(Q) = %.4f (target 2.0), k_inter = %.3f",
             achieved, dgm_cal$model_params$k_inter))

cat("\nALL GATES PASS -- Phase 4.5 count parity accepted.\n")
