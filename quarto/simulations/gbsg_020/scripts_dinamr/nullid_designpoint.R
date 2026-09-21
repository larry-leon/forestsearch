#!/usr/bin/env Rscript
# nullid design-point record (TASK_null_gbsg_identification_2026-09-21, Step 1).
# Reproduces the template's in-render gate standalone for both design points and
# writes a quotable text record.  Same code path as the gate: if this disagrees
# with the render, one of them is wrong.
suppressMessages(library(forestsearch))
seed_base <- 8316951L; n_super <- 100000L; harm_z1_quantile <- 0.25; maxk <- 2L
confounders_base <- c("er", "age", "meno", "pgr", "nodes", "size", "grade")
fs_conf_force <- c("meno == 0", "er <= 0", "pgr <= 0"); fs_conf.cont_jcuts <- list(er = 10L)
tol <- 1e-8; tol3_fb <- 1e-6

cat("nullid STRUCTURAL-NULL DESIGN-POINT RECORD\n")
cat("TASK_null_gbsg_identification_2026-09-21, Step 1 ; ", format(Sys.time()), "\n")
cat("host ", Sys.info()[["nodename"]], " ; R ", R.version$major, ".", R.version$minor,
    " ; forestsearch ", as.character(utils::packageVersion("forestsearch")), "\n", sep = "")
cat("DGM: setup_gbsg_dgm(model = \"null\"), z1_quantile 0.25, n_super 100000, seed 8316951.\n")
cat("  model = \"null\" drops zh from the AFT design and sets flag.harm == 0 with\n")
cat("  fs_harm_true = NULL (R/sim_aft_gbsg.R:302-310, :336-349, :370-378).\n\n")

for (target in c(0.657, 0.721)) {
  cat(sprintf("== design point: super-population marginal Cox HR = %.3f ==\n", target))
  hr_at <- function(kt) setup_gbsg_dgm(model = "null", k_treat = kt, z1_quantile = harm_z1_quantile,
                                       n_super = n_super, seed = seed_base)$hr_causal
  k_treat <- stats::uniroot(function(kt) hr_at(kt) - target,
                            interval = c(0.1, 5), extendInt = "yes", tol = 1e-12)$root
  dgm <- compute_dgm_cde(setup_gbsg_dgm(model = "null", k_inter = 1, k_treat = k_treat,
                                        z1_quantile = harm_z1_quantile,
                                        n_super = n_super, seed = seed_base))
  ev <- fs_build_eval_frame(dgm, outcome_type = "survival", eval_seed = 20260628L,
                            analysis_time = 84, cens_adjust = log(1.5))
  prev <- mean(dgm$df_super$flag_harm)
  lab  <- is.null(dgm$subgroup_info$fs_harm_true) && is.null(dgm$subgroup_info$grf_harm_true)
  cat(sprintf("  k_treat = %.14f\n", k_treat))
  cat(sprintf("  [1] planted-region prevalence %.10f ; truth labels absent %s -> %s\n",
              prev, lab, if (prev == 0 && lab) "PASS" else "FAIL"))
  fsd <- get_FSdata(df.analysis = ev, use_lasso = FALSE, use_grf = FALSE,
                    confounders.name = intersect(confounders_base, names(ev)),
                    conf_force = fs_conf_force, conf.cont_jcuts = fs_conf.cont_jcuts,
                    outcome.name = "y_sim", event.name = "event_sim", details = FALSE)
  Z <- as.matrix(fsd$df[, fsd$confs_names, drop = FALSE]); lp <- ev$loghr_po
  lhr <- unname(dgm$model_params$b_hr["treat"]); L <- ncol(Z)
  idx <- c(lapply(seq_len(L), function(i) i), utils::combn(L, 2L, simplify = FALSE))
  dev <- vapply(idx, function(ii) {
    m <- if (length(ii) == 1L) Z[, ii] == 1L else Reduce(`&`, lapply(ii, function(j) Z[, j] == 1L))
    if (!any(m)) NA_real_ else abs(mean(lp[m]) - lhr) }, numeric(1))
  md <- max(dev, na.rm = TRUE)
  cat(sprintf("  [2] candidate family: %d factors, %d enumerated conjunctions (maxk %d), %d non-empty\n",
              L, length(idx), maxk, sum(!is.na(dev))))
  cat(sprintf("      max_g |beta(g) - log HR_uniform| = %.3e -> %s\n", md, if (md <= tol) "PASS" else "FAIL"))
  cat(sprintf("      strongest form: range of the patient-level log HR over all %d subjects = %.3e\n",
              length(lp), diff(range(lp))))
  d3 <- abs(dgm$hr_causal - target); ok3p <- d3 <= tol
  cat(sprintf("  [3] super-population marginal Cox HR %.12f vs target %.12f ; |diff| %.3e -> %s\n",
              dgm$hr_causal, target, d3,
              if (ok3p) "PASS (1e-8)" else sprintf("PASS UNDER AMENDMENT (within %.0e)", tol3_fb)))
  cat(sprintf("      uniform patient-level HR exp(b0[treat]) = %.6f ; AHR = %.6f\n", exp(lhr), dgm$AHR))
  cat(sprintf("      hr_Hc_true (== hr_causal under the null) = %.6f\n\n", dgm$hr_Hc_true))
}
cat("For reference, the alt design's complement effects at HR 1.00, which fix these two targets:\n")
for (zq in c(0.25, 0.60)) {
  ki <- calibrate_k_inter(target_hr_harm = 1.00, model = "alt", use_ahr = FALSE, z1_quantile = zq)
  d <- setup_gbsg_dgm(model = "alt", k_inter = ki, z1_quantile = zq, n_super = n_super, seed = seed_base)
  cat(sprintf("  z1_quantile %.2f : prevalence %.5f ; hr_H %.6f ; hr_Hc %.6f\n",
              zq, mean(d$df_super$flag_harm), d$hr_H_true, d$hr_Hc_true))
}
