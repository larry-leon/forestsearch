# =============================================================================
# prerefactor_reference_2026-08-29.R -- the §1 pre-refactor reference for the
# fs_oc_grid task.  Run ONCE at 0.2.4 (HEAD 07aaab2a, before R/fs_oc_predict.R
# is refactored) to freeze fs_oc_predict()'s output at a fixed seed on both
# gates; the refactor guard (§5.1) re-runs the same calls and asserts
# identical().  Usage:
#   Rscript dev/glm-continuous-sims/prerefactor_reference_2026-08-29.R write
#   Rscript dev/glm-continuous-sims/prerefactor_reference_2026-08-29.R check
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
mode <- commandArgs(trailingOnly = TRUE)[1]
OUT <- "dev/glm-continuous-sims/prerefactor_reference_2026-08-29.rds"

make_oc_dgm <- function(N = 2000L, seed = 1L) {
  set.seed(seed)
  age <- round(stats::rnorm(N, 35, 9)); pre <- round(stats::rexp(N, 1 / 500))
  V   <- factor(stats::rbinom(N, 1L, 0.42), levels = 0:1)
  inQ <- as.integer(age > 34 & pre <= 745); mu0 <- 40 + 0.2 * age
  structure(list(df_super = data.frame(age = age, preanti = pre, V = V, mu0 = mu0,
                                       mu1 = mu0 - 26 - 14 * inQ, flag_harm = inQ),
                 outcome_type = "continuous", effect_measure = "MD",
                 model_params = list(sigma = 127.5)), class = c("glm_dgm", "list"))
}
piQ <- 0.34
hand <- structure(list(lab = c("Q", "P", "D"), Pg = c(piQ, .45, .31), PQg = c(1, .28/.45, 1),
  sens_g = c(1, .28/piQ, .31/piQ), spec_g = c(1, 1 - .17/(1 - piQ), 1),
  ovl = matrix(c(piQ, .28, .31, .28, .45, .28*.31/piQ, .31, .28*.31/piQ, .31), 3, 3),
  M = 3L, PQ = piQ), class = c("fs_oc_family", "list"))
hand$beta_g <- 26 + 14 * hand$PQg; hand$se_g <- 13.7 * sqrt(2) * sqrt(piQ / hand$Pg)

args <- list(confounders.name = c("age", "preanti", "V"),
             conf.cont_jcuts = list(age = 4, preanti = 4), n.min = 60,
             effect.threshold = 30, consistency.threshold = 10)
fam <- fs_oc_family_enumerate(make_oc_dgm(), args, n = 500)

strip <- function(p) { p$family <- NULL; p }
ref <- list(
  hand_resample = strip(fs_oc_predict(family = hand, n = 500, c1 = 30, c2 = 10,
                                      consistency_method = "resample", draws = 2e4, seed = 20260829)),
  hand_split    = strip(fs_oc_predict(family = hand, n = 500, c1 = 30, c2 = 10,
                                      consistency_method = "split", draws = 2e4, seed = 20260829)),
  fam_resample  = strip(fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
                                      consistency_method = "resample", draws = 2e4, seed = 20260829)),
  fam_split     = strip(fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
                                      consistency_method = "split", draws = 2e4, seed = 20260829)),
  fam_M = fam$M, pkg = as.character(utils::packageVersion("forestsearch")))

if (identical(mode, "write")) {
  saveRDS(ref, OUT); cat("reference written at", ref$pkg, ":", OUT, "(family M =", fam$M, ")\n")
} else {
  old <- readRDS(OUT)
  ok <- TRUE
  for (nm in c("hand_resample", "hand_split", "fam_resample", "fam_split")) {
    same <- identical(old[[nm]], ref[[nm]])
    cat(sprintf("%-14s identical: %s\n", nm, same)); ok <- ok && same
    if (!same) {
      for (q in c("det_rate", "EnH", "Esens", "Espec", "Eppv", "Enpv", "EbetaH", "Enaive_bias", "mass_below", "P1", "p_sel"))
        if (!identical(old[[nm]][[q]], ref[[nm]][[q]]))
          cat(sprintf("   %s differs: max|diff| = %.3e\n", q, max(abs(old[[nm]][[q]] - ref[[nm]][[q]]))))
      if (!identical(old[[nm]]$settings, ref[[nm]]$settings)) cat("   settings differ\n")
    }
  }
  cat("REFACTOR GUARD:", if (ok) "PASS (identical to the 0.2.4 reference)" else "FAIL", "\n")
  quit(status = if (ok) 0L else 1L)
}
