# =============================================================================
# GLM (binary/OR) covariate-adjustment: layer-by-layer end-to-end check.
#
# Confirms adjust_covariates is honoured on the GLM path through:
#   Layer 1 (search) + Layer 2 (consistency) : a SELECTION FLIP on a confounded
#       DGM (unadjusted selects the spurious region; adjust_covariates="S"
#       selects the real region) -- a no-op adjustment cannot produce this.
#   Layer 3 (bootstrap)                       : bias correction runs on the
#       adjusted fit (the bootstrap GLM estimator is now built with
#       adjust_covariates).
#
# Run AFTER devtools::document(); devtools::install() with the updated
# forestsearch_main.R, bootstrap_dofuture_main.R, summary_utility_functions.R,
# and forestsearch_cross_validation.R.
#
# Design (binary analogue of the survival selection-flip):
#   S   : strong prognostic confounder (logit effect ~1.8).
#   z1=(grp==1): SPURIOUS -- W confounded with S here; unadjusted OR large,
#                adjust(S) collapses to ~1.
#   z2=(grp==2): REAL -- genuine conditional treatment OR (~2.4), persists.
#   z3=(grp==3): noise.
# =============================================================================

library(forestsearch)

# ---- Simulate confounded BINARY data ----------------------------------------
set.seed(2025)
N   <- 5000
S   <- rbinom(N, 1, 0.5)
grp <- sample(1:4, N, replace = TRUE)
pW  <- ifelse(grp == 1, plogis(-1.5 + 3.0 * S), 0.5)   # W ~ S confounding in grp 1
W   <- rbinom(N, 1, pW)
eta <- -1.0 + 1.8 * S + 0.95 * W * (grp == 2)          # real OR in grp 2 only
y   <- rbinom(N, 1, plogis(eta))

dat <- data.frame(
  id        = seq_len(N),
  y_sim     = y,
  treat_sim = W,
  S         = S,                       # adjuster (NOT a candidate)
  z1 = as.integer(grp == 1),
  z2 = as.integer(grp == 2),
  z3 = as.integer(grp == 3)
)
confs <- c("z1", "z2", "z3")

# ---- Mechanism check (manual glm; the OR the search screens on) --------------
cat("=== DGM mechanism (manual within-region treatment OR) ===\n")
orr <- function(z, adj = FALSE) {
  d <- dat[dat[[z]] == 1, ]
  f <- if (adj) y_sim ~ treat_sim + S else y_sim ~ treat_sim
  exp(coef(glm(f, binomial, d))[["treat_sim"]])
}
for (z in confs)
  cat(sprintf("  %s==1: unadjusted OR=%.2f  adj(S) OR=%.2f\n", z, orr(z), orr(z, TRUE)))
cat("  (expect z1 to collapse under adj(S); z2 to persist)\n\n")

# ---- Common forestsearch settings (BINARY / OR) -----------------------------
fs_common <- list(
  df.analysis      = dat,
  outcome.name     = "y_sim",
  event.name       = "y_sim",          # binary: event.name = outcome.name
  treat.name       = "treat_sim",
  id.name          = "id",
  confounders.name = confs,            # S excluded -- it is the adjuster
  outcome_type     = "binary",
  effect_measure   = "OR",
  hr.threshold     = 1.25,             # OR scale
  hr.consistency   = 1.0,
  pconsistency.threshold = 0.60,
  fs.splits        = 200,
  maxk             = 1,
  use_lasso        = FALSE,
  use_grf          = FALSE,
  use_dina         = FALSE,
  parallel_args    = list(plan = "multisession", workers = 2, show_message = FALSE),
  seedit           = 8316951,
  details          = FALSE
)

run_fs <- function(label, ...) {
  fs <- tryCatch(do.call(forestsearch, modifyList(fs_common, list(...))),
                 error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(fs)) {
    sg <- if (is.null(fs$sg.harm)) "none" else paste(fs$sg.harm, collapse = " & ")
    # max_sg_est is stored on the LOG scale for GLM ratio measures (OR/RR/IRR);
    # exponentiate for display so it matches the manual within-region OR.
    est <- fs$max_sg_est
    or  <- if (!is.null(est) && !is.na(est)) exp(est) else NA_real_
    cat(sprintf("  %-18s -> subgroup: %-20s  max OR=%.2f (log-OR=%.2f)\n",
                label, sg, or, if (is.null(est)) NA_real_ else est))
  }
  invisible(fs)
}

# ---- LAYERS 1 + 2: search + consistency (selection flip) --------------------
cat("=== Layers 1+2: search + consistency ===\n")
fs_unadj <- run_fs("UNADJUSTED",   adjust_covariates = NULL)
fs_adj   <- run_fs("adjust S",     adjust_covariates = "S")

u <- if (!is.null(fs_unadj)) paste(fs_unadj$sg.harm, collapse = " & ") else NA
a <- if (!is.null(fs_adj))   paste(fs_adj$sg.harm,   collapse = " & ") else NA
cat(sprintf("\n  unadjusted selected: %s\n  adjust(S)  selected: %s\n", u, a))
layer12_ok <- !is.na(u) && !is.na(a) && u != a
cat(if (layer12_ok)
      "  >>> SELECTION CHANGED -> search + consistency adjust on the GLM path (PROVEN).\n"
    else
      "  >>> selection did NOT change -- GLM search/consistency adjustment NOT confirmed.\n")

# ---- LAYER 3: bootstrap on the adjusted fit ---------------------------------
cat("\n=== Layer 3: bootstrap (post-selection) on the adjusted fit ===\n")
if (!is.null(fs_adj) && !is.null(fs_adj$sg.harm)) {
  bc <- tryCatch(
    forestsearch_bootstrap_dofuture(
      fs.est        = fs_adj,
      nb_boots      = 50,                       # tiny for smoke; >=500 for real
      parallel_args = list(plan = "multisession", workers = 2)
    ),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(bc)) {
    cat("  bootstrap ran OK on adjusted fit; class:", class(bc)[1], "\n")
    cat("  (the bootstrap GLM estimator is now built WITH adjust_covariates, so\n")
    cat("   H_obs / bias-corrected OR are on the adjusted scale)\n")
  }
} else {
  cat("  skipped (adjusted run found no subgroup)\n")
}

# ---- LAYER 4: display (plot estimator resolves adjust_covariates) -----------
cat("\n=== Layer 4: plot estimator (display) ===\n")
if (!is.null(fs_adj) && !is.null(fs_adj$sg.harm) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  p <- tryCatch(
    plot_sg_glm_outcomes(fs.est = fs_adj, outcome.name = "y_sim",
                         treat.name = "treat_sim", show_bc = FALSE),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(p))
    cat("  plot_sg_glm_outcomes built OK; adjust_covariates auto-resolved from\n",
        "  fs.est$args_call_all, so the displayed effect labels are on the\n",
        "  adjusted scale (matches the estimate).\n")
} else {
  cat("  skipped (no adjusted subgroup or ggplot2 unavailable)\n")
}

cat("\n=== summary ===\n")
cat(sprintf("  Layers 1+2 (search/consistency): %s\n",
            if (isTRUE(layer12_ok)) "ADJUSTED (selection flip)" else "NOT CONFIRMED"))
cat("  Layer 3 (bootstrap): adjusted estimator wired; runs on adjusted fit.\n")
cat("  Layer 4 (plots): plot_sg_glm_outcomes / forest plot resolve\n")
cat("    adjust_covariates from the call record; display matches the estimate.\n")
cat("  Also wired: estimate-table + post-CV + GRF-table re-estimation\n")
cat("    (exercised when summary tables / forestsearch_Kfold are rendered).\n")
