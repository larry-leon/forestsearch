# =============================================================================
# Definitive end-to-end check: does adjust_covariates change WHICH subgroup
# forestsearch selects?  (Proves the adjusted model is actually used in the
# search/consistency engine -- not a silent no-op.)
#
# Design (Simpson's-paradox style confounding):
#   S   : a STRONG prognostic stratum factor (HR ~ 7).
#   grp : four disjoint regions (no overlap leak).
#     z1 = (grp==1): SPURIOUS harm -- W is confounded with S only here, so the
#          UNADJUSTED within-subgroup HR is large (~2.6) but it is an artefact
#          of S imbalance; stratifying by S collapses it to ~1.0.
#     z2 = (grp==2): REAL harm -- a genuine conditional treatment effect
#          (HR ~2.3) that PERSISTS under strata(S).
#     z3 = (grp==3): noise (HR ~1.0).
#     grp==4: implicit baseline (not a candidate; avoids collinearity).
#
# Expected result:
#   - UNADJUSTED forestsearch selects the z1 (spurious) region.
#   - adjust_covariates = "strata(S)" selects the z2 (real) region.
#   The SELECTION CHANGES -> the adjustment is provably flowing through the
#   search and consistency layers.
#
# Run AFTER devtools::document(); devtools::install()
# =============================================================================

library(forestsearch)
library(survival)

# ---- Simulate confounded survival data --------------------------------------
set.seed(2025)
N   <- 4000
S   <- rbinom(N, 1, 0.5)                     # strong prognostic stratum factor
grp <- sample(1:4, N, replace = TRUE)        # disjoint regions

pW <- ifelse(grp == 1, plogis(-1.5 + 3.0 * S), 0.5)  # W ~ S confounding in grp 1 only
W  <- rbinom(N, 1, pW)

lp   <- 2.0 * S + 0.90 * W * (grp == 2)      # real conditional harm in grp 2
Tlat <- rexp(N, rate = 0.04 * exp(lp))
Cens <- rexp(N, rate = 0.03)

dat <- data.frame(
  id    = seq_len(N),
  tte   = pmin(Tlat, Cens),
  event = as.integer(Tlat <= Cens),
  treat = W,
  S     = S,                                 # adjuster (NOT a candidate)
  z1    = as.integer(grp == 1),              # spurious region
  z2    = as.integer(grp == 2),              # real region
  z3    = as.integer(grp == 3)               # noise
)
confs <- c("z1", "z2", "z3")

# ---- Mechanism check (manual coxph -- the statistic the search screens on) ---
cat("=== DGM mechanism (manual within-region treatment HR) ===\n")
hr <- function(z, adj = FALSE) {
  d <- dat[dat[[z]] == 1, ]
  f <- if (adj) Surv(tte, event) ~ treat + strata(S) else Surv(tte, event) ~ treat
  exp(coef(coxph(f, data = d))[["treat"]])
}
for (z in confs)
  cat(sprintf("  %s==1: unadjusted HR=%.2f  strata(S) HR=%.2f\n",
              z, hr(z), hr(z, TRUE)))
cat("  (expect z1 to collapse under strata(S); z2 to persist)\n\n")

# ---- Common forestsearch settings -------------------------------------------
fs_common <- list(
  df.analysis      = dat,
  outcome.name     = "tte",
  event.name       = "event",
  treat.name       = "treat",
  id.name          = "id",
  confounders.name = confs,           # S deliberately excluded (it is the adjuster)
  hr.threshold     = 1.25,
  hr.consistency   = 1.0,
  pconsistency.threshold = 0.60,
  fs.splits        = 200,
  maxk             = 1,               # single-factor subgroups -> clean selection
  use_lasso        = FALSE,
  use_grf          = FALSE,
  parallel_args    = list(plan = "multisession", workers = 2, show_message = FALSE),
  seedit           = 8316951,
  details          = FALSE
)

run_fs <- function(label, ...) {
  fs <- tryCatch(do.call(forestsearch, modifyList(fs_common, list(...))),
                 error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(fs)) {
    sg <- if (is.null(fs$sg.harm)) "none" else paste(fs$sg.harm, collapse = " & ")
    cat(sprintf("  %s -> subgroup: %-28s  max HR=%.2f\n",
                label, sg,
                if (is.null(fs$max_sg_est)) NA_real_ else fs$max_sg_est))
  }
  invisible(fs)
}

cat("=== forestsearch selection ===\n")
fs_unadj  <- run_fs("UNADJUSTED        ")
fs_strata <- run_fs("strata(S)         ", adjust_covariates = "strata(S)")
fs_linear <- run_fs("linear S          ", adjust_covariates = "S")

cat("\n=== VERDICT ===\n")
u <- if (!is.null(fs_unadj))  paste(fs_unadj$sg.harm,  collapse = " & ") else NA
a <- if (!is.null(fs_strata)) paste(fs_strata$sg.harm, collapse = " & ") else NA
cat(sprintf("  unadjusted selected: %s\n  strata(S)  selected: %s\n", u, a))
if (!is.na(u) && !is.na(a) && u != a) {
  cat("  >>> SELECTION CHANGED -> adjustment is flowing end-to-end (PROVEN).\n")
} else {
  cat("  >>> selection did NOT change -- investigate whether adjustment reached\n")
  cat("      the scorer (try details = TRUE to see the adjustment confirmation).\n")
}
