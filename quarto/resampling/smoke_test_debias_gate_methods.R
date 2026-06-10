# =============================================================================
# smoke_test_debias_gate_methods.R
# -----------------------------------------------------------------------------
# Confirms the Tier-2 de-biased gate now populates out$debias_gate for ALL
# subgroup_method options (consistency / dina / grf-tree / grf-frontier), each
# de-biased against its OWN candidate family (approach (i)).
#
# Run from the package root after installing the modified source.  Uses
# devtools::install() (NOT load_all()): doFuture parallel workers spawn separate
# R processes that only see the *installed* package.
# =============================================================================

#devtools::install(quick = TRUE, upgrade = "never")
library(forestsearch)
suppressPackageStartupMessages(library(survival))

# -----------------------------------------------------------------------------
# 1. Minimal synthetic survival trial with a PLANTED harm subgroup
#    Harm region: er <= 30 & age >= 60  (treatment reverses to harm there).
#    n large enough that the ~17% harm region clears n.min for every method.
# -----------------------------------------------------------------------------
set.seed(8316951)
n         <- 1200
treat     <- rbinom(n, 1L, 0.5)
age       <- round(rnorm(n, 60, 10))
er        <- round(pmax(0, rnorm(n, 40, 25)))
nodes     <- rpois(n, 3)
biomarker <- rbinom(n, 1L, 0.4)

harm <- as.integer(er <= 30 & age >= 60)          # true subgroup
lp   <- 0.02 * (age - 60) - 0.01 * er + 0.15 * nodes +
        log(0.85) * treat +                       # overall benefit (HR < 1)
        log(2.4)  * treat * harm                  # reversed to harm (HR > 1)
t_event <- rexp(n, rate = 0.03 * exp(lp))
t_cens  <- rexp(n, rate = 0.02)
dat <- data.frame(
  id    = seq_len(n),
  tte   = pmin(t_event, t_cens),
  event = as.integer(t_event <= t_cens),
  treat = treat, age = age, er = er, nodes = nodes, biomarker = biomarker
)
conf <- c("age", "er", "nodes", "biomarker")
cat(sprintf("True harm subgroup size: %d (%.1f%%)\n",
            sum(harm), 100 * mean(harm)))

# -----------------------------------------------------------------------------
# 2. One forestsearch() call per method, gate ON.  Survival defaults apply:
#    outcome.name="tte", event.name="event", treat.name="treat", id.name="id".
#    adverse_outcome is left NULL: the gate forces harm-orientation (TRUE) for
#    survival internally, so it need not be set here.  sg_focus="hr" maps to the
#    gate re-selection rule "maxeff".
# -----------------------------------------------------------------------------
run_one <- function(method, grf_selection = "tree") {
  forestsearch(
    df.analysis            = dat,
    confounders.name       = conf,
    subgroup_method        = method,
    grf_selection          = grf_selection,
    grf_depth              = 2,
    sg_focus               = "hr",
    hr.threshold           = 1.25,
    pconsistency.threshold = 0.90,
    n.min                  = 50,
    fs.splits              = 200,                 # small: this is a smoke test
    debias_gate            = TRUE,
    debias_gate_args       = list(draws = 1000, gate = "point"),
    details                = FALSE
  )
}

fits <- list(
  consistency     = run_one("consistency"),
  dina            = run_one("dina"),
  `grf / tree`    = run_one("grf", "tree"),
  `grf / frontier`= run_one("grf", "frontier")
)

# -----------------------------------------------------------------------------
# 3. Inspect the gate.  Pass = out$debias_gate is non-NULL with a debiased CI.
#    (NULL is the legitimate "no subgroup found / gate skipped" outcome.)
# -----------------------------------------------------------------------------
show_gate <- function(fs, tag) {
  cat("\n== ", tag, " ==\n", sep = "")
  cat("  selected sg.harm: ",
      if (is.null(fs$sg.harm)) "none" else paste(fs$sg.harm, collapse = " & "),
      "\n", sep = "")
  g <- fs$debias_gate
  if (is.null(g)) { cat("  debias_gate = NULL\n"); return(invisible(FALSE)) }
  cat(sprintf("  family size (n_family) = %s   selected size (n_selected) = %s\n",
              g$n_family, g$n_selected))
  cat(sprintf("  measure = %s   selected_label = %s\n",
              g$measure, g$selected_label))
  if (!is.null(g$naive))
    cat(sprintf("  naive    est = %.3f  [%.3f, %.3f]\n",
                g$naive$est, g$naive$lower, g$naive$upper))
  if (!is.null(g$debiased))
    cat(sprintf("  debiased est = %.3f  [%.3f, %.3f]  (se_ij = %.3f, ij_source = %s)\n",
                g$debiased$est, g$debiased$lower, g$debiased$upper,
                g$debiased$se_ij, g$debiased$ij_source))
  cat(sprintf("  selection_bias = %.4f   selection_rate = %.3f   harm_flag = %s\n",
              g$selection_bias, g$selection_rate, g$harm_flag))
  invisible(!is.null(g$debiased))
}

cat("\n------------------------------------------------------------\n")
ok <- mapply(show_gate, fits, names(fits))
cat("\n------------------------------------------------------------\n")
cat("PASS (gate produced a debiased CI):\n")
print(ok)
