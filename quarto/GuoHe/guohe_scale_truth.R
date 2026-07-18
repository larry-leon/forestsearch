# guohe_scale_truth.R
#
# Post-hoc truth recovery for the enumerated-family (scale) cell of
# `test_guohe_algorithm3.qmd`. Run AFTER the full-tier render; nothing needs to
# be added to that document beforehand.
#
# Why this is recoverable after the fact: the scale cell's DGM is seeded
# (set.seed(9090)) and `print.guohe_a3` reports the selected candidate. Those
# two facts are sufficient to (a) regenerate the exact numeric cut points, and
# (b) rebuild the specific conjunction that won.
#
# What "truth" means here. The target is the least-false MARGINAL odds ratio of
# the *selected* rule -- beta(H-hat) -- not of the generating harm subgroup.
# The winner is a data-driven approximation to the true rule and will generally
# differ from it, so conditioning on the realized selection is the correct
# comparison, and is the estimand the MR/FB/Guo-He comparison is about.
# Because logistic regression is non-collapsible, this marginal OR also differs
# from the conditional exp(1.1) built into the DGM; that is expected, not error.
#
# Usage: edit the REPORTED block below with the numbers printed by the full-tier
# render, then `Rscript guohe_scale_truth.R`.

suppressWarnings(suppressMessages(library(survival)))

# ---- REPORTED: paste from the full-tier scale-cell output ------------------
SELECTED <- "z1.gt.q70_AND_z2.gt.q40"   # `selected` line
NAIVE    <- NA_real_                    # naive odds ratio
DEBIASED <- NA_real_                    # de-biased odds ratio
BOUND    <- NA_real_                    # one-sided lower bound  [CERTIFIED]
CI_LO    <- NA_real_                    # two-sided CI lower     [EXTENSION]
CI_HI    <- NA_real_                    # two-sided CI upper
# ---------------------------------------------------------------------------

N_TRUTH <- 500000L   # population size for the large-n truth
QGRID   <- c(.3, .4, .5, .6, .7, .8)

# ---- 1. regenerate the scale cell's sample to recover the CUT POINTS -------
# The rule's thresholds are sample quantiles of the original n = 1200 draw, so
# the numeric values must come from that exact sample, not from the population.
set.seed(9090)
n <- 1200
trt_s <- rbinom(n, 1, 0.5)
z_s   <- matrix(rnorm(n * 4), n, 4)
harm_s <- as.integer(z_s[, 1] > 0.5 & z_s[, 2] > 0)
lp_s   <- -0.4 + 0.3 * z_s[, 3] + trt_s * (-0.25 + 1.1 * harm_s)
y_s    <- rbinom(n, 1, plogis(lp_s))

cut_value <- function(j, qpct) {
  q <- QGRID[which(round(100 * QGRID) == qpct)]
  if (!length(q)) stop("unrecognised quantile tag: q", qpct)
  unname(stats::quantile(z_s[, j], q))
}

# ---- 2. parse the selected rule -------------------------------------------
parse_rule <- function(nm) {
  parts <- strsplit(nm, "_AND_", fixed = TRUE)[[1]]
  lapply(parts, function(p) {
    m <- regmatches(p, regexec("^z([0-9]+)\\.gt\\.q([0-9]+)$", p))[[1]]
    if (length(m) != 3L) stop("cannot parse component: ", p)
    list(j = as.integer(m[2]), qpct = as.integer(m[3]),
         thresh = cut_value(as.integer(m[2]), as.integer(m[3])))
  })
}
rule <- parse_rule(SELECTED)
cat("Selected rule:", SELECTED, "\n")
for (cmp in rule) {
  cat(sprintf("  z%d > %.4f   (sample q%02d of the n=1200 draw)\n",
              cmp$j, cmp$thresh, cmp$qpct))
}
apply_rule <- function(z) {
  keep <- rep(TRUE, nrow(z))
  for (cmp in rule) keep <- keep & (z[, cmp$j] > cmp$thresh)
  keep
}
cat(sprintf("  prevalence in the original sample: %.3f (%d of %d)\n\n",
            mean(apply_rule(z_s)), sum(apply_rule(z_s)), n))

# ---- 3. large-n truth for the SELECTED rule --------------------------------
set.seed(20260718)
trt <- rbinom(N_TRUTH, 1, 0.5)
z   <- matrix(rnorm(N_TRUTH * 4), N_TRUTH, 4)
harm <- as.integer(z[, 1] > 0.5 & z[, 2] > 0)
lp   <- -0.4 + 0.3 * z[, 3] + trt * (-0.25 + 1.1 * harm)
y    <- rbinom(N_TRUTH, 1, plogis(lp))

sel <- apply_rule(z)
fit <- stats::glm(y[sel] ~ trt[sel], family = stats::binomial())
truth_or <- unname(exp(stats::coef(fit)[2]))
mcse <- unname(summary(fit)$coefficients[2, 2])

# reference points for interpretation
fit_true <- stats::glm(y[harm == 1] ~ trt[harm == 1], family = stats::binomial())
truth_or_true_rule <- unname(exp(stats::coef(fit_true)[2]))

cat(sprintf("TRUE marginal OR of the SELECTED rule : %.4f  (log-scale SE %.4f at N=%d)\n",
            truth_or, mcse, N_TRUTH))
cat(sprintf("  -- of the GENERATING harm rule       : %.4f  (z1>0.5 & z2>0)\n",
            truth_or_true_rule))
cat(sprintf("  -- conditional OR built into the DGM : %.4f  (exp(1.1); non-collapsible, so larger)\n",
            exp(1.1)))
cat(sprintf("  selected-rule prevalence in population: %.3f\n\n", mean(sel)))

# ---- 4. compare against the reported run ----------------------------------
if (!is.na(NAIVE) && !is.na(DEBIASED)) {
  cat("Comparison against the full-tier run\n")
  cat(sprintf("  truth      %.4f\n", truth_or))
  cat(sprintf("  naive      %.4f   (bias %+.4f)\n", NAIVE, NAIVE - truth_or))
  cat(sprintf("  de-biased  %.4f   (bias %+.4f)\n", DEBIASED, DEBIASED - truth_or))
  red <- 100 * (1 - abs(DEBIASED - truth_or) / abs(NAIVE - truth_or))
  cat(sprintf("  bias reduction: %.1f%%%s\n", red,
              if (sign(DEBIASED - truth_or) != sign(NAIVE - truth_or))
                "  (NOTE: sign flipped -- over-corrected past truth)" else ""))
  if (!is.na(BOUND)) {
    cat(sprintf("  one-sided lower bound %.4f covers truth? %s\n",
                BOUND, if (BOUND <= truth_or) "YES" else "NO"))
  }
  if (!is.na(CI_LO) && !is.na(CI_HI)) {
    cat(sprintf("  two-sided CI (%.4f, %.4f) covers truth? %s\n",
                CI_LO, CI_HI,
                if (CI_LO <= truth_or && truth_or <= CI_HI) "YES" else "NO"))
  }
} else {
  cat("REPORTED block not filled in; truth computed only.\n")
}

# ---- 5. full per-candidate resample drop-rate ------------------------------
# `print.guohe_a3` shows only the count of non-zero rates and the worst one.
# Estimability depends solely on membership and treatment assignment, never on a
# model fit, so the whole vector is recoverable in seconds without re-running
# the expensive cell.
cat("\nPer-candidate resample drop-rate (recomputed; min_events = 5)\n")
base <- list(); nm <- character(0)
for (j in 1:4) for (q in QGRID) {
  base[[length(base) + 1L]] <- as.integer(z_s[, j] > stats::quantile(z_s[, j], q))
  nm <- c(nm, sprintf("z%d.gt.q%02d", j, round(100 * q)))
}
base <- as.data.frame(base); names(base) <- nm
prs <- list(); pnm <- character(0); ix <- utils::combn(ncol(base), 2)
for (k in seq_len(ncol(ix))) {
  v <- as.integer(base[[ix[1, k]]] == 1 & base[[ix[2, k]]] == 1)
  if (sum(v) >= 60) {
    prs[[length(prs) + 1L]] <- v
    pnm <- c(pnm, paste0(names(base)[ix[1, k]], "_AND_", names(base)[ix[2, k]]))
  }
}
fam <- cbind(base, as.data.frame(setNames(prs, pnm)))

set.seed(11)
NDRAW <- 400L
idx_list <- replicate(NDRAW, sample.int(n, n, replace = TRUE), simplify = FALSE)
drop_rate <- vapply(fam, function(memb) {
  mean(vapply(idx_list, function(ii) {
    s <- memb[ii] == 1L
    if (sum(s) < 10L) return(TRUE)
    tt <- table(trt_s[ii][s])
    length(tt) < 2L || min(tt) < 5L
  }, logical(1)))
}, numeric(1))

cat(sprintf("  candidates: %d | non-zero drop rate: %d | max %.3f\n",
            length(drop_rate), sum(drop_rate > 0), max(drop_rate)))
worst <- sort(drop_rate[drop_rate > 0], decreasing = TRUE)
if (length(worst)) {
  cat("  worst offenders (size, rate):\n")
  for (nmw in names(utils::head(worst, 8))) {
    cat(sprintf("    %-40s n=%4d  %.3f\n", nmw, sum(fam[[nmw]]), drop_rate[[nmw]]))
  }
  cat("\n  Non-zero rates mean the bootstrap maximum was taken over a SHRINKING\n",
      "  competitor set in those resamples, which UNDER-states the correction\n",
      "  (bias toward over-optimism -- the direction the method exists to prevent).\n", sep = "")
} else {
  cat("  no candidate dropped out; competitor set stable across resamples.\n")
}
