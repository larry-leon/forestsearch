# verify_conditional_convention.R
#
# V9 -- the option that was missing from the F13 analysis.
#
# F13 was framed as a choice between two conventions:
#
#   (i)   AS CODED    bias_sel over winner draws; bias_fix over B; IJ on winners
#   (iii) EQ. 9       both over B (D := 0 on no-winner draws); IJ over all B
#
# (i) is incoherent -- the two bias terms carry different denominators, so the
# residual mixes differently normalised quantities and mean(r) != 0.  That is
# the actual defect.  But (iii) is not the only coherent repair.  There is a
# third:
#
#   (iv)  CONDITIONAL both bias terms over WINNER draws; IJ on winner draws
#
# which is equally coherent, ALSO gives mean(r) = 0 exactly, and differs from
# the shipped code by ONE line: bias_fix averaged over the winner draws rather
# than over all B.
#
# The two repairs target different estimands:
#
#   (iii) unconditional -- averages the selection perturbation over all draws,
#         counting a no-winner draw as contributing zero selection bias.
#   (iv)  conditional on identification -- the correction is estimated on the
#         draws where a selection actually occurred, mirroring the reported
#         analysis, which exists only because a subgroup WAS identified.
#
# This script checks that (iv) is centered, and reports what each does to the
# variance.  .fs_mr_ij_var() is the package's own, bound directly from the
# installed package rather than copied -- so this script REQUIRES forestsearch
# installed, and is otherwise base R.  A snapshot would be free to drift from
# the shipped function; see the sibling verify_residual_centering.R.

.fs_mr_ij_var <- forestsearch:::.fs_mr_ij_var   # the SHIPPED function, not a copy

set.seed(20260805)

n <- 300L; draws <- 6000L; n_cand <- 3L; sel <- 1L
bh <- c(0.30, 0.26, 0.22); floor_t <- 0.42

db <- matrix(rnorm(n * n_cand), n, n_cand) / sqrt(n)
Xi <- matrix(rpois(n * draws, 1) - 1, n, draws)
P  <- crossprod(db, Xi)

beta_star <- bh + P
admitted  <- beta_star >= floor_t
sel_bias  <- rep(NA_real_, draws)
for (b in seq_len(draws)) {
  cand <- which(admitted[, b])
  if (!length(cand)) next
  sel_bias[b] <- P[cand[which.max(beta_star[cand, b])], b]
}
ok <- which(is.finite(sel_bias))
p  <- length(ok) / draws
sel_bias0 <- ifelse(is.finite(sel_bias), sel_bias, 0)

arms <- list(
  "(i)   as coded"   = list(bs = mean(sel_bias, na.rm = TRUE),
                            bf = mean(P[sel, ]),
                            sb = sel_bias,  idx = ok),
  "(iii) Eq. 9 lit." = list(bs = mean(sel_bias0),
                            bf = mean(P[sel, ]),
                            sb = sel_bias0, idx = seq_len(draws)),
  "(iv)  conditional" = list(bs = mean(sel_bias[ok]),
                             bf = mean(P[sel, ok]),
                             sb = sel_bias,  idx = ok)
)

cat(sprintf("selection_rate p = %.4f  (%d of %d draws admitted a candidate)\n\n",
            p, length(ok), draws))
cat(sprintf("  %-18s %12s %12s %13s %12s %12s\n",
            "arm", "bias_sel", "bias_fix", "mean(r)", "tilde_V", "se_ij"))

out <- list()
for (nm in names(arms)) {
  a  <- arms[[nm]]
  r  <- (a$bs + a$bf) - a$sb - P[sel, ]
  ij <- .fs_mr_ij_var(Xi, r, a$idx)
  rb <- r[a$idx]
  out[[nm]] <- c(bs = a$bs, rbar = mean(rb), tV = ij$tilde_V, hV = ij$hat_V)
  cat(sprintf("  %-18s %12.6f %12.6f %+13.3e %12.6f %12.6f\n",
              nm, a$bs, a$bf, mean(rb), ij$tilde_V,
              if (is.finite(ij$hat_V) && ij$hat_V > 0) sqrt(ij$hat_V) else NA))
}

cat("\n--- centering ---\n")
for (nm in names(arms)) {
  cat(sprintf("  %-18s mean(r) = %+.3e   %s\n", nm, out[[nm]][["rbar"]],
              if (abs(out[[nm]][["rbar"]]) < 1e-12) "CENTERED" else "not centered"))
}

cat("\n--- variance, relative to the shipped code ---\n")
base_tV <- out[["(i)   as coded"]][["tV"]]
for (nm in names(arms)) {
  cat(sprintf("  %-18s tilde_V %8.4f   %+7.1f%%\n", nm, out[[nm]][["tV"]],
              100 * (out[[nm]][["tV"]] / base_tV - 1)))
}

cat("\n--- point estimate: bias_sel relative to shipped ---\n")
base_bs <- out[["(i)   as coded"]][["bs"]]
for (nm in names(arms)) {
  cat(sprintf("  %-18s bias_sel %8.5f   ratio %.4f\n", nm,
              out[[nm]][["bs"]], out[[nm]][["bs"]] / base_bs))
}
cat(sprintf("\n  (iv) leaves bias_sel UNCHANGED from the shipped value;\n"))
cat(sprintf("  (iii) multiplies it by p = %.4f.\n", p))
cat("\n  So (iv) repairs the incoherence without moving the point estimate,\n")
cat("  and (iii) moves it. The two differ in estimand, not in correctness.\n")
