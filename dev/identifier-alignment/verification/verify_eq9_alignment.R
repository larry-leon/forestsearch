# verify_eq9_alignment.R
#
# V8 -- does aligning fully with Eq. 9 restore the centering condition by
# itself?
#
# Eq. 9 divides BOTH bias terms by B.  A draw admitting no candidate therefore
# contributes ZERO to bias_sel rather than being removed from the average.
# Read literally, that also fixes what the residual is on such a draw:
#
#     r_b = (bias_sel + bias_fix) - D_Hstar_b(b) - D_H(b),
#     with D_Hstar_b(b) := 0 when draw b admitted nobody.
#
# D_H(b) is always defined (H is the observed selection), so r_b is defined on
# every draw.  If the IJ is then evaluated over ALL B draws rather than the
# winner subset, the arithmetic gives
#
#     mean_b(r_b) = (bias_sel + bias_fix) - bias_sel - bias_fix = 0  EXACTLY,
#
# which is precisely the condition @sec-wager's bridge assumes.  If that holds,
# the manuscript needs no amendment at all: Eq. 13's uncentered rbar2 IS
# Wager's centered v-hat, and the gap was never in the theory -- only in an
# implementation that restricted the IJ to the winner subset while Eq. 9
# averaged over B.
#
# Three arms compared:
#   (i)   AS CODED       bias_sel over winners; IJ over winner draws
#   (ii)  HALF-ALIGNED   bias_sel over B;       IJ over winner draws
#   (iii) FULLY ALIGNED  bias_sel over B;       IJ over ALL B draws
#
# .fs_mr_ij_var() is the package's own, bound directly from the installed
# package rather than copied -- so this script REQUIRES forestsearch
# installed, and is otherwise base R.  A snapshot would be free to drift
# from the shipped function; see the sibling verify_residual_centering.R.

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

fixed_bias <- mean(P[sel, ])                       # over all B, both readings
bs_winners <- mean(sel_bias, na.rm = TRUE)         # as coded
bs_allB    <- sum(sel_bias, na.rm = TRUE) / draws  # Eq. 9

# sel_bias with no-winner draws set to 0, per Eq. 9's implicit convention
sel_bias0 <- ifelse(is.finite(sel_bias), sel_bias, 0)

arms <- list(
  "as coded"      = list(bs = bs_winners, sb = sel_bias,  idx = ok),
  "half-aligned"  = list(bs = bs_allB,    sb = sel_bias,  idx = ok),
  "fully aligned" = list(bs = bs_allB,    sb = sel_bias0, idx = seq_len(draws))
)

cat(sprintf("selection_rate p = %.4f   (%d of %d draws had a winner)\n\n",
            p, length(ok), draws))
cat(sprintf("  %-15s %12s %14s %14s %14s\n",
            "arm", "mean(r)", "tilde_V", "hat_V coded", "hat_V Wager"))

res <- list()
for (nm in names(arms)) {
  a  <- arms[[nm]]
  r  <- (a$bs + fixed_bias) - a$sb - P[sel, ]
  ij <- .fs_mr_ij_var(Xi, r, a$idx)
  rb <- r[a$idx]
  hatV_wager <- ij$tilde_V - (n / length(a$idx)) * mean((rb - mean(rb))^2)
  res[[nm]] <- c(rbar = mean(rb), tilde_V = ij$tilde_V,
                 hat_V = ij$hat_V, hat_V_w = hatV_wager)
  cat(sprintf("  %-15s %+12.3e %14.6f %14.6f %14.6f\n",
              nm, mean(rb), ij$tilde_V, ij$hat_V, hatV_wager))
}

cat("\n--- the claim under test ---\n")
rb_full <- res[["fully aligned"]][["rbar"]]
cat(sprintf("fully aligned mean(r)          %+.3e   %s\n", rb_full,
            if (abs(rb_full) < 1e-12) "ZERO to machine precision" else "NOT zero"))
cat(sprintf("fully aligned: coded == Wager? %s  (diff %.3e)\n",
            abs(res[["fully aligned"]][["hat_V"]] -
                res[["fully aligned"]][["hat_V_w"]]) < 1e-12,
            res[["fully aligned"]][["hat_V"]] - res[["fully aligned"]][["hat_V_w"]]))
cat(sprintf("as coded:      coded == Wager? %s  (diff %.3e)\n",
            abs(res[["as coded"]][["hat_V"]] -
                res[["as coded"]][["hat_V_w"]]) < 1e-12,
            res[["as coded"]][["hat_V"]] - res[["as coded"]][["hat_V_w"]]))

cat("\n--- point estimate consequences ---\n")
cat(sprintf("bias_sel  as coded (winners)   %.6f\n", bs_winners))
cat(sprintf("bias_sel  Eq. 9    (over B)    %.6f   ratio %.4f = p\n",
            bs_allB, bs_allB / bs_winners))
cat(sprintf("point estimate moves by        %+.6f\n", bs_winners - bs_allB))
