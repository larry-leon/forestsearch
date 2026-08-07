# The exact beta(Hhat) for the package's continuous DGM, and one correction.
#
# Read from R/generate_glm_dgm.R (continuous branch, lines 414-427) and
# R/simulate_from_glm_dgm.R (lines 97-120), the DGM is exactly:
#
#   mu0_i                                        (arbitrary prognostic function)
#   treat_effect = mu1 - mu0 = beta_treat        CONSTANT (base fit has no
#                                                treat:covariate interactions)
#   mu1_i = mu0_i + k_treat * beta_treat + beta_inter * 1{i in Q}
#   treat_sim ~ Bernoulli(p)                     independent of X
#   y_sim     = mu_obs + N(0, sigma)
#
# So the subject-level effect is TWO-VALUED:
#   tau_i = delta                 outside Q,   delta := k_treat * beta_treat
#   tau_i = delta + beta_inter    inside  Q
#
# THE CORRECTION.  The handoff states
#       beta(Hhat) = tau * P(H_true | Hhat),
# which sets the effect outside H_true to ZERO.  In this DGM it is delta, the
# fitted ACTG175 treatment effect, which is not zero -- build_actg175_glm_dgm.R
# says in as many words that "the COMPLEMENT inherits the fitted ACTG175
# treatment effect".  The correct target is
#       beta(Hhat) = delta + beta_inter * P(Q | Hhat).
#
# Consequences: the handoff's value is low by delta for EVERY rule, and its
# "disjoint from truth -> beta = 0 exactly" row is wrong -- a disjoint rule has
# beta = delta.
#
# Tests below: (A) the exact form matches a fitted MD; (B) the handoff's form
# fails, by exactly delta; (C) delta is already available as the DGM's own
# complement effect, so nothing new needs computing.

set.seed(20260806)

n_super <- 25000L
n_trial <- 400000L      # large so the fitted MD has small MC error

# --- Stand in for df_super: correlated covariates, arbitrary prognostic mu0 ---
z1 <- rnorm(n_super)
z2 <- 0.7 * z1 + sqrt(1 - 0.49) * rnorm(n_super)
x1 <- pnorm(z1)                       # ~ age
x2 <- pnorm(z2)                       # ~ preanti  (correlated with age)
mu0 <- 20 + 35 * x1^2 - 18 * x2 + 9 * sin(5 * x1 * x2)

# Package structure
k_treat    <- 1
beta_treat <- 6.5                     # fitted ACTG175 effect -> delta != 0
beta_inter <- 40                      # planted subgroup shift (template's +40)
delta      <- k_treat * beta_treat

flag_harm <- as.integer(x1 > 0.6 & x2 <= 0.55)     # Q = {z1=1 & z2=1} analogue
mu1 <- mu0 + delta + beta_inter * flag_harm

sigma <- 8

df_super <- data.frame(x1, x2, mu0, mu1, flag_harm)

# --- Draw a trial exactly as simulate_from_glm_dgm() does --------------------
idx       <- sample(n_super, n_trial, replace = TRUE)
d         <- df_super[idx, ]
treat_sim <- rbinom(n_trial, 1L, 0.5)
mu_obs    <- ifelse(treat_sim == 1L, d$mu1, d$mu0)
y_sim     <- mu_obs + rnorm(n_trial, sd = sigma)

# --- Realized rules, including shapes outside the handoff's stated coverage ---
rules <- list(
  "exact rule          {x1>.6} & {x2<=.55}" = with(d, x1 > 0.6  & x2 <= 0.55),
  "too narrow          {x1>.75} & {x2<=.4}" = with(d, x1 > 0.75 & x2 <= 0.40),
  "too wide            {x1>.4} & {x2<=.8}"  = with(d, x1 > 0.40 & x2 <= 0.80),
  "one cut only        {x1>.6}"             = with(d, x1 > 0.6),
  "DISJUNCTION         {x1>.8} | {x2<.15}"  = with(d, x1 > 0.8  | x2 <  0.15),
  "disjoint from truth {x1<.3}"             = with(d, x1 < 0.30)
)

cat("== A/B: exact target vs the handoff's form, against a fitted MD ==\n")
cat(sprintf("  %-42s %9s %9s %8s %9s %8s\n",
            "rule", "fitted", "exact", "diff", "handoff", "diff"))
for (nm in names(rules)) {
  R <- rules[[nm]]
  fitted <- mean(y_sim[R & treat_sim == 1L]) - mean(y_sim[R & treat_sim == 0L])

  # Exact target: the mean subject-level effect over the region. This is
  # .effect_one()'s MD branch evaluated at R instead of at flag_harm.
  exact <- mean(d$mu1[R]) - mean(d$mu0[R])

  # Handoff's form: tau * P(H_true | Hhat), with tau taken as beta_inter.
  handoff <- beta_inter * mean(d$flag_harm[R] == 1L)

  cat(sprintf("  %-42s %9.4f %9.4f %8.4f %9.4f %8.4f\n",
              nm, fitted, exact, abs(fitted - exact),
              handoff, abs(fitted - handoff)))
}
cat(sprintf("\n  delta = k_treat * beta_treat = %.4f\n", delta))
cat("  Every 'handoff' diff above should equal delta, and every 'exact' diff\n")
cat("  should be Monte Carlo error only.\n")

# --- C: delta needs no new computation -- it is the DGM's complement effect ---
cat("\n== C: delta is already reported by the DGM ==\n")
eff_Qc <- mean(df_super$mu1[df_super$flag_harm == 0]) -
          mean(df_super$mu0[df_super$flag_harm == 0])
cat(sprintf("  .effect_one() on Qc  = %.6f\n", eff_Qc))
cat(sprintf("  delta                = %.6f   diff %.2e\n",
            delta, abs(eff_Qc - delta)))

# --- D: the exact target carries NO Monte Carlo error ------------------------
# It is a finite mean over a fixed df_super, so repeated evaluation is identical.
cat("\n== D: exactness -- repeated evaluation of the target is bit-identical ==\n")
R <- rules[[1]]
v <- replicate(5, mean(d$mu1[R]) - mean(d$mu0[R]))
cat(sprintf("  5 evaluations, max pairwise difference = %.2e\n",
            max(v) - min(v)))
