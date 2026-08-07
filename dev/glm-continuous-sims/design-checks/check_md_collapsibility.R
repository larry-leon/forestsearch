# Why is MD special for beta(Hhat)?
#
# The handoff's answer: under INDEPENDENT covariates and CONJUNCTION rules,
#   beta(Hhat) = tau * P(H_true | Hhat)
# is a ratio of rectangle probabilities, hence a product of marginals. It then
# states the limits: correlated covariates need pmvnorm, GRF disjunctions need
# inclusion-exclusion.
#
# Alternative answer tested here: the operative property is COLLAPSIBILITY of the
# mean difference, not independence. Under randomization (A independent of X),
# for ANY region R,
#   beta(R) = E[Y | A=1, R] - E[Y | A=0, R] = E[tau(X) | R],
# because the prognostic part E[m(X) | R] is common to both arms and cancels.
# That is a within-region average of the subject-level effect. It needs no
# independence between covariates, no rectangle shape, and no step-shaped CATE.
#
# The ACTG175 continuous DGM resamples a FIXED finite super-population from real
# trial data, so if this holds, beta(Hhat) is an exact finite mean over df_super:
# no simulation, no model fit, no eval_seed, zero Monte Carlo error.
#
# Four tests. 1-3 stress the handoff's stated limits one at a time. Test 4 is the
# required negative control: the same construction on the OR scale must FAIL,
# since the odds ratio is non-collapsible.

set.seed(20260806)

n <- 3e6

# --- Population: CORRELATED covariates (violates the independence premise) ----
# X1, X2 correlated Gaussians pushed to uniform margins; X3 independent.
rho <- 0.7
Z1 <- rnorm(n)
Z2 <- rho * Z1 + sqrt(1 - rho^2) * rnorm(n)
X1 <- pnorm(Z1)
X2 <- pnorm(Z2)
X3 <- runif(n)
A  <- rbinom(n, 1, 0.5)               # randomized, independent of X

cat(sprintf("empirical cor(X1, X2) = %.4f  (independence premise violated)\n\n",
            cor(X1, X2)))

# --- GRADED CATE, not a step (violates the additive-step premise) -------------
# tau(x) is continuous in X1 and X2 and nonzero everywhere.
tau_fn <- function(x1, x2) 10 + 30 * x1 - 15 * x2^2

# Strong, nonlinear PROGNOSTIC signal that must cancel out of the contrast.
m_fn <- function(x1, x2, x3) 5 + 40 * x1^2 - 20 * x2 + 12 * sin(6 * x3)

tau_i <- tau_fn(X1, X2)
Y <- m_fn(X1, X2, X3) + tau_i * A + rnorm(n, sd = 8)

# --- Candidate regions, including shapes the handoff calls out as hard --------
regions <- list(
  "conjunction        {X1>.5} & {X2<.6}"    = (X1 > 0.5) & (X2 < 0.6),
  "DISJUNCTION        {X1>.7} | {X2<.2}"    = (X1 > 0.7) | (X2 < 0.2),
  "disjunction + conj  ({X1>.6}|{X2<.3}) & {X3>.5}" =
    ((X1 > 0.6) | (X2 < 0.3)) & (X3 > 0.5),
  "non-rectangle      {X1 + X2 > 1.2}"      = (X1 + X2 > 1.2),
  "negation           !{X1<=.4}"            = !(X1 <= 0.4)
)

cat("== Tests 1-3: MD, correlated covariates, graded CATE, arbitrary shapes ==\n")
cat(sprintf("  %-48s %10s %10s %9s\n",
            "region", "E[tau|R]", "fitted MD", "diff"))
for (nm in names(regions)) {
  R <- regions[[nm]]
  exact <- mean(tau_i[R])                                   # collapsibility claim
  fitted <- mean(Y[R & A == 1]) - mean(Y[R & A == 0])        # what the engine fits
  cat(sprintf("  %-48s %10.4f %10.4f %9.4f\n", nm, exact, fitted, abs(exact - fitted)))
}

# --- Test 1b: the handoff's product-of-marginals form under correlation -------
# Re-run the handoff's own setup but with X1, X2 CORRELATED and a step CATE, so
# the only changed premise is independence.
cat("\n== Test 1b: product-of-marginals vs truth when covariates are correlated ==\n")
tau0 <- 0.8
H_true <- (X1 > 0.5) & (X2 < 0.6)
Y2 <- 1 + tau0 * H_true * A + rnorm(n)
Hhat <- (X1 > 0.3) & (X2 < 0.8)

emp <- mean(Y2[Hhat & A == 1]) - mean(Y2[Hhat & A == 0])
collapse_form <- tau0 * mean(H_true[Hhat])          # tau * P(H_true | Hhat), counted
prod_marg <- tau0 * ((1 - max(0.5, 0.3)) * (min(0.6, 0.8) - 0)) /
                    ((1 - 0.3) * (0.8 - 0))          # handoff's closed form
cat(sprintf("  fitted MD                        %.5f\n", emp))
cat(sprintf("  counted  tau * P(H_true | Hhat)  %.5f   diff %.5f\n",
            collapse_form, abs(emp - collapse_form)))
cat(sprintf("  product-of-marginals form        %.5f   diff %.5f  <- FAILS\n",
            prod_marg, abs(emp - prod_marg)))

# --- Test 4: NEGATIVE CONTROL. The same move must fail for the odds ratio. ----
# Non-collapsibility: the marginal OR in a region is not the average of the
# subject-level ORs, so no analogue of the MD identity can hold.
cat("\n== Test 4 (negative control): OR is non-collapsible, so this must FAIL ==\n")
lp0 <- -0.5 + 1.8 * X1 - 1.2 * X2          # baseline logit, strong heterogeneity
delta <- 0.9                                # constant subject-level log-OR
p <- plogis(lp0 + delta * A)
Yb <- rbinom(n, 1, p)

for (nm in names(regions)[1:2]) {
  R <- regions[[nm]]
  fit <- glm(Yb[R] ~ A[R], family = binomial())
  marg_or <- exp(unname(coef(fit)[2]))
  mean_ind_or <- exp(delta)               # every subject's own OR
  cat(sprintf("  %-48s marginal OR %.4f   mean individual OR %.4f   diff %.4f\n",
              nm, marg_or, mean_ind_or, abs(marg_or - mean_ind_or)))
}
cat("\n  A nonzero diff here is the point: it confirms the MD result above is\n")
cat("  collapsibility and not an artifact that would pass for any effect measure.\n")
