# ============================================================================
# check_declaration_gaussian_md.R
#
# Design checks for detection-rate theory in the continuous/MD family:
# the candidate-effect process over a fixed enumerated family is (exactly,
# conditional on design, under Gaussian noise) multivariate normal with
#   mean_m  = collapsibility:  beta(g_m) = delta + beta_inter * P(Q | g_m)
#   cov     = overlap:         Cov = sigma^2 * C C^T  (exact two-arm formula)
# and declaration / selection events are MVN functionals computable without
# data-level simulation.
#
# Conventions (verified against package source, branch feature/glm-extension,
# HEAD 3c40467):
#   * ORIENTED scale throughout: positive effect = HARM. This is the scale on
#     which forestsearch's continuous estimator operates after negating Y when
#     adverse_outcome = FALSE (R/glm_effect_estimators.R:814-821).
#   * The estimator is the within-subgroup OLS treatment coefficient of
#     lm(y ~ treat) == difference in arm means (same file).
#   * The screening floor compares the POINT ESTIMATE to the threshold
#     (R/subgroup_search.R:626; R/subgroup_consistency_main.R:528-531), not a
#     test statistic. Estimated variance therefore enters the theory only via
#     plug-in Sigma-hat (Check 6c), never the event definition.
#   * Reference-cell shape mirrors the vignette calibration recorded in
#     dev/glm-continuous-sims/NOTE_target_is_collapsibility.md: two-valued
#     CATE, nonzero complement effect delta < 0, harm subgroup at +40
#     oriented. sigma = 70 is a synthetic ACTG175-magnitude stand-in (the
#     harness's residual sd is data-fitted, not a quotable constant).
#   * Thresholds: committed discriminating cell c1 = 30 (screening),
#     c2 = 10 (consistency), oriented scale.
#
# Standalone: base R + mvtnorm. Seeded. Runtime ~1-2 min.
# ============================================================================

suppressMessages({
  if (!requireNamespace("mvtnorm", quietly = TRUE))
    stop("mvtnorm is required")
})
library(mvtnorm)
library(stats)

set.seed(20260823)
t0 <- proc.time()[3]

`%+%` <- function(a, b) paste0(a, b)
hr <- function(x) cat("\n", strrep("=", 76), "\n", x, "\n", strrep("=", 76), "\n", sep = "")
chk <- function(label, value, tol, extra = "") {
  ok <- is.finite(value) && value <= tol
  cat(sprintf("  [%s] %-52s %10.3e (tol %.1e)%s\n",
              if (ok) "PASS" else "FAIL", label, value, tol,
              if (nzchar(extra)) paste0("  ", extra) else ""))
  invisible(ok)
}
note <- function(...) cat("  ", sprintf(...), "\n", sep = "")

# ----------------------------------------------------------------------------
# DGM over an enumerable covariate space
# ----------------------------------------------------------------------------
# 3 binary covariates; x1,x2 dependent (collapsibility needs no independence),
# x3 nearly ubiquitous (P = .95) to manufacture a near-duplicate candidate.
p_x1   <- 0.50
p_x2g1 <- 0.55   # P(x2=1 | x1=1)
p_x2g0 <- 0.30   # P(x2=1 | x1=0)
p_x3   <- 0.95

profiles <- expand.grid(x1 = 0:1, x2 = 0:1, x3 = 0:1)
p8 <- with(profiles,
  ifelse(x1 == 1, p_x1, 1 - p_x1) *
  ifelse(x1 == 1, ifelse(x2 == 1, p_x2g1, 1 - p_x2g1),
                  ifelse(x2 == 1, p_x2g0, 1 - p_x2g0)) *
  ifelse(x3 == 1, p_x3, 1 - p_x3))
stopifnot(abs(sum(p8) - 1) < 1e-12)

Qprof  <- with(profiles, x1 == 1 & x2 == 1)          # true harm rule Q
mu0    <- with(profiles, 20 * x1 - 15 * x2 + 30 * x3) # prognostic main effects

# Reference-cell CATE shape (oriented scale): tau = delta + beta_inter * 1{Q}
delta_ref <- -26.26
bint_ref  <-  66.26          # => tau_Q = +40.00 oriented
sigma     <-  70
c1 <- 30; c2 <- 10           # committed discriminating cell

# ----------------------------------------------------------------------------
# Candidate family: 1- and 2-factor conjunctions above a prevalence floor,
# plus an explicit near-duplicate of Q ({x1=1,x2=1,x3=1}, containment .95).
# ----------------------------------------------------------------------------
floor_prev <- 0.10
cand <- list()
for (v in 1:3) for (a in 0:1)
  cand[[length(cand) + 1]] <- list(vars = v, vals = a)
prs <- combn(3, 2)
for (j in seq_len(ncol(prs))) for (a in 0:1) for (b in 0:1)
  cand[[length(cand) + 1]] <- list(vars = prs[, j], vals = c(a, b))
cand[[length(cand) + 1]] <- list(vars = 1:3, vals = c(1, 1, 1))  # near-dup

memb_prof <- function(cd) {
  m <- rep(TRUE, 8)
  for (k in seq_along(cd$vars)) m <- m & (profiles[[cd$vars[k]]] == cd$vals[k])
  m
}
G8   <- vapply(cand, memb_prof, logical(8))          # 8 x n_cand
Pg   <- as.numeric(p8 %*% G8)
keep <- Pg >= floor_prev
G8   <- G8[, keep, drop = FALSE]
cand <- cand[keep]
M    <- ncol(G8)
lab  <- vapply(cand, function(cd)
  paste(sprintf("x%d=%d", cd$vars, cd$vals), collapse = "&"), character(1))
Pg   <- as.numeric(p8 %*% G8)
PQg  <- as.numeric((p8 * Qprof) %*% G8) / Pg          # P(Q | g_m), exact
iQ   <- which(lab == "x1=1&x2=1")
iDup <- which(lab == "x1=1&x2=1&x3=1")

hr("SETUP  family M = " %+% M %+% " candidates (floor " %+% floor_prev %+% ")")
note("P(Q) = %.4f ; near-dup containment P(dup)/P(Q) = %.3f",
     sum(p8[Qprof]), Pg[iDup] / Pg[iQ])
note("population corr(Q, dup) [overlap formula] = %.4f",
     Pg[iDup] / sqrt(Pg[iQ] * Pg[iDup]))

# ============================================================================
# CHECK 1 -- Family-wide collapsibility identity, machine precision
# beta(g) = delta + beta_inter * P(Q|g)  vs  direct enumeration of
# E[Y|A=1,g] - E[Y|A=0,g] (prognostic mu0 present and must cancel exactly).
# ============================================================================
hr("CHECK 1  family-wide collapsibility (exact enumeration)")
tau8     <- delta_ref + bint_ref * Qprof
beta_fml <- delta_ref + bint_ref * PQg
beta_dir <- vapply(seq_len(M), function(m) {
  w <- p8 * G8[, m]
  sum(w * (mu0 + tau8)) / sum(w) - sum(w * mu0) / sum(w)
}, numeric(1))
chk("max |formula - enumeration| over family", max(abs(beta_fml - beta_dir)), 1e-10)
note("range of beta(g) over family: [%.2f, %.2f]; beta(Q) = %.2f",
     min(beta_fml), max(beta_fml), beta_fml[iQ])

# ----------------------------------------------------------------------------
# One fixed design (the conditioning event for Checks 2-5)
# ----------------------------------------------------------------------------
n <- 600
prof_i <- sample.int(8, n, replace = TRUE, prob = p8)
A      <- rbinom(n, 1, 0.5)
Gn     <- G8[prof_i, , drop = FALSE]                  # n x M membership
n1g    <- colSums(Gn * A); n0g <- colSums(Gn * (1 - A))
Cmat   <- t(Gn * (A / rep(n1g, each = n) - (1 - A) / rep(n0g, each = n))) # M x n
EYfix  <- function(dl, bi) mu0[prof_i] + A * (dl + bi * Qprof[prof_i])
a0 <- as.numeric(Cmat %*% mu0[prof_i])
a1 <- as.numeric(Cmat %*% A)
a2 <- as.numeric(Cmat %*% (A * Qprof[prof_i]))
cond_mean <- function(dl, bi) a0 + dl * a1 + bi * a2
Sig  <- sigma^2 * tcrossprod(Cmat)                    # exact conditional cov
sdg  <- sqrt(diag(Sig))

# ============================================================================
# CHECK 2 -- Design-conditional joint law: mean, covariance, overlap formulas
# ============================================================================
hr("CHECK 2  design-conditional joint law (n = " %+% n %+% ")")
# (a) exact two-arm overlap formula == sigma^2 C C^T, machine precision
Sig_ovl <- matrix(0, M, M)
for (i in seq_len(M)) for (j in i:M) {
  int1 <- sum(Gn[, i] & Gn[, j] & A == 1)
  int0 <- sum(Gn[, i] & Gn[, j] & A == 0)
  Sig_ovl[i, j] <- Sig_ovl[j, i] <-
    sigma^2 * (int1 / (n1g[i] * n1g[j]) + int0 / (n0g[i] * n0g[j]))
}
chk("max |Sigma(CC^T) - two-arm overlap formula| (rel)",
    max(abs(Sig - Sig_ovl)) / max(abs(Sig)), 1e-12)
# (b) balanced approximation 4 sigma^2 n_int / (n_g n_g')
ng <- n1g + n0g
Sig_bal <- matrix(0, M, M)
for (i in seq_len(M)) for (j in i:M) {
  nint <- sum(Gn[, i] & Gn[, j])
  Sig_bal[i, j] <- Sig_bal[j, i] <- 4 * sigma^2 * nint / (ng[i] * ng[j])
}
note("balanced overlap approx: max rel dev from exact = %.3f (informational)",
     max(abs(Sig_bal - Sig) / pmax(abs(Sig), 1e-9)))
note("realized corr(Q, dup) = %.4f", Sig[iQ, iDup] / (sdg[iQ] * sdg[iDup]))
# (c) empirical mean & correlation vs exact, fixed design, Gaussian noise
nsim <- 4000
eps  <- matrix(rnorm(n * nsim, 0, sigma), n, nsim)
Nois <- Cmat %*% eps                                  # M x nsim, mean 0
bhat_ref <- cond_mean(delta_ref, bint_ref) + Nois     # reference scenario draws
emp_mean_z <- max(abs(rowMeans(Nois)) / (sdg / sqrt(nsim)))
chk("max |empirical - exact mean| / SE(mean)", emp_mean_z, 4)
emp_cor <- cor(t(bhat_ref))
the_cor <- Sig / tcrossprod(sdg)
chk("max |empirical - exact correlation|",
    max(abs(emp_cor - the_cor)), 4 / sqrt(nsim) + 0.02)

# ============================================================================
# CHECK 3 -- Rung A: family screening exceedance  P(max_m bhat_m >= c1)
# pmvnorm rectangle-complement vs data-level MC, three scenarios.
# ============================================================================
hr("CHECK 3  Rung A: P(max >= c1), pmvnorm vs data-level MC")
bi_bord <- uniroot(function(b) max(cond_mean(delta_ref, b)) - c1,
                   c(0, 200))$root
scens <- list(
  alt        = c(delta_ref, bint_ref),
  borderline = c(delta_ref, bi_bord),
  null       = c(0, 0)
)
note("borderline beta_inter (max conditional mean = c1): %.3f", bi_bord)
for (s in names(scens)) {
  m_s  <- cond_mean(scens[[s]][1], scens[[s]][2])
  pmv  <- pmvnorm(upper = rep(c1, M), mean = m_s, sigma = Sig)
  pA   <- 1 - as.numeric(pmv)
  rate <- mean(apply(m_s + Nois >= c1, 2, any))
  se   <- sqrt(rate * (1 - rate) / nsim) + attr(pmv, "error")
  chk(sprintf("%-10s computed %.4f vs MC %.4f", s, pA, rate),
      abs(pA - rate), 4 * se + 1e-4)
}

# ============================================================================
# CHECK 4 -- Rung B: split-half screening + consistency (idealized 2-split)
# pass_m = { W1+W2 >= 2 c1 , min(W1,W2) >= c2 } ; declare if any m passes.
# Gaussian-level MC on the exact 2M block law vs data-level MC; per-candidate
# P1 by 1-D integral (the April/SIM-2024 Section 2.1 object, per candidate).
# ============================================================================
hr("CHECK 4  Rung B: split-half consistency, family and per-candidate")
h1 <- sample.int(n, n %/% 2); h2 <- setdiff(seq_len(n), h1)
half_geom <- function(idx) {
  Gh <- Gn[idx, , drop = FALSE]; Ah <- A[idx]
  n1 <- colSums(Gh * Ah); n0 <- colSums(Gh * (1 - Ah))
  Ch <- t(Gh * (Ah / rep(n1, each = length(idx)) -
                (1 - Ah) / rep(n0, each = length(idx))))
  list(C = Ch, idx = idx,
       a0 = as.numeric(Ch %*% mu0[prof_i[idx]]),
       a1 = as.numeric(Ch %*% Ah),
       a2 = as.numeric(Ch %*% (Ah * Qprof[prof_i[idx]])),
       Sig = sigma^2 * tcrossprod(Ch))
}
H1 <- half_geom(h1); H2 <- half_geom(h2)
m1 <- H1$a0 + delta_ref * H1$a1 + bint_ref * H1$a2
m2 <- H2$a0 + delta_ref * H2$a1 + bint_ref * H2$a2
R  <- 2e5
W1 <- rmvnorm(R, mean = m1, sigma = H1$Sig)
W2 <- rmvnorm(R, mean = m2, sigma = H2$Sig)
passG   <- (W1 + W2 >= 2 * c1) & (W1 >= c2) & (W2 >= c2)
famB_G  <- mean(rowSums(passG) > 0)
W1d <- m1 + H1$C %*% eps[h1, ]; W2d <- m2 + H2$C %*% eps[h2, ]
passD  <- (W1d + W2d >= 2 * c1) & (W1d >= c2) & (W2d >= c2)
famB_D <- mean(colSums(passD) > 0)
seB <- sqrt(famB_G * (1 - famB_G) / R) + sqrt(famB_D * (1 - famB_D) / nsim)
chk(sprintf("family rate: Gaussian-level %.4f vs data %.4f", famB_G, famB_D),
    abs(famB_G - famB_D), 4 * seB)
p1_split <- function(mu1, s1, mu2, s2, k1, k2) {
  integrate(function(w) dnorm(w, mu1, s1) *
              pnorm(pmax(k2, 2 * k1 - w), mu2, s2, lower.tail = FALSE),
            lower = k2, upper = Inf, rel.tol = 1e-9)$value
}
p1_int <- vapply(seq_len(M), function(m)
  p1_split(m1[m], sqrt(H1$Sig[m, m]), m2[m], sqrt(H2$Sig[m, m]), c1, c2),
  numeric(1))
p1_G <- colMeans(passG)
chk("per-candidate P1: max |integral - Gaussian MC|",
    max(abs(p1_int - p1_G)), 4 * sqrt(max(p1_G) / R) + 1e-3)
note("P1(Q) = %.4f ; P1(dup) = %.4f ; family = %.4f  (competition inflates family over best single only mildly here)",
     p1_int[iQ], p1_int[iDup], famB_G)

# ============================================================================
# CHECK 5 -- Rung C: selection identity  P(argmax = g_m & max >= c1)
# per-m pmvnorm on the difference transform vs data-level selection freq.
# The near-duplicate makes this the sharpest competition test in the script.
# ============================================================================
hr("CHECK 5  Rung C: selection probabilities (argmax + screening)")
m_alt <- cond_mean(delta_ref, bint_ref)
sel_pmv <- vapply(seq_len(M), function(m) {
  Amat <- rbind(diag(M)[m, ],
                do.call(rbind, lapply(setdiff(seq_len(M), m),
                        function(j) diag(M)[m, ] - diag(M)[j, ])))
  mu_t <- as.numeric(Amat %*% m_alt)
  Sg_t <- Amat %*% Sig %*% t(Amat)
  as.numeric(pmvnorm(lower = c(c1, rep(0, M - 1)),
                     upper = rep(Inf, M),
                     mean = mu_t, sigma = Sg_t,
                     algorithm = GenzBretz(maxpts = 50000, abseps = 5e-4)))
}, numeric(1))
bh_alt  <- m_alt + Nois
win     <- max.col(t(bh_alt))
decl    <- apply(bh_alt >= c1, 2, any)
sel_mc  <- vapply(seq_len(M), function(m) mean(win == m & decl), numeric(1))
chk("max_m |pmvnorm - MC selection prob|",
    max(abs(sel_pmv - sel_mc)), 4 * sqrt(0.5 / nsim) + 2e-3)
chk("sum_m selection probs vs Rung A rate (consistency)",
    abs(sum(sel_pmv) - (1 - as.numeric(
      pmvnorm(upper = rep(c1, M), mean = m_alt, sigma = Sig)))), 5e-3)
p1_naive_Q <- pnorm(c1, m_alt[iQ], sdg[iQ], lower.tail = FALSE)
note("competition effect at Q:  naive per-candidate P(bhat_Q >= c1) = %.4f", p1_naive_Q)
note("  vs P(select Q) = %.4f, P(select dup) = %.4f, P(select Q or dup) = %.4f",
     sel_pmv[iQ], sel_pmv[iDup], sel_pmv[iQ] + sel_pmv[iDup])
note("  -> the near-duplicate cannibalizes selection; per-candidate P1 is not a selection rate.")

# ============================================================================
# CHECK 6 -- Edge-case stressors
# ============================================================================
hr("CHECK 6a  near-collinearity: pmvnorm stability at corr ~ .97")
pmv_a <- pmvnorm(upper = rep(c1, M), mean = m_alt, sigma = Sig)
chk("pmvnorm reported error at realized corr(Q,dup)",
    attr(pmv_a, "error"), 1e-3,
    sprintf("corr = %.4f", Sig[iQ, iDup] / (sdg[iQ] * sdg[iDup])))
note("exact duplicates (corr = 1) make Sigma singular: de-duplicate the family first.")

hr("CHECK 6b  small n + skewed noise: CLT quality of the Gaussian computation")
n_s <- 150
prof_s <- sample.int(8, n_s, replace = TRUE, prob = p8)
A_s    <- rbinom(n_s, 1, 0.5)
G_s    <- G8[prof_s, , drop = FALSE]
n1s <- colSums(G_s * A_s); n0s <- colSums(G_s * (1 - A_s))
if (min(n1s, n0s) < 3) note("WARNING: tiny arm count %d in some candidate", min(n1s, n0s))
C_s  <- t(G_s * (A_s / rep(n1s, each = n_s) - (1 - A_s) / rep(n0s, each = n_s)))
m_s  <- as.numeric(C_s %*% (mu0[prof_s] + A_s * (delta_ref + bint_ref * Qprof[prof_s])))
Sg_s <- sigma^2 * tcrossprod(C_s)
pA_s <- 1 - as.numeric(pmvnorm(upper = rep(c1, M), mean = m_s, sigma = Sg_s))
ln_s   <- 1.0                                  # lognormal shape: skewness ~ 6.2
ln_mu  <- exp(ln_s^2 / 2); ln_sd <- sqrt((exp(ln_s^2) - 1) * exp(ln_s^2))
mk_eps <- function(kind) {
  if (kind == "gauss") matrix(rnorm(n_s * nsim, 0, sigma), n_s, nsim)
  else matrix((rlnorm(n_s * nsim, 0, ln_s) - ln_mu) / ln_sd * sigma, n_s, nsim)
}
for (kind in c("gauss", "lognormal(skew~6.2)")) {
  e  <- mk_eps(if (kind == "gauss") "gauss" else "ln")
  rt <- mean(apply(m_s + C_s %*% e >= c1, 2, any))
  se <- sqrt(rt * (1 - rt) / nsim)
  tol <- if (kind == "gauss") 4 * se + 1e-3 else 6 * se + 5e-3
  chk(sprintf("n=%d %-20s computed %.4f vs MC %.4f", n_s, kind, pA_s, rt),
      abs(pA_s - rt), tol,
      if (kind != "gauss") "(diagnostic: CLT stress, looser tol)" else "")
}
note("min candidate arm size at n=%d: %d  (CLT floor for the smallest cell)",
     n_s, min(c(n1s, n0s)))

hr("CHECK 6c  plug-in Sigma-hat from ONE dataset (design known, sigma unknown)")
y1   <- EYfix(delta_ref, bint_ref) + eps[, 1]
sd_hat <- vapply(seq_len(M), function(m) {
  d <- data.frame(y = y1[Gn[, m]], tr = A[Gn[, m]])
  summary(lm(y ~ tr, data = d))$coefficients["tr", "Std. Error"]
}, numeric(1))
Sg_hat <- the_cor * tcrossprod(sd_hat)          # analyst corr from known design
pA_hat <- 1 - as.numeric(pmvnorm(upper = rep(c1, M), mean = m_alt, sigma = Sg_hat))
pA_tru <- 1 - as.numeric(pmvnorm(upper = rep(c1, M), mean = m_alt, sigma = Sig))
note("per-candidate lm SEs inflate over sigma-only truth by factor %.3f-%.3f",
     min(sd_hat / sdg), max(sd_hat / sdg))
note("  (lm residual variance absorbs prognostic mu0 spread + within-g CATE spread)")
chk("Rung A rate: plug-in Sigma-hat vs true-sigma computation",
    abs(pA_hat - pA_tru), 0.05,
    sprintf("plug-in %.4f vs true %.4f (diagnostic)", pA_hat, pA_tru))

hr("CHECK 6d  population level (random designs): naive sigma vs variance-corrected")
# Stressor: prognostic spread x3 and the borderline CATE (max pop mean = c1),
# where the declaration rate is steepest in the variance. The corrected
# variance is 4(sigma^2 + Var(mu0|g) + Var(tau|g)/2 + Cov(mu0,tau|g))/n_g:
# arm 1 carries Var(mu0+tau|g) = Vmu0 + Vtau + 2C, arm 0 carries Vmu0.
# Stress: prognostic spread x3. Corrected variance per candidate:
#   4(sigma^2 + Var(mu0|g) + Var(tau|g)/2 + Cov(mu0,tau|g)) / n_g
# (arm 1 carries Var(mu0+tau|g), arm 0 carries Var(mu0|g)).
# Correlations are kept at the sigma-overlap formula -- deliberately, to expose
# where that stops being enough (see borderline diagnostic below).
mu0d <- 3 * mu0
n_p  <- 600; nrep <- 4000
Emu0g <- as.numeric((p8 * mu0d) %*% G8) / Pg
Vmu0 <- vapply(seq_len(M), function(m) {
  w <- p8 * G8[, m] / Pg[m]; sum(w * mu0d^2) - sum(w * mu0d)^2 }, numeric(1))
EmQ  <- as.numeric((p8 * mu0d * Qprof) %*% G8) / Pg
cor_pop <- outer(seq_len(M), seq_len(M), Vectorize(function(i, j)
  sum(p8[G8[, i] & G8[, j]]) / sqrt(Pg[i] * Pg[j])))
idx16 <- function(pf, a) (pf - 1L) * 2L + a + 1L
run_pop <- function(bint_s, tag, pass_check) {
  beta_pop <- delta_ref + bint_s * PQg
  Cmt  <- bint_s * (EmQ - Emu0g * PQg)
  Vtau <- bint_s^2 * PQg * (1 - PQg)
  sd_naive <- sqrt(4 * sigma^2 / (n_p * Pg))
  sd_corr  <- sqrt(4 * (sigma^2 + Vmu0 + Vtau / 2 + Cmt) / (n_p * Pg))
  pA <- function(sdv) 1 - as.numeric(
    pmvnorm(upper = rep(c1, M), mean = beta_pop,
            sigma = cor_pop * tcrossprod(sdv),
            algorithm = GenzBretz(maxpts = 1e5, abseps = 1e-5)))
  decl <- vapply(seq_len(nrep), function(r) {
    pf <- sample.int(8, n_p, replace = TRUE, prob = p8)
    a  <- rbinom(n_p, 1, 0.5)
    y  <- mu0d[pf] + a * (delta_ref + bint_s * Qprof[pf]) + rnorm(n_p, 0, sigma)
    g16 <- idx16(pf, a)
    cnt <- tabulate(g16, 16); ys <- numeric(16)
    rs  <- rowsum(y, g16); ys[as.integer(rownames(rs))] <- rs
    cnt_m <- matrix(cnt, nrow = 2); ys_m <- matrix(ys, nrow = 2)
    n1v <- as.numeric(cnt_m[2, ] %*% G8); n0v <- as.numeric(cnt_m[1, ] %*% G8)
    b <- as.numeric(ys_m[2, ] %*% G8) / n1v - as.numeric(ys_m[1, ] %*% G8) / n0v
    any(b >= c1)
  }, logical(1))
  rate <- mean(decl); se <- sqrt(rate * (1 - rate) / nrep)
  cat(sprintf("  %-10s MC %.4f (SE %.4f) | naive %.4f | var-corrected %.4f\n",
              tag, rate, se, pA(sd_naive), pA(sd_corr)))
  if (pass_check)
    chk(sprintf("%s: var-corrected vs population MC", tag),
        abs(pA(sd_corr) - rate), 4 * se + 0.01)
  else
    note("%s gap %+.4f is a correlation effect, not a variance one: the imbalance",
         tag, pA(sd_corr) - rate)
  invisible(NULL)
}
run_pop(bint_ref, "alt", TRUE)
run_pop(c1 - delta_ref, "borderline", FALSE)
note("components of Q and its near-duplicate are much less correlated than their")
note("sigma components (the duplicate has no x3 variation inside it), so the")
note("sigma-overlap correlation overstates their joint dependence and understates")
note("P(max) exactly where both sit at the threshold. Population-level closed forms")
note("therefore need the full ADDITIVE covariance decomposition (all terms")
note("enumerable) -- identified deferred item; the design-conditional law (Checks")
note("2-5) is exact and unaffected.")

# ============================================================================
# CHECK 7 -- Implied L_eff as an OUTPUT (fixed family, null scenario)
# Bridge to theory_glm_leff_correction.qmd: there L_eff is fitted from a
# reference simulation; here the joint law yields the family rate directly and
# an implied count becomes a derived readout. Candidate variances are
# heterogeneous (prevalence-driven), so the implied count is normalized
# against max(P1): in the deep tail P_fam is dominated by the highest-variance
# candidates and a median-P1 normalization diverges.
# ============================================================================
hr("CHECK 7  implied L_eff from the joint law (null: delta = beta_inter = 0)")
cat(sprintf("  %6s %9s %9s %9s %11s %10s\n",
            "N", "P1_max", "P_indep", "P_fam", "Leff(max)", "fam/union"))
for (N in c(300, 600, 1200)) {
  pf <- sample.int(8, N, replace = TRUE, prob = p8)
  a  <- rbinom(N, 1, 0.5)
  Gp <- G8[pf, , drop = FALSE]
  n1 <- colSums(Gp * a); n0 <- colSums(Gp * (1 - a))
  Cp <- t(Gp * (a / rep(n1, each = N) - (1 - a) / rep(n0, each = N)))
  mN <- as.numeric(Cp %*% mu0[pf])          # null: imbalance-only means
  Sp <- sigma^2 * tcrossprod(Cp)
  p1 <- pnorm(c1, mN, sqrt(diag(Sp)), lower.tail = FALSE)
  pf_j <- 1 - as.numeric(pmvnorm(upper = rep(c1, M), mean = mN, sigma = Sp,
            algorithm = GenzBretz(maxpts = 2e5, abseps = 1e-6)))
  pind <- 1 - prod(1 - p1)
  leffx <- log(1 - pf_j) / log(1 - max(p1))
  cat(sprintf("  %6d %9.5f %9.5f %9.5f %11.2f %10.3f\n",
              N, max(p1), pind, pf_j, leffx, pf_j / sum(p1)))
}
note("Leff(max) = implied count of max-equivalent candidates: O(few), roughly stable;")
note("fam/union < 1 is the correlation discount. A single-P1 L_eff (April form) needs")
note("~homogeneous candidate variances, which prevalence heterogeneity breaks in the")
note("tail. The April growth-in-N came from family REGENERATION (n.min viability, new")
note("cuts), which the fixed-family computation deliberately isolates away.")

cat(sprintf("\nTotal runtime: %.1f s\n", proc.time()[3] - t0))
cat("check_declaration_gaussian_md.R done.\n")
