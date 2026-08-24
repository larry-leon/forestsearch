# Test the ADDITIVE population covariance decomposition against the borderline
# residual recorded in check_declaration_gaussian_md.R (CHECK 6d): MC 0.5545,
# naive/var-corrected 0.5359/0.5360 (gap -0.0185, diagnosed as a correlation
# effect). Claim under test: with the full additive decomposition
#   Sigma_pop(g,g') = (4/n) * (P_int/(P_g P_g')) *
#                     ( sigma^2 + E[ u_g u_g' + (u_g t_g' + t_g u_g' + t_g t_g')/2 | int ] )
# where u_g = mu0 - E[mu0|g], t_g = bint*(1Q - P(Q|g)),
# the computed rate matches the population MC. Diagonal must reproduce the
# verified variance correction sigma^2 + Var(mu0|g) + Var(tau|g)/2 + Cov(mu0,tau|g).
suppressMessages(library(mvtnorm))
set.seed(46)

# --- family setup: verbatim geometry from the committed check script ---------
p_x1 <- 0.50; p_x2g1 <- 0.55; p_x2g0 <- 0.30; p_x3 <- 0.95
profiles <- expand.grid(x1 = 0:1, x2 = 0:1, x3 = 0:1)
p8 <- with(profiles,
  ifelse(x1 == 1, p_x1, 1 - p_x1) *
  ifelse(x1 == 1, ifelse(x2 == 1, p_x2g1, 1 - p_x2g1),
                  ifelse(x2 == 1, p_x2g0, 1 - p_x2g0)) *
  ifelse(x3 == 1, p_x3, 1 - p_x3))
Qprof <- with(profiles, x1 == 1 & x2 == 1)
mu0   <- with(profiles, 20 * x1 - 15 * x2 + 30 * x3)
delta_ref <- -26.26; sigma <- 70; c1 <- 30
cand <- list()
for (v in 1:3) for (a in 0:1) cand[[length(cand) + 1]] <- list(vars = v, vals = a)
prs <- combn(3, 2)
for (j in seq_len(ncol(prs))) for (a in 0:1) for (b in 0:1)
  cand[[length(cand) + 1]] <- list(vars = prs[, j], vals = c(a, b))
cand[[length(cand) + 1]] <- list(vars = 1:3, vals = c(1, 1, 1))
memb_prof <- function(cd) {
  m <- rep(TRUE, 8)
  for (k in seq_along(cd$vars)) m <- m & (profiles[[cd$vars[k]]] == cd$vals[k])
  m
}
G8 <- vapply(cand, memb_prof, logical(8))
Pg <- as.numeric(p8 %*% G8); keep <- Pg >= 0.10
G8 <- G8[, keep, drop = FALSE]; M <- ncol(G8)
Pg  <- as.numeric(p8 %*% G8)
PQg <- as.numeric((p8 * Qprof) %*% G8) / Pg

# --- scenario: CHECK 6d stressor (mu0 x3, borderline CATE) -------------------
mu0d <- 3 * mu0
run_scenario <- function(bint_s, tag, nrep = 6000, n_p = 600) {
  beta_pop <- delta_ref + bint_s * PQg
  m_g  <- as.numeric((p8 * mu0d) %*% G8) / Pg          # E[mu0|g]
  # pairwise intersection moments, all enumerable on the 8 profiles
  Sig_add <- matrix(NA_real_, M, M)
  for (i in seq_len(M)) for (j in i:M) {
    w <- p8 * (G8[, i] & G8[, j])
    Pint <- sum(w)
    if (Pint < 1e-15) { Sig_add[i, j] <- Sig_add[j, i] <- 0; next }
    ui <- mu0d - m_g[i]; uj <- mu0d - m_g[j]
    ti <- bint_s * (Qprof - PQg[i]); tj <- bint_s * (Qprof - PQg[j])
    Emom <- sum(w * (ui * uj + (ui * tj + ti * uj + ti * tj) / 2)) / Pint
    Sig_add[i, j] <- Sig_add[j, i] <-
      (4 / n_p) * (Pint / (Pg[i] * Pg[j])) * (sigma^2 * 1 + 0) +   # sigma block
      (4 / n_p) * (Pint / (Pg[i] * Pg[j])) * Emom                  # imbalance block
  }
  # internal check: diagonal reproduces the verified variance correction
  Vmu0 <- vapply(seq_len(M), function(m) {
    w <- p8 * G8[, m] / Pg[m]; sum(w * mu0d^2) - sum(w * mu0d)^2 }, numeric(1))
  Cmt  <- bint_s * (as.numeric((p8 * mu0d * Qprof) %*% G8) / Pg - m_g * PQg)
  Vtau <- bint_s^2 * PQg * (1 - PQg)
  diag_ref <- 4 * (sigma^2 + Vmu0 + Vtau / 2 + Cmt) / (n_p * Pg)
  cat(sprintf("[%s] diag identity max rel dev: %.2e\n", tag,
              max(abs(diag(Sig_add) - diag_ref) / diag_ref)))
  # rates (eigen-repair: Sig_add is a Gram matrix, PSD in exact arithmetic,
  # but the near-duplicate pair puts it close enough to singular that floating
  # point can yield a tiny negative eigenvalue; floor and report)
  ev <- eigen(Sig_add, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf("[%s] Sig_add eigenvalues: min %.3e, max %.3e\n", tag,
              min(ev), max(ev)))
  repair <- function(S) {
    e <- eigen((S + t(S)) / 2, symmetric = TRUE)
    e$vectors %*% diag(pmax(e$values, 1e-8 * max(e$values))) %*% t(e$vectors)
  }
  pA <- function(S) 1 - as.numeric(pmvnorm(upper = rep(c1, M), mean = beta_pop,
          sigma = repair(S), algorithm = GenzBretz(maxpts = 2e5, abseps = 1e-5)))
  cor_ovl <- outer(seq_len(M), seq_len(M), Vectorize(function(i, j)
    sum(p8[G8[, i] & G8[, j]]) / sqrt(Pg[i] * Pg[j])))
  sd_naive <- sqrt(4 * sigma^2 / (n_p * Pg))
  sd_corr  <- sqrt(diag_ref)
  # fresh population MC
  idx16 <- function(pf, a) (pf - 1L) * 2L + a + 1L
  decl <- vapply(seq_len(nrep), function(r) {
    pf <- sample.int(8, n_p, replace = TRUE, prob = p8)
    a  <- rbinom(n_p, 1, 0.5)
    y  <- mu0d[pf] + a * (delta_ref + bint_s * Qprof[pf]) + rnorm(n_p, 0, sigma)
    g16 <- idx16(pf, a)
    cnt <- tabulate(g16, 16); ys <- numeric(16)
    rs <- rowsum(y, g16); ys[as.integer(rownames(rs))] <- rs
    cm <- matrix(cnt, nrow = 2); ym <- matrix(ys, nrow = 2)
    b <- as.numeric(ym[2, ] %*% G8) / as.numeric(cm[2, ] %*% G8) -
         as.numeric(ym[1, ] %*% G8) / as.numeric(cm[1, ] %*% G8)
    any(b >= c1)
  }, logical(1))
  rate <- mean(decl); se <- sqrt(rate * (1 - rate) / nrep)
  cat(sprintf("[%s] MC %.4f (SE %.4f) | naive %.4f | var-corr(mult) %.4f | ADDITIVE %.4f\n",
      tag, rate, se,
      pA(cor_ovl * tcrossprod(sd_naive)),
      pA(cor_ovl * tcrossprod(sd_corr)),
      pA(Sig_add)))
  cat(sprintf("[%s] additive gap %+.4f (%.1f MC SEs)\n\n", tag,
              pA(Sig_add) - rate, abs(pA(Sig_add) - rate) / se))
  # show the diagnosed pair: corr(Q, dup) under sigma-overlap vs additive
  iQ <- which(vapply(seq_len(M), function(m) all(G8[, m] == (profiles$x1==1 & profiles$x2==1)), logical(1)))
  iD <- which(vapply(seq_len(M), function(m) all(G8[, m] == (profiles$x1==1 & profiles$x2==1 & profiles$x3==1)), logical(1)))
  cat(sprintf("[%s] corr(Q,dup): sigma-overlap %.4f -> additive %.4f\n\n", tag,
      cor_ovl[iQ, iD],
      Sig_add[iQ, iD] / sqrt(Sig_add[iQ, iQ] * Sig_add[iD, iD])))
}
run_scenario(c1 - delta_ref, "borderline")
run_scenario(66.26,          "alt")
