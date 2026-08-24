# check_binary_rd_process.R
# ---------------------------------------------------------------------------
# B1 smoke test for extending the MD detection theory to BINARY outcomes on
# the risk-difference (RD) scale, and for the joint law with the package's
# Crump constant-CATE gate (test_hte_crump.R, sourced verbatim).
#
# Claims under test:
#  BR1  Means: RD is collapsible -> E[dphat(g)|D] = c(g)' p_W(X) exactly;
#       population value = E[tau_RD(X)|g] by enumeration.
#  BR2  Covariance: exact heteroskedastic overlap formula
#       Cov(g,g') = sum_a sum_{i in cap, W=a} p_a(1-p_a) / (n_{g,a} n_{g',a})
#       matches empirical covariance (Bernoulli CLT for the law itself).
#  BR3  Rung A on RD: pmvnorm family declaration vs data-level MC (alt+null).
#  BR4  Joint with Crump: the arm-difference OLS vector D = b1-b0 is linear
#       in Y; its exact conditional covariance and cross-covariance with the
#       candidate process follow from the same sum_i c c' Var(Y_i); the
#       gate-and-declare probability P(ConstCATE rejects & family declares)
#       computed at the Gaussian level matches data-level MC running the
#       ACTUAL package test (plug-in sandwich included).
# Base R + mvtnorm (+ MASS only via the sourced file's error fallbacks).
# ---------------------------------------------------------------------------
suppressMessages(library(mvtnorm))
source(file.path(rprojroot::find_root(rprojroot::is_r_package), "R", "test_hte_crump.R"))
set.seed(20260824)
t0 <- proc.time()[3]
hr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")
chk <- function(lab, v, tol, extra = "") cat(sprintf(
  "  [%s] %-52s %9.3e (tol %.1e)%s\n", if (is.finite(v) && v <= tol) "PASS" else "FAIL",
  lab, v, tol, if (nzchar(extra)) paste0("  ", extra) else ""))
note <- function(...) cat("  ", sprintf(...), "\n", sep = "")

# --- geometry: identical 8-profile family to the committed MD check --------
p_x1 <- .5; p_x2g1 <- .55; p_x2g0 <- .30; p_x3 <- .95
profiles <- expand.grid(x1 = 0:1, x2 = 0:1, x3 = 0:1)
p8 <- with(profiles, ifelse(x1==1,p_x1,1-p_x1) *
  ifelse(x1==1, ifelse(x2==1,p_x2g1,1-p_x2g1), ifelse(x2==1,p_x2g0,1-p_x2g0)) *
  ifelse(x3==1,p_x3,1-p_x3))
Qp <- with(profiles, x1==1 & x2==1)
cand <- list(); for (v in 1:3) for (a in 0:1) cand[[length(cand)+1]] <- list(v=v,a=a)
pr <- combn(3,2)
for (j in seq_len(ncol(pr))) for (a in 0:1) for (b in 0:1)
  cand[[length(cand)+1]] <- list(v=pr[,j], a=c(a,b))
cand[[length(cand)+1]] <- list(v=1:3, a=c(1,1,1))
G8 <- vapply(cand, function(cd){m<-rep(TRUE,8)
  for(k in seq_along(cd$v)) m <- m & (profiles[[cd$v[k]]]==cd$a[k]); m}, logical(8))
Pg <- as.numeric(p8 %*% G8); keep <- Pg >= .10
G8 <- G8[,keep,drop=FALSE]; M <- ncol(G8)
Pg <- as.numeric(p8 %*% G8)                      # recompute on the kept family
lab <- vapply(cand[keep], function(cd) paste(sprintf("x%d=%d",cd$v,cd$a),collapse="&"), character(1))

# --- binary DGM (adverse event; RD > 0 = harm) -----------------------------
eta0 <- with(profiles, -0.40 + 0.50*x1 - 0.40*x2 + 0.30*x3)
b_tr <- -0.35; b_int <- 1.25                    # benefit outside Q, harm inside
p0 <- plogis(eta0)
p1 <- plogis(eta0 + b_tr + b_int*Qp)
tauRD <- p1 - p0
c1 <- 0.10                                       # RD screening floor (harm)
betaRD <- as.numeric((p8*tauRD) %*% G8) / Pg     # exact E[tau_RD | g]
hr(sprintf("SETUP  M = %d; RD means range [%.3f, %.3f]; beta_RD(Q-rule) = %.3f; floor %.2f",
           M, min(betaRD), max(betaRD), betaRD[lab=="x1=1&x2=1"], c1))

# --- one fixed design -------------------------------------------------------
n <- 600
pf <- sample.int(8, n, TRUE, p8); W <- rbinom(n,1,.5)
Gn <- G8[pf,,drop=FALSE]
n1g <- colSums(Gn*W); n0g <- colSums(Gn*(1-W))
Cc  <- t(Gn * (W/rep(n1g,each=n) - (1-W)/rep(n0g,each=n)))   # M x n
pY  <- ifelse(W==1, p1[pf], p0[pf])              # E[Y_i | D]
vY  <- pY*(1-pY)                                 # Var[Y_i | D]
m_c <- as.numeric(Cc %*% pY)                     # exact conditional means
Sig <- Cc %*% (t(Cc)*vY)                         # exact heterosked. covariance

hr("CHECK BR1  RD means: conditional identity and collapsibility")
nsim <- 4000
Ymat <- matrix(rbinom(n*nsim, 1, rep(pY, nsim)), n, nsim)
Bh   <- Cc %*% Ymat                              # M x nsim candidate effects
chk("max |emp mean - c(g)'p| / SE", max(abs(rowMeans(Bh)-m_c)/sqrt(diag(Sig)/nsim)), 4)
note("max |cond mean - population E[tau_RD|g]| = %.4f (finite-n imbalance, informational)",
     max(abs(m_c - betaRD)))

hr("CHECK BR2  heteroskedastic overlap covariance, exact vs empirical")
Sig_ovl <- matrix(0, M, M)
for (i in 1:M) for (j in i:M) {
  s <- 0
  for (a in 0:1) {
    ii <- Gn[,i] & Gn[,j] & W==a
    na_i <- if (a==1) n1g[i] else n0g[i]; na_j <- if (a==1) n1g[j] else n0g[j]
    s <- s + sum(vY[ii])/(na_i*na_j)
  }
  Sig_ovl[i,j] <- Sig_ovl[j,i] <- s
}
chk("formula == C diag(v) C' (rel)", max(abs(Sig-Sig_ovl))/max(abs(Sig)), 1e-12)
empS <- cov(t(Bh)); D <- sqrt(diag(Sig))
chk("max |emp - exact correlation|", max(abs(empS/tcrossprod(sqrt(diag(empS))) -
    Sig/tcrossprod(D))), 4/sqrt(nsim)+.02)
Sig_hom <- outer(1:M,1:M,Vectorize(function(i,j){nint<-sum(Gn[,i]&Gn[,j]);
  4*mean(vY)*nint/((n1g[i]+n0g[i])*(n1g[j]+n0g[j]))}))
note("homoskedastic (mean-v) approx: max rel dev from exact = %.3f (informational)",
     max(abs(Sig_hom-Sig)/pmax(abs(Sig),1e-9)))

hr("CHECK BR3  Rung A on RD: pmvnorm vs data-level MC")
for (sc in c("alt","null")) {
  mm <- if (sc=="alt") m_c else as.numeric(Cc %*% ifelse(W==1, p0[pf], p0[pf]))
  YY <- if (sc=="alt") Ymat else matrix(rbinom(n*nsim,1,rep(p0[pf],nsim)),n,nsim)
  pA <- 1 - as.numeric(pmvnorm(upper=rep(c1,M), mean=mm,
        sigma=if (sc=="alt") Sig else Cc %*% (t(Cc)*(p0[pf]*(1-p0[pf]))),
        algorithm=GenzBretz(abseps=1e-5,maxpts=1e5)))
  rate <- mean(apply(Cc %*% YY >= c1, 2, any))
  chk(sprintf("%-4s computed %.4f vs MC %.4f", sc, pA, rate),
      abs(pA-rate), 4*sqrt(rate*(1-rate)/nsim)+3e-3)
}

hr("CHECK BR4  joint law with the ACTUAL Crump constant-CATE gate")
X3 <- as.matrix(profiles)[pf, ]                  # n x 3 numeric covariates
# arm-difference OLS vector D = b1 - b0 is linear in Y: build its coefficients
Praw <- build_basis(X3, poly_order = 1L)         # package basis (standardized)
K <- ncol(Praw)
Cd <- matrix(0, K, n)
for (a in 0:1) { ia <- which(W==a); Pa <- Praw[ia,,drop=FALSE]
  Cd[, ia] <- Cd[, ia] + (if (a==1) 1 else -1) * solve(crossprod(Pa), t(Pa)) }
mu_D  <- as.numeric(Cd %*% pY)
S_DD  <- Cd %*% (t(Cd)*vY)
S_BD  <- Cc %*% (t(Cd)*vY)                       # cross-covariance (M x K)
Demp  <- Cd %*% Ymat
chk("Crump D: max |emp mean - computed| / SE",
    max(abs(rowMeans(Demp)-mu_D)/sqrt(diag(S_DD)/nsim)), 4)
cx <- cor(t(Bh), t(Demp)) - S_BD/outer(sqrt(diag(Sig)), sqrt(diag(S_DD)))
chk("cross-corr (candidates x D): max |emp - computed|", max(abs(cx)),
    4/sqrt(nsim)+.02)
# gate & declare: Gaussian level (true V) vs data level (package test, plug-in V)
idxs <- 2:K; Vs <- S_DD[idxs, idxs]; Vsi <- solve(Vs); qc <- qchisq(.95, K-1)
Rg <- 2e5
Z  <- rmvnorm(Rg, mean = c(m_c, mu_D), sigma = rbind(cbind(Sig, S_BD),
                                                     cbind(t(S_BD), S_DD)))
decl_g <- apply(Z[,1:M] >= c1, 1, any)
dsl <- Z[, M+idxs, drop=FALSE]
gate_g <- rowSums((dsl %*% Vsi) * dsl) > qc
cat(sprintf("  Gaussian-level: P(gate)=%.4f  P(declare)=%.4f  P(gate & declare)=%.4f\n",
            mean(gate_g), mean(decl_g), mean(gate_g & decl_g)))
gate_d <- decl_d <- logical(nsim)
for (s in 1:nsim) {
  y <- Ymat[, s]
  tc <- test_constant_cate(y, W, X3, poly_order = 1L, regression = "ols")
  gate_d[s] <- !is.null(tc) && is.finite(tc$p_value_chi) && tc$p_value_chi < .05
  decl_d[s] <- any(Bh[, s] >= c1)
}
for (q in list(c("P(gate)", mean(gate_g), mean(gate_d)),
               c("P(gate & declare)", mean(gate_g & decl_g), mean(gate_d & decl_d)))) {
  g <- as.numeric(q[2]); d <- as.numeric(q[3])
  chk(sprintf("%-18s Gaussian %.4f vs data(pkg test) %.4f", q[1], g, d),
      abs(g-d), 4*sqrt(d*(1-d)/nsim)+4*sqrt(g*(1-g)/Rg)+.01,
      "(plug-in V-hat vs true V inside tol)")
}
note("conditional coupling: P(declare | gate) = %.3f vs P(declare) = %.3f (Gaussian)",
     mean(decl_g[gate_g]), mean(decl_g))
cat(sprintf("\nruntime %.1f s\n", proc.time()[3]-t0))
