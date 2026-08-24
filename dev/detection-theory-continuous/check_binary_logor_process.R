# check_binary_logor_process.R
# ---------------------------------------------------------------------------
# B2: the log-OR candidate-effect process, built as a TRANSFORM of the exact
# proportion core (B2 = B1 o logit).
#
# Executing-code anchors (R_23Aug2026 snapshot, verified):
#   * effect_measure "OR" fits logistic regression and RETURNS LOG-OR
#     (R/glm_effect_estimators.R:272); for the unadjusted within-region fit
#     the treatment coefficient equals logit(p1hat) - logit(p0hat) exactly
#     (saturated two-group model; verified numerically in BL7).
#   * Ratio-scale floors are mapped to the log scale
#     (R/forestsearch_main.R:1803-1804), and the screen compares the
#     ESTIMATE to the floor (R/subgroup_search.R:626; floor disabled only
#     under sg_focus = "maxeff").
#   * Binary default adverse_outcome = TRUE (event is bad), so positive
#     log-OR = harm with no negation (R/glm_effect_estimators.R:98).
#
# Structure of the theory under test:
#   CORE (exact, B1-class): the 2M-vector of region-arm event proportions
#     (phat_0(g), phat_1(g))_g is linear in Y with exact design-conditional
#     means pbar and exact heteroskedastic overlap covariance
#     Sigma_prop = C diag(p(1-p)) C', zero across arms.
#   TRANSFORM: theta_hat(g) = logit(phat_1) - logit(phat_0); the delta law
#     is N(h(pbar), J Sigma_prop J') with diagonal-per-(g,arm) Jacobian
#     1/(pbar(1-pbar)). The ONLY new approximation rung is the transform's
#     curvature, governed by the region-arm effective counts
#     n_{g,a} * pbar(1-pbar) (the d_eff halves of the April theory).
#   ERROR DECOMPOSITION (BL4/BL6): (i) delta-Gaussian pmvnorm vs
#     (ii) transform-exact Gaussian-level MC (sample the exact proportion
#     law, then logit) vs (iii) data-level MC. (iii)-(ii) = CLT of the
#     proportion core; (ii)-(i) = linearization.
#
# Base R + mvtnorm. Seed 20260825. Runtime ~15 s.
# ---------------------------------------------------------------------------
suppressMessages(library(mvtnorm))
set.seed(20260825)
t0 <- proc.time()[3]
hr <- function(x) cat("\n", strrep("=", 74), "\n", x, "\n", strrep("=", 74), "\n", sep = "")
chk <- function(lab, v, tol, extra = "") cat(sprintf(
  "  [%s] %-54s %9.3e (tol %.1e)%s\n", if (is.finite(v) && v <= tol) "PASS" else "FAIL",
  lab, v, tol, if (nzchar(extra)) paste0("  ", extra) else ""))
note <- function(...) cat("  ", sprintf(...), "\n", sep = "")

# --- geometry and DGM: identical to B1 (two transforms of one core) --------
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
Pg <- as.numeric(p8 %*% G8)
lab <- vapply(cand[keep], function(cd) paste(sprintf("x%d=%d",cd$v,cd$a),collapse="&"), character(1))
eta0 <- with(profiles, -0.40 + 0.50*x1 - 0.40*x2 + 0.30*x3)
b_tr <- -0.35; b_int <- 1.25
p0 <- plogis(eta0); p1 <- plogis(eta0 + b_tr + b_int*Qp)
c1 <- log(1.5); c2 <- log(1.15)               # floors on the log-OR scale

# Tier B exact reference: marginal region-arm risks and region log-OR
pa_g <- function(pa) as.numeric((p8*pa) %*% G8) / Pg
theta_pop <- qlogis(pa_g(p1)) - qlogis(pa_g(p0))
# The Tier-A-style WRONG answer (collapsed logit-CATE), for the gap readout
logor_cate <- qlogis(p1) - qlogis(p0)
theta_wrong <- as.numeric((p8*logor_cate) %*% G8) / Pg

hr(sprintf("SETUP  M = %d; region log-OR range [%.3f, %.3f]; theta(Q-rule) = %.3f; floors log(1.5)/log(1.15)",
           M, min(theta_pop), max(theta_pop), theta_pop[lab=="x1=1&x2=1"]))

hr("CHECK BL1  Tier B means and the non-collapsibility gap")
note("max over family |theta(g) - E[logOR-CATE|g]| = %.4f  (at Q-rule: %.4f)",
     max(abs(theta_pop - theta_wrong)),
     abs(theta_pop - theta_wrong)[lab=="x1=1&x2=1"])
note("-> the mixture identity fails on log-OR (Tier A -> Tier B), quantified.")

# --- one fixed design; the exact proportion core ---------------------------
n <- 600
pf <- sample.int(8, n, TRUE, p8); W <- rbinom(n,1,.5)
Gn <- G8[pf,,drop=FALSE]
n1g <- colSums(Gn*W); n0g <- colSums(Gn*(1-W))
Cp0 <- t(Gn*(1-W)) / matrix(n0g, M, n)
Cp1 <- t(Gn*W)     / matrix(n1g, M, n)
Cprop <- rbind(Cp0, Cp1)                       # 2M x n : rows (g,arm0) then (g,arm1)
pY <- ifelse(W==1, p1[pf], p0[pf]); vY <- pY*(1-pY)
pbar  <- as.numeric(Cprop %*% pY)              # exact conditional arm means
Sprop <- Cprop %*% (t(Cprop)*vY)               # exact heteroskedastic covariance
h_of  <- function(pv) qlogis(pv[(M+1):(2*M)]) - qlogis(pv[1:M])
theta_c <- h_of(pbar)                          # transform of exact means
Jw <- 1/(pbar*(1-pbar))
J  <- cbind(-diag(Jw[1:M]), diag(Jw[(M+1):(2*M)]))   # M x 2M
Sdel <- J %*% Sprop %*% t(J)

hr("CHECK BL2  proportion core: exactness and structure")
chk("cross-arm covariance block == 0 (machine)",
    max(abs(Sprop[1:M, (M+1):(2*M)])), 1e-15)
nsim <- 4000
Ymat <- matrix(rbinom(n*nsim, 1, rep(pY,nsim)), n, nsim)
Pemp <- Cprop %*% Ymat                          # 2M x nsim proportions
chk("core means: max |emp - exact| / SE",
    max(abs(rowMeans(Pemp)-pbar)/sqrt(diag(Sprop)/nsim)), 4)
chk("core corr: max |emp - exact|",
    max(abs(cor(t(Pemp)) - Sprop/tcrossprod(sqrt(diag(Sprop))))), 4/sqrt(nsim)+.02)

hr("CHECK BL3  delta law of theta-hat vs empirical")
Th <- qlogis(Pemp[(M+1):(2*M),]) - qlogis(Pemp[1:M,])
zerocell <- mean(!is.finite(Th))
note("non-finite theta-hat frequency (zero/full cells): %.5f", zerocell)
Thf <- Th; Thf[!is.finite(Thf)] <- NA
chk("theta means: max |emp - h(pbar)| (delta bias is O(1/n))",
    max(abs(rowMeans(Thf,na.rm=TRUE)-theta_c)), 4*max(sqrt(diag(Sdel)/nsim))+.01)
chk("theta corr: max |emp - delta|",
    max(abs(cor(t(Thf), use="pairwise") - Sdel/tcrossprod(sqrt(diag(Sdel))))),
    4/sqrt(nsim)+.03)
note("min region-arm effective count n_ga*pbar(1-pbar) = %.1f (governs the delta rung)",
     min(c(n0g,n1g)*pbar*(1-pbar)))

hr("CHECK BL4  Rung A on log-OR: three-way decomposition (alt and null)")
repair <- function(S){e<-eigen((S+t(S))/2,symmetric=TRUE)
  e$vectors %*% diag(pmax(e$values,1e-10*max(e$values))) %*% t(e$vectors)}
Rg <- 2e5
run_A <- function(pYs, tag) {
  pb <- as.numeric(Cprop %*% pYs); vs <- pYs*(1-pYs)
  Sp <- Cprop %*% (t(Cprop)*vs)
  th <- h_of(pb); Jl <- 1/(pb*(1-pb))
  Jm <- cbind(-diag(Jl[1:M]), diag(Jl[(M+1):(2*M)]))
  Sd <- Jm %*% Sp %*% t(Jm)
  pA_delta <- 1 - as.numeric(pmvnorm(upper=rep(c1,M), mean=th, sigma=repair(Sd),
                algorithm=GenzBretz(abseps=1e-5,maxpts=1e5)))
  Zp <- rmvnorm(Rg, mean=pb, sigma=repair(Sp))
  Zp <- pmin(pmax(Zp, 1e-9), 1-1e-9)
  Tg <- qlogis(Zp[,(M+1):(2*M),drop=FALSE]) - qlogis(Zp[,1:M,drop=FALSE])
  pA_texact <- mean(apply(Tg >= c1, 1, any))
  Ys <- matrix(rbinom(n*nsim,1,rep(pYs,nsim)), n, nsim)
  Pe <- Cprop %*% Ys
  Td <- qlogis(Pe[(M+1):(2*M),]) - qlogis(Pe[1:M,])
  pA_data <- mean(apply(Td >= c1, 2, any))
  se <- sqrt(pA_data*(1-pA_data)/nsim) + sqrt(pA_texact*(1-pA_texact)/Rg)
  cat(sprintf("  %-5s delta %.4f | transform-exact %.4f | data %.4f\n",
              tag, pA_delta, pA_texact, pA_data))
  chk(sprintf("%s: transform-exact vs data (core CLT)", tag),
      abs(pA_texact-pA_data), 4*se+3e-3)
  chk(sprintf("%s: delta vs transform-exact (linearization)", tag),
      abs(pA_delta-pA_texact), 4*sqrt(pA_texact*(1-pA_texact)/Rg)+8e-3)
  invisible(list(th=th, Sd=Sd, Td=Td, pA_delta=pA_delta))
}
alt <- run_A(pY, "alt")
nul <- run_A(p0[pf], "null")

hr("CHECK BL5  Rung C selection on log-OR (delta pmvnorm vs data)")
sel_pmv <- vapply(1:M, function(m){
  A <- rbind(diag(M)[m,], do.call(rbind, lapply(setdiff(1:M,m),
        function(j) diag(M)[m,]-diag(M)[j,])))
  as.numeric(pmvnorm(lower=c(c1,rep(0,M-1)), upper=rep(Inf,M),
    mean=as.numeric(A %*% alt$th), sigma=repair(A %*% alt$Sd %*% t(A)),
    algorithm=GenzBretz(abseps=5e-4,maxpts=5e4)))}, numeric(1))
Tda <- alt$Td; Tda[!is.finite(Tda)] <- -Inf
win <- max.col(t(Tda)); dec <- apply(Tda >= c1, 2, any)
sel_mc <- tabulate(win[dec], M)/nsim
chk("max_m |pmvnorm - MC selection prob|", max(abs(sel_pmv-sel_mc)),
    4*sqrt(.5/nsim)+8e-3)
chk("sum selection probs vs Rung A (delta) identity",
    abs(sum(sel_pmv)-alt$pA_delta), 8e-3)
top <- order(-sel_pmv)[1:3]
note("top selection: %s", paste(sprintf("%s %.3f", lab[top], sel_pmv[top]), collapse=" | "))

hr("CHECK BL6  small-n stress (n = 220): decomposition + boundary")
n_s <- 220
pf_s <- sample.int(8, n_s, TRUE, p8); W_s <- rbinom(n_s,1,.5)
G_s <- G8[pf_s,,drop=FALSE]
n1s <- colSums(G_s*W_s); n0s <- colSums(G_s*(1-W_s))
Cs <- rbind(t(G_s*(1-W_s))/matrix(n0s,M,n_s), t(G_s*W_s)/matrix(n1s,M,n_s))
pYs <- ifelse(W_s==1, p1[pf_s], p0[pf_s]); vs <- pYs*(1-pYs)
pb_s <- as.numeric(Cs %*% pYs); Sp_s <- Cs %*% (t(Cs)*vs)
th_s <- qlogis(pb_s[(M+1):(2*M)]) - qlogis(pb_s[1:M])
Jl <- 1/(pb_s*(1-pb_s))
Sd_s <- cbind(-diag(Jl[1:M]),diag(Jl[(M+1):(2*M)])) %*% Sp_s %*%
        t(cbind(-diag(Jl[1:M]),diag(Jl[(M+1):(2*M)])))
pA_d <- 1 - as.numeric(pmvnorm(upper=rep(c1,M), mean=th_s, sigma=repair(Sd_s),
          algorithm=GenzBretz(abseps=1e-5,maxpts=1e5)))
Zs <- pmin(pmax(rmvnorm(Rg, pb_s, repair(Sp_s)),1e-9),1-1e-9)
pA_t <- mean(apply(qlogis(Zs[,(M+1):(2*M)])-qlogis(Zs[,1:M]) >= c1, 1, any))
Ys <- matrix(rbinom(n_s*nsim,1,rep(pYs,nsim)), n_s, nsim)
Pe <- Cs %*% Ys
Ts <- qlogis(Pe[(M+1):(2*M),]) - qlogis(Pe[1:M,])
pA_dt <- mean(apply(Ts >= c1, 2, any))
note("min arm count %d; min effective count %.1f; nonfinite theta freq %.4f",
     min(n0s,n1s), min(c(n0s,n1s)*pb_s*(1-pb_s)), mean(!is.finite(Ts)))
cat(sprintf("  delta %.4f | transform-exact %.4f | data %.4f\n", pA_d, pA_t, pA_dt))
chk("n=220: transform-exact vs data (core CLT)", abs(pA_t-pA_dt),
    4*sqrt(pA_dt*(1-pA_dt)/nsim)+6e-3, "(diagnostic at small cells)")
chk("n=220: delta vs transform-exact (linearization)", abs(pA_d-pA_t),
    0.03, "(diagnostic: curvature rung at small effective counts)")

hr("CHECK BL7  estimator tie: glm treatment coefficient == logit difference")
g1 <- which(lab=="x1=1&x2=1")
idx <- Gn[,g1]; y1 <- Ymat[,1][idx]; w1 <- W[idx]
fit <- glm(y1 ~ w1, family=binomial())
chk("glm coef - (logit p1hat - logit p0hat), one region/dataset",
    abs(unname(coef(fit)[2]) - (qlogis(mean(y1[w1==1]))-qlogis(mean(y1[w1==0])))),
    1e-10, "(saturated two-group identity; R/glm_effect_estimators.R:272)")

cat(sprintf("\nruntime %.1f s\n", proc.time()[3]-t0))
