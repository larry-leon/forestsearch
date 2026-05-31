.normalize_sg_focus <- function(x) switch(x, eff="hr", effMaxSG="hrMaxSG",
                                          effMinSG="hrMinSG", x)
.validate_effect_neighborhood <- function(x) invisible(TRUE)
.validate_selection_rule <- function(rule, focus, en) invisible(TRUE)
.compute_inclusion_band <- function(hr_vec, n_vec, selection_rule,
                                     effect_neighborhood)
  as.integer(hr_vec >= (1 - effect_neighborhood) * max(hr_vec))
coef.dina <- function(object, ...) object$coef
vcov.dina <- function(object, ...) object$vcov
source("dina_subgroup.R")

med_time <- function(expr, reps = 5) {
  e <- substitute(expr); pf <- parent.frame()
  ts <- replicate(reps, system.time(eval(e, pf))[["elapsed"]])
  stats::median(ts)
}

## ---------------- Part A: gbsg, real covariate cardinalities ----------------
g <- survival::gbsg
covs <- c("age", "size", "nodes", "pgr", "er")
X <- as.matrix(g[covs]); storage.mode(X) <- "double"
Z <- scale(X)
# illustrative log-HR tau: harm rises with nodes, falls with pgr
beta_cov <- c(age = 0.10, size = 0.05, nodes = 0.25, pgr = -0.20, er = -0.05)
tau_hat <- as.numeric(Z %*% beta_cov)
# wrap as a "fit": coef on the RAW scale so dina_subgroup recomputes tau_hat
# from X identically. Recover raw-scale coef from standardized beta.
ctr <- attr(Z, "scaled:center"); scl <- attr(Z, "scaled:scale")
b_raw <- beta_cov / scl
b0    <- -sum(beta_cov * ctr / scl)
fit <- structure(list(family = "cox",
                      coef = c(`(Intercept)` = b0, b_raw),
                      vcov = diag(rep(1e-4, length(covs) + 1L))),
                 class = "dina")
stopifnot(max(abs(as.numeric(fit$coef[1] + X %*% fit$coef[-1]) - tau_hat)) < 1e-8)

m_diff <- log(1.3); n_min <- 60L
cat("=== Part A: gbsg (N=", nrow(g), ", K=", length(covs),
    " covariates), m_diff=log(1.3), n_min=", n_min, " ===\n", sep = "")

grids <- list(quartiles = c(.25, .5, .75),
              deciles   = seq(.1, .9, .1),
              ventiles  = seq(.05, .95, .05))

t1 <- med_time(dina_subgroup(fit, g, covs, m_diff = m_diff, n_min = n_min,
                             max_depth = 1L))
s1 <- dina_subgroup(fit, g, covs, m_diff = m_diff, n_min = n_min, max_depth = 1L)
cat(sprintf("\n depth-1            : %6.1f ms | searched %5d | qualifying %5d | %s\n",
            1000*t1, s1$n_candidates_searched, s1$n_candidates_qualifying, s1$label))

for (nm in names(grids)) {
  gp <- grids[[nm]]
  tt <- med_time(dina_subgroup(fit, g, covs, m_diff = m_diff, n_min = n_min,
                               max_depth = 2L, grid_probs = gp))
  ss <- dina_subgroup(fit, g, covs, m_diff = m_diff, n_min = n_min,
                      max_depth = 2L, grid_probs = gp)
  cat(sprintf(" depth-2 %-9s  : %6.1f ms | searched %5d | qualifying %5d | depth=%d %s (n=%d, HR=%.2f)\n",
              nm, 1000*tt, ss$n_candidates_searched, ss$n_candidates_qualifying,
              ss$depth, ss$label, ss$n_subgroup, exp(ss$mean_tau_hat)))
}

## ---------------- Part B: scaling in K (covariates) -------------------------
cat("\n=== Part B: scaling, deciles, direction='both', maxSG ===\n")
cat(sprintf(" %-5s %-6s %10s %12s %10s\n", "N", "K", "pairs_seen", "qualifying", "time_ms"))
set.seed(7)
for (N in c(686, 2000)) {
  for (K in c(5, 10, 15, 20)) {
    Xs <- matrix(runif(N * K), N, K)
    cs <- paste0("x", seq_len(K))
    colnames(Xs) <- cs
    dfx <- as.data.frame(Xs)
    bb  <- c(0, runif(K, -0.3, 0.3))
    fk  <- structure(list(family = "gaussian",
                          coef = setNames(bb, c("(Intercept)", cs)),
                          vcov = diag(rep(1e-4, K + 1L))), class = "dina")
    thh <- as.numeric(bb[1] + Xs %*% bb[-1])
    tt  <- med_time(.dina_collect_candidates_depth2(
              Xs, thh, cs, m_diff = -Inf, n_min = 60L, direction = "both",
              grid_probs = seq(.1, .9, .1)), reps = 3)
    hh  <- .dina_collect_candidates_depth2(Xs, thh, cs, m_diff = -Inf,
              n_min = 60L, direction = "both", grid_probs = seq(.1, .9, .1))
    cat(sprintf(" %-5d %-6d %10d %12d %10.1f\n",
                N, K, hh$n_searched, hh$n_qualifying, 1000*tt))
  }
}
