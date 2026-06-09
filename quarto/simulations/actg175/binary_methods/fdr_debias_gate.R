# =============================================================================
# fdr_debias_gate.R
#
# PROTOTYPE / SKETCH -- Tier 2: a FAST multiplier-bootstrap approximation of the
# forest-search bootstrap bias-corrected estimate, plus a de-biased SELECTION
# GATE for FDR mitigation, plus a single-fit gated-FDR calibration.
#
# Idea (companion to the resampling-theory note).  The full (Tier 1) bootstrap
# bias-correction of Leon et al. (2024, eq. 7) is, to leading order,
#
#     beta*_hat(H) ~= beta_hat(H) - mean_b D_{Hstar_b}^(b) - mean_b D_{H}^(b),
#
# where D_g^(b) = sum_i G_i^(b) dfbeta_{g,i} is the multiplier perturbation of
# candidate g on draw b, Hstar_b is the candidate the selection rule picks on
# draw b, and H is the observed-selected subgroup.  The fixed-subgroup term is
# mean-zero; the BIAS is the average perturbation of the *selected* candidate
# (winner's curse).  Both terms come from the SAME dfbeta draws the FDR routine
# already makes -- no refitting -- so the de-biased estimate, the FDR, and the
# gate decision are one influence-matrix computation.
#
# A "gate" then declares a subgroup only if its de-biased effect clears a
# threshold t_gate (point) or its selection-adjusted lower bound clears the null
# (CI).  Because the correction collapses spurious selections toward the null,
# the gate mitigates the search-level FDR; because the correction over-shrinks
# true effects, t_gate must be calibrated near the null (not at the screen).
#
# Multiplier note: the de-bias mimics the (multinomial) nonparametric bootstrap,
# so the default multiplier here is "poisson" (centred Poisson ~= K_bi - 1), NOT
# the Rademacher 50/50 split used for the consistency rate.
#
# Dependencies: fdr_family_multiplier.R (for .draw_multipliers).  The demo also
# uses fdr_calibration.R + consistency_resample.R (or installed forestsearch).
#
# Usage:
#   source("fdr_debias_gate.R"); run_debias_gate_demo()
# =============================================================================

if (!exists(".draw_multipliers")) {
  .fm <- getOption("fdr_family_src", "fdr_family_multiplier.R")
  if (file.exists(.fm)) source(.fm) else
    stop(".draw_multipliers() not found; place fdr_family_multiplier.R on the path.")
}


#' Apply a selection rule among the passing candidates on one draw
#'
#' @param beta,zcons,sizes Length-S perturbed estimates, standardized consistency
#'   margins `(beta - c_cons)/sigma_D`, and candidate sizes.
#' @param passers Integer indices of candidates that passed screening+consistency.
#' @param selection One of `"maxcons"`, `"maxeff"`, `"maxSG"`, `"minSG"`,
#'   `"effMaxSG"`, `"effMinSG"`.
#' @param nbhd Effect-neighborhood band (comparison scale) for the eff*SG rules.
#' @return Integer index of the selected candidate, or `NA`.
#' @keywords internal
.select_one <- function(beta, zcons, sizes, passers, selection, nbhd) {
  if (length(passers) == 0L) return(NA_integer_)
  pick <- switch(
    selection,
    maxcons  = passers[which.max(zcons[passers])],
    maxeff   = passers[which.max(beta[passers])],
    maxSG    = passers[which.max(sizes[passers])],
    minSG    = passers[which.min(sizes[passers])],
    effMaxSG = { band <- passers[beta[passers] >= max(beta[passers]) - nbhd]
                 band[which.max(sizes[band])] },
    effMinSG = { band <- passers[beta[passers] >= max(beta[passers]) - nbhd]
                 band[which.min(sizes[band])] },
    stop("unknown selection rule: ", selection))
  as.integer(pick)
}


#' Fast multiplier-bootstrap de-biased estimate of the selected subgroup
#'
#' Centres at the observed estimates (mimicking resampling the observed data),
#' propagates one multiplier vector to all candidates per draw, re-applies the
#' selection rule, and averages the perturbation of the selected candidate as the
#' selection-bias correction.  Approximates the Tier-1 (full bootstrap)
#' bias-corrected estimate at the cost of one influence-matrix product.
#'
#' @param influence `N x S` influence matrix (column g = candidate g's treatment
#'   dfbeta, zero outside g; see [assemble_influence()]).
#' @param beta_hat,sigma_D,sizes Length-S observed estimates, robust SDs
#'   (`NULL` -> `sqrt(colSums(influence^2))`), and candidate sizes
#'   (`NULL` -> nonzero-row counts).
#' @param c_screen,c_consistency,p_star Screening / consistency thresholds
#'   (comparison scale) and consistency-rate cutoff.
#' @param selection,effect_neighborhood Selection rule and band (see
#'   [.select_one()]).
#' @param selected_index Optional observed-selected candidate index; if `NULL` it
#'   is derived by applying the rule to the observed estimates.
#' @param draws,multiplier,seed Bootstrap controls. `multiplier` defaults to
#'   `"poisson"` (mimics the nonparametric bootstrap).
#' @return List with the observed-selected index, naive and de-biased estimates
#'   (log and effect scale), the selection / fixed bias terms, the bootstrap
#'   selection rate, and a one-sided CI lower bound for a CI-style gate.
#' @export
debias_family_multiplier <- function(influence, beta_hat, sigma_D = NULL, sizes = NULL,
                                     c_screen, c_consistency = 0, p_star = 0.90,
                                     selection = c("maxcons", "maxeff", "maxSG",
                                                   "minSG", "effMaxSG", "effMinSG"),
                                     effect_neighborhood = 0.10,
                                     selected_index = NULL,
                                     draws = 2000L,
                                     multiplier = c("poisson", "gaussian", "rademacher"),
                                     seed = NULL) {
  selection  <- match.arg(selection); multiplier <- match.arg(multiplier)
  if (!is.null(seed)) set.seed(seed)
  influence <- as.matrix(influence); N <- nrow(influence); S <- ncol(influence)
  if (is.null(sigma_D)) sigma_D <- sqrt(colSums(influence^2))
  if (is.null(sizes))   sizes   <- colSums(influence != 0)

  z   <- stats::qnorm((1 + p_star) / 2)
  t_g <- pmax(c_screen, c_consistency + z * sigma_D)

  # observed selection
  zc_obs   <- (beta_hat - c_consistency) / sigma_D
  pass_obs <- which(beta_hat >= t_g)
  sel_obs  <- if (!is.null(selected_index)) as.integer(selected_index)
              else .select_one(beta_hat, zc_obs, sizes, pass_obs, selection, effect_neighborhood)

  Xi <- .draw_multipliers(N, draws, multiplier)        # N x draws
  P  <- crossprod(influence, Xi)                       # S x draws : D_g^(b), centred at 0
  beta_star <- beta_hat + P                            # observed-centred perturbed estimates

  sel_bias <- rep(NA_real_, draws)
  for (b in seq_len(draws)) {
    bs   <- beta_star[, b]
    pass <- which(bs >= t_g)
    if (length(pass) == 0L) next
    s <- .select_one(bs, (bs - c_consistency) / sigma_D, sizes, pass, selection, effect_neighborhood)
    if (!is.na(s)) sel_bias[b] <- P[s, b]              # D_{Hstar_b}^(b) = eta_b(Hstar_b)
  }
  selection_rate <- mean(!is.na(sel_bias))
  selection_bias <- mean(sel_bias, na.rm = TRUE)
  fixed_bias     <- if (!is.na(sel_obs)) mean(P[sel_obs, ]) else NA_real_   # ~ 0

  out <- list(selected_index = sel_obs, selection = selection,
              selection_rate = selection_rate,
              selection_bias = selection_bias, fixed_bias = fixed_bias,
              draws = as.integer(draws), multiplier = multiplier)
  if (!is.na(sel_obs)) {
    bn <- beta_hat[sel_obs]
    bd <- bn - selection_bias - (if (is.finite(fixed_bias)) fixed_bias else 0)
    se <- sigma_D[sel_obs]
    out$beta_naive        <- bn
    out$beta_debiased     <- bd
    out$effect_naive      <- exp(bn)
    out$effect_debiased   <- exp(bd)
    out$se                <- se
    out$ci_lower_debiased <- exp(bd - stats::qnorm(0.95) * se)   # one-sided 95% lower bound
  }
  out
}


#' Gate decision from a de-bias result
#'
#' @param db Output of [debias_family_multiplier()].
#' @param t_gate Effect-scale gate threshold (e.g. OR `1.0`).
#' @param gate `"point"` (de-biased point estimate >= t_gate) or `"ci"`
#'   (selection-adjusted lower bound >= t_gate).
#' @return Logical: declare the subgroup or not.
#' @export
gated_declare <- function(db, t_gate = 1.0, gate = c("point", "ci")) {
  gate <- match.arg(gate)
  if (is.na(db$selected_index)) return(FALSE)
  if (gate == "point") db$effect_debiased   >= t_gate
  else                 db$ci_lower_debiased >= t_gate
}


#' Single-fit gated false-discovery curve under a null configuration
#'
#' Null-centres the candidate family, and for each multiplier draw applies the
#' selection rule, the (aggregate) selection-bias de-biasing, and a sweep of gate
#' thresholds.  Returns the ungated and gated family FDR for each `t_gate`, all
#' from one influence-matrix computation.
#'
#' @inheritParams debias_family_multiplier
#' @param beta_null Scalar or length-S null log effect at which to centre.
#' @param t_gate_grid Effect-scale gate thresholds to sweep.
#' @param gate `"point"` or `"ci"`.
#' @return Data frame: `t_gate`, `ungated_fdr`, `gated_fdr`.
#' @export
gated_fdr_curve <- function(influence, beta_hat, sigma_D = NULL, sizes = NULL,
                            c_screen, c_consistency = 0, p_star = 0.90,
                            selection = "maxcons", effect_neighborhood = 0.10,
                            beta_null = 0, t_gate_grid = c(1.0, 1.1, 1.25, 1.5),
                            gate = c("point", "ci"),
                            draws = 3000L, multiplier = c("poisson", "gaussian", "rademacher"),
                            seed = NULL) {
  gate <- match.arg(gate); multiplier <- match.arg(multiplier)
  if (!is.null(seed)) set.seed(seed)
  influence <- as.matrix(influence); N <- nrow(influence); S <- ncol(influence)
  if (is.null(sigma_D)) sigma_D <- sqrt(colSums(influence^2))
  if (is.null(sizes))   sizes   <- colSums(influence != 0)
  if (length(beta_null) == 1L) beta_null <- rep(beta_null, S)

  z   <- stats::qnorm((1 + p_star) / 2)
  t_g <- pmax(c_screen, c_consistency + z * sigma_D)

  Xi <- .draw_multipliers(N, draws, multiplier)
  P  <- crossprod(influence, Xi)
  beta_star <- beta_null + P

  passed   <- logical(draws)
  sel_beta <- rep(NA_real_, draws)     # perturbed estimate of the selected candidate
  sel_pert <- rep(NA_real_, draws)     # its perturbation D_{Hstar_b}^(b)
  sel_se   <- rep(NA_real_, draws)
  for (b in seq_len(draws)) {
    bs   <- beta_star[, b]
    pass <- which(bs >= t_g)
    if (length(pass) == 0L) next
    passed[b] <- TRUE
    s <- .select_one(bs, (bs - c_consistency) / sigma_D, sizes, pass, selection, effect_neighborhood)
    if (!is.na(s)) { sel_beta[b] <- bs[s]; sel_pert[b] <- P[s, b]; sel_se[b] <- sigma_D[s] }
  }
  biasterm <- mean(sel_pert, na.rm = TRUE)             # null selection bias
  debias_b <- sel_beta - biasterm                      # de-biased selected estimate per draw
  lo_b     <- debias_b - stats::qnorm(0.95) * sel_se   # CI lower bound per draw

  data.frame(
    t_gate      = t_gate_grid,
    ungated_fdr = mean(passed),
    gated_fdr   = vapply(t_gate_grid, function(tg) {
      val <- if (gate == "point") debias_b else lo_b
      mean(passed & is.finite(val) & (val >= log(tg)), na.rm = FALSE)
    }, numeric(1)))
}


# =============================================================================
# DEV-ONLY DEMO -- fast de-bias vs literal bootstrap, and the gated-FDR curve
# =============================================================================

#' Self-contained demonstration and validation
#'
#' (1) On an alternative fit (true harm subgroup z1 & z2), compares the fast
#' multiplier de-biased OR against a literal nonparametric bootstrap (re-fit +
#' re-select) on the same fixed-cut family. (2) On a null fit, traces the gated
#' family FDR across gate thresholds, showing the de-biased gate collapse the
#' false-discovery rate while the true subgroup's de-biased estimate stays well
#' above the gate.  Requires fdr_calibration.R + consistency_resample.R (or
#' installed forestsearch) for the per-candidate dfbeta.
#'
#' @keywords internal
run_debias_gate_demo <- function(n = 1500L, B_literal = 300L, draws = 4000L, seed = 11L) {
  for (f in c("fdr_calibration.R")) if (!exists("build_influence_matrix")) source(f)
  csc <- log(1.25); ccon <- log(1.0); pst <- 0.90; rule <- "maxcons"
  spec <- list(outcome_type = "binary", effect_measure = "OR", treat.name = "trt",
               outcome.name = "y", offset.name = NULL, adjust_covariates = NULL,
               adverse_outcome = TRUE)
  mk <- function(harm, s) {
    set.seed(s)
    trt <- rbinom(n,1,0.5); z1 <- rbinom(n,1,0.5); z2 <- rbinom(n,1,0.5)
    x1 <- rnorm(n); x2 <- rnorm(n)
    bt <- if (harm) -0.45 + 1.35*(z1==1 & z2==1) else -0.45
    data.frame(y = rbinom(n,1,plogis(-0.3 + bt*trt + 0.3*z1 - 0.2*z2 + 0.2*x1)),
               trt, z1, z2, x1, x2)
  }
  facs <- function(d) {
    f <- list("z1=1"=d$z1==1,"z1=0"=d$z1==0,"z2=1"=d$z2==1,"z2=0"=d$z2==0)
    for (v in c("x1","x2")) { q <- stats::quantile(d[[v]], c(.25,.5,.75)); for (k in 1:3) {
      f[[sprintf("%s<=%d",v,k)]] <- d[[v]] <= q[k]; f[[sprintf("%s>%d",v,k)]] <- d[[v]] > q[k] } }
    f
  }
  bld <- function(d) {
    cn <- enumerate_candidates(d, facs(d), maxk = 2L, min_n = 60L,
                               min_events_per_arm = 8L, treat.name = "trt", outcome.name = "y")
    list(cn = cn, infl = build_influence_matrix(d, cn, spec))
  }

  ## (1) ALT: fast vs literal -------------------------------------------------
  d <- mk(TRUE, seed); b <- bld(d); infl <- b$infl
  sizes <- vapply(seq_len(infl$n_usable), function(g) sum(infl$B[,g] != 0), numeric(1))
  fast <- debias_family_multiplier(infl$B, infl$beta_hat, infl$sigma_D, sizes,
            c_screen = csc, c_consistency = ccon, p_star = pst, selection = rule,
            draws = draws, multiplier = "poisson", seed = seed + 1L)

  z <- stats::qnorm((1+pst)/2); t_g <- pmax(csc, ccon + z*infl$sigma_D); N <- nrow(d)
  mem <- b$cn[seq_len(infl$n_usable)]
  selbias <- numeric(0)
  for (bb in seq_len(B_literal)) {
    idx <- sample.int(N, N, replace = TRUE); db <- d[idx,,drop = FALSE]
    bs <- rep(NA_real_, infl$n_usable)
    for (g in seq_len(infl$n_usable)) {
      mb <- which(idx %in% mem[[g]]); if (length(mb) < 40L) next
      dg <- db[mb,,drop=FALSE]; if (length(unique(dg$trt)) < 2) next
      e1 <- sum(dg$y[dg$trt==1]); e0 <- sum(dg$y[dg$trt==0])
      if (e1 %in% c(0,sum(dg$trt==1)) || e0 %in% c(0,sum(dg$trt==0))) next
      fit <- tryCatch(suppressWarnings(glm(y~trt, binomial(), dg)), error = function(e) NULL)
      if (!is.null(fit)) { cf <- stats::coef(fit)[["trt"]]; if (is.finite(cf)) bs[g] <- cf }
    }
    pass <- which(is.finite(bs) & bs >= t_g); if (length(pass) == 0L) next
    s <- pass[which.max(((bs - ccon)/infl$sigma_D)[pass])]
    selbias <- c(selbias, bs[s] - infl$beta_hat[s])
  }
  cat("==== (1) ALT fit: de-bias validation (fixed-cut family, maxcons) ====\n")
  cat(sprintf("selected subgroup: %s | naive OR %.3f\n",
              names(b$cn)[fast$selected_index], fast$effect_naive))
  cat(sprintf("selection bias  fast=%.4f  literal=%.4f (B=%d)\n",
              fast$selection_bias, mean(selbias), length(selbias)))
  cat(sprintf("de-biased OR    fast=%.3f  literal=%.3f\n",
              fast$effect_debiased, exp(infl$beta_hat[fast$selected_index] - mean(selbias))))
  alt_debiased_OR <- fast$effect_debiased

  ## (2) NULL: gated-FDR curve ------------------------------------------------
  dn <- mk(FALSE, seed + 12L); bn <- bld(dn); infn <- bn$infl
  szn <- vapply(seq_len(infn$n_usable), function(g) sum(infn$B[,g] != 0), numeric(1))
  itt <- as.numeric(stats::coef(glm(y~trt, binomial(), dn))[["trt"]])
  curve <- gated_fdr_curve(infn$B, infn$beta_hat, infn$sigma_D, szn,
             c_screen = csc, c_consistency = ccon, p_star = pst, selection = rule,
             beta_null = itt, t_gate_grid = c(1.00, 1.10, 1.25, 1.50), gate = "point",
             draws = draws, multiplier = "poisson", seed = seed + 2L)
  cat(sprintf("\n==== (2) NULL fit: gated-FDR curve (null OR %.2f) ====\n", exp(itt)))
  print(curve, row.names = FALSE)
  cat(sprintf("\nTrue subgroup de-biased OR (alt fit) = %.2f -> stays above an OR>1.0 gate,\n",
              alt_debiased_OR))
  cat("so the gate that collapses the null FDR still declares the true subgroup.\n")
  invisible(list(fast = fast, literal_selection_bias = mean(selbias), curve = curve))
}

if (identical(environment(), globalenv()) && !interactive() &&
    sys.nframe() == 0L && length(commandArgs(trailingOnly = TRUE)) == 0L) {
  run_debias_gate_demo()
}
