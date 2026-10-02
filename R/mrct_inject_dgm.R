## =============================================================================
## mrct_inject_dgm.R
##
## Inject the Zhang, Long, Bornkamp & Hou (2026, arXiv:2605.16885) regional
## heterogeneity structure into a real seed trial (cgd, actg175, ...) and
## return a forestsearch `aft_dgm_flex` object that simulate_from_dgm() and
## mrct_region_sims() consume unchanged.
##
## Covariate classes (Figure 1 of the paper):
##   X1  region-imbalanced, NOT an effect modifier     -> `x1_vars` (+ OR each)
##   X2  region-associated effect modifier             -> `x_pred` (spline TE + region logit)
##   X3  effect modifier, balanced across regions      -> `x3` (subgroup interaction)
##   U   unobserved factor; Region acts as proxy       -> analysis_covariates(case="unobserved")
##                                                        or `region_treat_loghr` (Region-U-TE)
##
## Region model (Section 4.2 of the paper, extended to several X1 terms):
##   logit P(Region = 1 | x) = a0 + log(OR_pred) * s(x_pred) + sum_j log(OR_j) * s(x1_j)
##   with s() in [0,1] (ECDF by default, or the paper's min-max) on the
##   super-population and a0 solved so that mean P(Region = 1) = `prevalence`.
##
## Expansion (any n): a Gaussian copula on the seed covariates (latent
## correlation via latentcor when available, else normal-score correlation)
## generates n_super novel covariate rows that preserve the seed's marginals
## and dependence.  This is what benchtm does for the paper's simulations.
##
## Pipeline
##   1. expand covariates -> df_pop (n_super rows)
##   2. draw Region on df_pop (and on the seed rows, needed for the AFT fit)
##   3. generate_aft_dgm_flex() on the seed: prognostic gamma from the data,
##      TE shape from spline_spec, Region prognostic log-HR via set_beta_spec,
##      X3 interaction via subgroup_vars + set_beta_spec("treat_harm")
##   4. swap dgm$df_super for df_pop with linear predictors recomputed through
##      the package's own calculate_linear_predictors()/prepare_censoring_model()
##
## Verified 2026-09-03 on survival::cgd (n = 128) and speff2trial::ACTG175
## (n = 2139) with n_super = 20000 and simulated trials of n = 500/700/1000.
## =============================================================================

## -----------------------------------------------------------------------------
## 0. Small utilities
## -----------------------------------------------------------------------------
.minmax <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

.is_binary <- function(x) {
  u <- unique(x[!is.na(x)])
  length(u) == 2L && all(u %in% c(0, 1))
}


## -----------------------------------------------------------------------------
## 1. Covariate expansion: Gaussian copula on the seed rows
## -----------------------------------------------------------------------------
#' Expand a seed covariate frame to n_super rows via a Gaussian copula.
#'
#' @param seed_df   data.frame of the seed trial (only `vars` are used)
#' @param vars      covariates to expand (numeric or 0/1 binary)
#' @param n_super   rows to generate
#' @param seed      RNG seed
#' @param integer_vars variables to round back to integers (auto-detected if NULL)
#' @param method    "copula" (default) or "bootstrap" (row resampling, as the
#'                  package does by default -- no novel covariate patterns)
#' @return data.frame with n_super rows and columns `vars`
#' @importFrom mvtnorm rmvnorm
#' @importFrom Matrix nearPD
#' @importFrom stats quantile pnorm complete.cases
#' @export
expand_covariates <- function(seed_df, vars, n_super, seed = 20260903L,
                              integer_vars = NULL,
                              method = c("copula", "bootstrap")) {
  method <- match.arg(method)
  X <- seed_df[, vars, drop = FALSE]
  X <- X[stats::complete.cases(X), , drop = FALSE]
  set.seed(seed)

  if (method == "bootstrap") {
    idx <- sample.int(nrow(X), n_super, replace = TRUE)
    out <- X[idx, , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }

  p      <- length(vars)
  binary <- vapply(X, .is_binary, logical(1))
  if (is.null(integer_vars))
    integer_vars <- vars[!binary & vapply(X, function(v) all(v == round(v)), logical(1))]

  ## latent correlation --------------------------------------------------------
  R <- NULL
  if (requireNamespace("latentcor", quietly = TRUE)) {
    R <- tryCatch({
      types <- ifelse(binary, "bin", "con")
      lc <- latentcor::latentcor(X, types = types, method = "original",
                                 use.nearPD = TRUE, showplot = FALSE)
      as.matrix(lc$R)
    }, error = function(e) NULL)
  }
  if (is.null(R)) {
    ## normal-score (van der Waerden) correlation as a fallback
    Zs <- sapply(X, function(v) stats::qnorm((rank(v, ties.method = "average") - 0.5) / length(v)))
    R  <- stats::cor(Zs)
    R  <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
  }
  dimnames(R) <- list(vars, vars)

  ## draw latent normals and map back through the seed marginals ------------
  Z <- mvtnorm::rmvnorm(n_super, sigma = R)
  colnames(Z) <- vars
  U <- stats::pnorm(Z)

  out <- as.data.frame(matrix(NA_real_, n_super, p, dimnames = list(NULL, vars)))
  for (j in seq_len(p)) {
    v <- vars[j]
    if (binary[j]) {
      pj <- mean(X[[v]])
      out[[v]] <- as.numeric(U[, j] > 1 - pj)
    } else {
      out[[v]] <- as.numeric(stats::quantile(X[[v]], probs = U[, j], type = 8,
                                             names = FALSE))
      if (v %in% integer_vars) out[[v]] <- round(out[[v]])
    }
  }
  attr(out, "latent_R") <- R
  out
}


## -----------------------------------------------------------------------------
## 2. Region model
## -----------------------------------------------------------------------------
#' Build the logistic region model and draw Region on a population.
#'
#' @param pop      data.frame holding x_pred and x1_vars (natural scale)
#' @param x_pred   name of the region-associated effect modifier (X2); NULL for
#'                 no X2 term (regional imbalance only through x1_vars)
#' @param or_pred  odds ratio for s(x_pred) -> Region  (OR = 1: none)
#' @param x1_vars  names of imbalanced-only covariates (X1); character(0) for none
#' @param or_x1    odds ratios for the X1 terms (recycled)
#' @param prevalence target marginal P(Region = 1)
#' @param scale_ref data.frame used to define the scaling (defaults to pop)
#' @param scale    "rank" (ECDF, default) or "minmax" (the paper's scaling)
#' @return list(alpha0, coefs, prob = function(df), draw = function(df, seed))
#' @importFrom stats plogis uniroot rbinom ecdf
#' @export
make_region_model <- function(pop, x_pred = NULL, or_pred = 1,
                              x1_vars = character(0), or_x1 = 1,
                              prevalence = 0.20, scale_ref = NULL,
                              scale = c("rank", "minmax")) {
  scale <- match.arg(scale)
  if (is.null(scale_ref)) scale_ref <- pop
  terms <- c(x_pred, x1_vars)
  ors   <- c(if (!is.null(x_pred)) or_pred else NULL,
             if (length(x1_vars)) rep_len(or_x1, length(x1_vars)) else NULL)
  names(ors) <- terms

  ## s(x) in [0,1], defined on the reference population so that the same
  ## transformation applies to seed rows and expanded rows.
  ##   "rank"  : empirical CDF of the reference population (default).  Robust
  ##             to skewed covariates such as CD4 counts; OR is then the odds
  ##             ratio between the top and bottom of the distribution.
  ##   "minmax": (x - min) / (max - min), the paper's parameterisation; with
  ##             heavy-tailed covariates most of the mass sits near 0 and the
  ##             induced imbalance is much weaker than the OR suggests.
  ref <- lapply(terms, function(v) scale_ref[[v]]); names(ref) <- terms
  scale_one <- function(x, v) {
    r <- ref[[v]]
    if (.is_binary(r)) return(x)
    if (scale == "minmax") {
      rg <- range(r, na.rm = TRUE)
      if (diff(rg) == 0) rep(0, length(x)) else (x - rg[1]) / diff(rg)
    } else {
      stats::ecdf(r)(x)
    }
  }
  scale_fn <- function(df) {
    if (!length(terms)) return(matrix(0, nrow(df), 0))
    S <- sapply(terms, function(v) scale_one(df[[v]], v))
    matrix(S, nrow = nrow(df), dimnames = list(NULL, terms))
  }

  eta_no_int <- function(df) as.vector(scale_fn(df) %*% log(ors))
  eta_pop    <- eta_no_int(pop)

  ## solve a0 so that mean(plogis(a0 + eta)) = prevalence
  f  <- function(a) mean(stats::plogis(a + eta_pop)) - prevalence
  a0 <- stats::uniroot(f, interval = c(-25, 25), tol = 1e-10)$root

  prob_fn <- function(df) stats::plogis(a0 + eta_no_int(df))

  list(
    alpha0     = a0,
    coefs      = log(ors),
    ors        = ors,
    prevalence = prevalence,
    scale      = scale,
    prob       = prob_fn,
    draw       = function(df, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      stats::rbinom(nrow(df), 1L, prob_fn(df))
    }
  )
}


## -----------------------------------------------------------------------------
## 3. Main: inject the MRCT structure and build the DGM
## -----------------------------------------------------------------------------
#' Build an aft_dgm_flex with synthetic regions and paper-style effect structure.
#'
#' @param seed_data     seed trial data.frame (one row per patient)
#' @param outcome_var,event_var,treatment_var  column names in seed_data
#' @param continuous_vars,factor_vars  covariates entering the AFT model
#'   (factor_vars must be 0/1 numeric or 2-level factors for clean expansion)
#' @param x_pred        region-associated effect modifier (X2), a continuous_var
#' @param spline_spec   list(knot, zeta, log_hrs) on the NATURAL scale of x_pred:
#'   true log-HR is piecewise linear through (0, `log_hrs[1]`), (knot, `log_hrs[2]`),
#'   (zeta, `log_hrs[3]`).  log_hrs constant => no heterogeneity (paper's beta1 = 0).
#' @param region        list(name = "region", prevalence = 0.2, or_pred = 10,
#'   x1_vars = character(0), or_x1 = 1, loghr = 0, scale = "rank") -- loghr is
#'   the prognostic log-HR of Region (no TE modification; 0 for the paper's
#'   DGM; -log(5) reproduces the Feb-2026 deck's strongly prognostic region);
#'   scale = "rank" (ECDF) or "minmax" (the paper's `[0,1]` scaling)
#' @param x3            NULL, or list(vars = c(...), cuts = list(...), loghr = ...)
#'   for a balanced effect modifier: cuts as in generate_aft_dgm_flex (use
#'   fixed values), loghr = log-HR change (treatment x subgroup).
#' @param region_treat_loghr  NULL, or a log-HR for a direct treat x Region
#'   interaction (the residual Region-U-TE pathway).  Mutually exclusive with x3
#'   (the package supports one interaction term).
#' @param n_super       size of the expanded super-population
#' @param expand        "copula" or "bootstrap"
#' @param continuous_vars_cens,factor_vars_cens,cens_type  passed through
#' @param seed          RNG seed
#' @param verbose       logical
#' @return aft_dgm_flex with df_super = expanded population (+ z_<region>), and
#'   a `$mrct` element documenting the injected structure.
#' @importFrom stats complete.cases setNames
#' @importFrom utils modifyList
#' @export
inject_mrct_structure <- function(seed_data,
                                  outcome_var, event_var, treatment_var,
                                  continuous_vars, factor_vars,
                                  x_pred,
                                  spline_spec,
                                  region = list(),
                                  x3 = NULL,
                                  region_treat_loghr = NULL,
                                  n_super = 5000L,
                                  expand = c("copula", "bootstrap"),
                                  continuous_vars_cens = NULL,
                                  factor_vars_cens = NULL,
                                  cens_type = "weibull",
                                  seed = 8316951L,
                                  verbose = FALSE) {
  expand <- match.arg(expand)
  region <- utils::modifyList(
    list(name = "region", prevalence = 0.20, or_pred = 10,
         x1_vars = character(0), or_x1 = 1, loghr = 0, scale = "rank"), region)
  if (!is.null(x3) && !is.null(region_treat_loghr))
    stop("x3 and region_treat_loghr are mutually exclusive: the AFT DGM ",
         "carries a single treatment interaction term.")
  stopifnot(x_pred %in% continuous_vars)

  ## keep only complete seed rows on the variables we use -------------------
  use_vars <- unique(c(outcome_var, event_var, treatment_var, continuous_vars,
                       factor_vars, continuous_vars_cens, factor_vars_cens,
                       if (!is.null(x3)) x3$vars))
  seed_df  <- seed_data[stats::complete.cases(seed_data[, use_vars]), use_vars]
  rownames(seed_df) <- NULL
  covs <- unique(c(continuous_vars, factor_vars, continuous_vars_cens,
                   factor_vars_cens, if (!is.null(x3)) x3$vars))
  ## factors -> 0/1
  for (v in covs) if (is.factor(seed_df[[v]])) {
    if (nlevels(seed_df[[v]]) != 2L)
      stop("Factor '", v, "' has ", nlevels(seed_df[[v]]),
           " levels; recode to 0/1 dummies before calling.")
    seed_df[[v]] <- as.numeric(seed_df[[v]] == levels(seed_df[[v]])[2])
  }

  ## 1. expand covariates ----------------------------------------------------
  df_pop <- expand_covariates(seed_df, covs, n_super, seed = seed, method = expand)

  ## 2. region -----------------------------------------------------------------
  rm_ <- make_region_model(df_pop, x_pred = x_pred, or_pred = region$or_pred,
                           x1_vars = region$x1_vars, or_x1 = region$or_x1,
                           prevalence = region$prevalence, scale_ref = df_pop,
                           scale = region$scale)
  df_pop[[region$name]]  <- rm_$draw(df_pop,  seed = seed + 1L)
  seed_df[[region$name]] <- rm_$draw(seed_df, seed = seed + 2L)
  ## guard: survreg needs both region levels in the seed
  if (length(unique(seed_df[[region$name]])) < 2L)
    seed_df[[region$name]][sample.int(nrow(seed_df), 2)] <- c(0, 1)

  ## 3. AFT DGM on the seed (prognostic structure from data) ------------------
  z_region <- paste0("z_", region$name)
  set_var  <- z_region
  beta_var <- region$loghr
  subgroup_vars <- NULL; subgroup_cuts <- NULL; model <- "alt"
  if (!is.null(x3)) {
    subgroup_vars <- x3$vars; subgroup_cuts <- x3$cuts
    set_var  <- c(set_var, "treat_harm"); beta_var <- c(beta_var, x3$loghr)
  } else if (!is.null(region_treat_loghr)) {
    subgroup_vars <- region$name
    subgroup_cuts <- stats::setNames(list(list(type = "multiple", values = 1)), region$name)
    set_var  <- c(set_var, "treat_harm"); beta_var <- c(beta_var, region_treat_loghr)
  }

  # generate_aft_dgm_flex() defaults *_cens to the outcome covariates when the
  # censoring model is fitted; replicate that so df_pop gets identical columns
  if (is.null(factor_vars_cens) && is.null(continuous_vars_cens)) factor_vars_cens <- factor_vars
  if (is.null(continuous_vars_cens)) continuous_vars_cens <- continuous_vars

  dgm <- generate_aft_dgm_flex(
    data                 = seed_df,
    continuous_vars      = continuous_vars,
    factor_vars          = c(factor_vars, region$name),
    continuous_vars_cens = continuous_vars_cens,
    factor_vars_cens     = factor_vars_cens,
    set_beta_spec        = list(set_var = set_var, beta_var = beta_var),
    outcome_var          = outcome_var,
    event_var            = event_var,
    treatment_var        = treatment_var,
    subgroup_vars        = subgroup_vars,
    subgroup_cuts        = subgroup_cuts,
    model                = model,
    k_inter              = 1,
    n_super              = nrow(seed_df),     # placeholder; replaced below
    cens_type            = cens_type,
    spline_spec          = list(var = paste0("z_", x_pred),
                                knot = spline_spec$knot, zeta = spline_spec$zeta,
                                log_hrs = spline_spec$log_hrs),
    seed                 = seed,
    verbose              = verbose,
    standardize          = FALSE
  )
  mp <- dgm$model_params

  ## 4. rebuild df_super on the expanded population --------------------------
  df_pop[[outcome_var]]   <- 1          # dummies: prepare_working_dataset needs them
  df_pop[[event_var]]     <- 1
  df_pop[[treatment_var]] <- rep_len(c(1, 0), nrow(df_pop))
  df_work <- prepare_working_dataset(df_pop, outcome_var, event_var, treatment_var,
                       continuous_vars, c(factor_vars, region$name), FALSE,
                       continuous_vars_cens, factor_vars_cens, verbose = FALSE)
  sg <- define_subgroups(df_work, df_pop, subgroup_vars, subgroup_cuts,
                continuous_vars, model, verbose = FALSE)
  df_work$flag_harm <- sg$flag_harm
  df_work <- create_spline_variables(df_work, mp$spline_info$var, mp$spline_info$knot)
  covariate_cols <- grep("^z_", names(df_work), value = TRUE)
  df_work <- calculate_linear_predictors(df_work, covariate_cols, mp$gamma, mp$b0, mp$spline_info)

  ## censoring linear predictors: re-run the package's censoring prep with the
  ## seed working frame (deterministic refit) and the new population
  seed_work <- mp$spline_info$df_work
  cr <- prepare_censoring_model(df_work = seed_work, cens_type = cens_type, cens_params = list(),
                  df_super = df_work, select_censoring = TRUE,
                  cens_intercept_only = FALSE, verbose = FALSE)
  df_super <- cr$df_super
  df_super$id <- seq_len(nrow(df_super))
  df_super[[outcome_var]] <- NULL; df_super[[event_var]] <- NULL
  df_super$y <- NULL; df_super$event <- NULL

  dgm$df_super <- df_super
  dgm$n_super  <- nrow(df_super)
  dgm$subgroup_info$size       <- sum(df_super$flag_harm)
  dgm$subgroup_info$proportion <- mean(df_super$flag_harm)
  set.seed(seed + 3L)
  dgm$hazard_ratios <- calculate_hazard_ratios(df_super, nrow(df_super), mp$mu, mp$tau, model, verbose)

  ## 5. document the injected structure ---------------------------------------
  dgm$mrct <- list(
    seed_n        = nrow(seed_df),
    expand        = expand,
    region        = list(name = region$name, z_name = z_region,
                         prevalence_target = region$prevalence,
                         prevalence_pop = mean(df_super[[z_region]]),
                         alpha0 = rm_$alpha0, coefs = rm_$coefs, scale = region$scale,
                         prognostic_loghr = region$loghr,
                         treat_interaction_loghr = region_treat_loghr),
    classes       = list(X1 = region$x1_vars, X2 = x_pred,
                         X3 = if (!is.null(x3)) x3$vars else character(0),
                         other = setdiff(c(continuous_vars, factor_vars),
                                         c(x_pred, region$x1_vars,
                                           if (!is.null(x3)) x3$vars))),
    spline        = spline_spec,
    region_model  = rm_,
    latent_R      = attr(df_pop, "latent_R")
  )
  dgm
}


## -----------------------------------------------------------------------------
## 4. Analysis-side helpers
## -----------------------------------------------------------------------------
#' Covariate set handed to forestsearch()/mrct_region_sims() for the paper's
#' observed vs unobserved analysis cases.
#' @export
analysis_covariates <- function(dgm, case = c("observed", "unobserved"),
                                drop = character(0)) {
  case <- match.arg(case)
  z <- grep("^z_", names(dgm$df_super), value = TRUE)
  z <- setdiff(z, c(dgm$mrct$region$z_name,
                    grep("_treat$|_k$", z, value = TRUE)))   # spline helper cols
  if (case == "unobserved") z <- setdiff(z, paste0("z_", dgm$mrct$classes$X2))
  setdiff(z, paste0("z_", drop))
}

#' Truth by region on the super-population: prevalence, x_pred shift,
#' causal AHR / CDE, and a large-sample Cox HR on UNCENSORED potential
#' outcomes.  Compare HR_cox_PO with the censored-trial HR from
#' sim_region_metrics(): a gap between them isolates the prognostic-shift /
#' follow-up-truncation mechanism (the Feb-2026 deck) from covariate shift
#' (AHR itself differing by region).
#' @importFrom stats rexp coef
#' @export
mrct_truth_by_region <- function(dgm, seed = 1L) {
  ds <- dgm$df_super
  zr <- dgm$mrct$region$z_name
  xp <- paste0("z_", dgm$mrct$classes$X2)
  mp <- dgm$model_params
  set.seed(seed)
  eps <- log(stats::rexp(nrow(ds)))
  T1 <- exp(mp$mu + mp$tau * eps + ds$lin_pred_1)
  T0 <- exp(mp$mu + mp$tau * eps + ds$lin_pred_0)
  po <- data.frame(time = c(T1, T0), event = 1,
                   treat = rep(c(1, 0), each = nrow(ds)), region = rep(ds[[zr]], 2))
  one <- function(idx, lab) {
    d <- ds[idx, ]; p <- po[po$region == (lab == "REGION") | lab == "ALL", ]
    hr <- exp(stats::coef(survival::coxph(survival::Surv(time, event) ~ treat, data = p)))
    data.frame(stratum = lab, n = nrow(d), prev = round(nrow(d) / nrow(ds), 3),
               x_pred_mean = mean(d[[xp]]), AHR = exp(mean(d$loghr_po)),
               CDE = mean(exp(d$theta_1)) / mean(exp(d$theta_0)), HR_cox_PO = unname(hr))
  }
  out <- rbind(one(seq_len(nrow(ds)), "ALL"),
               one(which(ds[[zr]] == 0), "NON-REGION"),
               one(which(ds[[zr]] == 1), "REGION"))
  smd <- (mean(ds[[xp]][ds[[zr]] == 1]) - mean(ds[[xp]][ds[[zr]] == 0])) /
    stats::sd(ds[[xp]])
  attr(out, "smd_x_pred") <- smd
  out
}

#' True log-HR psi0(x) along x_pred with region histograms (paper Fig. 3 analogue)
#' @importFrom grDevices png dev.off adjustcolor
#' @importFrom graphics par plot hist abline legend
#' @export
plot_region_modifier <- function(dgm, file = NULL, main = NULL) {
  ds <- dgm$df_super; zr <- dgm$mrct$region$z_name
  xp <- paste0("z_", dgm$mrct$classes$X2)
  o  <- order(ds[[xp]])
  if (!is.null(file)) grDevices::png(file, width = 900, height = 800, res = 120)
  op <- graphics::par(mfrow = c(2, 1), mar = c(4, 4.5, 3, 1))
  on.exit({ graphics::par(op); if (!is.null(file)) grDevices::dev.off() })
  graphics::plot(ds[[xp]][o], ds$loghr_po[o], type = "l", lwd = 3,
                 xlab = "", ylab = "true log(HR)  psi0(x)",
                 main = if (is.null(main)) paste("True causal log-HR along", xp) else main)
  graphics::abline(h = 0, lty = 3, col = "red")
  graphics::abline(h = log(dgm$hazard_ratios$AHR), lty = 2)
  brks <- pretty(ds[[xp]], 40)
  h0 <- graphics::hist(ds[[xp]][ds[[zr]] == 0], breaks = brks, plot = FALSE)
  h1 <- graphics::hist(ds[[xp]][ds[[zr]] == 1], breaks = brks, plot = FALSE)
  graphics::plot(h0, col = grDevices::adjustcolor("firebrick", 0.5), border = NA,
                 xlab = xp, main = "Distribution by region",
                 ylim = c(0, max(h0$counts, h1$counts)))
  graphics::plot(h1, col = grDevices::adjustcolor("steelblue", 0.7), border = NA, add = TRUE)
  graphics::legend("topright", fill = c(grDevices::adjustcolor("firebrick", 0.5),
                                        grDevices::adjustcolor("steelblue", 0.7)),
                   legend = c("Region = 0", "Region = 1"), bty = "n")
  invisible(NULL)
}

#' Region-consistency metrics on a simulated trial (ITT and optionally a subgroup)
#' @importFrom stats coef vcov
#' @export
sim_region_metrics <- function(sim, region_col, subgroup_index = NULL,
                               pi_preserve = 0.5) {
  hr <- function(d) {
    if (nrow(d) < 10 || sum(d$event_sim) < 5)
      return(c(hr = NA_real_, lo = NA_real_, hi = NA_real_, n = nrow(d), d = sum(d$event_sim)))
    f <- survival::coxph(survival::Surv(y_sim, event_sim) ~ treat_sim, data = d)
    b <- unname(stats::coef(f))[1]; se <- unname(sqrt(stats::vcov(f)[1, 1]))
    c(hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
      n = nrow(d), d = sum(d$event_sim))
  }
  r  <- sim[[region_col]] == 1
  ov <- hr(sim); nr <- hr(sim[!r, ]); rg <- hr(sim[r, ])
  ratio <- log(rg["hr"]) / log(ov["hr"])
  out <- list(overall = ov, non_region = nr, region = rg,
              consistency_ratio = unname(ratio),
              method1_pass = isTRUE(unname(ratio) >= pi_preserve),
              method2_pass = isTRUE(sign(log(rg["hr"])) == sign(log(ov["hr"]))),
              region_hr_below_1 = isTRUE(rg["hr"] < 1))
  if (!is.null(subgroup_index)) {
    out$region_subgroup     <- hr(sim[r & subgroup_index, ])
    out$non_region_subgroup <- hr(sim[!r & subgroup_index, ])
  }
  out
}
