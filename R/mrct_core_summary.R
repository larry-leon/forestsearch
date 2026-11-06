#' Core operating characteristics of an MRCT simulation campaign
#'
#' The rows of the regional-consistency core table, computed from one
#' [mrct_region_sims()] result run with `record = "full"`. Everything is read from the
#' result: its per-trial columns and its `"estimators"` attribute.
#'
#' The rows, with their denominators:
#' - `power`: P(ITT ub < 1.0), of all trials.
#' - `trigger`: trials reaching the trigger (ITT ub < 1.0 and region ITT > 1.0), of all
#'   trials.
#' - `declared`: a subgroup declared, of the triggered trials.
#' - `sens_nonregion`, `spec_nonregion`, `sens_region`, `spec_region`: the sensitivity
#'   and specificity of \eqn{\hat H} against \eqn{H^*}, trial by trial, averaged over
#'   the declared triggered trials. A trial where one is undefined does not enter its
#'   average.
#' - `declared_Ghat_lt1`: declared and region(\eqn{\hat G}) < 1.0, of the triggered
#'   trials.
#' - `Gstar_lt1`: region(\eqn{G^*}) < 1.0, of the triggered trials.
#' - `declared_Ghat_ge1`: declared and region(\eqn{\hat G}) not below 1.0, of the
#'   triggered trials.
#' - `Ghat_share`: the share of the region's patients kept in \eqn{\hat G}, trial by
#'   trial, averaged over the declared triggered trials.
#'
#' @param res A [mrct_region_sims()] result with `record = "full"`. Take the result
#'   itself: a `data.table` subset of it drops the `"estimators"` attribute.
#'
#' @return A plain data frame with one row per quantity above, in that order:
#'   - `quantity`, `label`: the row's key and its description.
#'   - `estimate`: the value, in percent.
#'   - `se`: its standard error, in percent, as the core table carries it: binomial for
#'     the proportions, `sd / sqrt(n)` over the trials for the Sens and Spec averages,
#'     `NA` for `trigger` and `Ghat_share`, which the core table gives without one.
#'   - `count`: the trials counted (`NA` for an average).
#'   - `denominator`: the trials the row is taken over.
#'   - `display`: the value as printed, at one decimal (for `trigger`, `"K of N"`).
#'
#'   A row with no trial in its denominator has `estimate` `NA` and `display` an en
#'   dash. Without a `flag_harm` column in the simulated trials, the Sens, Spec and
#'   `Gstar_lt1` rows are `NA`.
#'
#' @seealso [mrct_region_sims()], [mrct_select_case()].
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 400
#' d <- data.frame(age = rnorm(n, 50, 10), lab = rnorm(n, 300, 100),
#'                 f1 = rbinom(n, 1, 0.4), treat = rbinom(n, 1, 0.5))
#' d$tte <- rexp(n, 0.05 * exp(-0.5 * d$treat))
#' d$event <- as.integer(d$tte < 40)
#' d$tte <- pmin(d$tte, 40)
#' dgm <- suppressWarnings(inject_mrct_structure(
#'   d, outcome_var = "tte", event_var = "event", treatment_var = "treat",
#'   continuous_vars = c("age", "lab"), factor_vars = "f1", x_pred = "lab",
#'   spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
#'   region = list(prevalence = 0.25, or_pred = 1, x1_vars = "f1", or_x1 = 3,
#'                 loghr = -log(3), band = list(cut = 230, or = 1)),
#'   x3 = list(vars = "lab", cuts = list(lab = 230), loghr = log(3 / 0.55)),
#'   n_super = 2000, expand = "bootstrap", seed = 7L,
#'   continuous_vars_cens = c("age", "lab"), factor_vars_cens = "f1"))
#' dgm <- scale_dgm_hazard(dgm, m = 3)
#' res <- mrct_region_sims(
#'   dgm, n_sims = 8, n_sample = 800, region_var = "z_region", sg_focus = "hrMaxSG",
#'   confounders.name = analysis_covariates(dgm, "observed"), conf_force = character(0),
#'   fs_args = list(use_grf = FALSE, use_lasso = FALSE),
#'   sim_args = list(rand_ratio = 1, draw_treatment = TRUE, max_entry = 6, cens_adjust = 20),
#'   analysis_time = 30, cens_adjust = 20, parallel_args = list(plan = "sequential"),
#'   details = FALSE, verbose_n_sims = 0, seed = 5000L, record = "full")
#' mrct_core_summary(res)[, c("quantity", "estimate", "se", "denominator", "display")]
#' }
#'
#' @export
mrct_core_summary <- function(res) {
  est <- attr(res, "estimators")
  if (is.null(est))
    stop("res has no \"estimators\" attribute: mrct_core_summary() needs a ",
         "mrct_region_sims(record = \"full\") result (record = \"point\" has no per-set ",
         "estimates), and a data.table subset of it drops the attribute.", call. = FALSE)
  demo <- as.data.frame(res)
  needed <- c("sim", "any_found", "TP_nonregion", "FP_nonregion", "FN_nonregion", "TN_nonregion",
              "TP_region", "FP_region", "FN_region", "TN_region")
  absent <- setdiff(needed, names(demo))
  if (length(absent))
    stop("res lacks the columns of a record = \"full\" result: ", paste(absent, collapse = ", "),
         ".", call. = FALSE)
  est <- as.data.frame(est)

  # The definitions: the demo-summary chunk of the simulations article.
  g_est <- function(k) { y <- est[est$estimator == k, ]; y[match(demo$sim, y$sim), ] }
  itt <- g_est("itt_all"); rgn <- g_est("itt_region"); gh <- g_est("region_Ghat"); gs <- g_est("region_Gstar")
  lt1 <- function(v) !is.na(v) & v < 1; gt1 <- function(v) !is.na(v) & v > 1
  A <- lt1(itt$hi); K <- A & gt1(rgn$hr); dec <- demo$any_found == 1; DK <- dec & K
  ss <- function(tp, fn) tp / (tp + fn)
  sens_non <- ss(demo$TP_nonregion, demo$FN_nonregion); spec_non <- ss(demo$TN_nonregion, demo$FP_nonregion)
  sens_reg <- ss(demo$TP_region, demo$FN_region); spec_reg <- ss(demo$TN_region, demo$FP_region)
  has_truth <- any(est$estimator == "region_Gstar")
  N <- nrow(demo)

  dash <- "\u2013"
  pc <- function(x) formatC(100 * x, format = "f", digits = 1)
  out <- function(quantity, label, p, se, count, denominator, display = NULL) {
    defined <- is.finite(p)
    data.frame(quantity = quantity, label = label, estimate = if (defined) 100 * p else NA_real_,
               se = if (defined && !is.null(se)) 100 * se else NA_real_,
               count = as.numeric(count), denominator = as.numeric(denominator),
               display = if (!defined) dash else if (is.null(display)) pc(p) else display,
               stringsAsFactors = FALSE)
  }
  # a proportion over `d`; the binomial SE of the core table
  prop <- function(quantity, label, w, d = K) {
    k <- sum(w & d); m <- sum(d); p <- if (m) k / m else NA_real_
    out(quantity, label, p, if (m) sqrt(p * (1 - p) / m) else NULL, k, m)
  }
  # an average over the declared triggered trials; sd / sqrt(n) over its defined trials
  avg <- function(quantity, label, v, with_se = TRUE) {
    x <- v[DK]; x <- x[!is.na(x)]; m <- length(x)
    p <- if (m) mean(x) else NA_real_
    out(quantity, label, p, if (with_se && m > 1L) stats::sd(x) / sqrt(m) else NULL, NA_real_, m)
  }

  rows <- list(
    prop("power", "trial power (ITT ub < 1.0)", A, rep(TRUE, N)),
    {
      k <- sum(K)
      out("trigger", "trials reaching the trigger", k / N, NULL, k, N, sprintf("%d of %d", k, N))
    },
    prop("declared", "subgroup declared", dec),
    avg("sens_nonregion", "Sens of H-hat vs H*, non-region", sens_non),
    avg("spec_nonregion", "Spec of H-hat vs H*, non-region", spec_non),
    avg("sens_region", "Sens of H-hat vs H*, region", sens_reg),
    avg("spec_region", "Spec of H-hat vs H*, region", spec_reg),
    prop("declared_Ghat_lt1", "declared, region(G-hat) < 1.0", dec & lt1(gh$hr)),
    if (has_truth) prop("Gstar_lt1", "region(G*) < 1.0", lt1(gs$hr)) else
      out("Gstar_lt1", "region(G*) < 1.0", NA_real_, NULL, NA_real_, sum(K)),
    prop("declared_Ghat_ge1", "declared, region(G-hat) >= 1.0", dec & !lt1(gh$hr)),
    avg("Ghat_share", "G-hat, share of the region's patients", gh$n / rgn$n, with_se = FALSE))
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
