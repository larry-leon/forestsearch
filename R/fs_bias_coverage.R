# fs_bias_coverage.R
#
# Bias-vs-coverage display as a standard output
# (TASK_bias_coverage_display_2026-09-06; identity fixture
# dev/tasks/bias_coverage_points.csv).  Two add-only functions: a summariser
# turning one combined simulation bundle's per-replicate records into one
# bias/coverage row per estimator, and a three-panel plot over a stack of such
# rows.  The Gaussian reference is bookkeeping, not theory: with
# b = mean log error / empirical SD and r = mean reported SE / empirical SD,
# a Gaussian error model predicts one-sided coverage Phi(z_0.95 * r - b) and
# two-sided coverage Phi(z_0.975 * r - b) - Phi(-z_0.975 * r - b).
# `scale = "identity"` (add-only, TASK_continuous_field_mac_2026-09-07) runs
# the same bookkeeping on the natural scale for difference measures (MD/RD).

#' Bias and coverage summary for one simulation cell
#'
#' From a combined bundle's per-replicate `results` (the recorder columns of
#' the `fb_mr_field` template family: `nv_H_*`, `mr_H_*`, `fld_H_*`,
#' `betaHhat_H`, `detected`, and their `Hc` twins), computes, per estimator
#' over the detected replicates: the retained bias on the log scale, the
#' empirical SD of the log estimate, the mean reported SE, their ratios
#' `b = bias_log / sd_emp` and `r = se_mean / sd_emp`, the observed one- and
#' two-sided coverage of the target with Wilson intervals, and the Gaussian
#' reference coverages at (b, r).
#'
#' Conventions follow the committed template's diagnostics section: the
#' one-sided lower bound is the stored `fld_<block>_lo1s` for the field row
#' and `exp(log est - qnorm(level) * SE)` for naive and MR (IJ); two-sided
#' coverage uses each estimator's stored bounds; the field estimator is the
#' shrunk-field `est2` with the two-sided Lambda-quantile interval.  Under
#' `side = "upper"` (the complement's benefit-claim orientation,
#' TASK_mr_field_complement_2026-09-06) the one-sided coverage is
#' `truth <= upper bound`, with the field's stored `fld_<block>_up1s` and
#' `exp(log est + qnorm(level) * SE)` for naive and MR (IJ), and the Gaussian
#' reference is `Phi(z_level * r + b)` -- a positive retained bias helps an
#' upper bound where it hurts a lower one.
#'
#' @param results Data frame of per-replicate records from a combined bundle
#'   (`readRDS(...)$results`).
#' @param block `"H"` (harm, default) or `"Hc"` (complement).  The `"fld"`
#'   estimator is dropped with a message when the block's field columns are
#'   absent (`fld_Hc_*` exist only for bundles run with
#'   `field_complement = TRUE`).
#' @param estimators Subset of `c("naive", "mr", "fld", "mr_w", "mr_wf")`
#'   (default the first three).  `"mr_w"` / `"mr_wf"` are MR (IJ) with the
#'   winner-only / winner-floor SE (`mr_<block>_se_w` / `_se_wf` and their
#'   `_lo_w/_hi_w` / `_lo_wf/_hi_wf` bounds; TASK_complement_refinements_2026-09-06),
#'   available only on bundles that recorded them.
#' @param level Coverage level (default 0.95).
#' @param target `"betaHhat"` (default; the per-replicate conditional target
#'   `betaHhat_<block>`), `"oracle"` (the per-replicate oracle estimate
#'   `or_<block>_est`), or a numeric scalar on the HR scale (e.g. a CDE or
#'   marginal truth value from the bundle's `truth` element).
#' @param side `"lower"` (default; one-sided coverage of the lower bound,
#'   the harm-claim orientation) or `"upper"` (one-sided coverage of the
#'   upper bound, the benefit-claim orientation).  Two-sided quantities are
#'   unaffected.
#' @param scale `"log"` (default) or `"identity"`.  The working scale of the
#'   estimate columns: `"log"` for ratio measures (HR/OR/RR; bias, SD and the
#'   normal-based one-sided bound `exp(log est -/+ z * SE)` are computed on
#'   the log scale, and non-positive estimates are dropped), `"identity"`
#'   for difference measures (MD/RD; bias, SD and the bound `est -/+ z * SE`
#'   on the natural scale, no positivity guard).  Add-only
#'   (TASK_continuous_field_mac_2026-09-07): the default reproduces the
#'   previous output exactly.  Under `"identity"` the `bias_log` column
#'   holds the identity-scale bias (the column name is kept so stacked
#'   tables and [fs_plot_bias_coverage()] read either scale unchanged).
#' @return Data frame, one row per estimator: `estimator`, `n`, `bias_log`,
#'   `sd_emp`, `se_mean`, `b`, `r`, `cov1`, `cov1_wilson_lo`,
#'   `cov1_wilson_hi`, `cov2`, `cov2_wilson_lo`, `cov2_wilson_hi`,
#'   `cov1_ref`, `cov2_ref`.
#' @export
fs_sim_bias_coverage <- function(results,
                                 block = c("H", "Hc"),
                                 estimators = c("naive", "mr", "fld"),
                                 level = 0.95,
                                 target = "betaHhat",
                                 side = c("lower", "upper"),
                                 scale = c("log", "identity")) {
  block <- match.arg(block)
  side  <- match.arg(side)
  scale <- match.arg(scale)
  # Working-scale transform: log for ratio measures, identity for differences.
  .tr <- if (scale == "log") log else identity
  estimators <- match.arg(estimators, c("naive", "mr", "fld", "mr_w", "mr_wf"),
                          several.ok = TRUE)
  r <- results[results$detected %in% 1L, , drop = FALSE]
  if (!nrow(r)) stop("fs_sim_bias_coverage: no detected replicates")

  fld_pre <- paste0("fld_", block, "_")
  if ("fld" %in% estimators && is.null(r[[paste0(fld_pre, "est2")]])) {
    message("fs_sim_bias_coverage: no field columns (", fld_pre, "*) in ",
            "these results; dropping estimator 'fld' for block = \"",
            block, "\".")
    estimators <- setdiff(estimators, "fld")
  }
  for (k in intersect(c("mr_w", "mr_wf"), estimators)) {
    sfx <- if (k == "mr_w") "w" else "wf"
    if (is.null(r[[paste0("mr_", block, "_se_", sfx)]])) {
      message("fs_sim_bias_coverage: no mr_", block, "_se_", sfx,
              " column in these results; dropping estimator '", k, "'.")
      estimators <- setdiff(estimators, k)
    }
  }

  tgt <- if (is.numeric(target)) rep_len(target, nrow(r))
  else switch(target,
    betaHhat = r[[paste0("betaHhat_", block)]],
    oracle   = r[[paste0("or_", block, "_est")]],
    stop("unknown target: ", target))
  lt <- .tr(tgt)

  z1 <- stats::qnorm(level)
  z2 <- stats::qnorm(1 - (1 - level) / 2)

  wilson <- function(x, n, z = stats::qnorm(0.975)) {
    if (!is.finite(x) || n <= 0) return(c(NA_real_, NA_real_))
    ctr <- (x + z^2 / (2 * n)) / (1 + z^2 / n)
    hw  <- z * sqrt(x * (1 - x) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    c(ctr - hw, ctr + hw)
  }

  # `b1` is the one-sided bound on the requested side: the field's stored
  # lo1s / up1s; the Gaussian bound exp(log est -/+ z1 * SE) otherwise.
  cols <- function(k) {
    if (k == "fld") {
      f <- function(s) r[[paste0(fld_pre, s)]]
      list(e = f("est2"), lo = f("lo2s"), hi = f("hi2s"), se = f("se"),
           b1 = if (side == "lower") f("lo1s") else f("up1s"))
    } else {
      # mr_w / mr_wf share MR (IJ)'s point estimate; only the SE and bounds
      # come from the winner-only / winner-floor columns.
      pre <- if (k == "naive") "nv" else "mr"
      vs  <- switch(k, mr_w = "_w", mr_wf = "_wf", "")
      e  <- r[[paste0(pre, "_", block, "_est")]]
      se <- r[[if (k == "naive") paste0("nv_", block, "_se")
               else if (k == "mr") paste0("mr_", block, "_se_ij")
               else paste0("mr_", block, "_se", vs)]]
      sgn <- if (side == "lower") -1 else 1
      list(e = e,
           lo = r[[paste0(pre, "_", block, "_lo", vs)]],
           hi = r[[paste0(pre, "_", block, "_hi", vs)]],
           se = se,
           b1 = if (scale == "log")
                  ifelse(is.finite(e) & is.finite(se) & e > 0,
                         exp(log(e) + sgn * z1 * se), NA_real_)
                else ifelse(is.finite(e) & is.finite(se),
                            e + sgn * z1 * se, NA_real_))
    }
  }

  out <- lapply(estimators, function(k) {
    cc <- cols(k)
    le <- .tr(cc$e)
    okb <- is.finite(le) & is.finite(lt)
    bias_log <- mean(le[okb] - lt[okb])
    sd_emp   <- stats::sd(le[is.finite(le)])
    se_mean  <- mean(cc$se[is.finite(cc$se)])
    ok2 <- is.finite(lt) & is.finite(cc$lo) & is.finite(cc$hi)
    cov2 <- if (any(ok2)) mean(tgt[ok2] >= cc$lo[ok2] & tgt[ok2] <= cc$hi[ok2]) else NA_real_
    n2 <- sum(ok2)
    ok1 <- is.finite(lt) & is.finite(cc$b1)
    cov1 <- if (!any(ok1)) NA_real_
            else if (side == "lower") mean(tgt[ok1] >= cc$b1[ok1])
            else mean(tgt[ok1] <= cc$b1[ok1])
    n1 <- sum(ok1)
    w1 <- wilson(cov1, n1); w2 <- wilson(cov2, n2)
    b <- bias_log / sd_emp
    rr <- se_mean / sd_emp
    data.frame(
      estimator = k, n = as.integer(max(n1, n2)),
      bias_log = bias_log, sd_emp = sd_emp, se_mean = se_mean,
      b = b, r = rr,
      cov1 = cov1, cov1_wilson_lo = w1[1], cov1_wilson_hi = w1[2],
      cov2 = cov2, cov2_wilson_lo = w2[1], cov2_wilson_hi = w2[2],
      cov1_ref = if (side == "lower") stats::pnorm(z1 * rr - b)
                 else stats::pnorm(z1 * rr + b),
      cov2_ref = stats::pnorm(z2 * rr - b) - stats::pnorm(-z2 * rr - b),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}


#' Bias-vs-coverage display (three panels)
#'
#' Draws the standard bias-vs-coverage display from a stacked
#' [fs_sim_bias_coverage()] table carrying a `cell` column: (1) one-sided
#' coverage against `b` with Gaussian-reference curves `Phi(z1 * r - b)` for
#' the given `curves` values of `r`; (2) two-sided coverage against `b` with
#' curves `Phi(z2 * r - b) - Phi(-z2 * r - b)`; (3) observed coverage against
#' the Gaussian reference at each point's own `(b, r)`, with the diagonal.
#' Estimators are distinguished by marker shape, the nominal level by a
#' dashed line, and cells by text labels.
#'
#' @param tbl Stacked output of [fs_sim_bias_coverage()] with an added
#'   `cell` column (one factor level per simulation cell).
#' @param labels Draw cell labels next to the points (default `TRUE`).
#' @param curves Gaussian-reference `r` values for panels 1-2
#'   (default `c(1, 1.25, 1.5, 2)`).
#' @param level Nominal level for the reference lines (default 0.95).
#' @param side `"lower"` (default) or `"upper"`: which one-sided bound the
#'   table's `cov1` was computed for (pass the same value given to
#'   [fs_sim_bias_coverage()]).  Only panel 1's reference curves and title
#'   change: `Phi(z1 * r + b)` under `"upper"`.
#' @return A patchwork object of the three ggplot panels.
#' @export
fs_plot_bias_coverage <- function(tbl, labels = TRUE,
                                  curves = c(1, 1.25, 1.5, 2),
                                  level = 0.95,
                                  side = c("lower", "upper")) {
  stopifnot(is.data.frame(tbl), "cell" %in% names(tbl))
  side <- match.arg(side)
  z1 <- stats::qnorm(level)
  z2 <- stats::qnorm(1 - (1 - level) / 2)
  bg <- seq(min(tbl$b, -0.5) - 0.3, max(tbl$b, 0.5) + 0.3, length.out = 201)
  sgn1 <- if (side == "lower") -1 else 1
  cur1 <- do.call(rbind, lapply(curves, function(rv)
    data.frame(b = bg, cov = stats::pnorm(z1 * rv + sgn1 * bg), r = factor(rv))))
  cur2 <- do.call(rbind, lapply(curves, function(rv)
    data.frame(b = bg, cov = stats::pnorm(z2 * rv - bg) -
                             stats::pnorm(-z2 * rv - bg), r = factor(rv))))
  lab_layer <- function(mapping) if (isTRUE(labels))
    ggplot2::geom_text(mapping, size = 2.4, vjust = -0.9, na.rm = TRUE,
                       show.legend = FALSE) else NULL
  base_th <- ggplot2::theme_minimal(base_size = 11)

  p1 <- ggplot2::ggplot(tbl, ggplot2::aes(x = b, y = cov1)) +
    ggplot2::geom_line(data = cur1,
                       ggplot2::aes(y = cov, group = r, linewidth = NULL),
                       colour = "grey70") +
    ggplot2::geom_hline(yintercept = level, linetype = 2) +
    ggplot2::geom_point(ggplot2::aes(shape = estimator), size = 2.4) +
    lab_layer(ggplot2::aes(label = cell)) +
    ggplot2::labs(title = sprintf("One-sided %d%% %s-bound coverage",
                                  round(100 * level), side),
                  x = "residual bias b (SD units)", y = "coverage",
                  caption = sprintf("curves: Φ(%.3f·r %s b), r ∈ {%s}",
                                    z1, if (side == "lower") "−" else "+",
                                    paste(curves, collapse = ", "))) +
    base_th

  p2 <- ggplot2::ggplot(tbl, ggplot2::aes(x = b, y = cov2)) +
    ggplot2::geom_line(data = cur2,
                       ggplot2::aes(y = cov, group = r), colour = "grey70") +
    ggplot2::geom_hline(yintercept = level, linetype = 2) +
    ggplot2::geom_point(ggplot2::aes(shape = estimator), size = 2.4) +
    lab_layer(ggplot2::aes(label = cell)) +
    ggplot2::labs(title = sprintf("Two-sided %d%% coverage", round(100 * level)),
                  x = "residual bias b (SD units)", y = "coverage",
                  caption = sprintf("curves: Φ(%.2f·r − b) − Φ(−%.2f·r − b)",
                                    z2, z2)) +
    base_th

  long <- rbind(
    data.frame(cell = tbl$cell, estimator = tbl$estimator, side = "one-sided",
               ref = tbl$cov1_ref, obs = tbl$cov1),
    data.frame(cell = tbl$cell, estimator = tbl$estimator, side = "two-sided",
               ref = tbl$cov2_ref, obs = tbl$cov2))
  p3 <- ggplot2::ggplot(long, ggplot2::aes(x = ref, y = obs)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2,
                         colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(shape = estimator, colour = side),
                        size = 2.4) +
    lab_layer(ggplot2::aes(label = cell)) +
    ggplot2::labs(title = "Observed vs Gaussian reference at each point's own (b, r)",
                  x = "reference coverage", y = "observed coverage") +
    base_th

  patchwork::wrap_plots(p1, p2, p3, ncol = 3)
}
