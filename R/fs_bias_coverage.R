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
#' one-sided lower bound is the stored `fld_H_lo1s` for the field row and
#' `exp(log est - qnorm(level) * SE)` for naive and MR (IJ); two-sided
#' coverage uses each estimator's stored bounds; the field estimator is the
#' shrunk-field `est2` with the two-sided Lambda-quantile interval, and it
#' exists for the harm block only.
#'
#' @param results Data frame of per-replicate records from a combined bundle
#'   (`readRDS(...)$results`).
#' @param block `"H"` (harm, default) or `"Hc"` (complement; the `"fld"`
#'   estimator is dropped there with a message -- the gate computes no
#'   complement field).
#' @param estimators Subset of `c("naive", "mr", "fld")`.
#' @param level Coverage level (default 0.95).
#' @param target `"betaHhat"` (default; the per-replicate conditional target
#'   `betaHhat_<block>`), `"oracle"` (the per-replicate oracle estimate
#'   `or_<block>_est`), or a numeric scalar on the HR scale (e.g. a CDE or
#'   marginal truth value from the bundle's `truth` element).
#' @return Data frame, one row per estimator: `estimator`, `n`, `bias_log`,
#'   `sd_emp`, `se_mean`, `b`, `r`, `cov1`, `cov1_wilson_lo`,
#'   `cov1_wilson_hi`, `cov2`, `cov2_wilson_lo`, `cov2_wilson_hi`,
#'   `cov1_ref`, `cov2_ref`.
#' @export
fs_sim_bias_coverage <- function(results,
                                 block = c("H", "Hc"),
                                 estimators = c("naive", "mr", "fld"),
                                 level = 0.95,
                                 target = "betaHhat") {
  block <- match.arg(block)
  estimators <- match.arg(estimators, c("naive", "mr", "fld"),
                          several.ok = TRUE)
  r <- results[results$detected %in% 1L, , drop = FALSE]
  if (!nrow(r)) stop("fs_sim_bias_coverage: no detected replicates")

  if (identical(block, "Hc") && "fld" %in% estimators) {
    message("fs_sim_bias_coverage: the field block exists for the harm ",
            "block only; dropping estimator 'fld' for block = \"Hc\".")
    estimators <- setdiff(estimators, "fld")
  }

  tgt <- if (is.numeric(target)) rep_len(target, nrow(r))
  else switch(target,
    betaHhat = r[[paste0("betaHhat_", block)]],
    oracle   = r[[paste0("or_", block, "_est")]],
    stop("unknown target: ", target))
  lt <- log(tgt)

  z1 <- stats::qnorm(level)
  z2 <- stats::qnorm(1 - (1 - level) / 2)

  wilson <- function(x, n, z = stats::qnorm(0.975)) {
    if (!is.finite(x) || n <= 0) return(c(NA_real_, NA_real_))
    ctr <- (x + z^2 / (2 * n)) / (1 + z^2 / n)
    hw  <- z * sqrt(x * (1 - x) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    c(ctr - hw, ctr + hw)
  }

  cols <- function(k) {
    if (k == "fld") list(e = r$fld_H_est2, lo = r$fld_H_lo2s,
                         hi = r$fld_H_hi2s, se = r$fld_H_se,
                         lo1 = r$fld_H_lo1s)
    else {
      e  <- r[[paste0(if (k == "mr") "mr" else "nv", "_", block, "_est")]]
      se <- r[[if (k == "mr") paste0("mr_", block, "_se_ij")
               else paste0("nv_", block, "_se")]]
      list(e = e,
           lo = r[[paste0(if (k == "mr") "mr" else "nv", "_", block, "_lo")]],
           hi = r[[paste0(if (k == "mr") "mr" else "nv", "_", block, "_hi")]],
           se = se,
           lo1 = ifelse(is.finite(e) & is.finite(se) & e > 0,
                        exp(log(e) - z1 * se), NA_real_))
    }
  }

  out <- lapply(estimators, function(k) {
    cc <- cols(k)
    le <- log(cc$e)
    okb <- is.finite(le) & is.finite(lt)
    bias_log <- mean(le[okb] - lt[okb])
    sd_emp   <- stats::sd(le[is.finite(le)])
    se_mean  <- mean(cc$se[is.finite(cc$se)])
    ok2 <- is.finite(lt) & is.finite(cc$lo) & is.finite(cc$hi)
    cov2 <- if (any(ok2)) mean(tgt[ok2] >= cc$lo[ok2] & tgt[ok2] <= cc$hi[ok2]) else NA_real_
    n2 <- sum(ok2)
    ok1 <- is.finite(lt) & is.finite(cc$lo1)
    cov1 <- if (any(ok1)) mean(tgt[ok1] >= cc$lo1[ok1]) else NA_real_
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
      cov1_ref = stats::pnorm(z1 * rr - b),
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
#' @return A patchwork object of the three ggplot panels.
#' @export
fs_plot_bias_coverage <- function(tbl, labels = TRUE,
                                  curves = c(1, 1.25, 1.5, 2),
                                  level = 0.95) {
  stopifnot(is.data.frame(tbl), "cell" %in% names(tbl))
  z1 <- stats::qnorm(level)
  z2 <- stats::qnorm(1 - (1 - level) / 2)
  bg <- seq(min(tbl$b, -0.5) - 0.3, max(tbl$b, 0.5) + 0.3, length.out = 201)
  cur1 <- do.call(rbind, lapply(curves, function(rv)
    data.frame(b = bg, cov = stats::pnorm(z1 * rv - bg), r = factor(rv))))
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
    ggplot2::labs(title = sprintf("One-sided %d%% lower-bound coverage", round(100 * level)),
                  x = "residual bias b (SD units)", y = "coverage",
                  caption = sprintf("curves: Φ(%.3f·r − b), r ∈ {%s}",
                                    z1, paste(curves, collapse = ", "))) +
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
