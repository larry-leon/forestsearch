# =============================================================================
# fs_mr_oc_summary() -- operating-characteristics summary of an MR payload
# =============================================================================
#
# Transplanted from the inline `.build_block()` / `.build_block_med()` machinery
# duplicated across the MD40 batch documents.  That duplication is not cosmetic:
# the swapped-oracle-bounds defect fixed in commit 7c517531 sat corrected in one
# document for weeks without ever reaching its own sibling, because five copies
# of the same code drift independently.  One implementation, called by all of
# them, removes that failure mode.
#
# Everything the summary needs is already in the payload -- results, truth and
# meta -- so no DGM argument is required.  fs_dgm_scale() covers the separate
# question of the estimator's sampling scale.
# =============================================================================


#' Operating-Characteristics Summary of a Multiplier-Resampling Payload
#'
#' Summarises a simulation payload produced by the multiplier-resampling (MR)
#' harness into the estimation, coverage and classification quantities used by
#' the analytic documents, returning plain numbers rather than a formatted
#' table.
#'
#' @section Orientation:
#' Estimator columns (\code{or_}, \code{nv_}, \code{fb_}, \code{mr_}) are stored
#' on the ORIENTED scale, where positive means harm.  \code{betaHhat_H} and
#' \code{betaHhat_Hc} are stored on the RAW outcome scale.  The bridge is
#' \code{orient = if (adverse_outcome) 1 else -1}, applied here exactly as the
#' batch documents apply it.  Getting this backwards silently negates every
#' bias.
#'
#' @section Four targets, deliberately:
#' Every estimator is scored against every target in both blocks, so that
#' disagreement is visible rather than assumed away:
#' \describe{
#'   \item{\code{beta}}{Per-replicate exact \eqn{\beta(\hat H)} -- the
#'     conditional estimand at the realised region.}
#'   \item{\code{oracle}}{Per-replicate refit on the true region in that trial.}
#'   \item{\code{struct}}{Scalar structural potential-outcome effect in the true
#'     region; zero Monte Carlo error.}
#'   \item{\code{marg}}{Scalar fitted effect in the true region on one realised
#'     draw; carries sampling noise.}
#' }
#' Bias is ABSOLUTE, in outcome units. On an identity scale a complement effect
#' can sit near zero, so a relative bias would explode.
#'
#' @section Two classification conventions:
#' Both are returned, because the project contains both and they disagree
#' whenever detection is not near certain:
#' \describe{
#'   \item{\code{conditional}}{Averaged over DETECTED replicates only. This is
#'     what the batch documents report.}
#'   \item{\code{unconditional}}{Averaged over ALL replicates, following
#'     \code{build_classification_table()} and Leon et al. (2024, Table 1),
#'     with non-detection scored as sens = 0, ppv = 0, spec = 1, npv = 1.}
#' }
#' With detection at 999/1000 the two agree to three decimals; in a null cell
#' at 0.930 they do not. Quote which one you mean.
#'
#' @param payload A list with \code{results}, \code{truth} and \code{meta}, as
#'   written by the MR batch harness, or a path to such an \code{.rds} file.
#' @param estimators Character vector selecting estimators, from
#'   \code{"oracle"}, \code{"naive"}, \code{"FB"}, \code{"MR"}.  Default
#'   \code{NULL} takes all present, dropping \code{"FB"} when the run carried no
#'   full bootstrap (\code{meta$nb_boots} absent or below 1).
#' @param blocks Character vector of blocks, from \code{"H"} and \code{"Hc"}.
#'   Default both.
#' @param digits Integer or \code{NULL}.  \code{NULL} (default) returns
#'   unrounded values; a formatter should round, not the computation.
#'
#' @return An object of class \code{c("fs_mr_oc", "list")}:
#'   \describe{
#'     \item{\code{estimation}}{Data frame, one row per block and estimator:
#'       \code{n}, \code{avg}, \code{sd_emp}, \code{se_hat}, the four
#'       \code{bias_*}, the four \code{cov_*}, and \code{width}.}
#'     \item{\code{identification}}{Data frame with one row per convention:
#'       detection rate, mean subgroup size, and mean sens/spec/ppv/npv.}
#'     \item{\code{targets}}{The oriented target values used.}
#'     \item{\code{meta}}{Key run metadata, echoed.}
#'   }
#'
#' @seealso \code{\link{fs_dgm_scale}}
#'
#' @examples
#' \dontrun{
#' oc <- fs_mr_oc_summary("fs_maxeffCons_mr_md40_knoise0_n500_res_1_1000.rds")
#' oc$estimation[oc$estimation$block == "H", c("estimator", "avg", "bias_beta")]
#' oc$identification
#' }
#'
#' @export
#' @importFrom stats median sd
fs_mr_oc_summary <- function(payload,
                             estimators = NULL,
                             blocks = c("H", "Hc"),
                             digits = NULL) {

  if (is.character(payload) && length(payload) == 1L) {
    if (!file.exists(payload)) stop("no such payload: ", payload, call. = FALSE)
    payload <- readRDS(payload)
  }
  for (nm in c("results", "truth", "meta")) {
    if (is.null(payload[[nm]])) {
      stop("payload is missing its '", nm, "' block.", call. = FALSE)
    }
  }
  blocks <- match.arg(blocks, c("H", "Hc"), several.ok = TRUE)

  res   <- payload$results
  truth <- payload$truth
  meta  <- payload$meta

  orient <- if (isTRUE(meta$adverse_outcome)) 1 else -1
  rd <- res[res$status %in% "DETECTED", , drop = FALSE]

  keys <- c(oracle = "or", naive = "nv", FB = "fb", MR = "mr")
  run_fb <- !is.null(meta$nb_boots) && isTRUE(meta$nb_boots >= 1L)
  if (!run_fb) keys <- keys[names(keys) != "FB"]
  if (!is.null(estimators)) {
    unknown <- setdiff(estimators, names(keys))
    if (length(unknown)) {
      stop("unknown or unavailable estimator(s): ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    keys <- keys[estimators]
  }

  blk <- list(
    H  = list(beta = orient * rd$betaHhat_H,
              struct = orient * truth$effect_Q,
              marg   = orient * truth$marg_H,
              oracle = rd$or_H_est),
    Hc = list(beta = orient * rd$betaHhat_Hc,
              struct = orient * truth$effect_Qc,
              marg   = orient * truth$marg_Hc,
              oracle = rd$or_Hc_est)
  )

  rows <- list()
  for (sfx in blocks) {
    b <- blk[[sfx]]
    for (lab in names(keys)) {
      k  <- keys[[lab]]
      e  <- rd[[paste0(k, "_", sfx, "_est")]]
      if (is.null(e)) next
      lo <- rd[[paste0(k, "_", sfx, "_lo")]]
      hi <- rd[[paste0(k, "_", sfx, "_hi")]]
      se <- rd[[.fs_oc_se_col(k, sfx)]]
      if (is.null(se)) se <- rep(NA_real_, nrow(rd))

      rows[[length(rows) + 1L]] <- data.frame(
        block       = sfx,
        estimator   = lab,
        n           = sum(is.finite(e)),
        avg         = .fs_oc_mean(e),
        sd_emp      = .fs_oc_sd(e),
        se_hat      = .fs_oc_mean(se),
        bias_beta   = .fs_oc_mean(e - b$beta),
        bias_oracle = .fs_oc_mean(e - b$oracle),
        bias_struct = .fs_oc_mean(e - b$struct),
        bias_marg   = .fs_oc_mean(e - b$marg),
        cov_beta    = .fs_oc_cover(b$beta,   lo, hi),
        cov_oracle  = .fs_oc_cover(b$oracle, lo, hi),
        cov_struct  = .fs_oc_cover(b$struct, lo, hi),
        cov_marg    = .fs_oc_cover(b$marg,   lo, hi),
        width       = .fs_oc_mean(hi - lo),
        stringsAsFactors = FALSE
      )
    }
  }

  est_tab <- if (length(rows)) do.call(rbind, rows) else NULL
  if (!is.null(est_tab) && !is.null(digits)) {
    num <- vapply(est_tab, is.numeric, logical(1))
    est_tab[num] <- lapply(est_tab[num], round, digits = digits)
    rownames(est_tab) <- NULL
  }

  out <- list(
    estimation     = est_tab,
    identification = .fs_oc_identification(res, rd, digits),
    targets        = list(
      orient      = orient,
      struct_H    = orient * truth$effect_Q,
      struct_Hc   = orient * truth$effect_Qc,
      marg_H      = orient * truth$marg_H,
      marg_Hc     = orient * truth$marg_Hc,
      prevalence  = truth$prevalence_Q,
      beta_inter  = orient * truth$beta_inter
    ),
    meta = list(
      n_sample        = meta$n_sample,
      n_sims          = meta$n_sims,
      sg_focus        = meta$sg_focus,
      adverse_outcome = meta$adverse_outcome,
      nb_boots        = meta$nb_boots,
      mr_draws        = meta$mr_draws,
      run_fb          = run_fb,
      pkg_version     = meta$pkg_version
    )
  )
  class(out) <- c("fs_mr_oc", "list")
  out
}


#' Print method for \code{fs_mr_oc} objects
#'
#' @param x An object of class \code{"fs_mr_oc"} from
#'   \code{\link{fs_mr_oc_summary}}.
#' @param ... Unused; present for S3 compatibility.
#'
#' @return The input \code{x}, invisibly.
#'
#' @export
print.fs_mr_oc <- function(x, ...) {
  cat("MR operating-characteristics summary\n")
  cat(sprintf("  n = %s, sims = %s, focus = %s, FB = %s\n",
              x$meta$n_sample, x$meta$n_sims, x$meta$sg_focus, x$meta$run_fb))
  cat(sprintf("  oriented targets: struct_H %.4f  marg_H %.4f\n\n",
              x$targets$struct_H, x$targets$marg_H))
  if (!is.null(x$estimation)) {
    print(x$estimation[, c("block", "estimator", "n", "avg", "sd_emp",
                           "se_hat", "bias_beta", "cov_beta", "width")],
          digits = 5, row.names = FALSE)
  }
  cat("\n")
  print(x$identification, digits = 4, row.names = FALSE)
  invisible(x)
}


# -- internals ----------------------------------------------------------------

#' @keywords internal
#' @noRd
.fs_oc_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

#' @keywords internal
#' @noRd
.fs_oc_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1L) stats::sd(x) else NA_real_
}

#' @keywords internal
#' @noRd
.fs_oc_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

#' @keywords internal
#' @noRd
.fs_oc_cover <- function(target, lo, hi) {
  if (is.null(lo) || is.null(hi)) return(NA_real_)
  target <- rep_len(target, length(lo))
  ok <- is.finite(target) & is.finite(lo) & is.finite(hi)
  if (!any(ok)) return(NA_real_)
  mean(target[ok] >= lo[ok] & target[ok] <= hi[ok])
}

#' Standard-error column name for an estimator and block
#'
#' MR carries an infinitesimal-jackknife standard error under a different
#' suffix from the other estimators.
#' @keywords internal
#' @noRd
.fs_oc_se_col <- function(k, sfx) {
  if (identical(k, "mr")) paste0("mr_", sfx, "_se_ij") else paste0(k, "_", sfx, "_se")
}

#' Classification metrics under both conventions
#' @keywords internal
#' @noRd
.fs_oc_identification <- function(res, rd, digits = NULL) {
  n_all <- nrow(res)
  det   <- res$status %in% "DETECTED"

  # Unconditional: non-detection scores sens = 0, ppv = 0, spec = 1, npv = 1,
  # matching build_classification_table() and Leon et al. (2024, Table 1).
  fill <- function(col, miss) {
    v <- res[[col]]
    if (is.null(v)) return(rep(NA_real_, n_all))
    v[!det | is.na(v)] <- miss
    v
  }

  out <- rbind(
    data.frame(
      convention  = "conditional",
      n_used      = nrow(rd),
      detection   = mean(det),
      mean_size_H = .fs_oc_mean(rd$n_harm),
      sens = .fs_oc_mean(rd$sens), spec = .fs_oc_mean(rd$spec),
      ppv  = .fs_oc_mean(rd$ppv),  npv  = .fs_oc_mean(rd$npv),
      stringsAsFactors = FALSE
    ),
    data.frame(
      convention  = "unconditional",
      n_used      = n_all,
      detection   = mean(det),
      mean_size_H = .fs_oc_mean(rd$n_harm),
      sens = .fs_oc_mean(fill("sens", 0)), spec = .fs_oc_mean(fill("spec", 1)),
      ppv  = .fs_oc_mean(fill("ppv",  0)), npv  = .fs_oc_mean(fill("npv",  1)),
      stringsAsFactors = FALSE
    )
  )
  rownames(out) <- NULL
  if (!is.null(digits)) {
    num <- vapply(out, is.numeric, logical(1))
    out[num] <- lapply(out[num], round, digits = digits)
  }
  out
}
