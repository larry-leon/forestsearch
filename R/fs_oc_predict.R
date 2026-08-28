# =============================================================================
# fs_oc_predict() -- predicted operating characteristics of the maxeffCons
#                    search over a candidate family, by joint-normal draws
# -----------------------------------------------------------------------------
# The body is the `worked-predictions` chunk / predict_scn() of
# quarto/simulations/actg175/continuous/analytic_verification_and_prediction_
# md_harm.qmd with its family construction removed.  The operations are kept
# in the same order so that, given the same family, the same seed and the same
# draw count, the returned quantities are bit-identical to the chunk's:
#
#   Rho/Sg from ovl and se_g -> fs_sym_root(Sg, scale = 2) -> W1, W2
#   -> Bhat = (W1 + W2) / 2 -> the gate -> det / det_rate
#   -> Bmask / max.col argmax (maxeffCons) -> P1, p_sel, sel_c
#   -> EnH, Esens, Espec, Eppv, Enpv, EbetaH, Enaive_bias, mass_below.
#
# Gate:
#   "split"     (W1 + W2 >= 2*c1) & (W1 >= c2) & (W2 >= c2)   -- the document's
#               historical single-split gate: one half-sample pair.
#   "resample"  the package's production screen (R/consistency_resample.R):
#               rate = 2*Phi((beta_hat - c2) / sigma_D) - 1  >=  pconsistency
#               <=>  beta_hat >= c2 + qnorm((1 + pconsistency)/2) * sigma_D,
#               combined with the effect screen beta_hat >= c1.  sigma_D is
#               sqrt(sum(dfbeta[i, treat]^2)) (L86, L318), the sandwich SE of
#               the subgroup treatment coefficient; on one simulated MD40 trial
#               it equals the direct dfbeta sum of squares to 1e-10 and sits
#               within 4% of the model SE and within 7% of the wrapper's
#               population se_g with no prevalence trend
#               (dev/glm-continuous-sims/sigma_d_diagnostic_2026-08-29.R), so
#               the wrapper identifies sigma_D = se_g with no extra factor.
#               This branch draws only Bhat ~ N(beta_g, Sg): one matrix, not two.
# =============================================================================


#' Predicted operating characteristics of the search over a candidate family
#'
#' Monte-Carlo prediction of the maxeffCons search's operating characteristics
#' -- family declaration rate, per-candidate declaration, the selection
#' distribution, expected selected-subgroup size, sensitivity, specificity,
#' PPV, NPV, the oriented effect on the selected rule and the naive selection
#' bias -- from the joint normal law of the candidates' fitted effects.  The
#' family is a population enumeration from \code{\link{fs_oc_family_enumerate}}
#' or a supplied \code{fs_oc_family} object.
#'
#' @details
#' The full-sample fitted effects of the M candidates are taken as
#' \eqn{N(\beta_g, S_g)} with \eqn{S_g = \rho \circ (se_g se_g')},
#' \eqn{\rho_{ij} = P(g_i \cap g_j) / \sqrt{P(g_i) P(g_j)}}.  Two independent
#' half-sample draws \code{W1}, \code{W2}, each with covariance \code{2 * Sg},
#' are generated through \code{\link{fs_sym_root}}; the full-sample effect is
#' \code{Bhat = (W1 + W2) / 2}.
#'
#' \strong{The two gates.} Family declaration is any candidate declaring; the
#' selected rule is the effect maximiser among declaring candidates
#' (maxeffCons).
#' \describe{
#'   \item{\code{"resample"} -- the package's production screen}{The
#'     consistency stage of \code{\link{forestsearch}} (default
#'     \code{consistency_method = "resample"}) represents the random 50/50
#'     split's half effects as \eqn{\hat\beta \pm D} and computes the
#'     consistency rate in closed form as
#'     \eqn{2\Phi((\hat\beta - c_2)/\sigma_D) - 1}, a candidate passing when
#'     that rate is at least \code{pconsistency}
#'     (\code{R/consistency_resample.R}; \code{R/subgroup_consistency_helpers.R}
#'     drops a candidate when \code{p.consistency < pconsistency.threshold}).
#'     Inverting: \eqn{\hat\beta \ge c_2 + z_p \sigma_D} with
#'     \eqn{z_p = \Phi^{-1}((1+p)/2)}.  \eqn{\sigma_D} is
#'     \code{sqrt(sum(dfbeta[i, treat]^2))}, the sandwich standard error of the
#'     subgroup treatment coefficient; in the draw model above the analogous
#'     quantity is \code{se_g}, and on a simulated MD40 trial the two agree
#'     within a few percent with no prevalence trend
#'     (\code{dev/glm-continuous-sims/sigma_d_diagnostic_2026-08-29.R}), so the
#'     wrapper identifies \eqn{\sigma_D = se_g}.  The gate is then a single
#'     threshold on the full-sample draw,
#'     \code{(Bhat >= c1) & (Bhat - c2 >= z_p * se_g)}, and only
#'     \code{Bhat ~ N(beta_g, Sg)} is drawn -- one matrix, not two.}
#'   \item{\code{"split"} -- the analytic document's historical gate}{One
#'     half-sample pair \code{W1}, \code{W2} stands in for the consistency
#'     rate: \code{(W1 + W2 >= 2 * c1) & (W1 >= c2) & (W2 >= c2)}, the
#'     screening floor \code{c1} on the full-sample effect and the consistency
#'     floor \code{c2} on each half.  Kept bit-identical to the document's
#'     \code{worked-predictions} chunk.}
#' }
#'
#' All expectations are selection-weighted population functionals of the
#' family, conditional on declaration.  Monte-Carlo standard errors are given
#' for the proportions (\code{sqrt(p (1 - p) / draws)}).
#'
#' @param dgm An object of class \code{"glm_dgm"}; used only when
#'   \code{family = NULL}.
#' @param forestsearch_args Named list of \code{\link{forestsearch}} arguments.
#'   Supplies the default \code{c1} (\code{effect.threshold}) and \code{c2}
#'   (\code{consistency.threshold}) and, when \code{family = NULL}, the cut
#'   specification and floors for \code{\link{fs_oc_family_enumerate}}.
#' @param n Integer.  Trial size.  Overrides any size implied by the
#'   arguments; sets the size floor and the standard-error scale of an
#'   enumerated family and converts expected prevalence to expected subjects.
#' @param c1 Numeric.  Screening floor on the full-sample effect.  Default
#'   \code{forestsearch_args$effect.threshold}; an explicit value overrides.
#' @param c2 Numeric.  Consistency floor on each half-sample effect.  Default
#'   \code{forestsearch_args$consistency.threshold}; an explicit value
#'   overrides.
#' @param family \code{NULL} (enumerate with
#'   \code{\link{fs_oc_family_enumerate}}) or an \code{fs_oc_family} object
#'   used as-is.
#' @param consistency_method \code{"resample"} (the package's production
#'   closed-form screen, the default) or \code{"split"} (the analytic
#'   document's single-split gate).  See Details.
#' @param pconsistency Numeric in (0, 1).  Consistency-rate threshold for the
#'   \code{"resample"} gate.  Default
#'   \code{forestsearch_args$pconsistency.threshold}, itself defaulting to
#'   \code{forestsearch()}'s default (0.90); an explicit value overrides.
#'   Ignored by \code{"split"}.
#' @param draws Integer.  Number of Monte-Carlo draws.
#' @param seed Integer or \code{NULL}.  Passed to \code{set.seed()} before the
#'   draws when non-\code{NULL}.
#' @param ... Passed to \code{\link{fs_oc_family_enumerate}} when
#'   \code{family = NULL} (\code{max_M}, \code{verbose}).
#'
#' @return An object of class \code{c("fs_oc_predict", "list")}:
#'   \describe{
#'     \item{\code{det_rate}, \code{det_rate_se}}{Family declaration rate and
#'       its MC standard error.}
#'     \item{\code{P1}, \code{P1_se}}{Per-candidate declaration probability.}
#'     \item{\code{p_sel}, \code{p_sel_se}}{Per-candidate selection
#'       probability (unconditional).}
#'     \item{\code{sel_c}}{Selection distribution given declaration.}
#'     \item{\code{EnH}}{Expected selected-subgroup size in subjects,
#'       \code{n * sum(sel_c * Pg)}.}
#'     \item{\code{Esens}, \code{Espec}, \code{Eppv}, \code{Enpv}}{Expected
#'       classification metrics of the selected rule against Q.}
#'     \item{\code{EbetaH}}{Expected oriented true effect on the selected
#'       rule.}
#'     \item{\code{Enaive_bias}}{Expected naive minus true effect on the
#'       selected rule.}
#'     \item{\code{mass_below}}{Selection mass on rules whose true mean is
#'       below \code{c1}.}
#'     \item{\code{M}, \code{lab}}{Family size and labels.}
#'     \item{\code{settings}}{\code{n}, \code{c1}, \code{c2},
#'       \code{consistency_method}, \code{pconsistency} (\code{NA} for
#'       \code{"split"}), \code{draws}.}
#'     \item{\code{seed}}{The seed used.}
#'     \item{\code{family}}{The family object.}
#'   }
#'
#' @seealso \code{\link{fs_oc_family_enumerate}}, \code{\link{fs_sym_root}}
#'
#' @examples
#' # A hand-built three-candidate family
#' piQ <- 0.34
#' fam <- structure(list(
#'   lab = c("Q", "P", "D"), Pg = c(piQ, 0.45, 0.31),
#'   PQg = c(1, 0.28 / 0.45, 1), sens_g = c(1, 0.28 / piQ, 0.31 / piQ),
#'   spec_g = c(1, 1 - 0.17 / (1 - piQ), 1),
#'   ovl = matrix(c(piQ, 0.28, 0.31, 0.28, 0.45, 0.28 * 0.31 / piQ,
#'                  0.31, 0.28 * 0.31 / piQ, 0.31), 3, 3),
#'   M = 3L, PQ = piQ), class = c("fs_oc_family", "list"))
#' fam$beta_g <- 26 + 14 * fam$PQg
#' fam$se_g   <- 13.7 * sqrt(2) * sqrt(piQ / fam$Pg)
#' fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
#'               consistency_method = "resample", draws = 2e4, seed = 1)
#' fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
#'               consistency_method = "split", draws = 2e4, seed = 1)
#'
#' @export
fs_oc_predict <- function(dgm = NULL, forestsearch_args = list(), n,
                          c1 = NULL, c2 = NULL, family = NULL,
                          consistency_method = c("resample", "split"),
                          pconsistency = NULL,
                          draws = 2e5, seed = NULL, ...) {

  consistency_method <- match.arg(consistency_method)
  if (identical(consistency_method, "resample")) {
    if (is.null(pconsistency)) pconsistency <- forestsearch_args$pconsistency.threshold
    if (is.null(pconsistency)) pconsistency <- eval(formals(forestsearch)$pconsistency.threshold)
    if (!is.numeric(pconsistency) || length(pconsistency) != 1L ||
        is.na(pconsistency) || pconsistency <= 0 || pconsistency >= 1) {
      stop("'pconsistency' must be a single number in (0, 1).", call. = FALSE)
    }
  } else {
    pconsistency <- NA_real_
  }

  if (missing(n) || !is.numeric(n) || length(n) != 1L || is.na(n) || n <= 0) {
    stop("'n' must be a single positive number.", call. = FALSE)
  }
  if (is.null(c1)) c1 <- forestsearch_args$effect.threshold
  if (is.null(c2)) c2 <- forestsearch_args$consistency.threshold
  if (is.null(c1) || !is.numeric(c1) || length(c1) != 1L || is.na(c1)) {
    stop("'c1' is required: pass it explicitly or supply ",
         "forestsearch_args$effect.threshold.", call. = FALSE)
  }
  if (is.null(c2) || !is.numeric(c2) || length(c2) != 1L || is.na(c2)) {
    stop("'c2' is required: pass it explicitly or supply ",
         "forestsearch_args$consistency.threshold.", call. = FALSE)
  }
  if (!is.numeric(draws) || length(draws) != 1L || is.na(draws) || draws < 1) {
    stop("'draws' must be a single positive number.", call. = FALSE)
  }

  # ---- the family: enumerate, or use the supplied one as-is -----------------
  if (is.null(family)) {
    if (is.null(dgm)) {
      stop("either 'family' or 'dgm' (with 'forestsearch_args') is required.",
           call. = FALSE)
    }
    family <- fs_oc_family_enumerate(dgm, forestsearch_args, n, ...)
  } else if (!inherits(family, "fs_oc_family")) {
    stop("'family' must be an object of class 'fs_oc_family'.", call. = FALSE)
  }
  .fs_oc_check_family(family)

  lab    <- family$lab
  Pg     <- family$Pg
  PQg    <- family$PQg
  beta_g <- family$beta_g
  se_g   <- family$se_g
  sens_g <- family$sens_g
  spec_g <- family$spec_g
  ovl    <- family$ovl
  M      <- family$M
  PQ     <- family$PQ

  # ---- S2: covariance across candidates ------------------------------------
  Rho <- ovl / sqrt(outer(Pg, Pg))
  Sg  <- Rho * outer(se_g, se_g)

  if (!is.null(seed)) set.seed(seed)
  Rdraw <- draws
  if (identical(consistency_method, "resample")) {
    # ---- production gate: one draw of the full-sample effects, cov = Sg ----
    root_full <- fs_sym_root(Sg, scale = 1)
    Bhat <- matrix(stats::rnorm(Rdraw * M), Rdraw, M) %*% t(root_full) +
            rep(beta_g, each = Rdraw)
    # rate = 2*Phi((Bhat - c2)/sigma_D) - 1 >= p  <=>  Bhat - c2 >= z_p*sigma_D,
    # with sigma_D identified as se_g; combined with the effect screen at c1.
    z_p  <- stats::qnorm((1 + pconsistency) / 2)
    pass <- (Bhat >= c1) & (Bhat - c2 >= z_p * rep(se_g, each = Rdraw))
  } else {
    # ---- draws: two half-sample effects, cov = 2 * Sg each ----------------
    root_half <- fs_sym_root(Sg, scale = 2)
    W1 <- matrix(stats::rnorm(Rdraw * M), Rdraw, M) %*% t(root_half) +
          rep(beta_g, each = Rdraw)
    W2 <- matrix(stats::rnorm(Rdraw * M), Rdraw, M) %*% t(root_half) +
          rep(beta_g, each = Rdraw)
    Bhat <- (W1 + W2) / 2

    # ---- S3: the gate -------------------------------------------------------
    pass <- (W1 + W2 >= 2 * c1) & (W1 >= c2) & (W2 >= c2)
  }
  det  <- rowSums(pass) > 0

  # ---- S4: maxeffCons selection ---------------------------------------------
  Bmask <- Bhat; Bmask[!pass] <- -Inf
  winner <- max.col(Bmask, ties.method = "first"); winner[!det] <- NA_integer_
  P1 <- colMeans(pass)
  p_sel <- tabulate(winner[det], M) / Rdraw
  det_rate <- mean(det)
  sel_c <- p_sel / det_rate

  # ---- S5, S6: the functionals ----------------------------------------------
  EnH   <- n * sum(sel_c * Pg)
  Esens <- sum(sel_c * sens_g); Espec <- sum(sel_c * spec_g)
  Eppv  <- sum(sel_c * PQg)
  Enpv  <- sum(sel_c * vapply(seq_len(M), function(m) {
    pHc <- 1 - Pg[m]; (pHc - (PQ - PQg[m] * Pg[m])) / pHc }, numeric(1)))
  EbetaH <- sum(sel_c * beta_g)
  noise  <- Bhat - rep(beta_g, each = Rdraw)
  Enaive_bias <- mean(noise[cbind(which(det), winner[det])])
  mass_below <- sum(sel_c[beta_g < c1])

  out <- list(
    det_rate    = det_rate,
    det_rate_se = .fs_oc_mc_se_prop(det_rate, Rdraw),
    P1          = P1,
    P1_se       = .fs_oc_mc_se_prop(P1, Rdraw),
    p_sel       = p_sel,
    p_sel_se    = .fs_oc_mc_se_prop(p_sel, Rdraw),
    sel_c       = sel_c,
    EnH         = EnH,
    Esens       = Esens,
    Espec       = Espec,
    Eppv        = Eppv,
    Enpv        = Enpv,
    EbetaH      = EbetaH,
    Enaive_bias = Enaive_bias,
    mass_below  = mass_below,
    M           = M,
    lab         = lab,
    settings    = list(n = n, c1 = c1, c2 = c2,
                       consistency_method = consistency_method,
                       pconsistency = pconsistency,
                       draws = draws),
    seed        = seed,
    family      = family
  )
  class(out) <- c("fs_oc_predict", "list")
  out
}


#' @export
print.fs_oc_predict <- function(x, ...) {
  s <- x$settings
  cat("Predicted operating characteristics (fs_oc_predict)\n")
  cat(sprintf("  family      : M = %d;  n = %s;  gate = %s (c1 = %s, c2 = %s%s);  draws = %s%s\n",
              x$M, format(s$n), s$consistency_method, format(s$c1),
              format(s$c2),
              if (is.na(s$pconsistency)) "" else sprintf(", pcons = %s", format(s$pconsistency)),
              format(s$draws, big.mark = ","),
              if (is.null(x$seed)) "" else sprintf(";  seed = %s", format(x$seed))))
  cat(sprintf("  detection rate (DETECTED share)            : %.3f  (MC SE %.4f)\n",
              x$det_rate, x$det_rate_se))
  cat(sprintf("  E[|Hhat|] given detection                  : %.1f\n", x$EnH))
  cat(sprintf("  E[sens] / E[spec] / E[PPV] / E[NPV]        : %.3f / %.3f / %.3f / %.3f\n",
              x$Esens, x$Espec, x$Eppv, x$Enpv))
  cat(sprintf("  E[beta(Hhat)] oriented                     : %.2f\n", x$EbetaH))
  cat(sprintf("  E[naive - beta(Hhat)]                      : %.2f\n", x$Enaive_bias))
  cat(sprintf("  selection mass on rules with mean < c1     : %.3f\n", x$mass_below))
  top <- order(-x$p_sel)[seq_len(min(5L, x$M))]
  cat("  top selections (P_selected, sel | det):\n")
  for (m in top) {
    cat(sprintf("    %-40s %.3f  %.3f\n", substr(x$lab[m], 1, 40),
                x$p_sel[m], x$sel_c[m]))
  }
  invisible(x)
}


# -----------------------------------------------------------------------------
# internal helpers
# -----------------------------------------------------------------------------

# Monte-Carlo standard error of a proportion from R independent draws
# (the analytic document's .mc_se_prop).
.fs_oc_mc_se_prop <- function(p, R) sqrt(p * (1 - p) / R)

.fs_oc_check_family <- function(f) {
  need <- c("lab", "Pg", "PQg", "beta_g", "se_g", "sens_g", "spec_g",
            "ovl", "M", "PQ")
  miss <- setdiff(need, names(f))
  if (length(miss)) {
    stop("fs_oc_family is missing: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  M <- f$M
  vec <- c("lab", "Pg", "PQg", "beta_g", "se_g", "sens_g", "spec_g")
  bad <- vec[vapply(vec, function(nm) length(f[[nm]]) != M, logical(1))]
  if (length(bad)) {
    stop("fs_oc_family: length != M for ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  if (!is.matrix(f$ovl) || !all(dim(f$ovl) == c(M, M))) {
    stop("fs_oc_family: 'ovl' must be an M x M matrix.", call. = FALSE)
  }
  if (!is.numeric(f$PQ) || length(f$PQ) != 1L) {
    stop("fs_oc_family: 'PQ' must be a single number.", call. = FALSE)
  }
  invisible(TRUE)
}
