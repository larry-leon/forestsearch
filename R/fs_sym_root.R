# =============================================================================
# fs_sym_root() -- symmetric square root of a (possibly rank-deficient)
#                  covariance matrix, continuous in its input
# =============================================================================
#
# Promoted from the .sym_root() helper of the anchor chunk of
# quarto/simulations/actg175/continuous/analytic_verification_and_prediction_
# md_harm.qmd, where a document chunk could not carry the test asserting the
# property that motivated it.  The body is copied verbatim from that chunk.
#
# The failure this prevents: with the asymmetric root V D^{1/2}, an
# arbitrarily small change to the scale anchor moved every Monte-Carlo figure
# in 94fd4dad, because the eigenvector basis is not a continuous function of
# the input matrix.  The symmetric root V D^{1/2} V' is basis-free, so the
# sampled draws become a continuous function of the inputs.
# =============================================================================


#' Symmetric Square Root of a Covariance Matrix
#'
#' Computes the symmetric matrix square root of \code{scale * S} after
#' symmetrising \code{S} as \code{(S + t(S)) / 2}: with eigendecomposition
#' \eqn{S = V D V'}, the returned matrix is \eqn{R = V D^{1/2} V'} (eigenvalues
#' scaled by \code{scale} and clamped at zero before the square root), so that
#' \eqn{R R' = \mathrm{scale} \cdot S}.
#'
#' @section Why the symmetric root and not V D^(1/2):
#' Both roots satisfy \eqn{R R' = \mathrm{scale} \cdot S}, but the asymmetric
#' root \eqn{V D^{1/2}} depends on the eigenvector basis, which is not a
#' continuous function of its input: signs flip and vectors rotate freely
#' inside degenerate subspaces -- which a candidate covariance built from
#' complement pairs has exactly, since the complement constraints make it
#' rank-deficient by construction.  A pinned RNG seed fixes the standard-normal
#' draws Z; it cannot fix the basis Z is multiplied by, so the asymmetric root
#' is non-reproducible across arithmetically equivalent routes to the same
#' matrix.  The symmetric root \eqn{V D^{1/2} V'} is a matrix function of the
#' symmetrised input alone -- basis-free and continuous -- so samples built
#' from it are continuous, reproducible functions of the inputs.
#'
#' @param S Numeric square matrix, a covariance (or any matrix whose
#'   symmetrisation \code{(S + t(S)) / 2} is intended; the function symmetrises
#'   before decomposing).  May be rank-deficient; eigenvalues that are slightly
#'   negative from numerical noise are clamped to zero.
#' @param scale Numeric scalar multiplying \code{S} before the root is taken,
#'   so that the returned root \eqn{R} satisfies \eqn{R R' = \mathrm{scale}
#'   \cdot S}.  Default \code{2}, the half-sample covariance convention of the
#'   analytic documents.
#'
#' @return A symmetric numeric matrix \code{R} of the same dimension as
#'   \code{S}, satisfying \code{R %*% t(R)} equal to \code{scale} times the
#'   symmetrised \code{S} (up to the zero-clamp on negative eigenvalues).
#'
#' @examples
#' S <- crossprod(matrix(rnorm(25), 5))   # symmetric positive definite
#' R <- fs_sym_root(S, scale = 2)
#' max(abs(R %*% t(R) - 2 * S))           # ~ 1e-14
#' isSymmetric(R)
#'
#' # Rank-deficient input: the root exists and stays symmetric
#' X <- matrix(rnorm(15), 5, 3)
#' R2 <- fs_sym_root(tcrossprod(X), scale = 1)
#' max(abs(R2 %*% t(R2) - tcrossprod(X)))
#'
#' @export
fs_sym_root <- function(S, scale = 2) {
  eS <- eigen((S + t(S)) / 2, symmetric = TRUE)
  V  <- eS$vectors
  d  <- sqrt(pmax(scale * eS$values, 0))
  V %*% diag(d, nrow = length(d)) %*% t(V)
}
