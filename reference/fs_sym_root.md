# Symmetric Square Root of a Covariance Matrix

Computes the symmetric matrix square root of `scale * S` after
symmetrising `S` as `(S + t(S)) / 2`: with eigendecomposition \\S = V D
V'\\, the returned matrix is \\R = V D^{1/2} V'\\ (eigenvalues scaled by
`scale` and clamped at zero before the square root), so that \\R R' =
\mathrm{scale} \cdot S\\.

## Usage

``` r
fs_sym_root(S, scale = 2)
```

## Arguments

- S:

  Numeric square matrix, a covariance (or any matrix whose
  symmetrisation `(S + t(S)) / 2` is intended; the function symmetrises
  before decomposing). May be rank-deficient; eigenvalues that are
  slightly negative from numerical noise are clamped to zero.

- scale:

  Numeric scalar multiplying `S` before the root is taken, so that the
  returned root \\R\\ satisfies \\R R' = \mathrm{scale} \cdot S\\.
  Default `2`, the half-sample covariance convention of the analytic
  documents.

## Value

A symmetric numeric matrix `R` of the same dimension as `S`, satisfying
`R %*% t(R)` equal to `scale` times the symmetrised `S` (up to the
zero-clamp on negative eigenvalues).

## Why the symmetric root and not V D^(1/2)

Both roots satisfy \\R R' = \mathrm{scale} \cdot S\\, but the asymmetric
root \\V D^{1/2}\\ depends on the eigenvector basis, which is not a
continuous function of its input: signs flip and vectors rotate freely
inside degenerate subspaces – which a candidate covariance built from
complement pairs has exactly, since the complement constraints make it
rank-deficient by construction. A pinned RNG seed fixes the
standard-normal draws Z; it cannot fix the basis Z is multiplied by, so
the asymmetric root is non-reproducible across arithmetically equivalent
routes to the same matrix. The symmetric root \\V D^{1/2} V'\\ is a
matrix function of the symmetrised input alone – basis-free and
continuous – so samples built from it are continuous, reproducible
functions of the inputs.

## Examples

``` r
S <- crossprod(matrix(rnorm(25), 5))   # symmetric positive definite
R <- fs_sym_root(S, scale = 2)
max(abs(R %*% t(R) - 2 * S))           # ~ 1e-14
#> [1] 2.4869e-14
isSymmetric(R)
#> [1] TRUE

# Rank-deficient input: the root exists and stays symmetric
X <- matrix(rnorm(15), 5, 3)
R2 <- fs_sym_root(tcrossprod(X), scale = 1)
max(abs(R2 %*% t(R2) - tcrossprod(X)))
#> [1] 9.769963e-15
```
