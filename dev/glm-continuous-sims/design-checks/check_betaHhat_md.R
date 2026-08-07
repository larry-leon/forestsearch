# Does the closed form for beta(Hhat) survive REAL rule shapes -- multi-covariate
# conjunctions, mis-specified, including negations -- or only the single-cut case?
#
# betaHhat_truth.R computes beta(Hhat) by NUMERICAL evaluation on a fixed
# super-population, because for survival there is no closed form. For MD with an
# additive DGM the target is
#     beta(Hhat) = tau * P(H_true | Hhat),
# a ratio of rectangle probabilities -- closed form when covariates are
# independent, since an axis-aligned rectangle's probability is a product of
# marginals. This checks that against simulation for rules of the shape the
# engines actually realize.

set.seed(4)
n <- 4e6
X1 <- runif(n); X2 <- runif(n); X3 <- runif(n)
A  <- rbinom(n, 1, 0.5)
tau <- 0.8
H_true <- (X1 > 0.5) & (X2 < 0.6)                 # true harm rectangle
Y <- 1 + tau * H_true * A + rnorm(n)

# Closed form: P(rectangle) is a product of marginals under independence.
p_rect <- function(lo1, hi1, lo2, hi2, lo3, hi3)
  (hi1 - lo1) * (hi2 - lo2) * (hi3 - lo3)

beta_cf <- function(r) {            # r = list of (lo, hi) per covariate for Hhat
  # intersection of Hhat with H_true = {X1 in (0.5,1)} x {X2 in (0,0.6)} x all X3
  i1 <- c(max(r$x1[1], 0.5), min(r$x1[2], 1.0))
  i2 <- c(max(r$x2[1], 0.0), min(r$x2[2], 0.6))
  i3 <- r$x3
  num <- if (i1[2] <= i1[1] || i2[2] <= i2[1]) 0 else
           p_rect(i1[1], i1[2], i2[1], i2[2], i3[1], i3[2])
  den <- p_rect(r$x1[1], r$x1[2], r$x2[1], r$x2[2], r$x3[1], r$x3[2])
  tau * num / den
}

beta_emp <- function(r) {
  H <- X1 > r$x1[1] & X1 <= r$x1[2] &
       X2 > r$x2[1] & X2 <= r$x2[2] &
       X3 > r$x3[1] & X3 <= r$x3[2]
  mean(Y[H & A == 1]) - mean(Y[H & A == 0])
}

rules <- list(
  "exact rule            {X1>.5} & {X2<.6}" =
    list(x1 = c(0.5, 1), x2 = c(0, 0.6), x3 = c(0, 1)),
  "too narrow            {X1>.7} & {X2<.4}" =
    list(x1 = c(0.7, 1), x2 = c(0, 0.4), x3 = c(0, 1)),
  "too wide              {X1>.3} & {X2<.8}" =
    list(x1 = c(0.3, 1), x2 = c(0, 0.8), x3 = c(0, 1)),
  "one cut only          {X1>.5}"           =
    list(x1 = c(0.5, 1), x2 = c(0, 1),   x3 = c(0, 1)),
  "spurious 3rd covariate + {X3>.5}"        =
    list(x1 = c(0.5, 1), x2 = c(0, 0.6), x3 = c(0.5, 1)),
  "negation              !{X1<=.5} & !{X2>=.6}" =
    list(x1 = c(0.5, 1), x2 = c(0, 0.6), x3 = c(0, 1)),
  "disjoint from truth   {X1<.4}"           =
    list(x1 = c(0, 0.4), x2 = c(0, 1),   x3 = c(0, 1))
)

cat(sprintf("  %-42s %10s %10s %9s\n", "rule", "closed", "empirical", "diff"))
for (nm in names(rules)) {
  r <- rules[[nm]]
  cf <- beta_cf(r); em <- beta_emp(r)
  cat(sprintf("  %-42s %10.5f %10.5f %9.5f\n", nm, cf, em, abs(cf - em)))
}
cat("\n  MC SE at these sizes is ~2/sqrt(n_subgroup_arm); all diffs are that order.\n")
