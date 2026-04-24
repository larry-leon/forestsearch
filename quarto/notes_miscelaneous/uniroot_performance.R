# uniroot_performance.R
# ---------------------------------------------------------------
# Techniques for speeding up R's uniroot() plus alternatives.
# Run section by section; each benchmark is self-contained.
#
# Requires: microbenchmark (required), rootSolve/pracma (optional)
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(microbenchmark)
})

# ---------------------------------------------------------------
# A "moderately expensive" objective function we can reuse.
# Think of it as a stand-in for a likelihood or estimating equation.
# ---------------------------------------------------------------
make_objective <- function(target, work = 0L) {
  force(target); force(work)
  function(x) {
    # Optional fake CPU work to simulate a pricey evaluation
    if (work > 0L) for (i in seq_len(work)) sqrt(i)
    pnorm(x) - target
  }
}

f <- make_objective(0.975)

# ---------------------------------------------------------------
# 1. Bracket width matters
# ---------------------------------------------------------------
cat("\n[1] Bracket width ---------------------------------------\n")
print(microbenchmark(
  very_wide = uniroot(f, interval = c(-50, 50)),
  wide      = uniroot(f, interval = c(-10, 10)),
  tight     = uniroot(f, interval = c(1, 3)),
  times = 200L
))

# ---------------------------------------------------------------
# 2. Pre-computed endpoint values via f.lower / f.upper
#    Saves two function evaluations per call.
# ---------------------------------------------------------------
cat("\n[2] Passing f.lower / f.upper ---------------------------\n")
fl <- f(1); fu <- f(3)
print(microbenchmark(
  without_endpoints = uniroot(f, interval = c(1, 3)),
  with_endpoints    = uniroot(f, interval = c(1, 3),
                              f.lower = fl, f.upper = fu),
  times = 500L
))

# ---------------------------------------------------------------
# 3. Tolerance trade-off
#    Default tol is .Machine$double.eps^0.25 ~ 1.22e-4
# ---------------------------------------------------------------
cat("\n[3] Tolerance trade-off ---------------------------------\n")
print(microbenchmark(
  tol_1e_10 = uniroot(f, interval = c(1, 3), tol = 1e-10),
  tol_1e_8  = uniroot(f, interval = c(1, 3), tol = 1e-8),
  default   = uniroot(f, interval = c(1, 3)),
  tol_1e_3  = uniroot(f, interval = c(1, 3), tol = 1e-3),
  times = 500L
))

# ---------------------------------------------------------------
# 4. Expensive objective: compiled-function gain dominates
#    (illustrative — real wins come from Rcpp/C++ rewrites)
# ---------------------------------------------------------------
cat("\n[4] With artificial work in the objective ---------------\n")
f_slow <- make_objective(0.975, work = 1000L)
print(microbenchmark(
  wide_slow_f  = uniroot(f_slow, interval = c(-10, 10)),
  tight_slow_f = uniroot(f_slow, interval = c(1, 3)),
  times = 100L
))

# ---------------------------------------------------------------
# 5. Warm-start pattern for batched/bootstrap calls
#    Classic win when calling uniroot() many times with slowly
#    varying parameters (the common case inside a bootstrap).
# ---------------------------------------------------------------
cat("\n[5] Warm-start across many calls ------------------------\n")
targets <- seq(0.50, 0.999, length.out = 1000)

naive_batch <- function(targets) {
  vapply(targets, function(tg) {
    uniroot(make_objective(tg), interval = c(-10, 10))$root
  }, numeric(1))
}

warm_batch <- function(targets, half_width = 0.5) {
  out <- numeric(length(targets))
  prev <- 0
  for (i in seq_along(targets)) {
    # small window around previous solution; extendInt handles rare misses
    out[i] <- uniroot(
      make_objective(targets[i]),
      interval  = c(prev - half_width, prev + half_width),
      extendInt = "yes"
    )$root
    prev <- out[i]
  }
  out
}

# Sanity check (should be ~ 0)
stopifnot(max(abs(naive_batch(targets) - warm_batch(targets))) < 1e-3)

print(microbenchmark(
  naive = naive_batch(targets),
  warm  = warm_batch(targets),
  times = 5L
))

# ---------------------------------------------------------------
# 6. Endpoint caching inside a batched call
#    If many calls share bracket endpoints, cache endpoint values.
# ---------------------------------------------------------------
cat("\n[6] Shared endpoints across many calls ------------------\n")

# Suppose we're solving f_tg(x) = pnorm(x) - tg on a fixed bracket
# and the "expensive" part is independent of tg -- common in
# survival / estimating-equation problems after profiling out tg.

fixed_bracket_naive <- function(targets) {
  vapply(targets, function(tg) {
    g <- function(x) pnorm(x) - tg
    uniroot(g, interval = c(-5, 5))$root
  }, numeric(1))
}

fixed_bracket_cached <- function(targets) {
  lo <- -5; hi <- 5
  # Endpoint pieces that are target-independent
  plo <- pnorm(lo); phi <- pnorm(hi)
  vapply(targets, function(tg) {
    g  <- function(x) pnorm(x) - tg
    uniroot(g, interval = c(lo, hi),
            f.lower = plo - tg, f.upper = phi - tg)$root
  }, numeric(1))
}

print(microbenchmark(
  no_cache = fixed_bracket_naive(targets),
  cached   = fixed_bracket_cached(targets),
  times = 10L
))

# ---------------------------------------------------------------
# 7. Alternatives worth knowing
# ---------------------------------------------------------------
cat("\n[7] Alternative root-finders (if packages installed) ----\n")

# rootSolve::uniroot.all -- grids an interval, then polishes each sign change.
# Use when monotonicity is not guaranteed.
if (requireNamespace("rootSolve", quietly = TRUE)) {
  ff <- function(x) sin(x) - 0.3
  roots <- rootSolve::uniroot.all(ff, interval = c(0, 4 * pi))
  cat("rootSolve::uniroot.all roots on [0, 4pi] for sin(x)=0.3:\n")
  print(roots)
}

# pracma::brentDekker / pracma::fzero -- same family, sometimes more
# diagnostics available.
if (requireNamespace("pracma", quietly = TRUE)) {
  cat("\npracma::brentDekker vs base uniroot:\n")
  print(microbenchmark(
    base   = uniroot(f, interval = c(1, 3)),
    pracma = pracma::brentDekker(f, a = 1, b = 3),
    times = 200L
  ))
}

# ---------------------------------------------------------------
# 8. Newton's method for monotone C^1 objectives (hand-rolled)
#    When you have a cheap derivative, Newton typically beats Brent.
#    Here: quantile of Normal via f(x) = Phi(x) - p, f'(x) = phi(x).
# ---------------------------------------------------------------
cat("\n[8] Hand-rolled Newton vs uniroot (Normal quantile) -----\n")

newton_qnorm <- function(p, x0 = 0, tol = 1e-8, maxit = 50L) {
  x <- x0
  for (i in seq_len(maxit)) {
    fx  <- pnorm(x) - p
    if (abs(fx) < tol) return(x)
    dfx <- dnorm(x)
    if (dfx < 1e-300) break
    x <- x - fx / dfx
  }
  x
}

stopifnot(abs(newton_qnorm(0.975) - qnorm(0.975)) < 1e-6)

f975 <- make_objective(0.975)
print(microbenchmark(
  uniroot_call = uniroot(f975, interval = c(-5, 5))$root,
  newton_call  = newton_qnorm(0.975),
  qnorm_call   = qnorm(0.975),   # reference: the fully compiled option
  times = 500L
))

# ---------------------------------------------------------------
# Summary guidance
# ---------------------------------------------------------------
# Profile first:  profvis::profvis({ ... your code ... })
# Then, in order:
#   1. Narrow the bracket using domain knowledge
#   2. Pass f.lower / f.upper when available
#   3. Relax tol if 4-5 sig figs is enough
#   4. Avoid extendInt in hot loops
#   5. Warm-start across batched calls
#   6. Move f into Rcpp if it dominates the profile
#   7. Consider Newton/Halley for monotone C^1 objectives
# ---------------------------------------------------------------
