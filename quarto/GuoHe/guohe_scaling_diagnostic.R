# guohe_scaling_diagnostic.R
#
# How does the Guo & He (2021) de-biasing comparator behave when the candidate
# family grows from a handful to the hundreds, as in a `forestsearch`
# enumeration? Three diagnostics:
#
#   (1) family-size sweep on a NULL design (no effect modification anywhere),
#       so any apparent "best subgroup" is entirely winner's curse;
#   (2) correlation adaptation: independent vs near-duplicate candidate
#       families at matched prevalence;
#   (3) resample estimability: how often candidates drop out of the bootstrap
#       competition, which silently under-corrects.
#
# Run from the directory holding `guohe_comparator.R`. Single realizations
# throughout: read these for DIRECTION, not as variance estimates. For
# publication numbers, wrap in replicate() on the 127-core box.

suppressMessages(library(survival))
source("guohe_comparator.R")

N        <- 1000     # trial size
TRUE_MD  <- 0.1      # overall treatment effect; NO subgroup modifies it
B_DEMO   <- 200L     # bootstrap draws (raise for real runs)
K_GRID   <- c(8, 25, 50, 100, 200)

# ---- (1) family-size sweep, null design -----------------------------------
cat("=== (1) Family-size sweep (null design; true MD =", TRUE_MD, ") ===\n")
set.seed(2026)
trt  <- rbinom(N, 1, 0.5)
y    <- 1 + TRUE_MD * trt + rnorm(N)
xall <- as.data.frame(matrix(rbinom(N * max(K_GRID), 1, 0.4), N, max(K_GRID)))
names(xall) <- paste0("v", seq_len(max(K_GRID)))
dat <- cbind(data.frame(trt = trt, y = y), xall)

sweep <- do.call(rbind, lapply(K_GRID, function(k) {
  t0 <- Sys.time()
  r  <- guohe_debias(dat, "continuous", treatment = "trt",
                     candidates = paste0("v", seq_len(k)), y = "y",
                     orient = +1, B = B_DEMO, seed = 99L, symmetric = FALSE)
  data.frame(K = k, selected = r$selected,
             naive = round(r$naive, 3), debiased = round(r$debiased, 3),
             bound = round(r$bound, 3),
             covers_truth = r$bound <= TRUE_MD,
             secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
}))
print(sweep, row.names = FALSE)

# ---- (2) correlation adaptation, prevalence matched -----------------------
# forestsearch candidates are conjunctions on shared covariates: heavily
# overlapping near-duplicates. Prevalence is held equal across designs so that
# CORRELATION is isolated from subgroup SIZE (which independently drives the
# correction through candidate variance).
cat("\n=== (2) Correlation adaptation (K = 50, prevalence matched) ===\n")
set.seed(808)
K <- 50L; P <- 0.4
trt <- rbinom(N, 1, 0.5); y <- 1 + TRUE_MD * trt + rnorm(N)

x_ind <- as.data.frame(matrix(rbinom(N * K, 1, P), N, K))
names(x_ind) <- paste0("v", seq_len(K))

x_dup <- list(); nm <- character(0)
for (j in 1:5) {                          # 5 blocks x 10 near-duplicate clones
  base <- rbinom(N, 1, P)
  for (k in 1:10) {
    v <- base
    fl <- sample.int(N, round(0.04 * N))  # ~4% flips -> high within-block corr
    v[fl] <- 1L - v[fl]
    x_dup[[length(x_dup) + 1L]] <- v
    nm <- c(nm, sprintf("b%d_c%02d", j, k))
  }
}
x_dup <- as.data.frame(x_dup); names(x_dup) <- nm

run_family <- function(x, cand) {
  guohe_debias(cbind(data.frame(trt = trt, y = y), x), "continuous",
               treatment = "trt", candidates = cand, y = "y",
               orient = +1, B = B_DEMO, seed = 5L, symmetric = FALSE)
}
r_ind <- run_family(x_ind, names(x_ind))
r_dup <- run_family(x_dup, nm)
corr_of <- function(x) mean(abs(cor(x)[upper.tri(diag(ncol(x)))]))

cat(sprintf("  independent    mean|corr| %.3f | naive %.3f -> debiased %.3f | correction %.3f\n",
            corr_of(x_ind), r_ind$naive, r_ind$debiased, r_ind$naive - r_ind$debiased))
cat(sprintf("  near-duplicate mean|corr| %.3f | naive %.3f -> debiased %.3f | correction %.3f\n",
            corr_of(x_dup), r_dup$naive, r_dup$debiased, r_dup$naive - r_dup$debiased))
cat(sprintf("  correction ratio (dup / ind) = %.2f\n",
            (r_dup$naive - r_dup$debiased) / (r_ind$naive - r_ind$debiased)))

# ---- (3) resample estimability -------------------------------------------
# Candidates that fail to fit in a given resample are dropped from that
# resample's maximum. If drop rates are non-trivial, the bootstrap maximum is
# taken over a shrinking competitor set and the correction is UNDER-stated.
# Small conjunctions in a real forestsearch family are the risk case.
cat("\n=== (3) Resample estimability (small-subgroup stress) ===\n")
set.seed(31415)
trt <- rbinom(N, 1, 0.5); y <- 1 + TRUE_MD * trt + rnorm(N)
prev <- c(0.02, 0.05, 0.10, 0.25, 0.40)          # incl. very small subgroups
x_sm <- as.data.frame(lapply(prev, function(p) rbinom(N, 1, p)))
names(x_sm) <- sprintf("p%03d", round(1000 * prev))

drop_rate <- sapply(names(x_sm), function(cn) {
  memb <- x_sm[[cn]]
  mean(replicate(200, {
    idx <- sample.int(N, N, replace = TRUE)
    s   <- memb[idx] == 1L
    # non-estimable if either arm is absent within the resampled subgroup
    length(unique(trt[idx][s])) < 2L || sum(s) < 10L
  }))
})
print(data.frame(candidate = names(x_sm),
                 observed_size = as.integer(colSums(x_sm)),
                 resample_drop_rate = round(drop_rate, 3)),
      row.names = FALSE)
cat("\nInterpretation: non-zero drop rates mean the bootstrap maximum is taken\n",
    "over a varying competitor set. Raise `min_events`, or prune tiny\n",
    "candidates before supplying the family, so the competition is stable.\n", sep = "")
