# =============================================================================
# Contract test: MR's inclusion band is the shared band
#
# .fs_mr_select() carried a local reimplementation of the same three band
# rules the consistency engine and DINA get from .compute_inclusion_band().
# Two copies of one rule is how the GRF frontier came to band differently from
# MR in the first place, so MR now calls the shared helper.
#
# The two copies were NOT identical, and the differences needed deciding
# rather than reconciling silently:
#
#   NA handling  -- the local copy used a plain max() with no is.na() guard.
#                   One non-finite effect made the THRESHOLD NA, hence the
#                   whole band NA, hence `if (!any(ib))` an error.  Fixed by
#                   adopting the shared helper's na.rm/!is.na treatment.
#   Empty band   -- MR's "never empty" fallback is KEPT, but relocated to the
#                   call site.  In the consistency/DINA sort an all-zero band
#                   is harmless (the next sort key breaks the tie); in MR the
#                   band is a filter inside a draw, so returning nothing drops
#                   that draw from sel_bias -- and drops exactly the draws
#                   where the perturbed effects were extreme.
# =============================================================================

.old_fs_mr_select <- function(beta, zcons, sizes, passers, rule, nbhd,
                              selection_rule = "neighborhood", log_scale = TRUE) {
  if (!length(passers)) return(NA_integer_)
  .inband <- function() {
    eff <- if (log_scale) exp(beta[passers]) else beta[passers]
    sz  <- sizes[passers]
    nb  <- eff >= (1 - nbhd) * max(eff)
    ib <- switch(selection_rule,
      neighborhood = nb,
      pareto       = !forestsearch:::.pareto_dominated_xy(eff, sz),
      both         = nb & !forestsearch:::.pareto_dominated_xy(eff, sz),
      nb)
    if (!any(ib)) ib <- rep(TRUE, length(passers))
    passers[ib]
  }
  pick <- switch(rule,
    maxcons  = passers[which.max(zcons[passers])],
    maxeff   = passers[which.max(beta[passers])],
    maxSG    = passers[which.max(sizes[passers])],
    minSG    = passers[which.min(sizes[passers])],
    effMaxSG = { b <- .inband(); b[which.max(sizes[b])] },
    effMinSG = { b <- .inband(); b[which.min(sizes[b])] },
    stop("unknown reselection rule: ", rule))
  as.integer(pick)
}

test_that("the relocation is behaviour-preserving on finite inputs", {
  # Asserted, not claimed: the fallback moved from inside the band helper to
  # the call site, so every finite case must select the identical candidate.
  set.seed(21)
  nc <- 0L; mm <- 0L
  for (rep in 1:250) {
    S <- sample(2:40, 1)
    beta  <- rnorm(S, sd = 1.5)
    sizes <- sample(30:400, S, TRUE)
    zc    <- rnorm(S)
    passers <- sort(sample(seq_len(S), sample(1:S, 1)))
    for (sr in c("neighborhood", "pareto", "both")) {
      for (r in c("effMaxSG", "effMinSG")) {
        for (nb in c(0.05, 0.10, 0.30)) {
          for (ls in c(TRUE, FALSE)) {
            nc <- nc + 1L
            a <- .old_fs_mr_select(beta, zc, sizes, passers, r, nb, sr, ls)
            b <- forestsearch:::.fs_mr_select(beta, zc, sizes, passers, r, nb, sr, ls)
            if (!identical(a, b)) mm <- mm + 1L
          }
        }
      }
    }
  }
  expect_gt(nc, 8000L)
  expect_identical(mm, 0L)
})

test_that("the empty-band fallback is retained and still fires", {
  # The band can legitimately empty: with a NEGATIVE maximum effect,
  # (1 - nbhd) * max exceeds max, so even the maximum fails its own test.
  eff <- c(-1, -2)
  ib <- forestsearch:::.compute_inclusion_band(eff, c(100L, 200L),
                                               "neighborhood", 0.10)
  expect_true(all(ib == 0L))            # genuinely empty
  # ... and MR still returns a selection rather than dropping the draw.
  got <- forestsearch:::.fs_mr_select(
    beta = eff, zcons = c(0, 0), sizes = c(100L, 200L), passers = 1:2,
    rule = "effMaxSG", nbhd = 0.10, selection_rule = "neighborhood",
    log_scale = FALSE)
  expect_false(is.na(got))
  expect_identical(
    got,
    .old_fs_mr_select(eff, c(0, 0), c(100L, 200L), 1:2, "effMaxSG", 0.10,
                      "neighborhood", FALSE))
})

test_that("a non-finite effect no longer poisons the whole band", {
  eff <- c(2.0, 1.9, NA, 0.5); sz <- c(100L, 300L, 50L, 400L)

  # OLD: plain max() makes the THRESHOLD NA, so the entire band is NA -- not
  # merely the offending row -- and `if (!any(ib))` then errors.
  old_thresh <- (1 - 0.10) * max(eff)
  expect_true(is.na(old_thresh))
  expect_true(all(is.na(eff >= old_thresh)))
  expect_error(if (!any(eff >= old_thresh)) 1 else 2,
               "missing value where TRUE/FALSE needed")

  # NEW: only the non-finite candidate is excluded.
  ib <- forestsearch:::.compute_inclusion_band(eff, sz, "neighborhood", 0.10)
  expect_identical(ib, c(1L, 1L, 0L, 0L))
})

test_that("MR's band agrees with the consistency/DINA band", {
  # One rule, one implementation.  This is the property that stops the GRF
  # divergence recurring inside MR.
  set.seed(5)
  for (rep in 1:100) {
    S <- sample(3:30, 1)
    eff <- exp(rnorm(S, sd = 0.7)); sz <- sample(30:400, S, TRUE)
    for (sr in c("neighborhood", "pareto", "both")) {
      shared <- forestsearch:::.compute_inclusion_band(eff, sz, sr, 0.10)
      nb <- eff >= (1 - 0.10) * max(eff)
      manual <- switch(sr,
        neighborhood = nb,
        pareto       = !forestsearch:::.pareto_dominated_xy(eff, sz),
        both         = nb & !forestsearch:::.pareto_dominated_xy(eff, sz))
      expect_identical(shared == 1L, manual, info = paste(sr, rep))
    }
  }
})
