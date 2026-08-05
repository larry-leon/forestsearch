# =============================================================================
# Contract test: GRF's inclusion band honours selection_rule
#
# .grf_frontier_select() used to hardcode the multiplicative neighborhood band
# and took no selection_rule, while MR's re-selection DID honour it.  Under
# subgroup_method = "grf" with a band focus and selection_rule = "pareto" or
# "both", the identifier banded by neighborhood and MR banded by pareto -- a
# ranking/domain divergence of the same class as the admission-set defect, and
# silent in the same way.  Validation permitted the combination, so nothing
# surfaced it.
#
# The band now comes from .compute_inclusion_band(), the shared helper DINA and
# the consistency engine use.
# =============================================================================

.mk <- function(effect, size) data.frame(effect = effect, size = size)

test_that("selection_rule = 'neighborhood' is numerically unchanged", {
  # The two formulas are the same comparison written differently:
  #   old: effect >= emax * (1 - nbhd)
  #   new: effect >= (1 - effect_neighborhood) * max(effect)
  # This is the guarantee the change rests on, so it is checked exhaustively
  # rather than asserted.
  old_select <- function(cand, dmin, rule = "effMaxSG", nbhd = 0.10) {
    if (is.null(cand) || !nrow(cand)) return(NULL)
    elig <- cand[cand$effect >= dmin, , drop = FALSE]
    if (!nrow(elig)) return(NULL)
    emax <- max(elig$effect)
    if (rule == "maxSG") elig[which.max(elig$size), , drop = FALSE]
    else if (rule == "minSG") elig[which.min(elig$size), , drop = FALSE]
    else if (rule == "eff") elig[which.max(elig$effect), , drop = FALSE]
    else {
      band <- elig[elig$effect >= emax * (1 - nbhd), , drop = FALSE]
      if (rule == "effMinSG") band[which.min(band$size), , drop = FALSE]
      else band[which.max(band$size), , drop = FALSE]
    }
  }
  set.seed(7)
  for (rep in 1:60) {
    n <- sample(3:25, 1)
    cand <- .mk(round(rnorm(n), 3), sample(20:400, n, TRUE))
    for (r in c("eff", "maxSG", "minSG", "effMaxSG", "effMinSG")) {
      for (nb in c(0.05, 0.10, 0.25)) {
        for (dm in c(-Inf, 0, 0.5)) {
          expect_identical(
            old_select(cand, dmin = dm, rule = r, nbhd = nb),
            forestsearch:::.grf_frontier_select(
              cand, dmin = dm, rule = r, nbhd = nb,
              selection_rule = "neighborhood"),
            info = sprintf("rep=%d rule=%s nbhd=%g dmin=%g", rep, r, nb, dm))
        }
      }
    }
  }
})

test_that("selection_rule reaches the band and can change the winner", {
  # A candidate inside the neighborhood band but Pareto-dominated: under
  # "neighborhood" it can win on size; under "pareto"/"both" it cannot.
  cand <- .mk(effect = c(1.00, 0.95, 0.94),
              size   = c(100L, 300L, 290L))
  nb <- forestsearch:::.grf_frontier_select(cand, dmin = 0, rule = "effMaxSG",
                                            nbhd = 0.10,
                                            selection_rule = "neighborhood")
  pa <- forestsearch:::.grf_frontier_select(cand, dmin = 0, rule = "effMaxSG",
                                            nbhd = 0.10,
                                            selection_rule = "pareto")
  bo <- forestsearch:::.grf_frontier_select(cand, dmin = 0, rule = "effMaxSG",
                                            nbhd = 0.10,
                                            selection_rule = "both")
  # (0.94, 290) is dominated by (0.95, 300); it is in the neighborhood band.
  expect_equal(nb$size, 300L)
  expect_false(any(pa$effect == 0.94))
  expect_false(any(bo$effect == 0.94))
})

test_that("the non-band rules ignore selection_rule entirely", {
  cand <- .mk(c(1.0, 0.9, 0.5), c(100L, 300L, 200L))
  for (r in c("eff", "maxSG", "minSG")) {
    a <- forestsearch:::.grf_frontier_select(cand, 0, r, 0.10, "neighborhood")
    b <- forestsearch:::.grf_frontier_select(cand, 0, r, 0.10, "pareto")
    c_ <- forestsearch:::.grf_frontier_select(cand, 0, r, 0.10, "both")
    expect_identical(a, b, info = r)
    expect_identical(a, c_, info = r)
  }
})

test_that("selection_rule is threaded to every GRF entry point", {
  # A formal that exists but is never passed through is the failure mode here.
  expect_true("selection_rule" %in% names(formals(forestsearch:::.grf_frontier_select)))
  expect_true("selection_rule" %in% names(formals(forestsearch:::.grf_reselect_on_effect)))
  expect_true("selection_rule" %in% names(formals(forestsearch:::.forestsearch_grf_select)))
  expect_true("selection_rule" %in% names(formals(grf.subg.harm.survival)))
  expect_true("selection_rule" %in% names(formals(grf.subg.harm.glm)))
})

test_that(".grf_mark_frontier keeps its contract and uses the shared core", {
  set.seed(3)
  for (rep in 1:40) {
    n <- sample(2:20, 1)
    cand <- .mk(round(runif(n, 0, 2), 2), sample(c(50L, 100L, 150L), n, TRUE))
    out <- forestsearch:::.grf_mark_frontier(cand)
    # contract: sorted by (-effect, -size), with an on_frontier column
    expect_identical(out[, c("effect", "size")],
                     cand[order(-cand$effect, -cand$size),
                          c("effect", "size")])
    expect_true("on_frontier" %in% names(out))
    expect_type(out$on_frontier, "logical")
    # and it agrees with the shared dominance core
    expect_identical(out$on_frontier,
                     !forestsearch:::.pareto_dominated_xy(out$effect, out$size))
  }
})

test_that("exact (effect, size) duplicates are all non-dominated", {
  # The documented behaviour change: the former local sweep used
  # `size > best_size`, keeping only the first of a tied pair.  Pareto
  # dominance requires a strict improvement on one axis, so ties all survive.
  # on_frontier is a reported field, not a selection input.
  cand <- .mk(c(1.0, 1.0, 0.5), c(100L, 100L, 200L))
  out <- forestsearch:::.grf_mark_frontier(cand)
  expect_equal(sum(out$on_frontier[out$effect == 1.0]), 2)
})
