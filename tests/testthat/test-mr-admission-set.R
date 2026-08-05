# =============================================================================
# Contract test: the MR admission set equals the identifier's, by construction
#
# MR linearizes the identifier's selection MAP -- a ranking AND a domain.  It
# is valid only if the set it re-selects over is the set the identifier
# selected over.  Previously MR rebuilt that domain at its own call site from
# raw thresholds (t_g <- pmax(c_screen, c_consistency + z * sigma_D)), so the
# identifier decided its floors in one place and MR re-derived them in
# another.  Under sg_focus = "maxeff" the two disagreed: the identifier applied
# no floors, MR still applied the effect floor.  In the other direction, DINA
# and GRF were passed c_consistency = 0 with p_star still set, so t_g carried a
# z * sigma_D consistency term those engines never had.
#
# These tests pin the resolution table and the representation of absence.
# =============================================================================

test_that("every (focus, engine) pair resolves to the intended admission set", {
  expect_true(forestsearch:::.assert_admission_resolution_complete())
})

test_that("consistency + maxeff admits everything -- both floors NULL", {
  a <- forestsearch:::.fs_resolve_admission(
    "maxeff", "consistency",
    hr.threshold = log(1.25), hr.consistency = 0, pconsistency.threshold = 0.9)
  expect_null(a$effect_floor)
  expect_null(a$consistency)
})

test_that("consistency + every other focus keeps both floors", {
  for (f in setdiff(forestsearch:::.FS_SG_FOCUS_CANONICAL, "maxeff")) {
    a <- forestsearch:::.fs_resolve_admission(
      f, "consistency",
      hr.threshold = log(1.25), hr.consistency = 0, pconsistency.threshold = 0.9)
    expect_equal(a$effect_floor, log(1.25), info = f)
    expect_equal(a$consistency$c_cons, 0, info = f)
    expect_equal(a$consistency$p_star, 0.9, info = f)
  }
})

test_that("dina and grf get an effect floor and NO consistency term", {
  # The engines compute no Pcons and no sort key of theirs contains one, so a
  # consistency term in their admission set is a floor they never applied.
  for (m in c("dina", "grf")) {
    for (f in forestsearch:::.FS_SG_FOCUS_CANONICAL) {
      a <- forestsearch:::.fs_resolve_admission(f, m, hr.threshold = log(1.25))
      expect_equal(a$effect_floor, log(1.25), info = paste(m, f))
      expect_null(a$consistency)
    }
  }
})

test_that("absence is NULL, never a permissive sentinel", {
  # -Inf would make "no filter" arithmetically indistinguishable from "a filter
  # that admits everything", and only the second is testable.
  a <- forestsearch:::.fs_resolve_admission("maxeff", "consistency")
  expect_null(a$effect_floor)
  expect_false(identical(a$effect_floor, -Inf))
  expect_false(isTRUE(is.numeric(a$effect_floor)))
})

test_that("a floor that applies must have a finite value supplied", {
  # Silently defaulting a missing threshold would reintroduce reconstruction.
  expect_error(
    forestsearch:::.fs_resolve_admission("hr", "consistency"),
    "effect floor applies")
  expect_error(
    forestsearch:::.fs_resolve_admission("hr", "consistency",
                                         hr.threshold = log(1.25)),
    "consistency floor applies")
})

test_that("the identifier's presence pattern matches the resolver's", {
  # .fs_admission_applies() is what forestsearch()'s own disable_effect_floor
  # and rmin sites read.  If it drifts from .fs_resolve_admission(), the
  # identifier and MR disagree again -- the exact failure being fixed.
  for (m in c("consistency", "dina", "grf")) {
    for (f in forestsearch:::.FS_SG_FOCUS_CANONICAL) {
      ap <- forestsearch:::.fs_admission_applies(f, m)
      a  <- forestsearch:::.fs_resolve_admission(
        f, m, hr.threshold = log(1.25), hr.consistency = 0,
        pconsistency.threshold = 0.9)
      expect_equal(unname(ap[["effect"]]), !is.null(a$effect_floor),
                   info = paste(m, f))
      expect_equal(unname(ap[["consistency"]]), !is.null(a$consistency),
                   info = paste(m, f))
    }
  }
})

test_that("aliases resolve the same as their canonical spellings", {
  for (pair in list(c("eff", "hr"), c("maxcons", "hr"),
                    c("effMaxSG", "hrMaxSG"), c("effMinSG", "hrMinSG"))) {
    a <- forestsearch:::.fs_resolve_admission(
      pair[1], "consistency", hr.threshold = log(1.25),
      hr.consistency = 0, pconsistency.threshold = 0.9)
    b <- forestsearch:::.fs_resolve_admission(
      pair[2], "consistency", hr.threshold = log(1.25),
      hr.consistency = 0, pconsistency.threshold = 0.9)
    expect_equal(a, b, info = paste(pair, collapse = "/"))
  }
})

test_that("the resolution assertion fails when a pair is wrong", {
  # Negative control: without it, the assertion could pass vacuously.
  ns <- asNamespace("forestsearch")
  orig <- get(".fs_admission_applies", envir = ns)
  unlockBinding(".fs_admission_applies", ns)
  on.exit({
    assign(".fs_admission_applies", orig, envir = ns)
    lockBinding(".fs_admission_applies", ns)
  }, add = TRUE)

  # Give DINA a consistency floor it does not have.
  assign(".fs_admission_applies",
         function(sg_focus, subgroup_method = "consistency")
           c(effect = TRUE, consistency = TRUE),
         envir = ns)
  expect_error(forestsearch:::.assert_admission_resolution_complete(),
               "Admission resolution is wrong")
})

test_that("MR refuses to rank on maxcons with no consistency floor", {
  # Ranking on a consistency-standardized statistic the identifier never
  # computed is the same mismatch, in the ranking rather than the domain.
  expect_error(
    forestsearch:::fs_mr_inference(
      df = data.frame(y = rbinom(60, 1, .5), treat = rep(0:1, 30)),
      candidates = list(a = 1:40, b = 21:60),
      spec = list(outcome_type = "binary", effect_measure = "OR",
                  treat.name = "treat", outcome.name = "y",
                  event.name = "y", offset.name = NULL,
                  adjust_covariates = NULL, adverse_outcome = TRUE),
      selected_members = 1:40,
      admission = list(effect_floor = 0, consistency = NULL),
      reselection = "maxcons", draws = 10L),
    "maxcons")
})

test_that("the resolved admission is recorded on the object contract", {
  # $admission sits next to $family_status so a run's domain is auditable
  # rather than inferable from sg_focus and subgroup_method.
  a <- forestsearch:::.fs_resolve_admission(
    "hr", "consistency", hr.threshold = log(1.25),
    hr.consistency = 0, pconsistency.threshold = 0.9)
  s <- forestsearch:::.fs_format_admission(a)
  expect_type(s, "character")
  expect_match(s, "effect floor")
  expect_match(s, "consistency floor")
  expect_match(
    forestsearch:::.fs_format_admission(
      forestsearch:::.fs_resolve_admission("maxeff", "consistency")),
    "unrestricted")
  expect_null(forestsearch:::.fs_format_admission(NULL))
})
