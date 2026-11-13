# fs_attach_betaHhat(method = "exact") (TASK_fs-exact-betaHhat-mr-record_2026-11-01_r2,
# sections 1 and 3.1-3.2): the exact survival beta(Hhat) from the analytic
# layer against the frame route, the censoring guard, empty and unresolvable
# rules, and the identity with ahr_limit() on the replicate-5 membership.
# The DGM is the S2 build of helper-mrct-analytic.R.

skip_on_cran()

# Realised rules: the committed replicate-5 Hhat of the S2 case example
# (mrct-consistency results/s2_example/s2_example_definitions.md), the Hhat of
# replicates 1 and 2 (s2_example_trials.csv), and H* = {z_cd40 <= 230}.
ex_rules <- c("{z_cd40 <= 236} & !{z_preanti <= 196}",
              "{z_cd40 <= 241} & {z_preanti <= 243}",
              "{z_cd40 <= 212} & !{z_drugs}",
              "{z_cd40 <= 230}")
ex_rep5 <- ex_rules[1]
# Administrative censoring only: cens_adjust = 50, the campaigns' convention.
# max_entry = 24 is the entry maximum of fs_build_eval_frame()'s draw
# (simulate_from_dgm()'s default, which fs_build_eval_frame() does not set).
ex_at <- 40; ex_me <- 24; ex_ca <- 50
ex_seeds <- c(1L, 2L, 3L)

ex_memo <- new.env()
ex_get <- function(key, expr) {
  if (is.null(ex_memo[[key]])) ex_memo[[key]] <- expr
  ex_memo[[key]]
}
ex_exact <- function() ex_get("exact", fs_attach_betaHhat(
  data.frame(sg_def = ex_rules), NULL, focus = "harm", method = "exact",
  dgm = mrct_build_S2(), analysis_time = ex_at, max_entry = ex_me,
  cens_adjust = ex_ca))
ex_frames <- function() ex_get("frames", lapply(ex_seeds, function(sd) {
  f <- fs_build_eval_frame(mrct_build_S2(), eval_seed = sd,
                           analysis_time = ex_at, cens_adjust = ex_ca)
  fs_attach_betaHhat(data.frame(sg_def = ex_rules), f, focus = "harm")
}))

test_that("exact and frame agree within 3 SD of the frame route's Monte Carlo spread", {
  mrct_skip_if_no_build()
  ex <- ex_exact(); fr <- ex_frames()
  # same shape as the frame route
  expect_identical(names(ex), names(fr[[1]]))
  expect_identical(vapply(ex, class, ""), vapply(fr[[1]], class, ""))
  expect_identical(ex$betaHhat_status, rep("ok", length(ex_rules)))
  # same membership: the frame draw has every super-population row once
  expect_identical(ex$nH_eval, fr[[1]]$nH_eval)
  expect_identical(ex$nHc_eval, fr[[1]]$nHc_eval)
  for (blk in c("betaHhat_H", "betaHhat_Hc")) {
    L <- sapply(fr, function(r) log(r[[blk]]))
    sdv <- apply(L, 1, stats::sd)
    expect_true(all(is.finite(sdv) & sdv > 0))
    z <- (log(ex[[blk]]) - rowMeans(L)) / sdv
    expect_true(all(abs(z) <= 3), info = paste(blk, paste(round(z, 3), collapse = " ")))
  }
})

test_that("exact equals ahr_limit() on the same membership, on the whole super-population", {
  mrct_skip_if_no_build()
  dgm <- mrct_build_S2()
  ex <- ex_exact()
  pop <- ahr_population(dgm, analysis_time = ex_at, max_entry = ex_me)
  H <- get_dfpred(dgm$df_super, strsplit(ex_rep5, " & ", fixed = TRUE)[[1]])$treat.recommend == 0
  expect_identical(sum(H), ex$nH_eval[1])
  expect_identical(ex$betaHhat_H[1], ahr_limit(pop, set = H, psi = pop$psi, variance = FALSE)$hr)
  expect_identical(ex$betaHhat_Hc[1], ahr_limit(pop, set = !H, psi = pop$psi, variance = FALSE)$hr)
})

test_that("replicate 5: exact beta(Hhat) on the non-region patients' population is ahr_limit() on Hhat & non-region", {
  mrct_skip_if_no_build()
  dgm <- mrct_build_S2()
  zr <- dgm$mrct$region$z_name
  # the non-region patients' population: the DGM with df_super restricted
  dn <- dgm
  dn$df_super <- dgm$df_super[dgm$df_super[[zr]] == 0, , drop = FALSE]
  one <- fs_betaHhat_one(ex_rep5, NULL, focus = "harm", method = "exact", dgm = dn,
                         analysis_time = ex_at, max_entry = ex_me, cens_adjust = ex_ca)
  pop <- ahr_population(dgm, analysis_time = ex_at, max_entry = ex_me)
  H <- get_dfpred(dgm$df_super, strsplit(ex_rep5, " & ", fixed = TRUE)[[1]])$treat.recommend == 0
  nonreg <- dgm$df_super[[zr]] == 0
  expect_identical(one$nH_eval, sum(H & nonreg))
  expect_identical(one$nHc_eval, sum(!H & nonreg))
  expect_equal(one$betaHhat_H,
               ahr_limit(pop, set = H & nonreg, psi = pop$psi, variance = FALSE)$hr,
               tolerance = 1e-12)
  expect_equal(one$betaHhat_Hc,
               ahr_limit(pop, set = !H & nonreg, psi = pop$psi, variance = FALSE)$hr,
               tolerance = 1e-12)
  # The committed case example (mrct-consistency results/s2_example/) records
  # sample Cox HRs of the trial only, no population value of exp(beta(Hhat^c))
  # for the all-patient complement, so that comparison has nothing to check.
})

test_that("exact with random censoring stops with the documented message", {
  mrct_skip_if_no_build()
  dgm <- mrct_build_S2()
  rr <- data.frame(sg_def = ex_rep5)
  for (ca in c(0, log(1.5))) {
    expect_error(fs_attach_betaHhat(rr, NULL, focus = "harm", method = "exact", dgm = dgm,
                                    analysis_time = ex_at, max_entry = ex_me, cens_adjust = ca),
                 "covers administrative censoring only")
  }
  # the defaults of the exact route are the frame route's (cens_adjust = log(1.5))
  expect_error(fs_attach_betaHhat(rr, NULL, focus = "harm", method = "exact", dgm = dgm,
                                  analysis_time = ex_at),
               "covers administrative censoring only")
  # Inf removes the random censoring exactly
  a <- fs_attach_betaHhat(rr, NULL, focus = "harm", method = "exact", dgm = dgm,
                          analysis_time = ex_at, max_entry = ex_me, cens_adjust = Inf)
  expect_identical(a$betaHhat_H, ex_exact()$betaHhat_H[1])
  # survival only, and a DGM is required
  expect_error(fs_betaHhat_one(ex_rep5, NULL, focus = "harm", outcome_type = "binary",
                               method = "exact", dgm = dgm, cens_adjust = ex_ca),
               "survival only")
  expect_error(fs_betaHhat_one(ex_rep5, NULL, focus = "harm", method = "exact",
                               cens_adjust = ex_ca), "needs `dgm`")
})

test_that("a rule with no members, and an unresolvable rule, are NA on both routes", {
  mrct_skip_if_no_build()
  dgm <- mrct_build_S2()
  rr <- data.frame(sg_def = c("{z_cd40 <= -1}", "{z_nosuch <= 1}", NA))
  ex <- fs_attach_betaHhat(rr, NULL, focus = "harm", method = "exact", dgm = dgm,
                           analysis_time = ex_at, max_entry = ex_me, cens_adjust = ex_ca)
  f <- fs_build_eval_frame(dgm, eval_seed = ex_seeds[1], analysis_time = ex_at,
                           cens_adjust = ex_ca)
  fr <- fs_attach_betaHhat(rr, f, focus = "harm")
  for (r in list(ex, fr)) {
    # empty: H is NA, the complement is the whole population
    expect_true(is.na(r$betaHhat_H[1]))
    expect_identical(r$nH_eval[1], 0L)
    expect_true(is.finite(r$betaHhat_Hc[1]))
    # unresolvable: all NA
    expect_identical(r$betaHhat_status[2], "unresolved")
    expect_true(is.na(r$betaHhat_H[2]) && is.na(r$betaHhat_Hc[2]))
    # undetected: the no-subgroup record, the same as the empty rule
    expect_identical(r$betaHhat_Hc[3], r$betaHhat_Hc[1])
  }
  expect_identical(fs_betaHhat_counts(ex), fs_betaHhat_counts(fr))
  # the exact ITT is ahr_limit() on all patients
  pop <- ahr_population(dgm, analysis_time = ex_at, max_entry = ex_me)
  expect_identical(ex$betaHhat_Hc[1],
                   ahr_limit(pop, set = rep(TRUE, pop$N), psi = pop$psi, variance = FALSE)$hr)
})

test_that("method = 'frame' is the call without the argument", {
  mrct_skip_if_no_build()
  f <- fs_build_eval_frame(mrct_build_S2(), eval_seed = ex_seeds[1],
                           analysis_time = ex_at, cens_adjust = ex_ca)
  rr <- data.frame(sg_def = ex_rules)
  expect_identical(fs_attach_betaHhat(rr, f, focus = "harm", method = "frame"),
                   fs_attach_betaHhat(rr, f, focus = "harm"))
  expect_identical(fs_attach_betaHhat(rr, f, focus = "harm"), ex_frames()[[1]])
})
