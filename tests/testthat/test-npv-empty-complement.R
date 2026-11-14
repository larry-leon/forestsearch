# Rule F (TASK_cran_prep_2026-11-11, F4): when the identified set is
# everyone, the estimated complement is empty and npv = tn / (tn + fn) is
# 0/0.  It is scored 0, by the same 0/0 convention that gives ppv = 0 for an
# empty identified set.  fs_mr_oc_summary() applies it to a detected
# replicate with no stored npv whose identified set is everyone
# (n_harm equals meta$n_sample); every other missing npv keeps its fill.
# Hand-built inputs; no search, no simulation.

skip_on_cran()

.f4_data <- function(null = FALSE) {
  set.seed(1)
  N <- 200; nH <- 50
  data.frame(id = seq_len(N), y = rexp(N), ev = rbinom(N, 1, 0.7),
             tr = rbinom(N, 1, 0.5),
             flag_harm = if (null) rep(0L, N) else c(rep(1L, nH), rep(0L, N - nH)))
}

.f4_fs <- function(df, fs_res = NULL, fs_full = NULL) {
  .extract_fs_estimates_gen(
    df = df, fs_res = fs_res, fs_full = fs_full, dgm = NULL,
    cox_formula = NULL, cox_formula_adj = NULL,
    outcome_name = "y", event_name = "ev", treat_name = "tr",
    harm_col = "flag_harm", analysis = "FS", verbose = FALSE)
}

.f4_grf <- function(df, grf_est = NULL) {
  .extract_grf_estimates_gen(
    df = df, grf_est = grf_est, dgm = NULL,
    cox_formula = NULL, cox_formula_adj = NULL,
    outcome_name = "y", event_name = "ev", treat_name = "tr",
    harm_col = "flag_harm", analysis = "GRF", verbose = FALSE, debug = FALSE)
}

.f4_payload <- function(results, n_sample = 200L) {
  list(results = results,
       truth = list(effect_Q = 1, effect_Qc = 0, beta_inter = 1,
                    prevalence_Q = 0.25, marg_H = 1, marg_Hc = 0),
       meta = list(n_sample = n_sample, n_sims = nrow(results),
                   adverse_outcome = TRUE, mr_draws = 10L,
                   sg_focus = "test", pkg_version = "test"))
}

.f4_u <- function(pl) {
  id <- fs_mr_oc_summary(pl)$identification
  id[id$convention == "unconditional", ]
}

test_that("identified set is everyone: both extractors give npv = 0", {
  for (null in c(FALSE, TRUE)) {
    df <- .f4_data(null = null)
    N  <- nrow(df)
    outs <- list(
      `FS df.est` = .f4_fs(df,
        fs_res  = data.frame(M.1 = "{x <= 0}", stringsAsFactors = FALSE),
        fs_full = list(df.est = data.frame(treat.recommend = rep(0L, N)))),
      `GRF` = .f4_grf(df, grf_est = list(sg.harm.id = "x <= 0",
                                         data = cbind(df, treat.recommend = 0L))))
    for (nm in names(outs)) {
      out <- outs[[nm]]
      lab <- paste(nm, "everyone, null =", null)
      expect_identical(out$any.H, 1L, label = paste(lab, "any.H"))
      expect_identical(out$size.H, N, label = paste(lab, "size.H"))
      expect_identical(out$size.Hc, 0L, label = paste(lab, "size.Hc"))
      expect_equal(out$npv, 0, label = paste(lab, "npv"), tolerance = 1e-12)
      # the other three classification metrics keep the found-case formulas
      expect_equal(out$spec, 0, label = paste(lab, "spec"), tolerance = 1e-12)
      if (null) {
        expect_equal(out$ppv, 0, label = paste(lab, "ppv"), tolerance = 1e-12)
        expect_true(is.na(out$sens), label = paste(lab, "sens NA"))
      } else {
        expect_equal(out$ppv, 50 / 200, label = paste(lab, "ppv"), tolerance = 1e-12)
        expect_equal(out$sens, 1, label = paste(lab, "sens"), tolerance = 1e-12)
      }
    }
  }
})

test_that("non-empty complement: both extractors keep npv = tn / (tn + fn)", {
  df  <- .f4_data()
  N   <- nrow(df)
  rec <- rep(1L, N); rec[21:80] <- 0L   # identified: rows 21-80
  # truly in H: rows 1-50, so tp = 30, fp = 30, fn = 20, tn = 120
  fs  <- .f4_fs(df,
    fs_res  = data.frame(M.1 = "{x <= 0}", stringsAsFactors = FALSE),
    fs_full = list(df.est = data.frame(treat.recommend = rec)))
  grf <- .f4_grf(df, grf_est = list(sg.harm.id = "x <= 0",
                                    data = cbind(df, treat.recommend = rec)))
  for (out in list(fs, grf)) {
    expect_identical(out$size.Hc, 140L)
    expect_equal(out$npv, 120 / 140)
    expect_equal(out$ppv, 30 / 60)
    expect_equal(out$sens, 30 / 50)
    expect_equal(out$spec, 120 / 150)
  }
})

test_that("fs_mr_oc_summary(): a detected replicate with n_harm = n_sample and no stored npv contributes 0", {
  r <- data.frame(
    sim_id = 1:5,
    status = c("DETECTED", "DETECTED", "DETECTED", "DETECTED", "NO-DETECTION"),
    n_harm = c(200L, 60L, 60L, NA, NA),
    n_true = rep(50L, 5),
    sens = c(1, 0.6, NA, NA, NA), spec = c(0, 0.8, NA, NA, NA),
    ppv  = c(0.25, 0.5, NA, NA, NA),
    npv  = c(NA, 0.9, NA, NA, NA),
    stringsAsFactors = FALSE)

  # the replicate on its own: npv 0
  expect_equal(.f4_u(.f4_payload(r[1, , drop = FALSE]))$npv, 0, tolerance = 1e-12)

  # mixed: everyone -> 0; stored 0.9 kept; detected NA with n_harm < N or
  # n_harm NA -> 1 (unchanged); non-detected -> |Hc|/N = 150/200 (unchanged)
  u <- .f4_u(.f4_payload(r))
  expect_equal(u$npv, (0 + 0.9 + 1 + 1 + 150 / 200) / 5)

  # the detected-only (conditional) row is unchanged: stored values only
  id <- fs_mr_oc_summary(.f4_payload(r))$identification
  expect_equal(id$npv[id$convention == "conditional"], 0.9)

  # a stored npv is never overwritten, whatever n_harm is
  r2 <- r[1, , drop = FALSE]; r2$npv <- 0.3
  expect_equal(.f4_u(.f4_payload(r2))$npv, 0.3)
})

test_that("fs_mr_oc_summary(): without a usable n_sample or n_harm the detected fill stays 1", {
  r <- data.frame(sim_id = 1L, status = "DETECTED", n_harm = 200L,
                  n_true = 50L, sens = 1, spec = 0, ppv = 0.25, npv = NA_real_,
                  stringsAsFactors = FALSE)
  expect_equal(.f4_u(.f4_payload(r, n_sample = NULL))$npv, 1, tolerance = 1e-12)
  expect_equal(.f4_u(.f4_payload(r, n_sample = NA_integer_))$npv, 1, tolerance = 1e-12)
  expect_equal(.f4_u(.f4_payload(r, n_sample = "200"))$npv, 1, tolerance = 1e-12)
  r_no <- r; r_no$n_harm <- NULL
  expect_equal(.f4_u(.f4_payload(r_no))$npv, 1, tolerance = 1e-12)
  r_lt <- r; r_lt$n_harm <- 199L
  expect_equal(.f4_u(.f4_payload(r_lt))$npv, 1, tolerance = 1e-12)
})

test_that("fs_mr_oc_summary() agrees with the extractors on an everyone replicate", {
  df <- .f4_data()
  N  <- nrow(df)
  ext <- .f4_fs(df,
    fs_res  = data.frame(M.1 = "{x <= 0}", stringsAsFactors = FALSE),
    fs_full = list(df.est = data.frame(treat.recommend = rep(0L, N))))
  # The MR harness's record of the same replicate, npv not stored
  r <- data.frame(sim_id = 1L, status = "DETECTED", n_harm = ext$size.H,
                  n_true = sum(df$flag_harm), sens = ext$sens, spec = ext$spec,
                  ppv = ext$ppv, npv = NA_real_, stringsAsFactors = FALSE)
  expect_identical(.f4_u(.f4_payload(r, n_sample = N))$npv, ext$npv)
})
