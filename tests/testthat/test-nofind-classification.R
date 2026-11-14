# A-003: a replicate that identifies no subgroup has the whole sample as its
# estimated complement, so spec = sens(Hc) = 1 and npv = ppv(Hc) = |Hc|/N
# (Leon et al. 2024, Section 3), the values the found-case formulas
# tn / (tn + fp) and tn / (tn + fn) give for an empty identified set.
# Transplanted from the A-003 reproducer
# (dev/reports/REPORT_audit_revisit_2026-11-11.md, section A-003).
# No search, no simulation: the extractors are called on hand-built inputs.

skip_on_cran()

.a003_data <- function(null = FALSE) {
  set.seed(1)
  N <- 200; nH <- 50
  data.frame(id = seq_len(N), y = rexp(N), ev = rbinom(N, 1, 0.7),
             tr = rbinom(N, 1, 0.5),
             flag_harm = if (null) rep(0L, N) else c(rep(1L, nH), rep(0L, N - nH)))
}

.a003_fs <- function(df, fs_res = NULL, fs_full = NULL) {
  .extract_fs_estimates_gen(
    df = df, fs_res = fs_res, fs_full = fs_full, dgm = NULL,
    cox_formula = NULL, cox_formula_adj = NULL,
    outcome_name = "y", event_name = "ev", treat_name = "tr",
    harm_col = "flag_harm", analysis = "FS", verbose = FALSE)
}

.a003_grf <- function(df, grf_est = NULL) {
  .extract_grf_estimates_gen(
    df = df, grf_est = grf_est, dgm = NULL,
    cox_formula = NULL, cox_formula_adj = NULL,
    outcome_name = "y", event_name = "ev", treat_name = "tr",
    harm_col = "flag_harm", analysis = "GRF", verbose = FALSE, debug = FALSE)
}

# Found-case formulas of the two extractors, evaluated on an all-FALSE
# identified set.
.a003_formula <- function(df) {
  true_H <- df$flag_harm == 1L
  hat_H  <- rep(FALSE, nrow(df))
  tn <- sum(!true_H & !hat_H); fp <- sum(!true_H & hat_H)
  fn <- sum(true_H & !hat_H)
  c(spec = tn / (tn + fp), npv = tn / (tn + fn))
}

.a003_expect_alt <- function(out, df, label, check_any_H = TRUE) {
  if (check_any_H) expect_identical(out$any.H, 0L, label = paste(label, "any.H"))
  expect_equal(out$sens, 0,    label = paste(label, "sens"))
  expect_equal(out$ppv,  0,    label = paste(label, "ppv"))
  expect_equal(out$spec, 1,    label = paste(label, "spec"))
  expect_equal(out$npv,  0.75, label = paste(label, "npv"))
  f <- .a003_formula(df)
  expect_equal(out$spec, f[["spec"]], label = paste(label, "spec vs tn/(tn+fp)"))
  expect_equal(out$npv,  f[["npv"]],  label = paste(label, "npv vs tn/(tn+fn)"))
}

test_that("A-003 no-find, alternative: FS site 1 (no result table)", {
  df <- .a003_data()
  .a003_expect_alt(.a003_fs(df), df, "FS site 1")
})

test_that("A-003 no-find, alternative: GRF site 3 (no GRF estimate)", {
  df <- .a003_data()
  .a003_expect_alt(.a003_grf(df), df, "GRF site 3")
})

test_that("A-003 no-find, alternative: FS site 2 (rows, no recoverable M.* factor)", {
  df <- .a003_data()
  fs_res <- data.frame(M.1 = NA_character_, stringsAsFactors = FALSE)
  .a003_expect_alt(.a003_fs(df, fs_res = fs_res), df, "FS site 2")
})

test_that("A-003 no-find, alternative: GRF site 4 (no sg.harm.id / no treat.recommend)", {
  df <- .a003_data()
  no_sg <- list(sg.harm.id = NULL,
                data = data.frame(treat.recommend = rep(1L, nrow(df))))
  no_tr <- list(sg.harm.id = "x1 <= 0", data = df)
  .a003_expect_alt(.a003_grf(df, grf_est = no_sg), df, "GRF site 4, no sg.harm.id")
  .a003_expect_alt(.a003_grf(df, grf_est = no_tr), df, "GRF site 4, no treat.recommend")
})

test_that("A-003 no-find, null: spec and npv are 1 in both extractors", {
  df <- .a003_data(null = TRUE)
  outs <- list(
    `FS site 1`  = .a003_fs(df),
    `FS site 2`  = .a003_fs(df, fs_res = data.frame(M.1 = NA_character_)),
    `GRF site 3` = .a003_grf(df),
    `GRF site 4` = .a003_grf(df, grf_est = list(data = NULL)))
  for (nm in names(outs)) {
    expect_equal(outs[[nm]]$spec, 1, label = paste(nm, "spec"))
    expect_equal(outs[[nm]]$npv,  1, label = paste(nm, "npv"))
    expect_true(is.na(outs[[nm]]$sens), label = paste(nm, "sens NA"))
    expect_true(is.na(outs[[nm]]$ppv),  label = paste(nm, "ppv NA"))
  }
})

# Empty identified set (Rule E, TASK_empty_set_nofind_2026-11-11): whenever no
# subject is in the identified set, by whatever path, the replicate is scored
# exactly as the no-find path scores it for the same data.  Hand-built inputs;
# no search.

.es_same <- function(a, b, label) {
  expect_identical(as.list(a), as.list(b), label = paste(label, "columns"))
  expect_identical(attr(a, "sg_hat"), attr(b, "sg_hat"), label = paste(label, "sg_hat"))
}

test_that("empty set: FS site 2 (rows, no recoverable M.* factor) equals site 1", {
  for (null in c(FALSE, TRUE)) {
    df <- .a003_data(null = null)
    site1 <- .a003_fs(df)
    site2 <- .a003_fs(df, fs_res = data.frame(M.1 = NA_character_,
                                              stringsAsFactors = FALSE))
    expect_identical(site2$any.H, 0L)
    .es_same(site2, site1, paste("FS site 2, null =", null))
  }
})

test_that("empty set: GRF found path with no member equals the no-find output", {
  for (null in c(FALSE, TRUE)) {
    df <- .a003_data(null = null)
    out <- .a003_grf(df, grf_est = list(sg.harm.id = "x <= 0",
                                        data = cbind(df, treat.recommend = 1)))
    lab <- paste("GRF found, empty, null =", null)
    expect_identical(out$any.H, 0L, label = paste(lab, "any.H"))
    if (null) expect_true(is.na(out$ppv), label = paste(lab, "ppv NA"))
    else expect_identical(out$ppv, 0, label = paste(lab, "ppv 0"))
    .es_same(out, .a003_grf(df), lab)
  }
})

test_that("empty set: FS found path with no member equals the no-find output", {
  for (null in c(FALSE, TRUE)) {
    df <- .a003_data(null = null)
    nofind <- .a003_fs(df)
    fs_res <- data.frame(M.1 = "{x <= 0}", stringsAsFactors = FALSE)
    # treat.recommend all 1: no subject in the identified set
    via_df_est <- .a003_fs(df, fs_res = fs_res,
                           fs_full = list(df.est = data.frame(treat.recommend = rep(1, nrow(df)))))
    # M.* factor recovered but matching no subject (no df.est)
    via_m <- .a003_fs(df, fs_res = data.frame(M.1 = "flag_harm.7",
                                              stringsAsFactors = FALSE))
    for (nm in c("df.est", "M.*")) {
      out <- if (nm == "df.est") via_df_est else via_m
      lab <- paste("FS found, empty via", nm, "null =", null)
      expect_identical(out$any.H, 0L, label = paste(lab, "any.H"))
      if (null) expect_true(is.na(out$ppv), label = paste(lab, "ppv NA"))
      else expect_identical(out$ppv, 0, label = paste(lab, "ppv 0"))
      .es_same(out, nofind, lab)
    }
  }
})

test_that("non-empty set: FS and GRF found paths keep the found-case formulas", {
  df <- .a003_data()
  rec <- c(rep(0, 60), rep(1, nrow(df) - 60))      # 60 in the identified set
  true_H <- df$flag_harm == 1L; hat_H <- rec == 0
  tp <- sum(true_H & hat_H); fp <- sum(!true_H & hat_H)
  tn <- sum(!true_H & !hat_H); fn <- sum(true_H & !hat_H)
  want <- c(sens = tp / (tp + fn), spec = tn / (tn + fp),
            ppv = tp / (tp + fp), npv = tn / (tn + fn))
  fs <- .a003_fs(df, fs_res = data.frame(M.1 = "{x <= 0}", stringsAsFactors = FALSE),
                 fs_full = list(df.est = data.frame(treat.recommend = rec)))
  grf <- .a003_grf(df, grf_est = list(sg.harm.id = "x <= 0",
                                      data = cbind(df, treat.recommend = rec)))
  for (nm in c("FS", "GRF")) {
    out <- if (nm == "FS") fs else grf
    expect_identical(out$any.H, 1L, label = paste(nm, "any.H"))
    expect_identical(out$size.H, 60L, label = paste(nm, "size.H"))
    for (k in names(want))
      expect_equal(out[[k]], want[[k]], label = paste(nm, k))
  }

  # FS M.* path (no df.est): a recovered factor with members
  fs_m <- .a003_fs(df, fs_res = data.frame(M.1 = "flag_harm.1",
                                           stringsAsFactors = FALSE))
  hat_H <- df$flag_harm == 1L
  tp <- sum(true_H & hat_H); fp <- sum(!true_H & hat_H)
  tn <- sum(!true_H & !hat_H); fn <- sum(true_H & !hat_H)
  expect_identical(fs_m$any.H, 1L)
  expect_identical(fs_m$size.H, 50L)
  expect_equal(c(fs_m$sens, fs_m$spec, fs_m$ppv, fs_m$npv),
               c(tp / (tp + fn), tn / (tn + fp), tp / (tp + fp), tn / (tn + fn)))
})
