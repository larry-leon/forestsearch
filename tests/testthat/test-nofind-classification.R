# A-003: a replicate that identifies no subgroup has the whole sample as its
# estimated complement, so spec = sens(Hc) = 1 and npv = ppv(Hc) = |Hc|/N
# (Leon et al. 2024, Section 3), the values the found-case formulas
# tn / (tn + fp) and tn / (tn + fn) give for an empty identified set.
# Transplanted from the A-003 reproducer
# (dev/reports/REPORT_audit_revisit_2026-11-11.md, section A-003).
# No search, no simulation: the extractors are called on hand-built inputs.

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
  # Site 2 sets any.H before it finds no factor; that flag is outside A-003,
  # so it is not checked here.
  .a003_expect_alt(.a003_fs(df, fs_res = fs_res), df, "FS site 2",
                   check_any_H = FALSE)
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
