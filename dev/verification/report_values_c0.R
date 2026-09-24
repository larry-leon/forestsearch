suppressMessages(devtools::load_all(quiet = TRUE))
env <- new.env()
for (e in parse("tests/testthat/test-declaration-c0.R"))
  if (!(is.call(e) && identical(e[[1]], as.name("test_that")))) eval(e, envir = env)
with(env, {
  f <- function(x, d = 6) formatC(x, digits = d, format = "f")
  t05 <- fs_declaration_calibration(.fit_c0, alpha = 0.05, c0 = .c0_grid)$c0$table
  t10 <- fs_declaration_calibration(.fit_c0, alpha = 0.10, c0 = .c0_grid)$c0$table
  d0 <- fs_declaration_calibration(.fit_c0)
  cat("n_family", d0$n_family_prereduction, "admitted_current", d0$admitted_current, "p*", d0$p_star, "c_cons", d0$c_cons, "\n")
  cat("GBSG table: c0 | kappa05 | kappa10 | pstar_settable05 | fw(p*0.90, rounded rule) | n_adm_cal05 | q90 q95 q99\n")
  for (i in seq_len(nrow(t05))) cat(t05$c0[i], "|", f(t05$kappa_hat[i]), "|", f(t10$kappa_hat[i]), "|", f(t05$pstar_settable[i], 2),
    "|", f(t05$fw_size[i]), "|", t05$n_admitted_calibrated[i], "|", f(t05$Mstar_c0_q90[i],4), f(t05$Mstar_c0_q95[i],4), f(t05$Mstar_c0_q99[i],4), "\n")
  dl <- (0 - log(.c0_grid)) / d0$sigma_D
  cat("GBSG delta_g range per c0:\n"); print(t(apply(outer(1/d0$sigma_D, -log(.c0_grid)), 2, range)))
  cat("GBSG sigma_D range", f(range(d0$sigma_D)), "\n")
  cat("T5 rho", f(.rho5), "sig", f(.sig5, 10), "delta", f(.dl5), "fw target", f(.fw5_target), "kappa target", f(.kappa5_target), "\n")
  mr <- fs_mr_inference(.cfgB0$df, .cfgB0$cands, .spec_md0, selected_members = .cfgB0$cands$g1,
    admission = .adm_B0, reselection = "maxeff", draws = .B5, multiplier = "gaussian", ci_method = "ij", seed = 7L,
    keep_declaration_field = TRUE, declaration_c0 = .c0_5)
  d5 <- fs_declaration_calibration(mr, c0 = .c0_5); tb <- d5$c0$table
  cat("T5 exported realized fw", f(tb$fw_size), "kappa", f(tb$kappa_hat), "disc", f(tb$fw_size - .fw5_target), f(tb$kappa_hat - .kappa5_target), "\n")
  set.seed(808L); xi <- matrix(stats::rnorm(nrow(.db5) * .B5), nrow(.db5), .B5)
  m <- forestsearch:::.fs_decl_field(.db5, xi, keep_matrix = FALSE, shift = cbind(.dl5))$Mstar_shift[, 1]
  fw <- mean(m > d5$z_pstar)   # the calibration's rounded cutoff, as .fw5_target
  k <- quantile(m, 0.95, type = 1, names = FALSE)
  cat("T5 helper realized fw", f(fw), "kappa", f(k), "disc", f(fw - .fw5_target), f(k - .kappa5_target), "\n")
  cat("sdv vs field sigma identical (GBSG):", identical(unname(.fit_c0$mr_inference$declaration_field$sigma_D), unname(.fit_c0$mr_inference$declaration_field$meta$sigma_D_field)), "\n")
  gm <- .cont_fit$mr_inference$declaration_field
  cat("GLM: c_cons", gm$meta$c_cons, "cols", colnames(gm$Mstar_c0), "G", gm$meta$G, "identical c2 col", identical(unname(gm$Mstar_c0[, "0.25"]), gm$Mstar), "\n")
  print(fs_declaration_calibration(.fit_c0, c0 = .c0_grid))
})
