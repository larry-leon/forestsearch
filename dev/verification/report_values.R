suppressMessages(devtools::load_all(quiet = TRUE))
env <- new.env()
# reuse the test file's fixtures and helpers without running expectations
for (e in parse("tests/testthat/test-declaration-calibration.R"))
  if (!(is.call(e) && identical(e[[1]], as.name("test_that")))) eval(e, envir = env)
with(env, {
  fmt <- function(x, d = 6) formatC(x, digits = d, format = "f")
  pre <- fs_declaration_calibration(.fit_on); pre10 <- fs_declaration_calibration(.fit_on, alpha = 0.10)
  red <- fs_declaration_calibration(.fit_on, family = "reduced"); red10 <- fs_declaration_calibration(.fit_on, alpha = 0.10, family = "reduced")
  cat("GBSG: n_pre", pre$n_family_prereduction, "n_red", pre$n_family_reduced, "removed", pre$reduction$removed, "\n")
  cat("screened", length(pre$screened), "admitted_current", pre$admitted_current, "replay_check", pre$reduction$replay_check, "\n")
  cat("PRE kappa05", fmt(pre$kappa_hat), "kappa10", fmt(pre10$kappa_hat), "fw", fmt(pre$fw_size), "\n")
  cat("RED kappa05", fmt(red$kappa_hat), "kappa10", fmt(red10$kappa_hat), "fw", fmt(red$fw_size), "\n")
  cat("admitted_pstar(pre)", pre$admitted_pstar, "| admitted_calibrated(pre) n=", length(pre$admitted_calibrated), "\n")
  cat("T_hat selected", fmt(pre$T_hat[pre$admitted_current], 4), "z_pstar (rounded rule)", fmt(pre$z_pstar, 4), "\n")
  cat("Mstar quantiles:", fmt(quantile(pre$Mstar, c(.5,.9,.95,.99), type = 1), 4), "\n")
  cat("PC4 GBSG column_sd range", fmt(range(pre$column_sd), 4), "mean(Zstar)", fmt(pre$zstar_mean, 5), "tol", fmt(4/sqrt(pre$B), 4), "\n")
  d8 <- fs_declaration_calibration(.mr8)
  cat("PC4 OLS column_sd", fmt(d8$column_sd, 5), "mean", fmt(d8$zstar_mean, 6), "\n")
  ids <- .fit_on$grp.consistency$df_flag$id[.fit_on$grp.consistency$sg.harm.id == 1]
  sub <- .gb[.gb$id %in% ids, ]
  init <- survival::coxph(survival::Surv(time_months, status) ~ hormon, data = sub)$coefficients[1]
  rr <- consistency_resample(sub, method = "closed", tte.name = "time_months", event.name = "status", treat.name = "hormon", cox_init = init)
  cat("PC6 GBSG sigma_D field", fmt(pre$sigma_D[pre$admitted_current], 10), "screen", fmt(rr$sigma_D, 10), "diff", signif(pre$sigma_D[pre$admitted_current] - rr$sigma_D, 3), "\n")
  for (g in 1:2) { r8 <- consistency_resample(.cfgB$df[.cfgB$cands[[g]], ], method = "closed", outcome_type = "continuous", effect_measure = "MD", treat.name = "A", outcome.name = "Y", consistency_threshold = 0)
    cat("PC6 OLS g", g, "field", fmt(d8$sigma_D[g], 10), "screen", fmt(r8$sigma_D, 10), "\n") }
  fwt <- 1 - .phi2(d8$z_pstar, .rho8); kt <- .kappa_target(.rho8)   # fw target at the calibration's rounded cutoff
  cat("T8 rho_hat", fmt(.rho8), "targets: cor", fmt(.rho8), "fw", fmt(fwt), "kappa", fmt(kt), "\n")
  cat("T8 gaussian realized: cor", fmt(d8$field_cor[1,2]), "fw", fmt(d8$fw_size), "kappa", fmt(d8$kappa_hat), "\n")
  cat("T8 gaussian disc: cor", fmt(d8$field_cor[1,2]-.rho8), "fw", fmt(d8$fw_size-fwt), "kappa", fmt(d8$kappa_hat-kt), "\n")
  set.seed(808L); xi <- matrix(rnorm(nrow(.db8) * .B8), nrow(.db8), .B8)
  dh <- fs_declaration_calibration(.decl_fit(forestsearch:::.fs_decl_field(.db8, xi, keep_matrix = FALSE), bh = c(0, 0)))
  cat("T8 helper disc: cor", fmt(dh$field_cor[1,2]-.rho8), "fw", fmt(dh$fw_size-fwt), "kappa", fmt(dh$kappa_hat-kt), "\n")
  mp <- fs_mr_inference(.cfgB$df, .cfgB$cands, .spec_md, selected_members = .cfgB$cands$g1, admission = .adm_B,
    reselection = "maxeff", draws = .B8, multiplier = "poisson", ci_method = "ij", seed = 7L, keep_declaration_field = TRUE)
  dp <- fs_declaration_calibration(mp)
  cat("T8 poisson (production law) realized: cor", fmt(dp$field_cor[1,2]), "fw", fmt(dp$fw_size), "kappa", fmt(dp$kappa_hat), "\n")
  cat("T8 poisson disc: cor", fmt(dp$field_cor[1,2]-.rho8), "fw", fmt(dp$fw_size-fwt), "kappa", fmt(dp$kappa_hat-kt), "\n")
})
# D6: GLM (continuous) forestsearch path
source("tests/testthat/helper-synthetic-dgm.R")
cd <- .make_continuous_data(N = 400L, MD_harm = 2)
cargs <- .fs_args_for("continuous", confounders = c("age", "biomarker", "biomarker_hi", "sex"),
  extra = list(use_grf = FALSE, use_lasso = FALSE, quiet = TRUE, hr.threshold = 0.5, hr.consistency = 0.25,
               mr_inference = TRUE, mr_inference_args = list(draws = 300L, keep_declaration_field = TRUE, keep_field_matrix = TRUE)))
fc <- suppressWarnings(do.call(forestsearch, c(list(df.analysis = cd), cargs)))
dg <- fs_declaration_calibration(fc)
cat("GLM: n_pre", dg$n_family_prereduction, "n_red", dg$n_family_reduced, "replay_check", dg$reduction$replay_check,
    "unmatched", dg$reduction$n_unmatched, "admitted_current", dg$admitted_current, "kappa", dg$kappa_hat, "fw", dg$fw_size, "\n")
s <- dg$screened; rel <- intersect(s, dg$admitted_pstar)   # the package's rounded rule, not a local rebuild
cat("GLM relabel reproduces:", setequal(rel, dg$admitted_current), "\n")
cat("GLM hr.subgroups names:", paste(names(fc$find.grps$out.found$hr.subgroups)[1:10], collapse=","), "\n")
