# Builders for the MRCT analytic-layer tests (test-ahr-analytic.R,
# test-mrct-dgm-args.R). Functions only: nothing is built until a test asks,
# and each build is memoised for the session.
#
# The S2 cell of the mrct-consistency identification campaign at r = -log 4
# (sims/ident_rlog4_full_design.R, build_full("S2") = build_ident("S2")):
# ideal_set_m(build_maps_base(-log 4, b = 0.55, x3_loghr = log(2 / 0.55),
# or_band = 1, cut = 230), m = 3), its inject_mrct_structure() call replicated
# here with every argument as committed (sims/analytic_maps_design.R
# build_maps_base(); seed data of sims/cells.R and sims/ideal_design.R), and
# the hazard multiplier applied with scale_dgm_hazard().

mrct_fixture <- function(f) testthat::test_path("fixtures", "mrct_analytic", f)

mrct_seed_data <- function() {
  e <- new.env()
  utils::data("ACTG175", package = "speff2trial", envir = e)
  d <- e$ACTG175[e$ACTG175$arms %in% c(1, 3), ]
  d$treat_mrct <- as.integer(d$arms == 1)
  rownames(d) <- NULL
  d$time_months <- d$days / 30.4375
  d
}

mrct_memo <- new.env()

# The S2 base build (m = 1) at region prognostic log-HR r; extra arguments
# (e.g. latent_R) go to inject_mrct_structure().
mrct_build_S2_base <- function(r = -log(4), ...) {
  extra <- list(...)
  key <- if (length(extra)) NULL else sprintf("S2_base_r%.12g", r)
  if (!is.null(key) && !is.null(mrct_memo[[key]])) return(mrct_memo[[key]])
  b <- 0.55; h <- 2; cut <- 230
  w <- character(0)
  dgm <- withCallingHandlers(
    do.call(inject_mrct_structure, c(list(
      mrct_seed_data(),
      outcome_var = "time_months", event_var = "cens", treatment_var = "treat_mrct",
      continuous_vars = c("age", "wtkg", "karnof", "preanti", "cd40", "cd80"),
      factor_vars = c("hemo", "homo", "drugs", "race", "gender", "symptom"),
      x_pred = "cd40",
      spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(b), 3)),
      region = list(prevalence = 0.20, or_pred = 1, x1_vars = "symptom", or_x1 = 3, loghr = r,
                    band = list(cut = cut, or = 1)),
      x3 = list(vars = "cd40", cuts = list(cd40 = cut), loghr = log(h / b)),
      region_treat_loghr = NULL,
      n_super = 20000, expand = "copula", seed = 7L,
      continuous_vars_cens = c("age", "cd40"), factor_vars_cens = "symptom"), extra)),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  attr(dgm, "build_warnings") <- w
  if (!is.null(key)) mrct_memo[[key]] <- dgm
  dgm
}

# The S2 cell: hazards x 3.
mrct_build_S2 <- function(r = -log(4)) scale_dgm_hazard(mrct_build_S2_base(r), 3)

# The analytic population of S2 (analysis at month 40, entry U(0, 11)), memoised.
mrct_pop_S2 <- function(r = -log(4), entry_delay = NULL) {
  key <- sprintf("pop_r%.12g_%s", r, paste(format(unlist(entry_delay)), collapse = "_"))
  if (is.null(mrct_memo[[key]]))
    mrct_memo[[key]] <- ahr_population(mrct_build_S2(r), analysis_time = 40, max_entry = 11,
                                       entry_delay = entry_delay)
  mrct_memo[[key]]
}

mrct_skip_if_no_build <- function() {
  testthat::skip_if_not_installed("speff2trial")
  if (!requireNamespace("latentcor", quietly = TRUE))
    testthat::skip("latentcor is not installed: the S2 build of the records needs the latentcor fit")
}
