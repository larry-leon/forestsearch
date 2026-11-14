# run_simulation_analysis() effect scales (TASK_cran_prep_release_2026-11-12,
# Part 1).  In GLM mode the two estimate extractors, .extract_fs_estimates_gen()
# and .extract_grf_estimates_gen(), fill hr.itt, hr.H.hat, hr.Hc.hat, hr.H.true
# and hr.Hc.true from one two-arm glm(y ~ treat) per set.  "RD" and "RR" must
# be the difference and the ratio of the two fitted arm proportions, and
# "IRD" the difference of the two fitted arm rates (per unit of exposure when
# there is an offset, the convention under which the rate ratio is
# exp(coefficient)), treatment minus control, the direction the truth columns
# (compute_aor()) use.  (With an offset, the IRD truth is a difference in
# expected counts per patient, not per unit of exposure; see
# ?run_simulation_analysis.)
#
# Hand-built data: 160 patients in 8 cells of 20, one cell per
# (treat, H, Hhat), where H is the true subgroup (flag_harm) and Hhat the
# identified one.  A two-arm glm with an intercept and a 0/1 treatment term
# is saturated in the two arms, so its fitted arm proportions (rates) are the
# arm event shares (count / exposure); the expected values below are worked
# out from the cell totals by hand.  No search, no simulation.

skip_on_cran()

.es_cells <- data.frame(
  treat = c(1, 1, 1, 1, 0, 0, 0, 0),
  H     = c(1, 1, 0, 0, 1, 1, 0, 0),
  Hhat  = c(1, 0, 1, 0, 1, 0, 1, 0),
  # binary: events out of 20 in each cell
  ev    = c(12, 8, 6, 4, 4, 6, 5, 10),
  # count: total count in each cell; exposure per patient
  cnt   = c(48, 32, 24, 16, 15, 21, 12, 24),
  expo  = c(2, 2, 2, 2, 1.5, 1.5, 1.5, 1.5)
)

# One row per patient.  A cell's total `tot` is spread over its 20 patients
# as whole numbers (tot %/% 20 each, plus 1 for the first tot %% 20).
.es_spread <- function(tot) {
  c(rep(1L, tot %% 20), rep(0L, 20 - tot %% 20)) + tot %/% 20
}

.es_data <- function(kind) {
  rows <- lapply(seq_len(nrow(.es_cells)), function(i) {
    cl <- .es_cells[i, ]
    y <- if (kind == "binary") c(rep(1L, cl$ev), rep(0L, 20 - cl$ev))
         else .es_spread(cl$cnt)
    data.frame(y = y, treat = cl$treat, flag_harm = as.integer(cl$H),
               Hhat = as.integer(cl$Hhat), time = cl$expo)
  })
  df <- do.call(rbind, rows)
  df$id <- seq_len(nrow(df))
  # run_simulation_analysis() sets the event column so: the outcome itself
  # for a binary outcome, 1 for a count.
  df$event <- if (kind == "binary") df$y else 1L
  df
}

# The FS extractor reads the identified set from fs_full$df.est$treat.recommend
# (0 = in Hhat); the GRF extractor from grf_est$data$treat.recommend.
.es_fs <- function(df, em, offset_name = NULL, outcome_type = "binary") {
  .extract_fs_estimates_gen(
    df = df, fs_res = NULL,
    fs_full = list(df.est = data.frame(treat.recommend = 1L - df$Hhat),
                   sg.harm = "Hhat == 1"),
    dgm = NULL, cox_formula = NULL, cox_formula_adj = NULL,
    outcome_name = "y", event_name = "event", treat_name = "treat",
    harm_col = "flag_harm", analysis = "FS", verbose = FALSE,
    outcome_type = outcome_type, effect_measure = em,
    offset_name = offset_name)
}

.es_grf <- function(df, em, offset_name = NULL, outcome_type = "binary") {
  .extract_grf_estimates_gen(
    df = df,
    grf_est = list(sg.harm.id = "Hhat == 1",
                   data = cbind(df, treat.recommend = 1L - df$Hhat)),
    dgm = NULL, cox_formula = NULL, cox_formula_adj = NULL,
    outcome_name = "y", event_name = "event", treat_name = "treat",
    harm_col = "flag_harm", analysis = "GRF", verbose = FALSE,
    debug = FALSE, outcome_type = outcome_type, effect_measure = em,
    offset_name = offset_name, id_name = "id")
}

.es_cols <- c("hr.itt", "hr.H.hat", "hr.Hc.hat", "hr.H.true", "hr.Hc.true")

.es_check <- function(out, expected, lab) {
  for (cc in .es_cols) {
    expect_equal(out[[cc]], expected[[cc]], label = paste(lab, cc))
  }
}

test_that("binary RD and RR: difference and ratio of the fitted arm proportions", {
  # Arm proportions (treatment p1, control p0), from the cell events:
  #   ITT        p1 = (12+8+6+4)/80 = 30/80 = 0.375
  #              p0 = (4+6+5+10)/80 = 25/80 = 0.3125
  #              RD = 0.0625              RR = 0.375/0.3125 = 1.2
  #   Hhat = 1   p1 = (12+6)/40 = 0.45    p0 = (4+5)/40 = 0.225
  #              RD = 0.225               RR = 2
  #   Hhat = 0   p1 = (8+4)/40 = 0.3      p0 = (6+10)/40 = 0.4
  #              RD = -0.1                RR = 0.75
  #   H = 1      p1 = (12+8)/40 = 0.5     p0 = (4+6)/40 = 0.25
  #              RD = 0.25                RR = 2
  #   H = 0      p1 = (6+4)/40 = 0.25     p0 = (5+10)/40 = 0.375
  #              RD = -0.125              RR = 0.25/0.375 = 2/3
  # (Before the fix, "RD" gave the log odds ratio, e.g. ITT
  # log((0.375/0.625)/(0.3125/0.6875)) = log(1.32) = 0.2776, and "RR" the
  # odds ratio 1.32.)
  df <- .es_data("binary")
  rd <- list(hr.itt = 30/80 - 25/80, hr.H.hat = 0.45 - 0.225,
             hr.Hc.hat = 0.3 - 0.4, hr.H.true = 0.5 - 0.25,
             hr.Hc.true = 0.25 - 0.375)
  rr <- list(hr.itt = (30/80) / (25/80), hr.H.hat = 0.45 / 0.225,
             hr.Hc.hat = 0.3 / 0.4, hr.H.true = 0.5 / 0.25,
             hr.Hc.true = 0.25 / 0.375)
  .es_check(.es_fs(df, "RD"),  rd, "FS RD")
  .es_check(.es_grf(df, "RD"), rd, "GRF RD")
  .es_check(.es_fs(df, "RR"),  rr, "FS RR")
  .es_check(.es_grf(df, "RR"), rr, "GRF RR")
})

test_that("count IRD with an offset: difference of the fitted arm rates per unit of exposure", {
  # Exposure (time) is 2 per treated patient and 1.5 per control patient, so
  # each treated cell has 40 units and each control cell 30.  With the offset
  # log(time) the fitted arm rate is count / exposure:
  #   ITT        r1 = (48+32+24+16)/160 = 120/160 = 0.75
  #              r0 = (15+21+12+24)/120 = 72/120 = 0.6     IRD = 0.15
  #   Hhat = 1   r1 = (48+24)/80 = 0.9    r0 = (15+12)/60 = 0.45   IRD = 0.45
  #   Hhat = 0   r1 = (32+16)/80 = 0.6    r0 = (21+24)/60 = 0.75   IRD = -0.15
  #   H = 1      r1 = (48+32)/80 = 1      r0 = (15+21)/60 = 0.6    IRD = 0.4
  #   H = 0      r1 = (24+16)/80 = 0.5    r0 = (12+24)/60 = 0.6    IRD = -0.1
  # The same rates give the rate ratio r1/r0 that "IRR" returns
  # (ITT 0.75/0.6 = 1.25), so IRD and IRR share one exposure convention.
  # (Before the fix, "IRD" gave the log rate ratio, e.g. ITT log(1.25).)
  df <- .es_data("count")
  r1 <- list(hr.itt = 120/160, hr.H.hat = 72/80, hr.Hc.hat = 48/80,
             hr.H.true = 80/80, hr.Hc.true = 40/80)
  r0 <- list(hr.itt = 72/120, hr.H.hat = 27/60, hr.Hc.hat = 45/60,
             hr.H.true = 36/60, hr.Hc.true = 36/60)
  ird <- Map(`-`, r1, r0)
  irr <- Map(`/`, r1, r0)
  for (eng in c("FS", "GRF")) {
    f <- if (eng == "FS") .es_fs else .es_grf
    .es_check(f(df, "IRD", offset_name = "time", outcome_type = "count"),
              ird, paste(eng, "IRD offset"))
    # the rates above are on the exposure convention the "IRR" value uses
    .es_check(f(df, "IRR", offset_name = "time", outcome_type = "count"),
              irr, paste(eng, "IRR offset"))
  }
})

test_that("count IRD without an offset: difference of the fitted arm mean counts", {
  # No exposure: the fitted arm rate is the mean count per patient
  # (each cell has 20 patients, so a set's arm has 40 or 80 patients):
  #   ITT        r1 = 120/80 = 1.5       r0 = 72/80 = 0.9      IRD = 0.6
  #   Hhat = 1   r1 = 72/40 = 1.8        r0 = 27/40 = 0.675    IRD = 1.125
  #   Hhat = 0   r1 = 48/40 = 1.2        r0 = 45/40 = 1.125    IRD = 0.075
  #   H = 1      r1 = 80/40 = 2          r0 = 36/40 = 0.9      IRD = 1.1
  #   H = 0      r1 = 40/40 = 1          r0 = 36/40 = 0.9      IRD = 0.1
  df <- .es_data("count")
  ird <- list(hr.itt = 120/80 - 72/80, hr.H.hat = 72/40 - 27/40,
              hr.Hc.hat = 48/40 - 45/40, hr.H.true = 80/40 - 36/40,
              hr.Hc.true = 40/40 - 36/40)
  .es_check(.es_fs(df, "IRD", outcome_type = "count"),  ird, "FS IRD")
  .es_check(.es_grf(df, "IRD", outcome_type = "count"), ird, "GRF IRD")
})

test_that("a set with one arm only keeps the failure value NA", {
  # Treated patients only: the treatment coefficient is not estimable, so
  # no arm difference or ratio exists; the extractors return NA_real_.
  df <- .es_data("binary")
  df1 <- df[df$treat == 1, ]
  for (em in c("RD", "RR")) {
    expect_identical(.es_fs(df1, em)$hr.itt, NA_real_, label = paste("FS", em))
    expect_identical(.es_grf(df1, em)$hr.itt, NA_real_, label = paste("GRF", em))
  }
})
