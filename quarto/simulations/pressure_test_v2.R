# =============================================================================
# Pressure test for the three newly-updated vignettes + MD extension
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(codetools)
})

source("plot_sg_distribution.R")
source("plot_effect_distribution.R")

PASS <- 0L; FAIL <- 0L
ok <- function(label, expr) {
  res <- tryCatch(expr, error = function(e) e)
  if (inherits(res, "error")) {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s\n         %s\n", label, conditionMessage(res)))
    return(invisible(NULL))
  }
  if (isTRUE(res)) {
    PASS <<- PASS + 1L
    cat(sprintf("  [ OK ] %s\n", label))
  } else {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s (got: %s)\n", label,
                deparse(res, width.cutoff = 60)[1]))
  }
}
expect_error <- function(label, expr, pattern = NULL) {
  res <- tryCatch(expr, error = function(e) e)
  if (!inherits(res, "error")) {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s (expected error, got none)\n", label)); return()
  }
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(res))) {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s (wrong error: %s)\n", label,
                conditionMessage(res))); return()
  }
  PASS <<- PASS + 1L
  cat(sprintf("  [ OK ] %s\n", label))
}

set.seed(2026)

# =============================================================================
# 1. MD extension regression tests
# =============================================================================
cat("\n=== MD extension ===\n\n")

n <- 80L
det <- rbinom(n, 1L, 0.55)
# ITT MD near -5 (control better on switched scale), H(switched) near +40
md_alt <- data.frame(
  analysis = "FS", any.H = det,
  hr.itt   = rnorm(n, mean = -5, sd = 2),
  hr.H.hat  = ifelse(det == 1L, rnorm(n, mean = 40, sd = 6), NA_real_),
  hr.Hc.hat = ifelse(det == 1L, rnorm(n, mean =  0, sd = 3), NA_real_)
)

# --- 1a. MD harm notation (no flip) ---------------------------------------
p_md_h <- plot_effect_distribution(md_alt, effect_measure = "MD",
                                   subgroup_notation = "harm")
ok("1a.1 MD harm returns ggplot", inherits(p_md_h, "ggplot"))
ok("1a.2 MD harm y-label",
   identical(p_md_h$labels$y, "Mean Difference"))
pd_h <- attr(p_md_h, "panel_data")
ok("1a.3 MD harm: ITT not flipped (mean near -5)",
   abs(mean(pd_h$est[pd_h$group == "ITT"], na.rm = TRUE) - (-5)) < 1)

# --- 1b. MD benefit notation (sign flip for additive measures) -------------
p_md_b <- plot_effect_distribution(md_alt, effect_measure = "MD",
                                   subgroup_notation = "benefit",
                                   reference_effect = -40)
ok("1b.1 MD benefit returns ggplot", inherits(p_md_b, "ggplot"))
ok("1b.2 MD benefit y-label",
   identical(p_md_b$labels$y, "Mean Difference (benefit scale)"))
pd_b <- attr(p_md_b, "panel_data")
ok("1b.3 MD benefit: ITT flipped (mean near +5)",
   abs(mean(pd_b$est[pd_b$group == "ITT"], na.rm = TRUE) - 5) < 1)
# Flip of H(switched = +40) -> benefit H-hat = -40, matches reference
H_mean_benefit <- mean(
  pd_b$est[grepl("^Identified H-hat", as.character(pd_b$group))],
  na.rm = TRUE)
ok("1b.4 MD benefit: H-hat flipped (mean near -40)",
   abs(H_mean_benefit - (-40)) < 3)

# --- 1c. match.arg still rejects unknown measures ---------------------------
expect_error("1c.1 unknown effect_measure errors",
             plot_effect_distribution(md_alt, effect_measure = "MAD"),
             pattern = "'arg' should be one of")

# =============================================================================
# 2. Vignette-parity smoke tests — all three vignettes, full chunks
# =============================================================================
cat("\n=== Vignette chunks: survival / continuous / poisson ===\n\n")

# Fake run_simulation_analysis() output factory for the three effect scales
build_fake_results <- function(nsim, p_detect, sg_labels, effect_type) {
  n_methods <- 2L  # FS + GRF
  n_rows <- nsim * n_methods
  det <- rbinom(n_rows, 1L, p_detect)
  scale_params <- switch(effect_type,
    "HR"  = list(itt = c(log(0.85), 0.15), H = c(log(2.0),  0.25),
                 Hc = c(log(0.60), 0.15), logscale = TRUE),
    "MD"  = list(itt = c(-5, 2),  H = c(40, 6),  Hc = c(0, 3),
                 logscale = FALSE),
    "IRR" = list(itt = c(log(1.0), 0.15), H = c(log(1.8), 0.2),
                 Hc = c(log(0.9), 0.15), logscale = TRUE)
  )
  draw <- function(pars, detected) {
    n <- length(detected)
    x <- if (scale_params$logscale)
      rlnorm(n, pars[1], pars[2]) else rnorm(n, pars[1], pars[2])
    ifelse(detected == 1L, x, NA_real_)
  }
  data.frame(
    analysis  = rep(c("FS", "GRF"), each = nsim),
    any.H     = det,
    hr.itt    = if (scale_params$logscale)
      rlnorm(n_rows, scale_params$itt[1], scale_params$itt[2])
      else rnorm(n_rows, scale_params$itt[1], scale_params$itt[2]),
    hr.H.hat  = draw(scale_params$H,  det),
    hr.Hc.hat = draw(scale_params$Hc, det),
    sg.def    = ifelse(det == 1L,
      sample(sg_labels, n_rows, replace = TRUE,
             prob = c(rep(0.2, min(3, length(sg_labels))),
                      rep(0.02 / max(1L, length(sg_labels) - 3L),
                          max(0L, length(sg_labels) - 3L)))[seq_along(sg_labels)]),
      ""),
    stringsAsFactors = FALSE
  )
}

# Common setup-chunk variables
sgdist_min_pct_alt     <- 5
sgdist_min_pct_null    <- 5
sgdist_top_k           <- 15L
sgdist_show_other_alt  <- FALSE
sgdist_show_other_null <- FALSE
sim_n_sample           <- 1000L

# ---------------------------------------------------------------------------
# 2A. SURVIVAL VIGNETTE (HR, benefit search)
# ---------------------------------------------------------------------------
cat("\n-- survival vignette --\n\n")

surv_labels <- c("{z1 = 1}", "{z2 = 1}", "{z1 = 1 & z2 = 1}",
                 paste0("{noise_z", 3:12, "}"))
results_alt  <- build_fake_results(50L, 0.60, surv_labels, "HR")
results_null <- build_fake_results(100L, 0.04,
                                   c(surv_labels, paste0("{rare", 1:10, "}")),
                                   "HR")
dgm_calibrated <- list(hazard_ratios = list(harm_subgroup = 2.0))
nsims_alt <- 50L; nsims_null <- 100L
sg_notation <- "benefit"

# Chunk: oc-sgdist-build
plot_sgdist_alt <- plot_sg_distribution(
  results_alt[results_alt$analysis == "FS", ],
  any_col = "any.H", label_col = "sg.def",
  min_pct = sgdist_min_pct_alt, top_k = sgdist_top_k,
  show_other = sgdist_show_other_alt, max_halvings = 2L,
  bar_label_inside = TRUE, placeholder_on_empty = TRUE,
  title = "HTE: Subgroups Identified (H-hat)", wrap_width = 30)
plot_sgdist_null <- plot_sg_distribution(
  results_null[results_null$analysis == "FS", ],
  any_col = "any.H", label_col = "sg.def",
  min_pct = sgdist_min_pct_null, top_k = sgdist_top_k,
  show_other = sgdist_show_other_null, max_halvings = 2L,
  bar_label_inside = TRUE, placeholder_on_empty = TRUE,
  title = "Null: Subgroups Identified (H-hat)", wrap_width = 30)
h_alt  <- sgdist_fig_height(attr(plot_sgdist_alt,  "n_bars"))
h_null <- sgdist_fig_height(attr(plot_sgdist_null, "n_bars"))

ok("2A.1 survival sgdist H1 is ggplot",  inherits(plot_sgdist_alt,  "ggplot"))
ok("2A.2 survival sgdist H0 is ggplot",  inherits(plot_sgdist_null, "ggplot"))
ok("2A.3 survival h_alt in [4, 10]",     h_alt  >= 4 && h_alt  <= 10)
ok("2A.4 survival h_null in [4, 10]",    h_null >= 4 && h_null <= 10)

# Chunk: oc-alt-hrdist + oc-null-hrdist
hr_Q_true <- 1 / dgm_calibrated$hazard_ratios$harm_subgroup
p_alt_hr  <- plot_effect_distribution(
  results_alt, "FS", effect_measure = "HR",
  subgroup_notation = sg_notation, reference_effect = hr_Q_true,
  title = "ACTG175 Survival (HTE)",
  subtitle = sprintf("n=%d, %d sims | HR(Q) = %.2f",
                     sim_n_sample, nsims_alt, hr_Q_true))
p_null_hr <- plot_effect_distribution(
  results_null, "FS", effect_measure = "HR",
  subgroup_notation = sg_notation,
  title = "ACTG175 Survival (Null)",
  subtitle = sprintf("n=%d, %d sims | null", sim_n_sample, nsims_null))

ok("2A.5 survival hrdist H1 is ggplot", inherits(p_alt_hr,  "ggplot"))
ok("2A.6 survival hrdist H0 is ggplot", inherits(p_null_hr, "ggplot"))

# ---------------------------------------------------------------------------
# 2B. CONTINUOUS VIGNETTE (MD, benefit search)
# ---------------------------------------------------------------------------
cat("\n-- continuous vignette --\n\n")

cont_labels <- c("{z1 = 1}", "{z2 = 1}", "{z1 = 1 & z2 = 1}",
                 paste0("{cd4_q", 1:4, "}"))
results_alt_c  <- build_fake_results(50L, 0.55, cont_labels, "MD")
results_null_c <- build_fake_results(100L, 0.04,
                                     c(cont_labels, paste0("{n", 1:8, "}")),
                                     "MD")
dgm_calibrated_c <- list(hazard_ratios = list(harm_subgroup = 40))  # switched MD

# sgdist chunks (identical structure)
pA <- plot_sg_distribution(
  results_alt_c[results_alt_c$analysis == "FS", ],
  any_col = "any.H", label_col = "sg.def",
  min_pct = sgdist_min_pct_alt, top_k = sgdist_top_k,
  show_other = sgdist_show_other_alt, max_halvings = 2L,
  bar_label_inside = TRUE, placeholder_on_empty = TRUE,
  title = "HTE: Subgroups Identified (H-hat)", wrap_width = 30)
pN <- plot_sg_distribution(
  results_null_c[results_null_c$analysis == "FS", ],
  any_col = "any.H", label_col = "sg.def",
  min_pct = sgdist_min_pct_null, top_k = sgdist_top_k,
  show_other = sgdist_show_other_null, max_halvings = 2L,
  bar_label_inside = TRUE, placeholder_on_empty = TRUE,
  title = "Null: Subgroups Identified (H-hat)", wrap_width = 30)
ok("2B.1 continuous sgdist H1 is ggplot", inherits(pA, "ggplot"))
ok("2B.2 continuous sgdist H0 is ggplot", inherits(pN, "ggplot"))

# md chunks — reference effect = -dgm$hazard_ratios$harm_subgroup
md_Q_true <- -dgm_calibrated_c$hazard_ratios$harm_subgroup
p_alt_md <- plot_effect_distribution(
  results_alt_c, "FS", effect_measure = "MD",
  subgroup_notation = "benefit", reference_effect = md_Q_true,
  title = "ACTG175 Continuous (HTE)",
  subtitle = sprintf("MD(Q) = %.1f", md_Q_true))
p_null_md <- plot_effect_distribution(
  results_null_c, "FS", effect_measure = "MD",
  subgroup_notation = "benefit",
  title = "ACTG175 Continuous (Null)")

ok("2B.3 continuous mddist H1 is ggplot", inherits(p_alt_md,  "ggplot"))
ok("2B.4 continuous mddist H0 is ggplot", inherits(p_null_md, "ggplot"))
# Critical check: ref line overlay actually shows MD(Q) on benefit scale
# (i.e. -40 from flipped +40); the line is drawn via geom_hline
pd_alt_md <- attr(p_alt_md, "panel_data")
H_md_benefit <- mean(
  pd_alt_md$est[grepl("^Identified H-hat", as.character(pd_alt_md$group))],
  na.rm = TRUE)
ok("2B.5 continuous MD benefit: H-hat flipped (mean near -40)",
   abs(H_md_benefit - (-40)) < 5)

# ---------------------------------------------------------------------------
# 2C. POISSON VIGNETTE (IRR, harm search — no inversion)
# ---------------------------------------------------------------------------
cat("\n-- poisson vignette --\n\n")

pois_labels <- c("{er_low = 1}", "{age_high = 1}", "{size_large = 1}",
                 paste0("{covar", 1:8, "}"))
results_alt_p  <- build_fake_results(100L, 0.50, pois_labels, "IRR")
results_null_p <- build_fake_results(200L, 0.04,
                                     c(pois_labels, paste0("{n", 1:12, "}")),
                                     "IRR")
dgm_calibrated_p <- list(hazard_ratios = list(harm_subgroup = 1.8))

pA_p <- plot_sg_distribution(
  results_alt_p[results_alt_p$analysis == "FS", ],
  any_col = "any.H", label_col = "sg.def",
  min_pct = sgdist_min_pct_alt, top_k = sgdist_top_k,
  show_other = sgdist_show_other_alt, max_halvings = 2L,
  bar_label_inside = TRUE, placeholder_on_empty = TRUE,
  title = "HTE: Subgroups Identified (H-hat)", wrap_width = 30)
pN_p <- plot_sg_distribution(
  results_null_p[results_null_p$analysis == "FS", ],
  any_col = "any.H", label_col = "sg.def",
  min_pct = sgdist_min_pct_null, top_k = sgdist_top_k,
  show_other = sgdist_show_other_null, max_halvings = 2L,
  bar_label_inside = TRUE, placeholder_on_empty = TRUE,
  title = "Null: Subgroups Identified (H-hat)", wrap_width = 30)
ok("2C.1 poisson sgdist H1 is ggplot", inherits(pA_p, "ggplot"))
ok("2C.2 poisson sgdist H0 is ggplot", inherits(pN_p, "ggplot"))

# IRR chunks — harm search, no inversion
irr_H_true <- dgm_calibrated_p$hazard_ratios$harm_subgroup
p_alt_irr <- plot_effect_distribution(
  results_alt_p, "FS", effect_measure = "IRR",
  subgroup_notation = "harm", reference_effect = irr_H_true,
  title = "GBSG Poisson (HTE)",
  subtitle = sprintf("IRR(H) = %.2f", irr_H_true))
p_null_irr <- plot_effect_distribution(
  results_null_p, "FS", effect_measure = "IRR",
  subgroup_notation = "harm",
  title = "GBSG Poisson (Null)")

ok("2C.3 poisson irrdist H1 is ggplot", inherits(p_alt_irr,  "ggplot"))
ok("2C.4 poisson irrdist H0 is ggplot", inherits(p_null_irr, "ggplot"))
# Harm search: H-hat should NOT be inverted — mean should match fake DGM (1.8)
pd_alt_irr <- attr(p_alt_irr, "panel_data")
H_irr_harm <- mean(
  pd_alt_irr$est[grepl("^Identified H-hat", as.character(pd_alt_irr$group))],
  na.rm = TRUE)
ok("2C.5 poisson IRR harm: H-hat NOT inverted (mean near 1.8)",
   abs(H_irr_harm - 1.8) < 0.3)
ok("2C.6 poisson axis label (no benefit scale)",
   identical(p_alt_irr$labels$y, "Incidence Rate Ratio"))

# =============================================================================
# 3. Grob materialisation — catches ggplot-layer issues the above miss
# =============================================================================
cat("\n=== Grob materialisation ===\n\n")

all_plots <- list(
  surv_sgdist_alt  = plot_sgdist_alt,
  surv_sgdist_null = plot_sgdist_null,
  surv_hr_alt      = p_alt_hr,
  surv_hr_null     = p_null_hr,
  cont_sgdist_alt  = pA,
  cont_sgdist_null = pN,
  cont_md_alt      = p_alt_md,
  cont_md_null     = p_null_md,
  pois_sgdist_alt  = pA_p,
  pois_sgdist_null = pN_p,
  pois_irr_alt     = p_alt_irr,
  pois_irr_null    = p_null_irr
)
for (nm in names(all_plots)) {
  ok(sprintf("3.x grob: %s", nm),
     inherits(ggplot_gtable(ggplot_build(all_plots[[nm]])), "gtable"))
}

# =============================================================================
# 4. codetools::checkUsage() on updated plot_effect_distribution
# =============================================================================
cat("\n=== codetools::checkUsage() on updated plot_effect_distribution ===\n\n")
checkUsage(plot_effect_distribution, all = TRUE)

# =============================================================================
cat(sprintf("\n=== TOTAL: %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0L) quit(save = "no", status = 1L)
