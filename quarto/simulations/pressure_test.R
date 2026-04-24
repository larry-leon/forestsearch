# =============================================================================
# Pressure test for the two new/modified package files
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

#source("plot_sg_distribution.R")
#source("plot_effect_distribution.R")

PASS <- 0L
FAIL <- 0L

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
    cat(sprintf("  [FAIL] %s (got: %s)\n", label, deparse(res, width.cutoff = 60)[1]))
  }
}

expect_error <- function(label, expr, pattern = NULL) {
  res <- tryCatch(expr, error = function(e) e, warning = function(w) w)
  if (!inherits(res, "error")) {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s (expected error, got none)\n", label))
    return(invisible(NULL))
  }
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(res))) {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s (wrong error: %s)\n", label, conditionMessage(res)))
    return(invisible(NULL))
  }
  PASS <<- PASS + 1L
  cat(sprintf("  [ OK ] %s\n", label))
}

set.seed(42)

# =============================================================================
# 1. plot_sg_distribution() — schema flexibility + core behaviours
# =============================================================================
cat("\n=== plot_sg_distribution() ===\n\n")

# --- 1a. MRCT schema (legacy, any_found / sg_found) -------------------------
# 100 reps, 70 found, three dominant labels + a long tail
sg_labels <- c(rep("{er <= 0}",          30),
               rep("{size <= 35}",       20),
               rep("{er <= 0 & size <= 35}", 10),
               paste0("{noise", 1:10, "}"))      # 10 distinct singletons
sg_labels <- c(sg_labels, rep("", 30))            # 30 non-detections
stopifnot(length(sg_labels) == 100)
mrct_df <- data.frame(
  any_found = c(rep(1L, 70), rep(0L, 30)),
  sg_found  = sg_labels,
  stringsAsFactors = FALSE
)

p1 <- plot_sg_distribution(mrct_df, min_pct = 5)
ok("1a.1 legacy schema returns ggplot",         inherits(p1, "ggplot"))
ok("1a.2 n_bars attribute set",                 attr(p1, "n_bars") == 3L)
ok("1a.3 effective_min_pct preserved",          attr(p1, "effective_min_pct") == 5)

# --- 1b. run_simulation_analysis schema (any.H / sg.def) -------------------
fs_rows <- data.frame(
  any.H  = c(rep(1L, 70), rep(0L, 30)),
  sg.def = c(rep("{age <= 34}",     25),
             rep("{preanti <= 745}",20),
             rep("{age <= 34 & preanti <= 745}", 10),
             paste0("{noise", 1:15, "}"),
             rep(NA_character_, 30)),
  stringsAsFactors = FALSE
)
p2 <- plot_sg_distribution(fs_rows,
                           any_col   = "any.H",
                           label_col = "sg.def",
                           min_pct   = 5)
ok("1b.1 run_sim schema returns ggplot",        inherits(p2, "ggplot"))
ok("1b.2 n_bars = 3 (dominants only)",          attr(p2, "n_bars") == 3L)

# --- 1c. top-K cap + show_other -------------------------------------------
# Force a scenario with 10 labels each at >= min_pct
many_labels <- c(rep(paste0("sg", 1:10), each = 5))      # 50 finds, 10 labels
df_many <- data.frame(
  any.H  = rep(1L, 50),
  sg.def = many_labels,
  stringsAsFactors = FALSE
)
p3 <- plot_sg_distribution(df_many, any_col = "any.H", label_col = "sg.def",
                           min_pct = 5, top_k = 4L, show_other = TRUE)
ok("1c.1 top_k cap: 4 top + 1 pooled = 5 bars", attr(p3, "n_bars") == 5L)

p4 <- plot_sg_distribution(df_many, any_col = "any.H", label_col = "sg.def",
                           min_pct = 5, top_k = 4L, show_other = FALSE)
ok("1c.2 show_other = FALSE: only 4 bars",      attr(p4, "n_bars") == 4L)

# --- 1d. auto-halving -----------------------------------------------------
# Force a scenario where no label clears the initial min_pct but one clears
# the halved threshold.  70 finds, max single-label count = 5 (= 7.1%)
sparse_labels <- c(rep("{rare_a}", 5),
                   rep("{rare_b}", 5),
                   paste0("uniq", 1:60))
sparse_df <- data.frame(
  any.H  = rep(1L, 70),
  sg.def = sparse_labels,
  stringsAsFactors = FALSE
)
# min_pct = 20% initial — 7.1% < 20; halve to 10% — 7.1% < 10; halve to 5% — passes
p5 <- plot_sg_distribution(sparse_df, any_col = "any.H", label_col = "sg.def",
                           min_pct = 20, max_halvings = 2L,
                           placeholder_on_empty = TRUE)
ok("1d.1 auto-halving returns ggplot",          inherits(p5, "ggplot"))
ok("1d.2 auto-halved threshold = 5",            attr(p5, "effective_min_pct") == 5)

# Same data but max_halvings = 0 (original behaviour) — should hit placeholder
p5b <- plot_sg_distribution(sparse_df, any_col = "any.H", label_col = "sg.def",
                            min_pct = 20, max_halvings = 0L,
                            placeholder_on_empty = TRUE)
ok("1d.3 max_halvings=0 + empty -> placeholder", attr(p5b, "n_bars") == 0L)

# --- 1e. empty input behaviour --------------------------------------------
empty_df <- data.frame(any.H = rep(0L, 50), sg.def = rep("", 50),
                       stringsAsFactors = FALSE)
p6 <- plot_sg_distribution(empty_df, any_col = "any.H", label_col = "sg.def",
                           placeholder_on_empty = TRUE)
ok("1e.1 empty + placeholder=T returns ggplot", inherits(p6, "ggplot"))
ok("1e.2 empty placeholder has n_bars = 0",     attr(p6, "n_bars") == 0L)

p7 <- suppressMessages(
  plot_sg_distribution(empty_df, any_col = "any.H", label_col = "sg.def",
                       placeholder_on_empty = FALSE))
ok("1e.3 empty + placeholder=F returns NULL",   is.null(p7))

# --- 1f. missing-columns errors -------------------------------------------
expect_error("1f.1 missing any_col errors",
             plot_sg_distribution(fs_rows, any_col = "nope"),
             pattern = "'nope' not found")
expect_error("1f.2 missing label_col errors",
             plot_sg_distribution(fs_rows, any_col = "any.H", label_col = "nope"),
             pattern = "'nope' not found")

# --- 1g. backward-compat invocation (no new args) -------------------------
p8 <- plot_sg_distribution(mrct_df, min_pct = 5)
ok("1g.1 legacy call still works",              inherits(p8, "ggplot"))
ok("1g.2 legacy call has n_bars attr",          !is.null(attr(p8, "n_bars")))

# =============================================================================
# 2. sgdist_fig_height()
# =============================================================================
cat("\n=== sgdist_fig_height() ===\n\n")

ok("2.1 n_bars = 0 -> floor_in",                sgdist_fig_height(0L)  == 4.0)
ok("2.2 n_bars = 1 -> floor_in",                sgdist_fig_height(1L)  == 4.0)  # 3 + 0.35 = 3.35 < 4
ok("2.3 n_bars = 5 -> 4.75",                    abs(sgdist_fig_height(5L) - 4.75) < 1e-9)
ok("2.4 n_bars = 50 -> ceil_in",                sgdist_fig_height(50L) == 10.0)
ok("2.5 n_bars = NA defaults safely",           sgdist_fig_height(NA_integer_) == 4.0)
ok("2.6 custom per_bar works",
   abs(sgdist_fig_height(10L, per_bar = 0.5, base_in = 2.0) - 7.0) < 1e-9)

# =============================================================================
# 3. plot_effect_distribution()
# =============================================================================
cat("\n=== plot_effect_distribution() ===\n\n")

# --- 3a. Build a realistic run_simulation_analysis() output ---------------
# 100 sims, detection rate 60%, OR(ITT) near 0.85, OR(H) near 2.0 switched
n <- 100L
detected <- c(rep(1L, 60), rep(0L, 40))
results_alt <- data.frame(
  analysis  = "FS",
  any.H     = detected,
  hr.itt    = rlnorm(n, meanlog = log(0.85), sdlog = 0.15),
  hr.H.hat  = ifelse(detected == 1L,
                     rlnorm(n, meanlog = log(2.0), sdlog = 0.25),
                     NA_real_),
  hr.Hc.hat = ifelse(detected == 1L,
                     rlnorm(n, meanlog = log(0.60), sdlog = 0.15),
                     NA_real_),
  stringsAsFactors = FALSE
)

# Also a GRF-labelled row set (to test analysis_method slicing)
results_grf <- results_alt
results_grf$analysis <- "GRF"
results_all <- rbind(results_alt, results_grf)

# --- 3b. Benefit-scale OR inversion ---------------------------------------
p_or <- plot_effect_distribution(
  results           = results_alt,
  analysis_method   = "FS",
  effect_measure    = "OR",
  subgroup_notation = "benefit",
  reference_effect  = 2.0,   # benefit-scale truth
  title             = "Test: OR benefit"
)
ok("3b.1 OR benefit returns ggplot",            inherits(p_or, "ggplot"))
panel_data <- attr(p_or, "panel_data")
ok("3b.2 panel_data has 3 groups",              length(unique(panel_data$group)) == 3L)
# On benefit scale, ITT estimates 1/hr.itt should be > 1 (since hr.itt < 1)
itt_vals <- panel_data$est[panel_data$group == "ITT"]
ok("3b.3 ITT inverted correctly (mean > 1)",    mean(itt_vals, na.rm = TRUE) > 1)
# H-hat estimates 1/hr.H.hat where hr.H.hat ~ 2.0 switched => ~0.5 benefit
H_vals <- panel_data$est[grepl("^Identified H-hat", as.character(panel_data$group))]
ok("3b.4 H-hat inverted (mean ~ 0.5)",
   abs(mean(H_vals, na.rm = TRUE) - 0.5) < 0.1)

# --- 3c. Harm-scale (no inversion) -----------------------------------------
p_hr <- plot_effect_distribution(
  results           = results_alt,
  effect_measure    = "HR",
  subgroup_notation = "harm"
)
ok("3c.1 HR harm returns ggplot",               inherits(p_hr, "ggplot"))
panel_hr <- attr(p_hr, "panel_data")
itt_hr <- panel_hr$est[panel_hr$group == "ITT"]
# harm notation = identity transform, so should match original hr.itt values
ok("3c.2 HR harm identity (mean < 1)",          mean(itt_hr, na.rm = TRUE) < 1)

# --- 3d. RD (additive) flip -----------------------------------------------
results_rd <- results_alt
results_rd$hr.itt    <- rnorm(n, mean = -0.05, sd = 0.02)   # RD: control better
results_rd$hr.H.hat  <- ifelse(detected == 1L, rnorm(n, mean = -0.15, sd = 0.03), NA)
results_rd$hr.Hc.hat <- ifelse(detected == 1L, rnorm(n, mean =  0.00, sd = 0.02), NA)

p_rd_benefit <- plot_effect_distribution(
  results           = results_rd,
  effect_measure    = "RD",
  subgroup_notation = "benefit"
)
panel_rd <- attr(p_rd_benefit, "panel_data")
itt_rd <- panel_rd$est[panel_rd$group == "ITT"]
# Original ITT is -0.05, benefit flip -> +0.05
ok("3d.1 RD benefit: sign flipped (ITT > 0)",
   abs(mean(itt_rd, na.rm = TRUE) - 0.05) < 0.02)

# Same data, harm notation: no flip
p_rd_harm <- plot_effect_distribution(results_rd, effect_measure = "RD",
                                      subgroup_notation = "harm")
panel_rd_h <- attr(p_rd_harm, "panel_data")
itt_rd_h <- panel_rd_h$est[panel_rd_h$group == "ITT"]
ok("3d.2 RD harm: no flip (ITT < 0)",
   abs(mean(itt_rd_h, na.rm = TRUE) - (-0.05)) < 0.02)

# --- 3e. analysis_method slicing + error handling -------------------------
p_grf <- plot_effect_distribution(results_all, analysis_method = "GRF")
ok("3e.1 slicing GRF works",                    inherits(p_grf, "ggplot"))

expect_error("3e.2 unknown method errors",
             plot_effect_distribution(results_all, analysis_method = "NOPE"),
             pattern = "No rows with analysis")

# Missing required columns
results_bad <- results_alt[, c("analysis", "any.H")]  # drop effect cols
expect_error("3e.3 missing hr.* cols errors",
             plot_effect_distribution(results_bad),
             pattern = "missing.*hr")

# --- 3f. drop_undetected behaviour ----------------------------------------
# With drop_undetected = FALSE the H/Hc panels should include NAs for the
# non-detecting replicates; the group counts should reflect all 100 rows.
p_keep <- plot_effect_distribution(results_alt, drop_undetected = FALSE)
panel_keep <- attr(p_keep, "panel_data")
H_count_keep <- sum(panel_keep$group == "Identified H-hat")
ok("3f.1 drop_undetected = FALSE keeps all rows", H_count_keep == 100L)

p_drop <- plot_effect_distribution(results_alt, drop_undetected = TRUE)
panel_drop <- attr(p_drop, "panel_data")
H_count_drop <- sum(panel_drop$group == "Identified H-hat")
ok("3f.2 drop_undetected = TRUE drops 40 rows",  H_count_drop == 60L)
# ITT panel always keeps all 100
ok("3f.3 ITT panel always has all N",
   sum(panel_drop$group == "ITT") == 100L)

# --- 3g. IRR (ratio-family sanity) ----------------------------------------
p_irr <- plot_effect_distribution(results_alt, effect_measure = "IRR",
                                  subgroup_notation = "benefit")
ok("3g.1 IRR benefit returns ggplot",           inherits(p_irr, "ggplot"))

# =============================================================================
# 4. End-to-end vignette smoke test
# =============================================================================
cat("\n=== Vignette-chunk smoke test ===\n\n")

# Simulate what the vignette sets up in its setup chunk
sgdist_min_pct_alt     <- 0
sgdist_min_pct_null    <- 0
sgdist_top_k           <- 15L
sgdist_show_other_alt  <- FALSE
sgdist_show_other_null <- FALSE
sim_n_sample           <- 1000L
nsims_alt              <- 50L
nsims_null             <- 200L
sg_notation            <- "benefit"

# Fake results_alt and results_null exactly matching run_simulation_analysis
# schema (including leading FS + GRF rows and the full column superset)
build_fake_results <- function(nsim, p_detect, true_H_hr = 2.0) {
  n_methods <- 2L  # FS + GRF
  n_rows <- nsim * n_methods
  det <- rbinom(n_rows, 1L, p_detect)
  data.frame(
    analysis  = rep(c("FS", "GRF"), each = nsim),
    any.H     = det,
    hr.itt    = rlnorm(n_rows, log(0.85), 0.15),
    hr.H.hat  = ifelse(det == 1L, rlnorm(n_rows, log(true_H_hr), 0.25), NA_real_),
    hr.Hc.hat = ifelse(det == 1L, rlnorm(n_rows, log(0.60), 0.15), NA_real_),
    sg.def    = ifelse(det == 1L,
                       sample(c("{age <= 34}", "{preanti <= 745}",
                                "{age <= 34 & preanti <= 745}",
                                paste0("{noise", 1:10, "}")),
                              n_rows, replace = TRUE,
                              prob = c(0.25, 0.20, 0.10, rep(0.045, 10))),
                       ""),
    stringsAsFactors = FALSE
  )
}
results_alt  <- build_fake_results(nsims_alt,  p_detect = 0.6, true_H_hr = 2.0)
results_null <- build_fake_results(nsims_null, p_detect = 0.04, true_H_hr = 1.0)

# Fake calibrated DGM for the reference_effect argument
dgm_calibrated <- list(hazard_ratios = list(harm_subgroup = 2.0))

# === Chunk: oc-sgdist-build (from the simplified vignette) ================
plot_sgdist_alt <- plot_sg_distribution(
  results              = results_alt[results_alt$analysis == "FS", ],
  any_col              = "any.H",
  label_col            = "sg.def",
  min_pct              = sgdist_min_pct_alt,
  top_k                = sgdist_top_k,
  show_other           = sgdist_show_other_alt,
  max_halvings         = 2L,
  bar_label_inside     = TRUE,
  placeholder_on_empty = TRUE,
  title                = "HTE: Subgroups Identified (H-hat)",
  wrap_width           = 30
)

plot_sgdist_null <- plot_sg_distribution(
  results              = results_null[results_null$analysis == "FS", ],
  any_col              = "any.H",
  label_col            = "sg.def",
  min_pct              = sgdist_min_pct_null,
  top_k                = sgdist_top_k,
  show_other           = sgdist_show_other_null,
  max_halvings         = 2L,
  bar_label_inside     = TRUE,
  placeholder_on_empty = TRUE,
  title                = "Null: Subgroups Identified (H-hat)",
  wrap_width           = 30
)

h_alt  <- sgdist_fig_height(attr(plot_sgdist_alt,  "n_bars"))
h_null <- sgdist_fig_height(attr(plot_sgdist_null, "n_bars"))

ok("4.1 vignette chunk oc-sgdist-build: alt is ggplot",  inherits(plot_sgdist_alt,  "ggplot"))
ok("4.2 vignette chunk oc-sgdist-build: null is ggplot", inherits(plot_sgdist_null, "ggplot"))
ok("4.3 h_alt in [4, 10]",  h_alt >= 4 && h_alt <= 10)
ok("4.4 h_null in [4, 10]", h_null >= 4 && h_null <= 10)

# === Chunk: oc-alt-ordist + oc-null-ordist (from the simplified vignette) ==
or_Q_true <- 1 / dgm_calibrated$hazard_ratios$harm_subgroup

p_alt_or <- plot_effect_distribution(
  results           = results_alt,
  analysis_method   = "FS",
  effect_measure    = "OR",
  subgroup_notation = sg_notation,
  reference_effect  = or_Q_true,
  title             = "ACTG175 Binary (HTE): OR Estimates Across Simulations",
  subtitle          = sprintf(
    "n = %d, %d sims | dotted green = true OR(Q) on benefit scale = %.2f",
    sim_n_sample, nsims_alt, or_Q_true)
)

p_null_or <- plot_effect_distribution(
  results           = results_null,
  analysis_method   = "FS",
  effect_measure    = "OR",
  subgroup_notation = sg_notation,
  title             = "ACTG175 Binary (Null): OR Estimates Across Simulations",
  subtitle          = sprintf(
    "n = %d, %d sims | no true harm/benefit subgroup; dashed red = OR = 1",
    sim_n_sample, nsims_null)
)

ok("4.5 vignette chunk oc-alt-ordist:  returns ggplot",  inherits(p_alt_or,  "ggplot"))
ok("4.6 vignette chunk oc-null-ordist: returns ggplot",  inherits(p_null_or, "ggplot"))

# Actually *build* the plot grobs (ggplot errors sometimes surface only on build)
grob_tests <- list(
  "alt  sgdist"   = plot_sgdist_alt,
  "null sgdist"   = plot_sgdist_null,
  "alt  ordist"   = p_alt_or,
  "null ordist"   = p_null_or
)
for (nm in names(grob_tests)) {
  ok(sprintf("4.X grob builds: %s", nm),
     inherits(ggplot2::ggplot_gtable(ggplot2::ggplot_build(grob_tests[[nm]])),
              "gtable"))
}

# =============================================================================
cat(sprintf("\n=== TOTAL: %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0L) quit(save = "no", status = 1L)
