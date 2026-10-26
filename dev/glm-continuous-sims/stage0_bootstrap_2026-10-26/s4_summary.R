# s4_summary.R -- assembles results/stage0_summary.md from the key = value
# record written by s0-s3.  Pure arithmetic on measured numbers; every figure
# in the summary is traceable to a line in results/keyvalues.txt.

source("stage0_common.R")
lines <- readLines(out_path("keyvalues.txt"))
kvs <- strsplit(lines, " = ", fixed = TRUE)
vals <- setNames(lapply(kvs, function(z) paste(z[-1L], collapse = " = ")), vapply(kvs, `[[`, "", 1L))
# later writes of the same key win
vals <- vals[!duplicated(names(vals), fromLast = TRUE)]
num <- function(k) { v <- vals[[k]]; if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(strsplit(v, ", ")[[1L]][1L])) }
chr <- function(k) { v <- vals[[k]]; if (is.null(v)) "<missing>" else v }
pct <- function(x) if (is.na(x)) "n/a" else sprintf("%.1f%%", 100 * x)
s   <- function(x, d = 2) if (is.na(x)) "n/a" else formatC(x, format = "f", digits = d)

BAR <- 0.30

# ---- A ------------------------------------------------------------------------
rep_s     <- num("boot_S_mean_rep_s")          # sequential replicate seconds, measured now
call_s    <- num("call_s_median")              # single sequential call, measured now
mach_s    <- num("machinery_s_per_call")       # direct timers (s2)
mach_rp   <- num("rprof_formula_machinery_s")  # Rprof corroboration (s1)
A_share_rep  <- mach_s / rep_s
A_share_call <- mach_s / call_s
A_clears  <- isTRUE(A_share_rep >= BAR)

# ---- C / D --------------------------------------------------------------------
W <- num("boot_workers_W"); NB <- num("boot_nb_main")
g <- function(tag, f) num(sprintf("boot_%s_%s", tag, f))
startup_M <- g("M", "spawn_s") + g("M", "overhead_s")
startup_C <- g("C", "spawn_s") + g("C", "overhead_s")
# Projected application-shape wall at B replicates:
#   startup (plan-specific) + B * rep / W + static-chunk imbalance (plan-specific,
#   as measured at nb = 4W).  The replicate term is COMMON to both plans (mean of
#   the two loaded per-replicate means) so that the C comparison isolates startup
#   and is not swamped by run-to-run noise in the replicate mean.
rep_common <- mean(c(g("M", "mean_rep_s"), g("C", "mean_rep_s")))
proj <- function(tag, B, rep_mult = 1) {
  (g(tag, "spawn_s") + g(tag, "overhead_s")) + B * rep_common * rep_mult / W +
    (g(tag, "contig_s") - g(tag, "balanced_s"))
}
B <- 1000
wall_M_1000 <- proj("M", B); wall_C_1000 <- proj("C", B)
C_share_1000 <- 1 - wall_C_1000 / wall_M_1000
C_share_meas <- 1 - g("C", "outer_s") / g("M", "outer_s")
D_share_M <- (g("M", "contig_s") - g("M", "balanced_s")) / g("M", "outer_s")
# A inside the application shape: replicate seconds shrink by the machinery share
A_rep_mult <- 1 - A_share_rep
wall_M_1000_A  <- proj("M", B, A_rep_mult)
wall_C_1000_A  <- proj("C", B, A_rep_mult)

md <- c(
"# Stage 0 -- bootstrap efficiency measurements (CRAN forestsearch 0.4.0)",
"",
sprintf("Installed: forestsearch %s (%s; repository %s) · survival %s · future %s · %s · %s/%s · cores %s · multicore %s",
        chr("installed_version"), chr("installed_packaged"), chr("installed_repository"),
        chr("survival_version"), chr("future_version"), chr("r_version"), chr("sysname"), chr("machine"),
        chr("available_cores"), chr("supports_multicore")),
sprintf("Repo: %s @ %s (%s); files matching CRAN 0.4.0: %s, differing: %s",
        chr("repo_root"), chr("repo_head"), chr("repo_branch"),
        chr("repo_files_matching_cran"), chr("repo_files_differing")),
sprintf("Fixture: %s (n = %s)", chr("fixture"), chr("fixture_n")),
"",
"## 1. The single sequential call (fixture V)",
"",
"| quantity | value |",
"|---|---|",
sprintf("| call wall-clock, median of 3 (s) | %s (reps %s) |", s(call_s, 3), chr("call_s_rep")),
sprintf("| candidates evaluated in the search | %s |", chr("search_n_evaluated")),
sprintf("| candidates reaching the Cox fit | %s |", chr("search_n_reached_cox")),
sprintf("| candidates passing the effect screen | %s |", chr("search_n_passed_effect")),
sprintf("| consistency stage: candidates total / evaluated / passed | %s / %s / %s (early stop: %s) |",
        chr("consistency_n_total"), chr("consistency_n_evaluated"), chr("consistency_n_passed"), chr("consistency_early_stop")),
sprintf("| identified subgroup | %s |", chr("identified_subgroup")),
"",
"Rprof, median of 3 profiles (seconds of total time; see rprof_watched_median.csv for the full list):",
"",
"| function | total s |",
"|---|---|",
sprintf("| coxph (formula interface, all callers) | %s |", s(num("rprof_coxph_total_s"))),
sprintf("| coxph.fit (+ agreg.fit) -- the fit proper | %s |", s(num("rprof_coxph_fit_total_s"))),
sprintf("| formula-interface machinery = coxph - coxph.fit | %s (%s of coxph) |", s(mach_rp), pct(num("rprof_machinery_share_of_coxph"))),
sprintf("| of which concordancefit | %s |", s(num("rprof_concordancefit_s"))),
sprintf("| search_combinations_parallel (whole search stage) | %s |", s(num("rprof_search_stage_s"))),
sprintf("| fit_cox_for_subgroup (search-stage fits) | %s |", s(num("rprof_fit_cox_for_subgroup_s"))),
sprintf("| survfit (medians on effect-screen passers) | %s |", s(num("rprof_survfit_s"))),
sprintf("| evaluate_consistency_twostage (whole consistency evaluation) | %s |", s(num("rprof_consistency_stage_s"))),
sprintf("| .consistency_cox_pieces / residuals.coxph | %s / %s |", s(num("rprof_consistency_cox_pieces_s")), s(num("rprof_residuals_coxph_s"))),
sprintf("| coxph time outside the search stage (Rprof coxph total - direct-timed search fits; candidate B's territory) | %s |",
        s(num("rprof_coxph_total_s") - num("time_formula_s_median"))),
"",
"## 2. Candidate A -- lean search-stage Cox fit (coxph.fit instead of formula coxph)",
"",
"| quantity | value |",
"|---|---|",
sprintf("| captured search-stage fits (= candidates reaching the fit) | %s (matches filter counts: %s; adjusted: %s) |",
        chr("capture_n_fits"), chr("capture_matches_filter_counts"), chr("capture_n_adjusted")),
sprintf("| CI rows identical, lean vs formula (identical()) | %s of %s |", chr("ident_ci_row_identical"), chr("ident_n_candidates")),
sprintf("| CI values identical ignoring dimnames | %s |", chr("ident_ci_values_identical")),
sprintf("| failure sets agree (try-errors) | %s (formula %s, lean %s) |", chr("ident_error_sets_agree"), chr("ident_errors_formula"), chr("ident_errors_lean")),
sprintf("| warning sets agree | %s (candidates warned: %s; identical among them: %s) |", chr("ident_warn_sets_agree"), chr("ident_n_warned"), chr("ident_identical_among_warned")),
sprintf("| formula path, all captured fits, median of 5 (s) | %s (%s ms/fit) |", s(num("time_formula_s_median"), 3), s(num("time_ms_per_fit_formula"), 3)),
sprintf("| lean path, all captured fits, median of 5 (s) | %s (%s ms/fit) |", s(num("time_lean_s_median"), 3), s(num("time_ms_per_fit_lean"), 3)),
sprintf("| **machinery removed per call (s)** | **%s** (direct timers) vs %s (Rprof coxph - coxph.fit, all callers) |", s(mach_s, 3), s(mach_rp, 3)),
sprintf("| machinery share of the formula path | %s |", pct(num("machinery_share_of_formula"))),
sprintf("| sequential replicate now, mean of 5 (s) | %s |", s(rep_s, 3)),
sprintf("| **A: share of a sequential replicate** | **%s** (of the single call: %s) |", pct(A_share_rep), pct(A_share_call)),
sprintf("| **A clears the %s bar on the per-replicate (campaign) shape?** | **%s** |", pct(BAR), if (A_clears) "YES" else "NO"),
"",
"Captured coxph.fit() conventions (the lean recipe, from the installed survival):",
"",
sprintf("- coxph.fit formals: `%s`", chr("conv_coxph_fit_formals")),
sprintf("- x: storage `%s`, colnames `%s`, values = treatment column: %s", chr("conv_x_storage"), chr("conv_x_colnames"), chr("conv_x_values_equal_treat")),
sprintf("- y: type `%s`; equals Surv(Y,E): %s; equals aeqSurv(Surv(Y,E)): %s; control$timefix = %s",
        chr("conv_y_type"), chr("conv_y_values_equal_raw"), chr("conv_y_values_equal_aeqSurv"), chr("conv_timefix")),
sprintf("- strata NULL: %s · offset: %s · init: %s · control = coxph.control(): %s · weights NULL: %s · method: %s · rownames: %s · resid: %s · nocenter: %s",
        chr("conv_strata_is_null"), chr("conv_offset_desc"), chr("conv_init_desc"), chr("conv_control_is_default"),
        chr("conv_weights_is_null"), chr("conv_method"), chr("conv_rownames_desc"), chr("conv_resid"), chr("conv_nocenter")),
"",
"## 3. Candidates C and D -- application shape (bootstrap at W workers)",
"",
sprintf("W = %s workers, nb = %s for the M/C runs; identity pair at nb = 20.", chr("boot_workers_W"), chr("boot_nb_main")),
"",
"| run | plan | nb | outer s | inner s | spawn s | overhead s | sum rep s | mean rep s | contiguous makespan s | balanced s | imbalance s |",
"|---|---|---|---|---|---|---|---|---|---|---|---|",
vapply(c("S", "M", "C", "Mi", "Ci"), function(t) sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
        t, chr(sprintf("boot_%s_plan", t)), chr(sprintf("boot_%s_nb", t)), s(g(t, "outer_s"), 1), s(g(t, "inner_s"), 1),
        s(g(t, "spawn_s"), 1), s(g(t, "overhead_s"), 1), s(g(t, "sum_rep_s"), 1), s(g(t, "mean_rep_s"), 2),
        s(g(t, "contig_s"), 1), s(g(t, "balanced_s"), 1), s(g(t, "contig_s") - g(t, "balanced_s"), 1)), character(1)),
"",
"| quantity | value |",
"|---|---|",
sprintf("| identity across plans, nb = 20 (all non-volatile components) | %s (results table only: %s) |", chr("identity_multicore_vs_multisession_nb20"), chr("identity_results_table_only")),
sprintf("| per-replicate slowdown under load (M / sequential, C / sequential) | %s / %s |", s(num("rep_slowdown_M"), 2), s(num("rep_slowdown_C"), 2)),
sprintf("| startup (spawn + load): multisession vs multicore (s) | %s vs %s; removed: %s |", s(startup_M, 1), s(startup_C, 1), s(num("startup_removed_by_multicore_s"), 1)),
sprintf("| **C: measured wall reduction at nb = %s** | **%s** |", chr("boot_nb_main"), pct(C_share_meas)),
sprintf("| **C: projected wall reduction at B = 1000** (startup + B·rep/W + imbalance) | **%s** (M %s s -> C %s s) |", pct(C_share_1000), s(wall_M_1000, 0), s(wall_C_1000, 0)),
sprintf("| **D: static-chunk imbalance, share of the M run's wall** | **%s** (%s s of %s s) |", pct(D_share_M), s(g("M", "contig_s") - g("M", "balanced_s"), 1), s(g("M", "outer_s"), 1)),
sprintf("| A inside the application shape at B = 1000: multisession %s s -> %s s; multicore %s s -> %s s | A+C combined: %s |",
        s(wall_M_1000, 0), s(wall_M_1000_A, 0), s(wall_C_1000, 0), s(wall_C_1000_A, 0), pct(1 - wall_C_1000_A / wall_M_1000)),
"",
"## 4. Reading",
"",
sprintf("- Per-replicate (campaign) shape: A removes %s s of a %s s replicate = %s; bar %s: %s.",
        s(mach_s, 2), s(rep_s, 2), pct(A_share_rep), pct(BAR), if (A_clears) "cleared" else "not cleared"),
sprintf("- Bit-identity of the lean fit: %s of %s CI rows identical; failure sets agree: %s.",
        chr("ident_ci_row_identical"), chr("ident_n_candidates"), chr("ident_error_sets_agree")),
sprintf("- Application shape at B = 1000 (projection from measured components): C alone %s; D (imbalance) %s of the measured M wall; A+C %s.",
        pct(C_share_1000), pct(D_share_M), pct(1 - wall_C_1000_A / wall_M_1000)),
"- Projections carry the per-replicate slowdown under load as measured and the imbalance as measured at nb = 4W; they are not a substitute for a timed B = 1000 run.",
""
)
writeLines(md, out_path("stage0_summary.md"))
cat(md, sep = "\n")
cat("\ns4 complete\n")
