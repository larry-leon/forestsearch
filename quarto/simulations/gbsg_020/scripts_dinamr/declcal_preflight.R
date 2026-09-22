# declcal_preflight.R -- Block C design-feasibility firewall
# (TASK_declcal_CAMPAIGN_2026-09-22_v2, section 4).  Run BEFORE any Block C
# replicate.  For each Block C cell: the planted region's super-population
# prevalence under the template's alt DGM (z1_quantile 0.25, the template
# default), the size floor actually in force (the template passes n.min = NULL,
# which forestsearch() resolves to max(60, ceiling(0.10 n)) --
# R/forestsearch_main.R:1680-1689 -- so the rule is max, not the min the task
# text quotes "where that is the rule"), and the share of replicates whose
# realized region falls below that floor: exact binomial, and realized on the
# campaign's own draws (seed_base + rep, rep 1..2000).
#
# Materiality, fixed here before any number is computed: a cell is EXCLUDED
# when the realized share below the floor exceeds 0.05.
# Writes logs/declcal_preflight.txt; exit status 0 always (an exclusion is the
# firewall working); the excluded cells are listed on the line "EXCLUDE:".
suppressPackageStartupMessages(library(forestsearch))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
seed_base <- 8316951L; n_super <- 100000L; z1q <- 0.25
analysis_time <- 84; cens_adjust <- log(1.5)
materiality <- 0.05
cells <- data.frame(cell = c("C1", "C2", "C3", "C4"), hr = c(1.5, 1.5, 2.0, 2.0),
                    n = c(1000L, 1500L, 1000L, 1500L))
out <- character(0)
say <- function(...) { s <- sprintf(...); cat(s, "\n"); out <<- c(out, s) }
say("declcal Block C pre-flight -- %s -- forestsearch %s", format(Sys.time()),
    as.character(packageVersion("forestsearch")))
say("materiality: exclude a cell when the realized share of replicates with |region| < floor exceeds %.2f", materiality)
excl <- character(0)
for (h in unique(cells$hr)) {
  k_inter <- calibrate_k_inter(target_hr_harm = h, model = "alt", use_ahr = FALSE, z1_quantile = z1q)
  dgm <- setup_gbsg_dgm(model = "alt", k_inter = k_inter, k_treat = 1, z1_quantile = z1q,
                        n_super = n_super, seed = seed_base)
  p <- mean(dgm$df_super$flag_harm)
  say("HR %.2f: k_inter %.5f ; planted-region prevalence (super-population) %.4f ; hr_H_true %.4f",
      h, k_inter, p, dgm$hr_H_true %||% NA_real_)
  for (i in which(cells$hr == h)) {
    n <- cells$n[i]; fl <- max(60L, as.integer(ceiling(0.10 * n)))
    sz <- vapply(1:2000, function(r) sum(simulate_from_dgm(dgm, n = n, analysis_time = analysis_time,
                 cens_adjust = cens_adjust, seed = seed_base + r)$flag_harm == 1L), integer(1))
    p_bin <- stats::pbinom(fl - 1L, n, p)
    p_real <- mean(sz < fl)
    dec <- if (p_real > materiality) "EXCLUDE" else "RUN"
    if (dec == "EXCLUDE") excl <- c(excl, cells$cell[i])
    say("  %s (HR %.2f, n %d): floor %d ; E|region| %.1f (sd %.1f) ; P(|region| < floor) binomial %.4f, realized on reps 1-2000 %.4f (%d of 2000; min %d, median %d) -> %s",
        cells$cell[i], h, n, fl, n * p, sqrt(n * p * (1 - p)), p_bin, p_real, sum(sz < fl),
        min(sz), as.integer(median(sz)), dec)
  }
}
say("EXCLUDE: %s", if (length(excl)) paste(excl, collapse = " ") else "none")
dir.create("logs", showWarnings = FALSE)
writeLines(out, "logs/declcal_preflight.txt")
