# b3_summary.R -- assembles results/stage0b_summary.md from results/keyvalues.txt:
# per fixture the call and replicate times, the Rprof buckets, each candidate's
# measured machinery seconds and share of a replicate against the 30% bar, the
# identity counts, and the modification map (where package functions would
# change) with the measured status of each entry.  Pure arithmetic on
# recorded numbers.

source("stage0b_common.R")
lines <- readLines(out_path("keyvalues.txt"))
kvs <- strsplit(lines, " = ", fixed = TRUE)
vals <- setNames(lapply(kvs, function(z) paste(z[-1L], collapse = " = ")), vapply(kvs, `[[`, "", 1L))
vals <- vals[!duplicated(names(vals), fromLast = TRUE)]
num <- function(k) { v <- vals[[k]]; if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(strsplit(v, ", ")[[1L]][1L])) }
chr <- function(k) { v <- vals[[k]]; if (is.null(v)) "<missing>" else v }
pct <- function(x) if (is.na(x)) "n/a" else sprintf("%.1f%%", 100 * x)
s   <- function(x, d = 2) if (is.na(x)) "n/a" else formatC(x, format = "f", digits = d)
BAR <- 0.30

# Stage-0 (committed) survival A figures, quoted for the survival row
A_mach_s <- 1.509; A_rep_s <- 3.2112

md <- c("# Stage 0b -- GLM paths and the consistency stage (CRAN forestsearch 0.4.0)", "",
  sprintf("Installed: forestsearch %s (%s; repository %s) · survival %s · data.table %s · %s · cores %s",
          chr("installed_version"), chr("installed_packaged"), chr("installed_repository"), chr("survival_version"),
          chr("data_table_version"), chr("r_version"), chr("available_cores")),
  sprintf("Repo: %s @ %s (%s); files matching CRAN 0.4.0: %s, differing: %s",
          chr("repo_root"), chr("repo_head"), chr("repo_branch"), chr("repo_files_matching_cran"), chr("repo_files_differing")),
  "")

# ---- per-fixture profile table ----------------------------------------------------
md <- c(md, "## 1. The single sequential call and the replicate, per fixture", "",
  "| fixture | source | n | call s (median of 3) | sequential replicate s (nb = 5) | search: evaluated / reached fit / passed effect | consistency: total / evaluated / passed | early stop | identified subgroup |",
  "|---|---|---|---|---|---|---|---|---|")
for (t in names(FIXTURES)) md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s / %s / %s | %s / %s / %s | %s | %s |",
  t, chr(paste0("fixture_", t, "_source")), chr(paste0("fixture_", t, "_n")), s(num(paste0(t, "_call_s_median")), 3),
  s(num(paste0(t, "_rep_s_mean")), 3), chr(paste0(t, "_n_evaluated_search")), chr(paste0(t, "_n_reached_fit")),
  chr(paste0(t, "_n_passed_effect")), chr(paste0(t, "_n_consistency_total")), chr(paste0(t, "_n_consistency_evaluated")),
  chr(paste0(t, "_n_consistency_passed")), chr(paste0(t, "_early_stop_triggered")), chr(paste0(t, "_sg_harm"))))

md <- c(md, "", "Rprof buckets, median of 3 profiles, seconds of total time (a function reached through do.call() is not attributable by name; its children are):", "",
  "| fixture | call s | search stage | of which per-candidate fit | consistency stage | of which pieces (fit + dfbeta) | consistency stage minus pieces (extraction, bookkeeping) | dedup / sort | GRF | LASSO | survfit medians |",
  "|---|---|---|---|---|---|---|---|---|---|---|")
for (t in names(FIXTURES)) md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
  t, s(num(paste0(t, "_rprof_sampling_s"))), s(num(paste0(t, "_rprof_search_stage_s"))), s(num(paste0(t, "_rprof_search_fit_s"))),
  s(num(paste0(t, "_rprof_consistency_stage_s"))), s(num(paste0(t, "_rprof_consistency_pieces_s"))),
  s(num(paste0(t, "_rprof_consistency_stage_s")) - num(paste0(t, "_rprof_consistency_pieces_s"))),
  s(num(paste0(t, "_rprof_dedup_sort_s"))), s(num(paste0(t, "_rprof_grf_s"))), s(num(paste0(t, "_rprof_lasso_s"))), s(num(paste0(t, "_rprof_survfit_s")))))

# ---- candidates --------------------------------------------------------------------
md <- c(md, "", "## 2. Candidates: machinery removed per call, identity, and the 30% bar (per-replicate shape)", "",
  "Machinery = (package function on the recorded inputs) - (lean equivalent on the same inputs), direct timers, median of 3 alternating reps. Share = machinery / measured sequential replicate. Identity = identical() on what the package consumes.", "",
  "| candidate | fixture | what | n calls | package s | lean s | machinery s | share of replicate | identity | clears 30%? |",
  "|---|---|---|---|---|---|---|---|---|---|")
row <- function(letter, t, what, n, a, b, ident) {
  rep_s <- num(paste0(t, "_rep_s_mean")); mach <- a - b; sh <- mach / rep_s
  sprintf("| **%s** | %s | %s | %s | %s | %s | **%s** | **%s** | %s | %s |", letter, t, what, n, s(a, 3), s(b, 3), s(mach, 3), pct(sh), ident,
          if (is.na(sh)) "n/a" else if (sh >= BAR) "YES" else "no")
}
# survival A (Stage 0, quoted)
md <- c(md, sprintf("| **A** | survival | search-stage Cox fit (Stage 0, committed) | 1410 | 1.887 | 0.378 | **1.509** | **%s** | 1410/1410 | YES |", pct(A_mach_s / A_rep_s)))
for (t in names(FIXTURES)) {
  L <- chr(paste0(t, "_pieces_letter"))
  md <- c(md, row(L, t, "consistency pieces (fit + dfbeta), narrow frame",
                  chr(paste0(t, "_pieces_n")), num(paste0(t, "_pieces_formula_s")), num(paste0(t, "_pieces_lean_s")),
                  sprintf("rate identical %s/%s; NULL sets agree %s", chr(paste0(t, "_pieces_rate_identical")), chr(paste0(t, "_pieces_n")), chr(paste0(t, "_pieces_null_sets_agree")))))
  if (t == "binary") md <- c(md, row("F", t, "search-stage OR fit (glm.fit on vectors, SE from QR)",
                  chr("binary_F_n"), num("binary_F_formula_s"), num("binary_F_lean_s"),
                  sprintf("identical %s/%s; kinds agree %s", chr("binary_F_identical"), chr("binary_F_n"), chr("binary_F_kinds_agree"))))
  md <- c(md, row(if (t == "survival") "H (+B)" else if (t == "continuous") "H+E" else "H+G", t,
                  "whole consistency evaluator: lean extraction + lean pieces",
                  chr(paste0(t, "_stage_n")), num(paste0(t, "_stage_formula_s")), num(paste0(t, "_stage_lean_s")),
                  sprintf("Pcons identical: %s (passed %s vs %s; lean fallbacks %s)", chr(paste0(t, "_stage_pcons_identical_all")),
                          chr(paste0(t, "_stage_n_passed_pkg")), chr(paste0(t, "_stage_n_passed_lean")), chr(paste0(t, "_stage_n_lean_fallback")))))
  md <- c(md, row("H alone", t, if (t == "survival") "lean extraction + lean warm start, package pieces" else "lean extraction, package pieces",
                  chr(paste0(t, "_stage_n")), num(paste0(t, "_stage_formula_s")), num(paste0(t, "_stage_lean_Honly_s")),
                  sprintf("Pcons identical: %s", chr(paste0(t, "_stage_pcons_identical_H_only")))))
}
md <- c(md, "", sprintf("Continuous search-stage fit (already lean since 0.3.3): %s s per call over %s candidates -- no candidate.",
                        s(num("continuous_searchfit_fast_s"), 3), chr("continuous_searchfit_fast_n")))

# ---- combined per fixture ----------------------------------------------------------
md <- c(md, "", "## 3. Combined per fixture (per-replicate shape; machinery seconds are additive across stages)", "",
  "| fixture | replicate s | search-stage machinery s | consistency-stage machinery s (whole evaluator) | total s | share | clears 30%? |", "|---|---|---|---|---|---|---|")
for (t in names(FIXTURES)) {
  rep_s <- num(paste0(t, "_rep_s_mean"))
  srch <- if (t == "survival") A_mach_s else if (t == "binary") num("binary_F_machinery_s") else 0
  cons <- num(paste0(t, "_stage_machinery_s"))
  tot <- srch + cons; sh <- tot / rep_s
  md <- c(md, sprintf("| %s | %s | %s | %s | **%s** | **%s** | %s |", t, s(rep_s, 3), s(srch, 3), s(cons, 3), s(tot, 3), pct(sh),
                      if (is.na(sh)) "n/a" else if (sh >= BAR) "YES" else "no"))
}

md <- c(md, "", "The survival search-stage figure (A) is Stage 0's committed measurement on this machine (1.509 s of a 3.2112 s replicate); the survival replicate measured here is the denominator of the combined row.", "")

# ---- modification map --------------------------------------------------------------
st <- function(k) chr(k)
md <- c(md, "", "## 4. Modification map (CRAN 0.4.0 line references; nothing here was changed)", "",
  "| candidate | function(s) that would change | file:lines (0.4.0) | what replaces what | new helper | measured status |", "|---|---|---|---|---|---|",
  "| A | `fit_cox_for_subgroup()` unadjusted branch | R/subgroup_search.R:809-845 | formula `coxph()` -> `survival::coxph.fit()` with captured conventions; fallback to the formula block on error | `.fs_cox_lean_ci()` | Stage 0: 1410/1410 identical; 1.509 s/call |",
  sprintf("| B | `.consistency_cox_pieces()` (and the warm-start fits in `evaluate_subgroup_consistency()` / `evaluate_consistency_twostage()`) | R/consistency_resample.R:55-91; R/subgroup_consistency_helpers.R:1457-1467, 1827-1837 | two formula `coxph()` + `residuals(type = 'dfbeta')` with model-frame re-evaluation -> `coxph.fit()` x2 + `residuals.coxph()` on a constructed object | shares `.fs_cox_lean_ci()`'s fit | pieces rate identical %s/%s; machinery %s s/call |",
          st("survival_pieces_rate_identical"), st("survival_pieces_n"), s(num("survival_pieces_machinery_s"), 3)),
  sprintf("| E | `.consistency_glm_pieces()` MD branch + `.dfbeta_glm()` | R/consistency_resample.R:292, 168-176 | `lm()` + `lm.influence()` -> `lm.fit()` + hat values from the fit's QR | `.fs_glm_pieces_lean()` | pieces rate identical %s/%s; machinery %s s/call |",
          st("continuous_pieces_rate_identical"), st("continuous_pieces_n"), s(num("continuous_pieces_machinery_s"), 3)),
  sprintf("| F | routing predicate in `forestsearch()` + a fast OR closure + `fit_glm_for_subgroup_fast()` dispatch | R/forestsearch_main.R:2085-2104; R/glm_effect_estimators.R (new closure beside `.make_lm_estimator_fast`, :1058); R/subgroup_search.R:975-1003 | data-frame slice + `glm()` + `vcov()` -> vectors + `glm.fit()` + SE from the QR; existence check via the package's own `.fs_existence_reason()` | `.make_or_estimator_fast()`, `.fs_glm_two_group()` | search fits identical %s/%s; machinery %s s/call |",
          st("binary_F_identical"), st("binary_F_n"), s(num("binary_F_machinery_s"), 3)),
  sprintf("| G | `.consistency_glm_pieces()` OR branch + `.dfbeta_glm()` | R/consistency_resample.R:271-272, 168-176 | `glm()` + `lm.influence()` -> `glm.fit()` + hat values from the fit's QR | shares `.fs_glm_pieces_lean()` | pieces rate identical %s/%s; machinery %s s/call |",
          st("binary_pieces_rate_identical"), st("binary_pieces_n"), s(num("binary_pieces_machinery_s"), 3)),
  sprintf("| H | member extraction in `evaluate_consistency_twostage()` and `evaluate_subgroup_consistency()` | R/subgroup_consistency_helpers.R:1803-1811 and 1431-1439 (+ index row extraction :1778, :1391) | `as.data.table(df)` + sequential row filters per candidate -> one precomputed logical column cache and an index matrix, subset only the columns the pieces need | none (hoisting into `subgroup.consistency()`) | Pcons identical, H alone: survival %s, continuous %s, binary %s; machinery (H alone) %s / %s / %s s/call |",
          st("survival_stage_pcons_identical_H_only"), st("continuous_stage_pcons_identical_H_only"), st("binary_stage_pcons_identical_H_only"),
          s(num("survival_stage_H_only_machinery_s"), 3), s(num("continuous_stage_H_only_machinery_s"), 3), s(num("binary_stage_H_only_machinery_s"), 3)),
  "| C, D | none | -- | closed by Stage 0 measurement | -- | -- |",
  "", "Count/IRR is not measured: its search fit and pieces are the binary code paths with `poisson(log)` and an offset; the binary figures are the proxy. RR/RD (multi-tier fallbacks), adjusted and propensity-weighted fits, DINA and GRF identifiers are outside every candidate above.",
  "", "## 5. Reading", "")
for (t in names(FIXTURES)) {
  rep_s <- num(paste0(t, "_rep_s_mean"))
  srch <- if (t == "survival") A_mach_s else if (t == "binary") num("binary_F_machinery_s") else 0
  cons <- num(paste0(t, "_stage_machinery_s"))
  md <- c(md, sprintf("- %s: replicate %s s; search-stage machinery %s s, consistency-stage machinery %s s; combined %s (bar %s: %s). Consistency-stage identity: %s.",
    t, s(rep_s, 2), s(srch, 2), s(cons, 2), pct((srch + cons) / rep_s), pct(BAR), if ((srch + cons) / rep_s >= BAR) "cleared" else "not cleared",
    chr(paste0(t, "_stage_pcons_identical_all"))))
}
md <- c(md, "- Shares depend on each fixture's candidate flow (the consistency stage's cost is proportional to the number of candidates evaluated under full evaluation); the fixture sources above say which configuration each row describes.", "")
writeLines(md, out_path("stage0b_summary.md"))
cat(md, sep = "\n")
cat("\nb3 complete\n")
