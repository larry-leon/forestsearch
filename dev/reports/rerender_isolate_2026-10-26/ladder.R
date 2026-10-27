# ladder.R -- build-ladder comparison (task cc_task_isolate_rerender_2026-10-26, §6)
# Run from dev/reports/rerender_isolate_2026-10-26/:  Rscript ladder.R
# Reuses ../rerender_2026-10-26/extract_lib.R and the C / B key-quantity specs of ../rerender_2026-10-26/compare.R
# (evaluated from that file's parse tree, not copied); adds G's specs.
# Outputs: ladder_<D>.csv, ladder_<D>.md, outputs_<D>_<rung>.txt, rungdiff_counts.csv
prev <- "../rerender_2026-10-26"
source(file.path(prev, "extract_lib.R"))
ex <- parse(file.path(prev, "compare.R"))
keep <- vapply(ex, function(e) {
  d <- deparse(e)[1]
  grepl("^(extract_prose|S|specs|grab|clauses|decimals|verdict) <- ", d)
}, TRUE)
for (e in ex[keep]) eval(e, globalenv())

specs$G <- rbind(
  S("combinations evaluated", "^Combinations evaluated: (\\d+)"),
  S("passed effect threshold", "^Passed effect threshold .*: (\\d+)"),
  S("candidates found", "^Found (\\d+) subgroup"),
  S("to evaluate / on frontier / in band", "^To evaluate: (.*)$"),
  S("consistency summary line", "^Evaluated: (\\d+ +Passed: .*)$"),
  S("Pcons (selected)", "^% consistency criteria met= (\\S+)"),
  S("rank-1 candidate row", "^1 +([0-9.]+ .*)$"),
  S("selected subgroup", "^Subgroup identified: (.*)$"),
  S("MR harm confirmation (MR point estimate)", "^MR harm confirmation: (.*)$"),
  S("max HR estimate", "Maximum HR estimate = (\\S+)"),
  S("H row (naive)", "^Questionable1 \\| (.*)$"),
  S("Hc row (naive)", "^Recommend \\| (.*)$", 1L),
  S("Pareto selected row", "^★ \\| (.*)$"),
  S("Pareto frontier size", "^(\\d+) non-dominated candidate"),
  S("FB success", "^Successful subgroup ID: +(.*)$"),
  S("MR candidate family size", "^Candidate family size: (\\d+)"),
  S("scale cross-check", "^Scale cross-check.*: (.*)$"),
  S("10-fold CV Sensitivity (H)", "^Sensitivity \\(H\\) \\| Agreement.* \\| (\\S+)$"),
  S("10-fold CV PPV (H)", "^PPV \\(H\\) \\| Positive.* \\| (\\S+)$"),
  S("10-fold CV Exact Match", "^Exact Match \\| Exact.* \\| (\\S+)$"),
  S("10-fold CV partial match", "^Partial match \\(>=1 shared factor\\) \\| (.*)$"),
  S("LOO line", "^LOO \\[(.*)$"),
  S("LOO agreement", "^agreement: (.*)$"),
  S("LOO KfA_0 row", "^KfA_0 \\| (.*)$"),
  S("LOO KfA_1 row", "^KfA_1 \\| (.*)$"),
  S("summary: Naive", "^Naive \\| ([0-9.]+)$"),
  S("summary: Full bootstrap", "^Full bootstrap \\| ([0-9.]+ \\(.*)$"),
  S("summary: Multiplier resampling (headline)", "^Multiplier resampling \\| (.*)$"),
  S("summary: Leave-one-out", "^Leave-one-out \\(out-of-sample\\) \\| (.*)$"))

docs <- list(
  C = list(app = "actg175", stem = "analysis_actg175_continuous_intervals", era = "era_C"),
  G = list(app = "gbsg",    stem = "analysis_gbsg_survival_effMaxSG",        era = "era_G"),
  B = list(app = "actg175", stem = "analysis_actg175_binary_multimethod_fixed_family", era = "era_B"))
root <- "/home/larryleon/Documents/GitHub/forestsearch/quarto/applications"
cause <- c(era = "not reproducible under the era build on this machine (platform/environment or stale record)",
           cran040 = "package history between the era commit and 0.4.0",
           depr = "deprecation (candidate family)",
           lean = "lean fits")

dcounts <- list()
for (D in names(docs)) {
  d <- docs[[D]]
  rungs <- c(legacy = file.path(root, d$app, paste0(d$stem, ".html")),
             era     = file.path(root, d$app, "rerun_2026-10-26", paste0(D, "_", d$era), paste0(d$stem, ".html")),
             cran040 = file.path(root, d$app, "rerun_2026-10-26", paste0(D, "_cran040"), paste0(d$stem, ".html")),
             depr    = file.path(root, d$app, "rerun_2026-10-26", paste0(D, "_depr"), paste0(d$stem, ".html")),
             lean    = file.path(root, d$app, "rerun_2026-10-26", paste0(D, "_lean"), paste0(d$stem, ".html")))
  have <- setNames(file.exists(rungs), names(rungs))
  O <- lapply(names(rungs), function(r) if (have[[r]]) extract_outputs(rungs[[r]]) else NULL); names(O) <- names(rungs)
  P <- lapply(names(rungs), function(r) if (have[[r]]) extract_prose(rungs[[r]]) else NULL);   names(P) <- names(rungs)
  for (r in names(rungs)) if (have[[r]]) writeLines(O[[r]], sprintf("outputs_%s_%s.txt", D, r))
  sp <- specs[[D]]
  tab <- do.call(rbind, lapply(seq_len(nrow(sp)), function(i) {
    v <- vapply(names(rungs), function(r) {
      if (!have[[r]]) return("(not rendered)")
      src <- if (sp$src[i] == "prose") P[[r]] else O[[r]]
      grab(src, sp$re[i], sp$k[i])
    }, "")
    path <- vapply(names(rungs)[-1], function(r) {
      if (!have[[r]]) return("-")
      switch(verdict(v[["legacy"]], v[[r]]), identical = "=", "same clauses, reordered" = "~", "absent in both" = "=", "x")
    }, "")
    first <- names(path)[path %in% c("x", "~")][1]
    rn <- names(rungs)
    moves <- rn[-1][vapply(2:length(rn), function(j) have[[rn[j]]] && have[[rn[j - 1]]] && !identical(v[[rn[j]]], v[[rn[j - 1]]]), TRUE)]
    data.frame(quantity = sp$quantity[i], legacy = v[["legacy"]], era = v[["era"]], cran040 = v[["cran040"]],
               depr = v[["depr"]], lean = v[["lean"]], path = paste(path, collapse = ""),
               first_departs = if (is.na(first)) "none" else paste0(first, if (path[[first]] == "~") " (clause order only)" else ""),
               cause = if (is.na(first)) "" else cause[[first]],
               moves_at = paste(moves, collapse = ","))
  }))
  write.csv(tab, sprintf("ladder_%s.csv", D), row.names = FALSE)
  esc <- function(x) gsub("|", "\\|", ifelse(is.na(x), "NA", x), fixed = TRUE)
  md <- c(sprintf("| quantity | legacy | %s | cran040 | depr | lean | path (era,cran040,depr,lean vs legacy) | first departs | changes at |", d$era),
          "|---|---|---|---|---|---|---|---|---|",
          sprintf("| %s | %s | %s | %s | %s | %s | `%s` | %s | %s |", esc(tab$quantity), esc(tab$legacy), esc(tab$era), esc(tab$cran040),
                  esc(tab$depr), esc(tab$lean), tab$path, tab$first_departs, tab$moves_at))
  writeLines(md, sprintf("ladder_%s.md", D))
  # line-diff counts between adjacent rungs
  pairs <- list(c("legacy", "era"), c("era", "cran040"), c("cran040", "depr"), c("depr", "lean"), c("legacy", "lean"))
  for (p in pairs) {
    if (!all(have[p])) { dcounts[[length(dcounts) + 1]] <- data.frame(doc = D, from = p[1], to = p[2], lines_from = NA, lines_to = NA, removed = NA, added = NA); next }
    f1 <- sprintf("outputs_%s_%s.txt", D, p[1]); f2 <- sprintf("outputs_%s_%s.txt", D, p[2])
    dd <- suppressWarnings(system2("diff", c(f1, f2), stdout = TRUE))
    dcounts[[length(dcounts) + 1]] <- data.frame(doc = D, from = p[1], to = p[2], lines_from = length(O[[p[1]]]), lines_to = length(O[[p[2]]]),
                                                 removed = sum(grepl("^< ", dd)), added = sum(grepl("^> ", dd)))
  }
  cat(sprintf("%s: %d quantities; first departs: %s\n", D, nrow(tab), paste(names(table(tab$first_departs)), table(tab$first_departs), sep = "=", collapse = ", ")))
}
dc <- do.call(rbind, dcounts)
write.csv(dc, "rungdiff_counts.csv", row.names = FALSE)
print(dc)
