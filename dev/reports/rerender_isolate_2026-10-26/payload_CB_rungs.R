# payload_CB_rungs.R -- supplementary: C and B payload leaves compared across adjacent rungs (read-only)
A <- "/home/larryleon/Documents/GitHub/forestsearch/quarto/applications/actg175"
flat <- function(x) unlist(rapply(x, function(z) if (is.atomic(z)) z else NULL, how = "list"))
drop <- "built_at|timestamp|timing|seconds|minutes|elapsed|forestsearch_version|machine|wall|version"
P <- list(C = c("C_era_C", "analysis_actg175_continuous_oc_intervals/analysis_actg175_continuous_oc_intervals_effmaxsg_payload.rds"),
          B = c("B_era_B", "analysis_actg175_binary_multimethod_fixed_family/analysis_actg175_binary_multimethod_fixed_family_payload.rds"))
leg <- file.path(A, "_payloads/analysis_actg175_continuous_oc_intervals/analysis_actg175_continuous_oc_intervals_effmaxsg_payload.rds")
for (D in names(P)) {
  runs <- c(P[[D]][1], paste0(D, c("_cran040", "_depr", "_lean")))
  X <- lapply(runs, function(r) { f <- flat(readRDS(file.path(A, "rerun_2026-10-26", r, "_payloads", P[[D]][2]))); f[!grepl(drop, names(f), ignore.case = TRUE)] })
  names(X) <- runs
  if (D == "C") { L <- flat(readRDS(leg)); L <- L[!grepl(drop, names(L), ignore.case = TRUE)]; X <- c(list(legacy = L), X) }
  for (j in 2:length(X)) {
    a <- X[[j - 1]]; b <- X[[j]]; cm <- intersect(names(a), names(b))
    neq <- cm[vapply(cm, function(k) !identical(as.character(a[[k]]), as.character(b[[k]])), TRUE)]
    cat(sprintf("== %s: %s -> %s: leaves %d/%d, differing %d, only-from %d, only-to %d\n", D, names(X)[j - 1], names(X)[j],
                length(a), length(b), length(neq), length(setdiff(names(a), names(b))), length(setdiff(names(b), names(a)))))
    for (k in head(neq, 25)) cat(sprintf("   %-55s %s -> %s\n", k, format(a[[k]], digits = 7), format(b[[k]], digits = 7)))
  }
}
