# payload_G_compare.R -- §6 G payload comparison (read-only)
A <- "/home/larryleon/Documents/GitHub/forestsearch/quarto/applications/gbsg"
leg <- file.path(A, "_payloads/analysis_gbsg_survival_effMaxSG/analysis_gbsg_survival_effMaxSG_payload.rds")
f <- function(r) file.path(A, "rerun_2026-10-26", r, "_payloads/analysis_gbsg_survival_effMaxSG/analysis_gbsg_survival_effMaxSG_payload.rds")
flat <- function(x) { u <- unlist(rapply(x, function(z) if (is.atomic(z)) z else NULL, how = "list")); u }
cat("Legacy payload path:", leg, "\nexists:", file.exists(leg), "\n")
cat("No legacy G payload exists under quarto/applications/gbsg/_payloads/ (nor _payloads_2026-08-31/09-01*): the legacy\n",
    "record is the HTML only. The leaf-by-leaf comparison below is therefore supplementary: the fresh lean payload vs the\n",
    "era_G payload (the legacy-reproducing build), and lean vs cran040.\n\n", sep = "")
L <- readRDS(f("G_lean")); cat("G_lean payload top-level names:", paste(names(L), collapse = ", "), "\n")
cat("G_lean forestsearch_version:", format(L$forestsearch_version), "\n\n")
for (r in c("G_era_G", "G_cran040", "G_depr")) {
  X <- readRDS(f(r)); a <- flat(X); b <- flat(L)
  drop <- "built_at|timestamp|timing|seconds|minutes|elapsed|forestsearch_version|machine|wall"
  a <- a[!grepl(drop, names(a), ignore.case = TRUE)]; b <- b[!grepl(drop, names(b), ignore.case = TRUE)]
  common <- intersect(names(a), names(b))
  neq <- common[vapply(common, function(k) !identical(as.character(a[[k]]), as.character(b[[k]])), TRUE)]
  cat(sprintf("== %s vs G_lean: leaves %d / %d, common %d, differing %d, only in %s %d, only in lean %d\n", r, length(a), length(b),
              length(common), length(neq), r, length(setdiff(names(a), names(b))), length(setdiff(names(b), names(a)))))
  for (k in neq) cat(sprintf("   %-60s %s -> %s\n", k, format(a[[k]], digits = 7), format(b[[k]], digits = 7)))
  if (length(setdiff(names(b), names(a)))) cat("   only in lean:", paste(head(setdiff(names(b), names(a)), 40), collapse = ", "), "\n")
  if (length(setdiff(names(a), names(b)))) cat("   only in", r, ":", paste(head(setdiff(names(a), names(b)), 40), collapse = ", "), "\n")
}
