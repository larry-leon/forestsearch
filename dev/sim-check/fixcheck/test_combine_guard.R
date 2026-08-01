# Exercise the NA-tolerant focus_tag guard added by fix C, in isolation.
`%||%` <- function(a,b) if (is.null(a)||length(a)==0||all(is.na(a))) b else a
.meta_get <- function(b, k) b$meta[[k]] %||% NA

guard <- function(bundles, focus_tag_setup = "maxcons") {
  ftags <- unlist(lapply(bundles, .meta_get, k = "focus_tag"))
  ftags <- unique(ftags[!is.na(ftags)])
  if (length(ftags) > 1L)
    stop(sprintf("combine mode: batches disagree on meta$focus_tag (%s).",
                 paste(ftags, collapse=", ")))
  if (!length(ftags)) cat("    note: no batch records focus_tag; NOT verified\n")
  if (length(ftags) && !identical(ftags, focus_tag_setup))
    warning(sprintf("batches record '%s' but setup has '%s'", ftags, focus_tag_setup),
            call. = FALSE)
  invisible(ftags)
}
mk <- function(ft) if (is.null(ft)) list(meta=list(n_sample=500L)) else
                     list(meta=list(n_sample=500L, focus_tag=ft))
run <- function(label, bundles, setup="maxcons") {
  cat(sprintf("  %-52s ", label))
  r <- tryCatch(withCallingHandlers({guard(bundles, setup); "PASS"},
        warning=function(w){cat("[warn: ", conditionMessage(w), "] ", sep=""); invokeRestart("muffleWarning")}),
        error=function(e) paste("REFUSED ->", conditionMessage(e)))
  cat(r, "\n")
}
cat("### fix C guard behaviour\n")
run("legacy only (no focus_tag) -- must still pool",  list(mk(NULL), mk(NULL)))
run("legacy + new, same focus -- must still pool",    list(mk(NULL), mk("maxcons")))
run("all new, same focus -- must pool",               list(mk("maxcons"), mk("maxcons")))
run("all new, DIFFERENT focus -- must refuse",        list(mk("maxcons"), mk("effMaxSG")))
run("legacy + two different foci -- must refuse",     list(mk(NULL), mk("maxcons"), mk("minSG")))
run("new batches vs mismatched setup knob -- warn",   list(mk("effMaxSG"), mk("effMaxSG")))
cat("\n### the real on-disk batch (pre-fix, no focus_tag)\n")
b <- readRDS("../../../quarto/simulations/gbsg_redux/results/legacy_fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_100.rds")
cat("  meta keys:", paste(names(b$meta), collapse=", "), "\n")
run("existing res_1_100.rds alone", list(b))
run("existing res_1_100.rds + a future post-fix batch", list(b, mk("maxcons")))
