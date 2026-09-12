# Gate 3 -- GRF alignment, per batch, stop-on-failure
# (TASK_grfmr_campaign_2026-09-11, PART A).
#
# The three values the gate names are TEMPLATE LITERALS, not FS_S7_* knobs:
# grf_selection, grf_select_statistic, grf_depth and dmin.grf are set at
# sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:503-506 and are not overridable
# from the environment.  So the gate resolves them three ways, and all three
# must agree:
#
#   1. FROM SOURCE -- the literal at its line in the template, quoted.
#   2. FROM THE RENDERED HTML -- the template's own audit line prints
#      "Run config: method=grf/<grf_selection>", which is the resolved value
#      that batch actually ran under.  It goes to the DOCUMENT, not to the
#      render log (render.sh captures quarto's progress output only), so the
#      gate reads the .html and not the .log.  The grep skips the format
#      string itself, which appears earlier in the echoed source.
#   3. FROM THE BUNDLE -- `admitted_n` is written ONLY by
#      .grf_reselect_on_effect() (R/forestsearch_helpers.R:1651, :1654), which
#      .forestsearch_grf_select() reaches only when grf_select_statistic ==
#      "effect" AND grf_selection == "frontier" (R/forestsearch_helpers.R:1772-
#      1775).  A finite admitted_n on the GRF path is therefore direct
#      per-batch evidence that BOTH resolved as required -- stronger than the
#      audit line, because it is produced by the code path itself.
#
# dmin.grf has no read-back: it is the DR-score PRE-FILTER inside grf
# (grf_main.R:291), consumed before anything that reaches the bundle.  Source
# is the only resolution available for it, and the gate says so rather than
# implying a check it did not make.
#
# usage: Rscript gate3.R <bundle.rds> [<render.html>]
args <- commandArgs(trailingOnly = TRUE)
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
TPL <- file.path(QMD_DIR, "sim_fs_maxeffCons_fb_mr_field_m1_template.qmd")
fail <- 0L; pass <- 0L
say <- function(ok, msg) { if (ok) pass <<- pass + 1L else fail <<- fail + 1L
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", msg)) }

## 1. source
src <- readLines(TPL, warn = FALSE)
grab <- function(nm) {
  i <- grep(sprintf("^%s\\s*<-", gsub(".", "[.]", nm, fixed = TRUE)), src)
  if (!length(i)) return(list(line = NA_integer_, txt = NA_character_, val = NULL))
  i <- i[1]
  list(line = i, txt = trimws(src[i]),
       val = eval(parse(text = sub("#.*$", "", src[i]))))
}
cat("=== 1. RESOLVED FROM SOURCE ===\n")
res <- list()
for (nm in c("grf_selection", "grf_select_statistic", "dmin.grf", "grf_depth")) {
  g <- grab(nm); res[[nm]] <- g$val
  cat(sprintf("  %s:%d  %s\n", basename(TPL), g$line, g$txt))
}
say(identical(res$grf_select_statistic, "effect"),
    sprintf("grf_select_statistic resolves to \"%s\" (required \"effect\")",
            res$grf_select_statistic))
say(identical(res$grf_selection, "frontier"),
    sprintf("grf_selection resolves to \"%s\" (required \"frontier\")",
            res$grf_selection))
say(isTRUE(all.equal(res$dmin.grf, 0.0)),
    sprintf("dmin.grf resolves to %s (required 0.0)", format(res$dmin.grf)))

## 2. render log
if (length(args) >= 2L && file.exists(args[2])) {
  ht <- paste(readLines(args[2], warn = FALSE), collapse = "\n")
  m  <- regmatches(ht, gregexpr("Run config: method=[^<\n]*", ht))[[1]]
  m  <- m[!grepl("%s", m, fixed = TRUE)]   # drop the echoed format string
  cat("\n=== 2. RESOLVED IN THE RENDERED DOCUMENT (the batch's own audit line) ===\n")
  if (length(m)) cat(paste0("  ", m[1], "\n"))
  say(length(m) >= 1L && grepl("method=grf/frontier", m[1]),
      "the batch's audit line reads method=grf/frontier")
} else cat("\n=== 2. rendered document not supplied; skipped ===\n")

## 3. bundle
cat("\n=== 3. CORROBORATED BY THE BUNDLE (the effect/frontier path's own output) ===\n")
b <- readRDS(args[1]); r <- b$results
det <- r$status %in% "DETECTED"
say(identical(as.character(b$meta$subgroup_method), "grf"),
    sprintf("meta$subgroup_method = %s", b$meta$subgroup_method))
say("admitted_n" %in% names(r), "admitted_n column present (Part T2)")
say(sum(is.finite(r$admitted_n)) > 0L,
    sprintf("admitted_n finite on %d of %d rows -- the effect/frontier re-selection ran",
            sum(is.finite(r$admitted_n)), nrow(r)))
cat(sprintf("  admitted_n: median %s, range %s-%s ; n_family: median %s, range %s-%s\n",
            format(median(r$admitted_n, na.rm = TRUE)),
            format(min(r$admitted_n, na.rm = TRUE)), format(max(r$admitted_n, na.rm = TRUE)),
            format(median(r$n_family, na.rm = TRUE)),
            format(min(r$n_family, na.rm = TRUE)), format(max(r$n_family, na.rm = TRUE))))
cat("\nNOTE: dmin.grf is resolved from SOURCE ONLY.  It is the DR-score pre-filter\n")
cat("(grf_main.R:291) and is consumed before any quantity that reaches the bundle,\n")
cat("so there is no read-back for it and the gate does not claim one.\n")
cat(sprintf("\nGATE 3: %d passes, %d failures -- %s\n", pass, fail,
            if (fail) "FAIL: STOP" else "PASS"))
if (fail) quit(status = 1L)
