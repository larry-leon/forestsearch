# Main-analysis invariance guard — payload-to-payload, field by field, full precision
# (task: cc_task_applications_complete_render_2026-09-01.md §6)
# Transplant of dev/tasks/compare_payloads_applications_0.3.5_2026-09-01.R (committed form
# 7fef766b); changed lines only: pair set (the two FLAGGED gbsg documents, reference =
# the verified _payloads_2026-09-01, new = _payloads_2026-09-01_complete), output default,
# the NEW-CONTENT classification (fields present only on the complete side: absent-in-
# reference names, NULL-in-reference leaves that gained content, and the flag-off
# placeholder cells whose reference value literally encodes the disabled flag), and the
# named-lines footer restricted to the gbsg documents.
#
# Pairs each _payloads_2026-09-01_complete/ file with its _payloads_2026-09-01/ counterpart
# by name and classifies every leaf:
#   EXCLUDED-VOLATILE  — version strings, built_at/timestamps, wall-clock/timing fields,
#                        session-time echoes; enumerated by path with both values printed.
#   STRUCTURAL-VERSION — fields present in one vintage only; listed by path with the side.
#   COMPARED           — everything else; identical() required. Subtrees that pass
#                        identical() are pruned (everything under them is compared-equal);
#                        descent happens only where identical() fails, to localize.
#   NONCOMP-TYPE       — environments (and pointers whose addresses differ): not comparable
#                        across sessions by construction; enumerated, never silently skipped.
# Closures and language objects (ggplot quosures, stored calls) are compared by
# deparse() after removeSource() — their captured environments are session artifacts;
# a deparse-equal pair is recorded under NORMALIZED-EQUAL, a deparse-unequal pair is a
# COMPARED difference.
# Character blobs that fail identical() (captured console text) are diffed line by line;
# a leaf whose differing lines ALL match the volatile-line pattern (timings, worker
# counts, batch topology, paths, version stamps, dates) is EXCLUDED-VOLATILE with every
# differing line printed; any other differing line makes the leaf a COMPARED difference.

OUT <- Sys.getenv("FS_COMPARE_OUT", "/tmp/fs_compare_guard_complete_2026-09-01.md")
REPO <- "."

vol_base <- c("built_at", "forestsearch_version", "time_search", "minutes_all",
              "tmins_search", "tmins_iteration", "boot_sec", "mr_sec", "speedup",
              "timing_seconds")
vol_timing_leaves <- c("loo", "kfold", "fit", "fb", "mr", "gh")
vol_line_rx <- paste0(
  "(second|minute|elapsed|Elapsed|completed in|workers|Workers|Batch |batch_size|",
  "maxRSS|/home/|/Users/|UTC|forestsearch overall|forestsearch[ _]0\\.3\\.|",
  "Sys\\.time|20[0-9]{2}-[0-9]{2}-[0-9]{2} [0-9]{2}:|built_at)")

is_volatile_path <- function(path) {
  base <- sub(".*/", "", path)
  if (base %in% vol_base) return(TRUE)
  if (grepl("/meta/timings/", path, fixed = TRUE) && base %in% vol_timing_leaves) return(TRUE)
  if (grepl("/machine/timestamp$", path)) return(TRUE)
  FALSE
}

fmt <- function(x, max_chars = 4000) {
  s <- if (is.double(x) && is.null(dim(x)) && length(x) <= 20 && !inherits(x, "POSIXt")) {
    paste(sprintf("%.17g", x), collapse = ", ")
  } else {
    paste(utils::capture.output(print(x)), collapse = "\n")
  }
  if (nchar(s) > max_chars) s <- paste0(substr(s, 1, max_chars), " …[truncated]")
  s
}

REC <- new.env(parent = emptyenv())
reset_rec <- function() {
  REC$volatile <- list(); REC$structural <- list(); REC$diff <- list()
  REC$noncomp <- list(); REC$normeq <- list(); REC$newcontent <- list()
}
rec <- function(kind, path, txt) REC[[kind]][[length(REC[[kind]]) + 1L]] <- c(path = path, txt = txt)

deparse_norm <- function(x) {
  x <- tryCatch(utils::removeSource(x), error = function(e) x)
  attributes(x) <- NULL
  paste(deparse(x, width.cutoff = 500L), collapse = "\n")
}

diff_char_blob <- function(a, b, path) {
  fa <- tempfile(); fb <- tempfile()
  writeLines(unlist(strsplit(paste(a, collapse = "\n"), "\n", fixed = TRUE)), fa)
  writeLines(unlist(strsplit(paste(b, collapse = "\n"), "\n", fixed = TRUE)), fb)
  d <- suppressWarnings(system2("diff", c(shQuote(fa), shQuote(fb)), stdout = TRUE))
  unlink(c(fa, fb))
  chg <- grep("^[<>] ", d, value = TRUE)
  bad <- chg[!grepl(vol_line_rx, chg)]
  if (length(chg) == 0L) {                       # differs only in attributes/encoding
    rec("diff", path, paste0("character leaf not identical() but no line differs; ",
                             "attr a: ", fmt(attributes(a)), " | attr b: ", fmt(attributes(b))))
  } else if (length(bad) == 0L) {
    rec("volatile", path, paste0("console/text blob; ", length(chg),
        " differing line(s), all volatile-pattern:\n", paste(chg, collapse = "\n")))
  } else {
    rec("diff", path, paste0("character blob with NON-volatile differing lines (",
        length(bad), " of ", length(chg), "):\n", paste(bad, collapse = "\n"),
        "\n--- all differing lines:\n", paste(chg, collapse = "\n")))
  }
}

atomic_diff <- function(a, b, path) {
  if (length(a) != length(b)) {
    rec("diff", path, paste0("length ", length(a), " vs ", length(b),
                             "\nA: ", fmt(a, 1500), "\nB: ", fmt(b, 1500)))
    return(invisible())
  }
  if (is.numeric(a) && is.numeric(b)) {
    idx <- which(xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b))
    show <- utils::head(idx, 50L)
    rec("diff", path, paste0(length(idx), " differing element(s) of ", length(a),
        "; shown (index: A | B):\n",
        paste(sprintf("  [%d]: %.17g | %.17g", show, as.double(a[show]), as.double(b[show])),
              collapse = "\n"),
        if (length(idx) > 50L) sprintf("\n  … %d more", length(idx) - 50L) else ""))
  } else if (is.character(a)) {
    idx <- which(xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b))
    if (length(idx) == 0L) {
      rec("diff", path, paste0("equal values, differing attributes\nA attr: ",
                               fmt(attributes(a)), "\nB attr: ", fmt(attributes(b))))
    } else if (any(grepl("\n", a[idx])) || any(nchar(a[idx]) > 300)) {
      for (i in utils::head(idx, 10L)) diff_char_blob(a[i], b[i], sprintf("%s[%d]", path, i))
    } else {
      show <- utils::head(idx, 50L)
      lines <- paste(sprintf("  [%d]: %s | %s", show, a[show], b[show]), collapse = "\n")
      # console text captured as a character vector: same triage as blob leaves —
      # if every differing element (both sides) matches the volatile-line pattern,
      # the leaf is a timing/session echo, enumerated with the lines printed
      if (all(grepl(vol_line_rx, a[idx])) && all(grepl(vol_line_rx, b[idx]))) {
        rec("volatile", path, paste0("console-text vector; ", length(idx),
            " differing element(s), all volatile-pattern:\n", lines))
      } else if (all(grepl("run_(loo|cv|cv_kfold)\\s*=\\s*FALSE", a[idx]))) {
        # flag-off placeholder cells: the reference value literally encodes the
        # disabled flag ("— (run_loo = FALSE)"); the complete side carries the value
        rec("newcontent", path, paste0("flag-off placeholder -> complete-side value:\n", lines))
      } else {
        rec("diff", path, paste0(length(idx), " differing element(s):\n", lines))
      }
    }
  } else {
    rec("diff", path, paste0("A: ", fmt(a, 1500), "\nB: ", fmt(b, 1500)))
  }
}

cmp <- function(a, b, path, depth = 0L) {
  if (is_volatile_path(path)) {
    rec("volatile", path, paste0("A(reference): ", fmt(a, 300), "\nB(complete): ", fmt(b, 300)))
    return(invisible())
  }
  if (identical(a, b)) return(invisible())
  if (is.null(a) && !is.null(b)) {   # NEW-CONTENT: NULL in the reference, populated complete-side
    rec("newcontent", path, paste0("reference NULL -> complete-side content:\n", fmt(b, 1500)))
    return(invisible())
  }
  if (!is.null(a) && is.null(b)) {   # content vanished on the complete side — alarming, a real diff
    rec("diff", path, paste0("reference content -> complete-side NULL\nA: ", fmt(a, 1500)))
    return(invisible())
  }
  if (depth > 80L) { rec("diff", path, "max depth exceeded with non-identical subtree"); return(invisible()) }
  if (is.environment(a) && is.environment(b)) {
    rec("noncomp", path, "environment pair, not identical() — session artifact, not compared by content")
    return(invisible())
  }
  if (is.function(a) && is.function(b)) {
    da <- deparse_norm(a); db <- deparse_norm(b)
    if (identical(da, db)) rec("normeq", path, "closure: deparse-equal (differs only by environment/srcref)")
    else rec("diff", path, paste0("closure bodies differ:\nA: ", substr(da, 1, 2000), "\nB: ", substr(db, 1, 2000)))
    return(invisible())
  }
  if (is.language(a) && is.language(b)) {
    da <- deparse_norm(a); db <- deparse_norm(b)
    if (identical(da, db)) rec("normeq", path, "language object: deparse-equal (attrs/env stripped)")
    else rec("diff", path, paste0("language objects differ:\nA: ", da, "\nB: ", db))
    return(invisible())
  }
  if (is.list(a) && is.list(b) && !is.function(a)) {
    na <- names(a); nb <- names(b)
    if (!is.null(na) && !is.null(nb) && !anyDuplicated(na) && !anyDuplicated(nb) &&
        all(nzchar(na)) && all(nzchar(nb))) {
      for (n in setdiff(na, nb)) rec("structural", paste0(path, "/", n), "present in A(reference) only")
      for (n in setdiff(nb, na)) rec("newcontent", paste0(path, "/", n),
        paste0("present only on the complete side:\n", fmt(b[[n]], 1500)))
      for (n in intersect(na, nb)) cmp(a[[n]], b[[n]], paste0(path, "/", n), depth + 1L)
    } else {
      la <- length(a); lb <- length(b)
      if (la != lb) rec("structural", path,
        sprintf("unnamed list length %d in A vs %d in B; extra tail entries uncompared", la, lb))
      for (i in seq_len(min(la, lb))) cmp(a[[i]], b[[i]], sprintf("%s[[%d]]", path, i), depth + 1L)
    }
    aa <- attributes(a); ab <- attributes(b)
    aa$names <- NULL; ab$names <- NULL
    if (!identical(aa, ab)) cmp(aa, ab, paste0(path, "@attr"), depth + 1L)
    return(invisible())
  }
  if (is.atomic(a) && is.atomic(b) && typeof(a) == typeof(b)) {
    ua <- a; ub <- b; attributes(ua) <- NULL; attributes(ub) <- NULL
    if (identical(ua, ub)) {
      aa <- attributes(a); ab <- attributes(b)
      if (!identical(aa, ab)) cmp(aa, ab, paste0(path, "@attr"), depth + 1L)
      return(invisible())
    }
    if (is.character(a) && length(a) == 1L && (grepl("\n", a) || nchar(a) > 300)) {
      diff_char_blob(a, b, path); return(invisible())
    }
    atomic_diff(ua, ub, path); return(invisible())
  }
  if ((typeof(a) == "object" && typeof(b) == "object") || (isS4(a) && isS4(b))) {
    # S7 (ggplot2 >= 4) and S4 objects: properties/slots live in attributes —
    # decompose and compare those recursively rather than falling through
    cmp(attributes(a), attributes(b), paste0(path, "@obj"), depth + 1L)
    return(invisible())
  }
  rec("diff", path, paste0("type/class mismatch: ", typeof(a), "/", paste(class(a), collapse = ","),
                           " vs ", typeof(b), "/", paste(class(b), collapse = ","),
                           "\nA: ", fmt(a, 800), "\nB: ", fmt(b, 800)))
}

# Pass 1 — enumerate every volatile-named field in the NEW payload (deliverable exclusion
# list), with both sides' values, whether or not they differ.
enum_volatile <- function(a, b, path, acc, depth = 0L) {
  if (is_volatile_path(path)) {
    acc[[length(acc) + 1L]] <- c(path = path, txt = paste0("A(reference): ", fmt(a, 300),
                                                          "\nB(complete): ", fmt(b, 300)))
    return(acc)
  }
  if (depth > 80L) return(acc)
  if (is.list(b) && !is.function(b) && !is.null(names(b)) && all(nzchar(names(b)))) {
    for (n in names(b)) {
      av <- if (is.list(a) && n %in% names(a)) a[[n]] else "<absent in reference>"
      acc <- enum_volatile(av, b[[n]], paste0(path, "/", n), acc, depth + 1L)
    }
  }
  acc
}

pairs <- list(
  list(name = "gbsg frozen_family payload (complete vs reference)",
       ref = "quarto/applications/gbsg/_payloads_2026-09-01/analysis_gbsg_survival_frozen_family/analysis_gbsg_survival_frozen_family_payload.rds",
       new = "quarto/applications/gbsg/_payloads_2026-09-01_complete/analysis_gbsg_survival_frozen_family/analysis_gbsg_survival_frozen_family_payload.rds"),
  list(name = "gbsg multimethod payload (complete vs reference)",
       ref = "quarto/applications/gbsg/_payloads_2026-09-01/analysis_gbsg_survival_multimethod/analysis_gbsg_survival_multimethod_payload.rds",
       new = "quarto/applications/gbsg/_payloads_2026-09-01_complete/analysis_gbsg_survival_multimethod/analysis_gbsg_survival_multimethod_payload.rds")
)

sink(OUT)
cat("# Invariance guard — complete (flags ON) vs the verified _payloads_2026-09-01 reference\n\n")
overall_clean <- TRUE
for (p in pairs) {
  cat("\n## ", p$name, "\n\nref: `", p$ref, "`\nnew: `", p$new, "`\n\n", sep = "")
  if (!file.exists(p$ref) || !file.exists(p$new)) {
    cat("**MISSING FILE** — ref exists: ", file.exists(p$ref),
        ", new exists: ", file.exists(p$new), "\n", sep = ""); overall_clean <- FALSE; next
  }
  a <- readRDS(p$ref); b <- readRDS(p$new)
  reset_rec()
  vols <- enum_volatile(a, b, "", list())
  cmp(a, b, "")
  seen <- vapply(vols, `[[`, "", "path")
  for (v in REC$volatile) if (!(v[["path"]] %in% seen)) vols[[length(vols) + 1L]] <- v
  cat("### EXCLUDED-VOLATILE (", length(vols), ")\n\n", sep = "")
  for (v in vols) cat("- `", v[["path"]], "`\n```\n", v[["txt"]], "\n```\n", sep = "")
  cat("\n### STRUCTURAL-VERSION (", length(REC$structural), ")\n\n", sep = "")
  for (v in REC$structural) cat("- `", v[["path"]], "` — ", v[["txt"]], "\n", sep = "")
  cat("\n### NONCOMP-TYPE (", length(REC$noncomp), ")\n\n", sep = "")
  for (v in REC$noncomp) cat("- `", v[["path"]], "` — ", v[["txt"]], "\n", sep = "")
  cat("\n### NORMALIZED-EQUAL (", length(REC$normeq), ")\n\n", sep = "")
  for (v in REC$normeq) cat("- `", v[["path"]], "` — ", v[["txt"]], "\n", sep = "")
  cat("\n### NEW-CONTENT (", length(REC$newcontent), ")\n\n", sep = "")
  for (v in REC$newcontent) cat("- `", v[["path"]], "`\n```\n", v[["txt"]], "\n```\n", sep = "")
  cat("\n### COMPARED differences (", length(REC$diff), ")\n\n", sep = "")
  for (v in REC$diff) cat("- `", v[["path"]], "`\n```\n", v[["txt"]], "\n```\n", sep = "")
  verdict <- if (length(REC$diff) == 0L) "IDENTICAL on every COMPARED leaf" else "NOT CLEAN"
  if (length(REC$diff) > 0L) overall_clean <- FALSE
  cat("\n**Verdict:** ", verdict, "\n", sep = "")
}

cat("\n\n# Named selection lines (gbsg documents, reference vs complete)\n\n")
nm <- function(path, get) { x <- readRDS(path); get(x) }
for (side in c("2026-09-01", "2026-09-01_complete")) {
  ff  <- sprintf("quarto/applications/gbsg/_payloads_%s/analysis_gbsg_survival_frozen_family/analysis_gbsg_survival_frozen_family_payload.rds", side)
  mmp <- sprintf("quarto/applications/gbsg/_payloads_%s/analysis_gbsg_survival_multimethod/analysis_gbsg_survival_multimethod_payload.rds", side)
  cat("## ", side, "\n", sep = "")
  cat("- gbsg frozen_family sg_harm: ", deparse(nm(ff, function(x) x$labels$sg_harm)), "\n", sep = "")
  cat("- gbsg multimethod labels: ", deparse(nm(mmp, function(x) x$labels[c("fs","dina","grf")])), "\n", sep = "")
}
cat("\nOVERALL: ", if (overall_clean) "CLEAN — no COMPARED differences in any pair"
                   else "NOT CLEAN — COMPARED differences present", "\n", sep = "")
sink()
cat("wrote ", OUT, "\n", sep = "")
