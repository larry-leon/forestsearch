# REPORT_classification_outside_package_2026-11-11_resummarise.R
#
# Part B of dev/tasks/TASK_classification_outside_package_2026-11-11_v2.md.
# Re-summarise the stored $oc (fs_mr_oc_summary() output) of the tracked MR payloads
# with the branch's code, and correct the stored unconditional NPV in place.
#
#   * forestsearch is loaded from this tree with pkgload::load_all(); sequential.
#   * Payloads: every tracked .rds whose object holds an fs_mr_oc object at $oc.
#   * Recompute: fs_mr_oc_summary(payload), all defaults -- exactly as the batch documents
#     call it (sim_fs_mr_field_or_template.qmd:1220-1222, :1316-1317;
#     sim_fs_maxeffCons_mr_field_md_template.qmd:1135, :1221;
#     sim_fs_maxeffCons_mr_md{40,120}_knoise0_n{500,700}_batch_1_1000.qmd:875/893, :968/986).
#     The function reads only $results, $truth and $meta; the stored $oc is not an input.
#   * Per-file gate (Amendment 1, 2026-11-11): the recomputed $oc is compared with the stored
#     $oc using a numeric tolerance of 1e-8 (absolute), not exact equality. A file is BLOCKED
#     (not written) only if some value other than the unconditional NPV,
#     $oc$identification$npv[convention == "unconditional"], differs by more than 1e-8, or
#     the structure (names, classes, dimensions, non-numeric values) differs.
#   * Write (Amendment 1): only where |new - old| unconditional NPV > 1e-8. Only that one cell
#     is written into the stored payload; every other stored value is kept bit for bit. The
#     file is saved in place with the original's serialization version and compression
#     (checked from the file header); after writing, the re-read object is identical() to the
#     original except that one cell (and identical() to the original once it is restored).
#   * Output: <stem>_payloads.csv (file, old value, new value) for the files written;
#     the full per-file log is printed.
#
# Run from the forestsearch repo root:
#   Rscript dev/reports/REPORT_classification_outside_package_2026-11-11_resummarise.R
# Optional: DRY_RUN=1 writes nothing (neither payloads nor the CSV).

t_start <- Sys.time()
dry_run <- identical(Sys.getenv("DRY_RUN"), "1")
repo <- normalizePath(".")
stopifnot(file.exists(file.path(repo, "DESCRIPTION")),
          grepl("^Package: forestsearch$", readLines(file.path(repo, "DESCRIPTION"), 1)))
suppressMessages(pkgload::load_all(repo, quiet = TRUE, export_all = FALSE))
cat("forestsearch", as.character(utils::packageVersion("forestsearch")),
    "loaded from", getNamespaceInfo("forestsearch", "path"), "\n")
stopifnot(normalizePath(getNamespaceInfo("forestsearch", "path")) == repo)

csv_out <- file.path(repo, "dev/reports/REPORT_classification_outside_package_2026-11-11_payloads.csv")

# ---- enumerate tracked payloads holding an fs_mr_oc at $oc ----------------------------
rds <- system2("git", c("-c", "core.quotepath=off", "ls-files", "--", shQuote("*.rds")), stdout = TRUE)
rds <- rds[nzchar(rds)]
stopifnot(all(file.exists(rds)))
cat("tracked .rds:", length(rds), "\n")
has_oc <- vapply(rds, function(f) {
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  is.list(x) && !is.null(names(x)) && "oc" %in% names(x) && inherits(x[["oc"]], "fs_mr_oc")
}, logical(1))
files <- rds[has_oc]
cat("payloads with $oc of class fs_mr_oc:", length(files), "\n")

# ---- serialization header of a file ---------------------------------------------------
rds_header <- function(f) {
  con <- file(f, "rb"); magic <- readBin(con, "raw", 6); close(con)
  comp <- if (identical(magic[1:2], as.raw(c(0x1f, 0x8b)))) "gzip"
          else if (identical(magic[1:3], charToRaw("BZh"))) "bzip2"
          else if (identical(magic[1:6], as.raw(c(0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00)))) "xz"
          else "none"
  con <- if (comp == "none") file(f, "rb") else gzcon(file(f, "rb"))
  if (comp %in% c("bzip2", "xz")) { close(con); con <- (if (comp == "bzip2") bzfile else xzfile)(f, "rb") }
  hdr <- readBin(con, "raw", 6); close(con)
  fmt <- rawToChar(hdr[1:2])
  ver <- if (fmt == "X\n") sum(as.integer(hdr[3:6]) * 256^(3:0)) else NA
  list(compress = comp, format = fmt, version = ver,
       gz10 = if (comp == "gzip") { con <- file(f, "rb"); g <- readBin(con, "raw", 10); close(con); g } else NULL)
}


TOL <- 1e-8
# Differences between two objects beyond TOL: numeric values compared by absolute difference
# (NA pattern must match), everything else (names, classes, dims, non-numeric values) exactly.
tol_diff <- function(a, b, path) {
  out <- character(0); mx <- 0
  rec <- function(a, b, path) {
    if (!identical(class(a), class(b))) { out <<- c(out, sprintf("%s: class %s vs %s", path,
      paste(class(a), collapse = "/"), paste(class(b), collapse = "/"))); return() }
    if (!identical(names(a), names(b))) { out <<- c(out, sprintf("%s: names differ", path)); return() }
    if (!identical(dim(a), dim(b))) { out <<- c(out, sprintf("%s: dims differ", path)); return() }
    if (is.list(a)) {
      if (length(a) != length(b)) { out <<- c(out, sprintf("%s: length differs", path)); return() }
      for (j in seq_along(a)) rec(a[[j]], b[[j]], sprintf("%s$%s", path,
        if (is.null(names(a)) || !nzchar(names(a)[j])) j else names(a)[j]))
      return()
    }
    if (is.numeric(a) && is.numeric(b)) {
      if (length(a) != length(b) || !identical(is.na(a), is.na(b))) {
        out <<- c(out, sprintf("%s: length or NA pattern differs", path)); return() }
      d <- abs(a - b); d <- d[!is.na(d)]
      if (length(d)) { mx <<- max(mx, d); if (any(d > TOL))
        out <<- c(out, sprintf("%s: max abs diff %.3g > 1e-8", path, max(d))) }
      return()
    }
    if (!identical(unname(a), unname(b))) out <<- c(out, sprintf("%s: values differ", path))
  }
  rec(a, b, path)
  structure(out, max_abs = mx)
}

log <- vector("list", length(files)); written <- list()
for (i in seq_along(files)) {
  f <- files[i]
  p <- readRDS(f); old <- p$oc
  new <- tryCatch(fs_mr_oc_summary(p), error = function(e) e)
  if (inherits(new, "error")) {
    log[[i]] <- data.frame(file = f, outcome = "BLOCKED", old = NA, new = NA, max_other = NA,
                           detail = paste("recompute error:", conditionMessage(new))); next
  }
  ou <- which(old$identification$convention == "unconditional")
  nu <- which(new$identification$convention == "unconditional")
  stopifnot(length(ou) == 1L, length(nu) == 1L)
  old_v <- old$identification$npv[ou]; new_v <- new$identification$npv[nu]
  oc_repl <- old; oc_repl$identification$npv[ou] <- new_v
  gate <- tol_diff(oc_repl, new, "$oc")             # per-file gate, tolerance 1e-8
  max_other <- attr(gate, "max_abs")
  if (length(gate)) {
    log[[i]] <- data.frame(file = f, outcome = "BLOCKED", old = old_v, new = new_v,
                           max_other = max_other, detail = paste(gate, collapse = " || ")); next
  }
  if (isTRUE(abs(new_v - old_v) <= TOL) || identical(old_v, new_v)) {
    log[[i]] <- data.frame(file = f, outcome = "unchanged (|change| <= 1e-8)", old = old_v, new = new_v,
                           max_other = max_other, detail = ""); next
  }
  # replace the one cell; every other element must be identical before and after
  p2 <- p; p2$oc$identification$npv[ou] <- new_v
  stopifnot(identical(names(p), names(p2)), identical(attributes(p), attributes(p2)),
            all(vapply(seq_along(p), function(j) names(p)[j] %in% "oc" || identical(p[[j]], p2[[j]]), TRUE)),
            identical(attributes(p$oc), attributes(p2$oc)),
            all(vapply(setdiff(names(p$oc), "identification"), function(k) identical(p$oc[[k]], p2$oc[[k]]), TRUE)),
            identical(p$oc$identification[-which(names(p$oc$identification) == "npv")],
                      p2$oc$identification[-which(names(p2$oc$identification) == "npv")]),
            identical(p$oc$identification$npv[-ou], p2$oc$identification$npv[-ou]))
  h0 <- rds_header(f)
  stopifnot(h0$format == "X\n", h0$version %in% 2:3, h0$compress %in% c("gzip", "bzip2", "xz", "none"))
  if (!dry_run) {
    tmp <- paste0(f, ".tmp_resummarise")
    saveRDS(p2, tmp, version = h0$version, compress = if (h0$compress == "none") FALSE else h0$compress)
    h1 <- rds_header(tmp)
    stopifnot(identical(h1[c("compress", "format", "version")], h0[c("compress", "format", "version")]),
              identical(h1$gz10, h0$gz10), identical(readRDS(tmp), p2))
    stopifnot(file.rename(tmp, f))
    back <- readRDS(f)
    stopifnot(identical(back, p2))                      # the one cell, nothing else
    back$oc$identification$npv[ou] <- old_v
    stopifnot(identical(back, p))                       # restored: the original, bit for bit
  }
  log[[i]] <- data.frame(file = f, outcome = if (dry_run) "would write" else "written",
                         old = old_v, new = new_v, max_other = max_other,
                         detail = sprintf("%s, serialization v%d", h0$compress, h0$version))
  written[[length(written) + 1L]] <- data.frame(file = f, old_npv_unconditional = old_v,
                                                new_npv_unconditional = new_v)
}
log <- do.call(rbind, log)
options(width = 250)
cat("\n---- per-file outcome ----\n")
print(log[, c("file", "outcome", "old", "new", "max_other")], digits = 15, row.names = FALSE)
cat("\n---- BLOCKED (not written) ----\n")
b <- log[log$outcome == "BLOCKED", ]
if (nrow(b)) for (j in seq_len(nrow(b))) cat(sprintf("%s\n    old %.15g new %.15g\n    %s\n", b$file[j], b$old[j], b$new[j], b$detail[j]))
cat("\n---- counts ----\n"); print(table(log$outcome))
if (length(written)) {
  w <- do.call(rbind, written)
  if (!dry_run) write.csv(w, csv_out, row.names = FALSE)
  cat(sprintf("\n%s %d rows to %s\n", if (dry_run) "would write" else "wrote", nrow(w), csv_out))
}
cat(sprintf("wall clock: %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
