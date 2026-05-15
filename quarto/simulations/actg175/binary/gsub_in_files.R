# gsub_in_files.R
#
# Find-and-replace across a vector of files.  Dry-run by default.
#
# Usage:
#   source("gsub_in_files.R")
#   files <- Sys.glob("forestsearch/quarto/simulations/*_fs2.qmd")
#   # or recursively:
#   # files <- list.files("forestsearch/quarto/simulations",
#   #                     pattern = "_fs2\\.qmd$",
#   #                     recursive = TRUE, full.names = TRUE)
#
#   # Preview
#   gsub_in_files(files, 'sim_mode <- "demo"', 'sim_mode <- "full"')
#
#   # Commit
#   gsub_in_files(files, 'sim_mode <- "demo"', 'sim_mode <- "full"',
#                 dry_run = FALSE, backup = TRUE)

gsub_in_files <- function(files, find, replace,
                          dry_run = TRUE,
                          backup  = FALSE,
                          fixed   = TRUE) {
  stopifnot(
    is.character(files),   length(files)   >= 1L,
    is.character(find),    length(find)    == 1L, !is.na(find),
    is.character(replace), length(replace) == 1L, !is.na(replace),
    is.logical(dry_run),   length(dry_run) == 1L,
    is.logical(backup),    length(backup)  == 1L,
    is.logical(fixed),     length(fixed)   == 1L
  )

  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0L) {
    stop("File(s) not found:\n  ", paste(missing_files, collapse = "\n  "),
         call. = FALSE)
  }

  mode_str <- if (dry_run) "[DRY RUN]" else "[COMMIT]"
  cat(sprintf("%s gsub_in_files: %d file(s)\n", mode_str, length(files)))
  cat(sprintf("  find:    %s\n", encodeString(find,    quote = "\"")))
  cat(sprintf("  replace: %s\n", encodeString(replace, quote = "\"")))
  cat(strrep("-", 60), "\n", sep = "")

  summary_rows <- vector("list", length(files))

  for (i in seq_along(files)) {
    f <- files[[i]]
    lines     <- readLines(f, warn = FALSE, encoding = "UTF-8")
    new_lines <- gsub(find, replace, lines, fixed = fixed)
    changed   <- which(lines != new_lines)
    n_match   <- sum(vapply(lines[changed], function(ln) {
      length(gregexpr(find, ln, fixed = fixed)[[1]])
    }, integer(1)))

    if (length(changed) == 0L) {
      cat(sprintf("  %s\n    (no matches)\n", f))
    } else {
      cat(sprintf("  %s\n    %d match(es) on line(s): %s\n",
                  f, n_match, paste(changed, collapse = ", ")))
      if (!dry_run) {
        if (backup) file.copy(f, paste0(f, ".bak"), overwrite = TRUE)
        writeLines(new_lines, f, useBytes = FALSE)
        cat("    -> written\n")
      }
    }

    summary_rows[[i]] <- data.frame(
      file      = f,
      n_matches = n_match,
      modified  = !dry_run && length(changed) > 0L,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, summary_rows)
  cat(strrep("-", 60), "\n", sep = "")
  cat(sprintf("Summary: %d total match(es) across %d file(s); modified: %d\n",
              sum(out$n_matches),
              sum(out$n_matches > 0L),
              sum(out$modified)))
  if (dry_run && any(out$n_matches > 0L)) {
    cat("Re-run with dry_run = FALSE to commit.\n")
  }

  invisible(out)
}
