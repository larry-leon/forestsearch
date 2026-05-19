# gsub_in_files.R
#
# Find-and-replace across a vector of files.  Dry-run by default.
# Three modes:
#   (1) Edit in place:  default behavior.
#   (2) rename = TRUE:  move + edit (original file gone after commit).
#   (3) copy   = TRUE:  duplicate + edit (original UNTOUCHED; new file
#                       created alongside it).
# rename and copy are mutually exclusive.
#
# Usage:
#   source("gsub_in_files.R")
#   files <- Sys.glob("forestsearch/quarto/simulations/*_fs2.qmd")
#
#   # (1) Content edit only (existing behavior)
#   gsub_in_files(files, 'sim_mode <- "demo"', 'sim_mode <- "full"',
#                 dry_run = FALSE, backup = TRUE)
#
#   # (2) Rename + content edit (fs2 -> fs1, fs2 files are GONE)
#   gsub_in_files(files, "fs2", "fs1", rename = TRUE,
#                 dry_run = FALSE, backup = TRUE)
#
#   # (3) Copy + content edit (fs2 -> fs1, BOTH versions exist)
#   gsub_in_files(files, "fs2", "fs1", copy = TRUE)          # PREVIEW
#   gsub_in_files(files, "fs2", "fs1", copy = TRUE,
#                 dry_run = FALSE)                            # COMMIT

gsub_in_files <- function(files, find, replace,
                          dry_run = TRUE,
                          backup  = FALSE,
                          fixed   = TRUE,
                          rename  = FALSE,
                          copy    = FALSE) {
  stopifnot(
    is.character(files),   length(files)   >= 1L,
    is.character(find),    length(find)    == 1L, !is.na(find),
    is.character(replace), length(replace) == 1L, !is.na(replace),
    is.logical(dry_run),   length(dry_run) == 1L,
    is.logical(backup),    length(backup)  == 1L,
    is.logical(fixed),     length(fixed)   == 1L,
    is.logical(rename),    length(rename)  == 1L,
    is.logical(copy),      length(copy)    == 1L
  )
  if (rename && copy) {
    stop("rename = TRUE and copy = TRUE are mutually exclusive.  ",
         "rename moves originals; copy keeps originals AND writes new files.",
         call. = FALSE)
  }

  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0L) {
    stop("File(s) not found:\n  ", paste(missing_files, collapse = "\n  "),
         call. = FALSE)
  }

  # Compute target paths (basenames only -- don't touch directories).
  # Both rename and copy use the same target-path rule.
  apply_to_basename <- rename || copy
  new_paths <- if (apply_to_basename) {
    file.path(dirname(files), gsub(find, replace, basename(files), fixed = fixed))
  } else {
    files
  }

  # Pre-flight checks for rename/copy safety
  if (apply_to_basename) {
    dup <- duplicated(new_paths)
    if (any(dup)) {
      stop("Target name collision -- multiple inputs would map to:\n  ",
           paste(unique(new_paths[dup]), collapse = "\n  "), call. = FALSE)
    }
    needs_new_path <- new_paths != files
    clobber <- needs_new_path & file.exists(new_paths) & !(new_paths %in% files)
    if (any(clobber)) {
      stop("Target file(s) already exist (would overwrite):\n  ",
           paste(new_paths[clobber], collapse = "\n  "), call. = FALSE)
    }
  }

  mode_str  <- if (dry_run) "[DRY RUN]" else "[COMMIT]"
  op_str    <- if (rename)      " (with rename)"
  else if (copy)   " (with copy)"
  else             ""
  cat(sprintf("%s gsub_in_files%s: %d file(s)\n", mode_str, op_str, length(files)))
  cat(sprintf("  find:    %s\n", encodeString(find,    quote = "\"")))
  cat(sprintf("  replace: %s\n", encodeString(replace, quote = "\"")))
  cat(strrep("-", 60), "\n", sep = "")

  summary_rows <- vector("list", length(files))

  for (i in seq_along(files)) {
    f        <- files[[i]]
    new_path <- new_paths[[i]]

    lines     <- readLines(f, warn = FALSE, encoding = "UTF-8")
    new_lines <- gsub(find, replace, lines, fixed = fixed)
    changed   <- which(lines != new_lines)
    n_match   <- sum(vapply(lines[changed], function(ln) {
      length(gregexpr(find, ln, fixed = fixed)[[1]])
    }, integer(1)))

    contents_changed <- length(changed) > 0L
    will_rename      <- rename && (new_path != f)
    will_copy        <- copy   && (new_path != f)

    # Per-file diagnostic
    if (will_rename) {
      cat(sprintf("  %s -> %s  (rename)\n", f, new_path))
    } else if (will_copy) {
      cat(sprintf("  %s => %s  (copy)\n", f, new_path))
    } else {
      cat(sprintf("  %s\n", f))
    }
    if (contents_changed) {
      cat(sprintf("    %d content match(es) on line(s): %s\n",
                  n_match, paste(changed, collapse = ", ")))
    } else {
      cat("    (no content matches)\n")
    }

    if (!dry_run) {
      if (will_copy) {
        # Original file is left untouched; new file holds edited content.
        # backup = TRUE is meaningless under copy (original IS the backup);
        # ignore silently.
        writeLines(new_lines, new_path, useBytes = FALSE)
        cat("    -> new file written\n")
      } else {
        # rename or in-place edit -- both may mutate `f`
        if (contents_changed) {
          if (backup) file.copy(f, paste0(f, ".bak"), overwrite = TRUE)
          writeLines(new_lines, f, useBytes = FALSE)
          cat("    -> contents written\n")
        }
        if (will_rename) {
          file.rename(f, new_path)
          cat("    -> renamed\n")
        }
      }
    }

    summary_rows[[i]] <- data.frame(
      file       = f,
      new_file   = new_path,
      n_matches  = n_match,
      modified   = !dry_run && contents_changed && !will_copy,
      renamed    = !dry_run && will_rename,
      copied     = !dry_run && will_copy,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, summary_rows)
  cat(strrep("-", 60), "\n", sep = "")
  cat(sprintf(
    "Summary: %d total match(es) across %d file(s); modified: %d, renamed: %d, copied: %d\n",
    sum(out$n_matches),
    sum(out$n_matches > 0L),
    sum(out$modified),
    sum(out$renamed),
    sum(out$copied)
  ))
  if (dry_run && (any(out$n_matches > 0L) || any(out$file != out$new_file))) {
    cat("Re-run with dry_run = FALSE to commit.\n")
  }

  invisible(out)
}
