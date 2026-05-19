# =============================================================================
# restructure_bootstrap_chunks.R --- Apply the single-trial-{alt,null}-analysis
#                                     chunk-split restructure across qmd files
# =============================================================================
# Location: quarto/_utils/restructure_bootstrap_chunks.R
# Dependencies: base R only.
#
# WHAT THIS DOES
# Replaces the monolithic `single-trial-alt-analysis` and `single-trial-null-
# analysis` chunks in a per-config article qmd with the six-chunk variants
# (one compute chunk + five render chunks).  See the canonical example in
# actg175_binary_m1_harm_effMinSG-both_fs3.qmd for the structure being
# applied.  Designed as a template for the fs3 family of files; copy and
# update the templates if/when applying to fs1, fs4, fs5, fs6 etc.
#
# WHY YOU NEED THIS
# Bootstrap-results gt tables (summaries$table and pareto_frontier_table())
# do not render correctly from inside an if(){} block when other expressions
# follow them, because knitr's auto-print only fires on the if-block's last
# expression.  Splitting into separate render-only chunks makes each gt the
# sole expression of its own chunk -- the path knit_print.gt_tbl is built for.
#
# SAFETY FEATURES
# * Dry-run by default; nothing is written unless dry_run = FALSE.
# * Fingerprint check before any substitution: a file is touched only when
#   its alt-analysis and null-analysis chunks match the expected pre-fix
#   shape, AND have not already been restructured.
# * Per-file disposition report (transformed / already-done / shape-mismatch /
#   not-found) at the end.
# * Optional .bak files for paranoia (off by default).
# * `.knit_show()` helper -- if present in the setup chunk from a prior
#   attempted fix -- is removed.
# * Removed-only matching: the script never adds chunks to a file that
#   doesn't have the corresponding old chunk in the expected shape.
#
# USAGE
#
#   source("quarto/_utils/restructure_bootstrap_chunks.R")
#
#   # 1. Preview (dry-run, default):
#   files <- Sys.glob("quarto/simulations/actg175/binary/*_fs3.qmd")
#   restructure_bootstrap_chunks(files)
#
#   # 2. Apply for real:
#   restructure_bootstrap_chunks(files, dry_run = FALSE)
#
#   # 3. With backups (.qmd.bak per file):
#   restructure_bootstrap_chunks(files, dry_run = FALSE, backup = TRUE)
#
# This file is NOT part of the forestsearch package; it is an
# exploratory utility for one-off qmd-restructure operations.
# =============================================================================


# ---------------------------------------------------------------------------
# TEMPLATES -- copied verbatim from the canonical fixed file
# (actg175_binary_m1_harm_effMinSG-both_fs3.qmd).
#
# If you adapt this script for a different fs configuration, regenerate
# these templates from the corresponding canonical file: extract the lines
# spanning the opening fence of the compute chunk through the closing
# fence of the *-render-explain chunk, and paste them in as the new
# value of the relevant template.
# ---------------------------------------------------------------------------

.ALT_TEMPLATE <- c(
  "```{r single-trial-alt-analysis}",
  "#| code-fold: true",
  "#| code-summary: \"Single simulated trial + FS analysis alt (click to expand)\"",
  "",
  "itt_fit <- glm(y_sim ~ treat_sim, family = \"binomial\", data = df_sim)",
  "itt_OR <- exp(coef(itt_fit)[\"treat_sim\"])",
  "itt_event <- with(df_sim,mean(y_sim))",
  "",
  "fs_params_single <- fs_params",
  "fs_params_single$show_candidate_summary <- TRUE ",
  "",
  "fs_single <- do.call(forestsearch, c(",
  "  list(df.analysis = df_sim, confounders.name = confounders),",
  "  fs_params_single",
  "))",
  "",
  "if (!is.null(fs_single$sg.harm)) {",
  "  fs_single_tabs <- sg_tables(fs_single, ndecimals = 3)",
  "}",
  "",
  "t0 <- proc.time()",
  "fs <- fs_single",
  "",
  "if (!is.null(fs$sg.harm)) {",
  "  fs_bc <- forestsearch_bootstrap_dofuture(",
  "    fs.est = fs,",
  "    nb_boots = 300,",
  "    show_three = FALSE,",
  "    details = FALSE,",
  "    parallel_args = list(",
  "      plan = \"multisession\",",
  "      workers = n_workers,",
  "      show_message = TRUE",
  "    )",
  "  )",
  "",
  "  plan(\"sequential\")",
  "  timings <- (proc.time() - t0)[\"elapsed\"]",
  "",
  "  cat(\"\\nBootstrap completed in\",",
  "      round(timings / 60, 1), \"minutes\\n\")",
  "",
  "  # Comprehensive summary with diagnostics.  All variables built here are",
  "  # rendered by the small downstream chunks (single-trial-alt-render-*),",
  "  # which exist solely so each gt table is the last (and only) expression",
  "  # of its chunk -- the only way knitr's auto-print reliably emits raw",
  "  # HTML for a gt_tbl object.",
  "  summaries <- summarize_bootstrap_results(",
  "    sgharm = fs$sg.harm,",
  "    boot_results = fs_bc,",
  "    create_plots = TRUE,",
  "    est.scale = \"hr\"",
  "  )",
  "",
  "  ci_tab <- compute_frontier_cis(",
  "    fs,",
  "    n_splits = frontier_params$n_splits,",
  "    seed     = frontier_params$seed",
  "  )",
  "",
  "  # Plot rendered inline here -- the graphics-device path is independent",
  "  # of knitr's stdout capture, so print(p) works inside the if-block.",
  "  p <- plot_pareto_frontier(",
  "    fs,",
  "    ci_table  = ci_tab,",
  "    show_band = TRUE,",
  "    xlim      = frontier_params$xlim",
  "  )",
  "  print(p)",
  "}",
  "```",
  "",
  "```{r single-trial-alt-render-sg10}",
  "#| echo: false",
  "# Render sg10 candidate-summary table.  If fs_single$sg.harm is NULL the",
  "# if-block returns invisible(NULL) and nothing is emitted.",
  "if (!is.null(fs_single$sg.harm)) fs_single_tabs$sg10_out",
  "```",
  "",
  "```{r single-trial-alt-render-tab-estimates}",
  "#| echo: false",
  "# Render the formatted estimates table (gt).",
  "if (!is.null(fs_single$sg.harm)) fs_single_tabs$tab_estimates",
  "```",
  "",
  "```{r single-trial-alt-render-bootstrap-table}",
  "#| echo: false",
  "# Render the bias-corrected bootstrap-results table (gt).",
  "if (!is.null(fs$sg.harm)) summaries$table",
  "```",
  "",
  "```{r single-trial-alt-render-pareto-table}",
  "#| echo: false",
  "# Render the Pareto frontier table.  pareto_frontier_table() may return",
  "# a gt_tbl (normal case), a data.table (fallback when gt is missing), or",
  "# NULL (no consistency output).  All three behave correctly via auto-print.",
  "if (!is.null(fs$sg.harm)) {",
  "  pareto_frontier_table(",
  "    fs,",
  "    ci_table      = ci_tab,",
  "    digits_effect = frontier_params$digits_effect,",
  "    digits_pcons  = frontier_params$digits_pcons,",
  "    digits_ci     = frontier_params$digits_ci",
  "  )",
  "}",
  "```",
  "",
  "```{r single-trial-alt-render-explain}",
  "#| echo: false",
  "#| results: asis",
  "# Markdown explanation -- results: asis allows cat() output to flow through",
  "# as markdown rather than being captured into a <pre> block.  Safe here",
  "# because the explanation is intentional markdown.",
  "if (!is.null(fs$sg.harm)) cat(explain_pareto_selection(fs, ci_table = ci_tab))",
  "```"
)

.NULL_TEMPLATE <- c(
  "```{r single-trial-null-analysis}",
  "#| code-fold: true",
  "#| code-summary: \"Single simulated trial + FS analysis Null (click to expand)\"",
  "",
  "frontier_params <- list(",
  "  n_splits      = 1000L,",
  "  seed          = 1L,",
  "  xlim          = c(1.25, 3.0),",
  "  digits_effect = 2L,",
  "  digits_pcons  = 2L,",
  "  digits_ci     = 2L",
  ")",
  "",
  "itt_fit <- glm(y_sim ~ treat_sim, family = \"binomial\", data = df_sim)",
  "itt_OR <- exp(coef(itt_fit)[\"treat_sim\"])",
  "itt_event <- with(df_sim,mean(y_sim))",
  "",
  "fs_params_single <- fs_params",
  "fs_params_single$show_candidate_summary <- TRUE ",
  "",
  "",
  "fs_single <- do.call(forestsearch, c(",
  "  list(df.analysis = df_sim, confounders.name = confounders),",
  "  fs_params_single",
  "))",
  "",
  "if (!is.null(fs_single$sg.harm)) {",
  "  fs_single_tabs <- sg_tables(fs_single, ndecimals = 3)",
  "}",
  "",
  "t0 <- proc.time()",
  "fs <- fs_single",
  "",
  "# Under the null DGM, the FS analysis typically produces no harm subgroup,",
  "# so the entire bootstrap-and-rendering pathway below is a no-op.  The",
  "# six-chunk structure (compute + 5 render) is kept symmetric with the",
  "# alt-analysis section so the document scaffolding is uniform.",
  "if (!is.null(fs$sg.harm)) {",
  "  fs_bc <- forestsearch_bootstrap_dofuture(",
  "    fs.est = fs,",
  "    nb_boots = 300,",
  "    show_three = FALSE,",
  "    details = FALSE,",
  "    parallel_args = list(",
  "      plan = \"multisession\",",
  "      workers = n_workers,",
  "      show_message = TRUE",
  "    )",
  "  )",
  "",
  "  plan(\"sequential\")",
  "  timings <- (proc.time() - t0)[\"elapsed\"]",
  "",
  "  cat(\"\\nBootstrap completed in\",",
  "      round(timings / 60, 1), \"minutes\\n\")",
  "",
  "  summaries <- summarize_bootstrap_results(",
  "    sgharm = fs$sg.harm,",
  "    boot_results = fs_bc,",
  "    create_plots = TRUE,",
  "    est.scale = \"hr\"",
  "  )",
  "",
  "  ci_tab <- compute_frontier_cis(",
  "    fs,",
  "    n_splits = frontier_params$n_splits,",
  "    seed     = frontier_params$seed",
  "  )",
  "",
  "  p <- plot_pareto_frontier(",
  "    fs,",
  "    ci_table  = ci_tab,",
  "    show_band = TRUE,",
  "    xlim      = frontier_params$xlim",
  "  )",
  "  print(p)",
  "}",
  "```",
  "",
  "```{r single-trial-null-render-sg10}",
  "#| echo: false",
  "if (!is.null(fs_single$sg.harm)) fs_single_tabs$sg10_out",
  "```",
  "",
  "```{r single-trial-null-render-tab-estimates}",
  "#| echo: false",
  "if (!is.null(fs_single$sg.harm)) fs_single_tabs$tab_estimates",
  "```",
  "",
  "```{r single-trial-null-render-bootstrap-table}",
  "#| echo: false",
  "if (!is.null(fs$sg.harm)) summaries$table",
  "```",
  "",
  "```{r single-trial-null-render-pareto-table}",
  "#| echo: false",
  "if (!is.null(fs$sg.harm)) {",
  "  pareto_frontier_table(",
  "    fs,",
  "    ci_table      = ci_tab,",
  "    digits_effect = frontier_params$digits_effect,",
  "    digits_pcons  = frontier_params$digits_pcons,",
  "    digits_ci     = frontier_params$digits_ci",
  "  )",
  "}",
  "```",
  "",
  "```{r single-trial-null-render-explain}",
  "#| echo: false",
  "#| results: asis",
  "if (!is.null(fs$sg.harm)) cat(explain_pareto_selection(fs, ci_table = ci_tab))",
  "```"
)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Locate a single labelled chunk by line range (open fence to next close).
# Returns c(open_line, close_line) inclusive, or NULL when not unique / not
# found / no closing fence.
.find_chunk_lines <- function(lines, chunk_label) {
  open_pattern <- sprintf("^```\\{r %s[},]", chunk_label)
  open_idx     <- grep(open_pattern, lines, perl = TRUE)
  if (length(open_idx) != 1L) return(NULL)
  close_search <- grep("^```$", lines, perl = TRUE)
  close_idx    <- close_search[close_search > open_idx][1]
  if (is.na(close_idx)) return(NULL)
  c(open_idx, close_idx)
}

# Fingerprint check: does the chunk body have the "pre-fix" shape?
# Required: summarize_bootstrap_results + pareto_frontier_table.
# Disqualifying: presence of the post-fix render chunks (e.g. "-render-").
.chunk_is_pre_fix <- function(chunk_body) {
  required <- c("summarize_bootstrap_results",
                "pareto_frontier_table",
                "forestsearch_bootstrap_dofuture")
  ok_required <- all(vapply(required,
                            function(p) any(grepl(p, chunk_body, fixed = TRUE)),
                            logical(1)))
  if (!ok_required) return(FALSE)
  # If the file has already been restructured, the original chunk no longer
  # contains all four render-helpers in its body -- they live in sibling
  # chunks.  But the original chunk body would, post-fix, lack the
  # render-helper invocations.  The single-source check: does the chunk body
  # still contain the bare .knit_show() / summaries$table / pareto_frontier_
  # table line patterns that we replaced?
  has_old_render <- any(grepl("summaries\\$table|\\.knit_show", chunk_body)) ||
                    any(grepl("^\\s*pareto_frontier_table\\(", chunk_body)) ||
                    any(grepl("^\\s*\\.knit_show\\(", chunk_body))
  has_old_render
}

# Detect whether the file has already been restructured: presence of any of
# the new render chunk labels means we should skip.
.file_already_restructured <- function(lines) {
  any(grepl("^```\\{r single-trial-(alt|null)-render-", lines, perl = TRUE))
}

# Remove the .knit_show() helper block from the setup chunk if present.
# Returns modified lines.  No-op if the helper is not present.
.strip_knit_show <- function(lines) {
  start_pattern <- "^# -- Type-safe knit rendering helper"
  start_idx     <- grep(start_pattern, lines)
  if (length(start_idx) == 0L) return(lines)
  if (length(start_idx) > 1L) {
    warning(".knit_show() block found in multiple locations; ",
            "skipping removal to avoid surprises.")
    return(lines)
  }
  # The block ends at the closing brace of the function definition.
  # We look forward for the literal "^}$" that closes the function body.
  end_idx <- grep("^\\}$", lines)
  end_idx <- end_idx[end_idx > start_idx][1]
  if (is.na(end_idx)) {
    warning(".knit_show() block has no clean end; skipping removal.")
    return(lines)
  }
  # Also strip the two blank lines immediately preceding the comment
  # header (they were inserted as a separator from the previous content).
  drop_start <- start_idx
  while (drop_start > 1L && lines[drop_start - 1L] == "") {
    drop_start <- drop_start - 1L
  }
  lines[-seq.int(drop_start, end_idx)]
}

# Replace a chunk's [open_line, close_line] range with the supplied
# replacement-lines vector.  Returns the modified lines.
.replace_chunk <- function(lines, range, replacement) {
  c(
    if (range[1L] > 1L) lines[seq_len(range[1L] - 1L)] else character(0),
    replacement,
    if (range[2L] < length(lines)) lines[seq.int(range[2L] + 1L, length(lines))] else character(0)
  )
}


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

#' Apply the bootstrap chunk-split restructure to a set of qmd files.
#'
#' Arguments
#'   files     character vector of qmd file paths (e.g. from Sys.glob()).
#'   dry_run   logical.  When TRUE (default), no files are written; the
#'             function reports what WOULD be done.  Pass FALSE to apply.
#'   backup    logical.  When TRUE and dry_run = FALSE, writes a .qmd.bak
#'             alongside each modified file before overwriting.  Default
#'             FALSE.
#'   verbose   logical.  When TRUE (default), prints per-file disposition.
#'
#' Returns: invisibly, a data.frame with per-file disposition:
#'   path, alt_status, null_status, knit_show_removed, written.
#'
#' Dispositions (alt_status / null_status):
#'   "transformed"      -- chunk was in pre-fix shape; replacement queued/applied
#'   "already_done"     -- file already has *-render-* chunks; no action
#'   "shape_mismatch"   -- chunk found but body doesn't match expected
#'                         fingerprints; SKIPPED with warning
#'   "not_found"        -- chunk label not found in file; SKIPPED
#'   "ambiguous"        -- chunk label appears more than once; SKIPPED
restructure_bootstrap_chunks <- function(files,
                                         dry_run = TRUE,
                                         backup  = FALSE,
                                         verbose = TRUE) {

  if (!is.character(files) || length(files) == 0L) {
    stop("'files' must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.logical(dry_run) || length(dry_run) != 1L || is.na(dry_run)) {
    stop("'dry_run' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(backup) || length(backup) != 1L || is.na(backup)) {
    stop("'backup' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  }

  results <- data.frame(
    path              = files,
    alt_status        = NA_character_,
    null_status       = NA_character_,
    knit_show_removed = FALSE,
    written           = FALSE,
    stringsAsFactors  = FALSE
  )

  for (i in seq_along(files)) {
    f <- files[i]

    if (!file.exists(f)) {
      results$alt_status[i]  <- "file_missing"
      results$null_status[i] <- "file_missing"
      if (verbose) cat(sprintf("[SKIP] %s -- file not found\n", f))
      next
    }

    lines <- readLines(f, warn = FALSE, encoding = "UTF-8")

    # Short-circuit: file already restructured?
    if (.file_already_restructured(lines)) {
      results$alt_status[i]  <- "already_done"
      results$null_status[i] <- "already_done"
      if (verbose) cat(sprintf("[SKIP] %s -- already restructured\n", f))
      next
    }

    # --- Process alt chunk ---
    alt_range <- .find_chunk_lines(lines, "single-trial-alt-analysis")
    if (is.null(alt_range)) {
      occurrences <- length(grep("^```\\{r single-trial-alt-analysis[},]",
                                 lines, perl = TRUE))
      results$alt_status[i] <- if (occurrences > 1L) "ambiguous" else "not_found"
    } else {
      alt_body <- lines[seq.int(alt_range[1L], alt_range[2L])]
      if (.chunk_is_pre_fix(alt_body)) {
        results$alt_status[i] <- "transformed"
      } else {
        results$alt_status[i] <- "shape_mismatch"
      }
    }

    # --- Process null chunk (look it up fresh, alt_range still valid for now) ---
    null_range <- .find_chunk_lines(lines, "single-trial-null-analysis")
    if (is.null(null_range)) {
      occurrences <- length(grep("^```\\{r single-trial-null-analysis[},]",
                                 lines, perl = TRUE))
      results$null_status[i] <- if (occurrences > 1L) "ambiguous" else "not_found"
    } else {
      null_body <- lines[seq.int(null_range[1L], null_range[2L])]
      if (.chunk_is_pre_fix(null_body)) {
        results$null_status[i] <- "transformed"
      } else {
        results$null_status[i] <- "shape_mismatch"
      }
    }

    # Decide whether to actually rewrite this file
    do_alt  <- results$alt_status[i]  == "transformed"
    do_null <- results$null_status[i] == "transformed"

    if (!do_alt && !do_null) {
      if (verbose) cat(sprintf("[SKIP] %s -- alt:%s null:%s\n",
                               f, results$alt_status[i], results$null_status[i]))
      next
    }

    # Apply substitutions.  IMPORTANT: do null FIRST (it sits later in
    # the file), then alt, so line-number references stay valid.
    new_lines <- lines

    if (do_null) {
      null_range_now <- .find_chunk_lines(new_lines, "single-trial-null-analysis")
      new_lines      <- .replace_chunk(new_lines, null_range_now, .NULL_TEMPLATE)
    }
    if (do_alt) {
      alt_range_now <- .find_chunk_lines(new_lines, "single-trial-alt-analysis")
      new_lines     <- .replace_chunk(new_lines, alt_range_now, .ALT_TEMPLATE)
    }

    # Strip .knit_show() block from setup if present.
    had_knit_show <- any(grepl("knit_show", lines, fixed = TRUE))
    new_lines     <- .strip_knit_show(new_lines)
    results$knit_show_removed[i] <- had_knit_show &&
                                    !any(grepl("knit_show", new_lines, fixed = TRUE))

    if (verbose) {
      cat(sprintf("[%s] %s -- alt:%s null:%s%s\n",
                  if (dry_run) "DRY"  else "WRITE",
                  f,
                  results$alt_status[i],
                  results$null_status[i],
                  if (results$knit_show_removed[i]) " (.knit_show removed)" else ""))
    }

    if (!dry_run) {
      if (backup) {
        bak <- paste0(f, ".bak")
        if (file.exists(bak)) {
          if (verbose) cat(sprintf("  -- backup '%s' already exists, will overwrite\n", bak))
        }
        ok <- file.copy(f, bak, overwrite = TRUE)
        if (!ok) {
          warning(sprintf("Backup copy failed for '%s'; skipping write.", f))
          next
        }
      }
      writeLines(new_lines, f, useBytes = FALSE)
      results$written[i] <- TRUE
    }
  }

  # Summary report
  if (verbose) {
    cat("\n", strrep("-", 60), "\n", sep = "")
    cat(sprintf("Summary: %d file(s) examined\n", nrow(results)))
    tab <- table(alt = results$alt_status, null = results$null_status,
                 useNA = "ifany")
    print(tab)
    if (!dry_run) {
      cat(sprintf("Files written: %d\n", sum(results$written)))
      cat(sprintf(".knit_show removed: %d\n", sum(results$knit_show_removed)))
    } else {
      cat("(dry-run; nothing written.  Re-run with dry_run = FALSE to apply.)\n")
    }
  }

  invisible(results)
}
