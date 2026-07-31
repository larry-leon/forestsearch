# =============================================================================
# Extract comparable quantities from a rendered GBSG notebook HTML
# =============================================================================
# The payload (.rds) is the source of truth for the FS / DINA / GRF rows and
# the timings.  Two things the payload does NOT carry and that the brief asks
# for must come from the HTML:
#   * analysis (A), `fs_dina_screen` -- never written to the payload;
#   * the MR harm-confirmation flag (`mr_harm_confirmed` / legacy `harm flag`);
#   * the Computational Timing table (search time).
# Read-only; writes nothing outside dev/replication-check/.

html_text <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  x <- paste(x, collapse = "\n")
  x
}

# strip tags -> plain text, decode the few entities the console output uses
to_text <- function(h) {
  h <- gsub("<script[^>]*>.*?</script>", " ", h, perl = TRUE)
  h <- gsub("<style[^>]*>.*?</style>", " ", h, perl = TRUE)
  h <- gsub("<[^>]+>", "", h, perl = TRUE)
  h <- gsub("&lt;", "<", h, fixed = TRUE)
  h <- gsub("&gt;", ">", h, fixed = TRUE)
  h <- gsub("&amp;", "&", h, fixed = TRUE)
  h <- gsub("&quot;", '"', h, fixed = TRUE)
  h <- gsub("&#39;", "'", h, fixed = TRUE)
  h <- gsub("&nbsp;", " ", h, fixed = TRUE)
  h
}

grab <- function(txt, pat, n = Inf) {
  m <- regmatches(txt, gregexpr(pat, txt, perl = TRUE))[[1]]
  if (is.finite(n)) utils::head(m, n) else m
}

extract_html <- function(path) {
  h <- html_text(path)
  tx <- to_text(h)

  out <- list(file = basename(path))

  # --- selected subgroups from the consistency-search console output --------
  # occurrence 1 = main FS fit; occurrence 2 = (A) fs_dina_screen
  sg <- grab(tx, "Subgroup identified: \\{[^\\n]*")
  sg <- trimws(sub("Subgroup identified: ", "", sg))
  out$subgroup_identified <- sg

  # --- MR / gate harm confirmation (rendered, not the format string) --------
  # pre-rename renders say "Gate harm flag = TRUE"; post-rename may differ.
  gh <- grab(tx, "(Gate harm flag|MR harm confirm[a-z]*|harm_flag_debiased|mr_harm_confirmed)[^\\n]{0,200}")
  gh <- gh[!grepl("%s|%\\.[0-9]f|%d", gh)]          # drop sprintf templates
  out$harm_flag_lines <- trimws(gh)

  # de-biased HR + selection bias, in document order
  out$mr_lines <- trimws(grab(tx, "de-biased HR [0-9.]+[^\\n]{0,160}"))
  out$mr_est   <- as.numeric(sub("^de-biased HR ([0-9.]+).*$", "\\1", out$mr_lines))
  out$sel_bias <- as.numeric(sub(".*selection bias ([0-9.]+).*", "\\1", out$mr_lines))

  # --- de-biased CI lines ---------------------------------------------------
  out$debias_ci <- trimws(grab(tx, "De-biased CI[^\\n]{0,200}"))

  # --- subgroup sizes: "n = 61 (8.9%)" style, plus H/Hc counts --------------
  out$sg_size_lines <- trimws(grab(tx, "[Ss]ubgroup (size|n)[^\\n]{0,120}"))

  # --- consistency / search timing ------------------------------------------
  out$fs_overall_sec <- trimws(grab(tx,
    "Seconds and minutes forestsearch overall = [0-9.]+ [0-9.]+"))
  out$fs_completed <- trimws(grab(tx, "ForestSearch completed in [^\\n]{0,60}"))
  out$boot_completed <- trimws(grab(tx, "Bootstrap completed in [^\\n]{0,60}"))
  out$consistency_method_used <- trimws(grab(tx, "Consistency algorithm used: [^\\n]{0,40}"))

  # --- DINA standalone frontier selection ----------------------------------
  out$dina_sel <- trimws(grab(tx, "(grade3|pgr|er|nodes|size|age|meno) >=?\\s*[0-9.]+[^\\n]{0,40}", 40))

  # --- vocabulary fingerprint (rename attribution) -------------------------
  out$n_debias_gate  <- length(grab(tx, "debias_gate"))
  out$n_run_debias   <- length(grab(tx, "run_debias_gate"))
  out$n_gate_draws   <- length(grab(tx, "gate_draws"))
  out$n_mr_inference <- length(grab(tx, "mr_inference"))
  out$n_mr_draws     <- length(grab(tx, "mr_draws"))
  out$n_mr_harm_conf <- length(grab(tx, "mr_harm_confirmed"))
  out$n_harm_flag_db <- length(grab(tx, "harm_flag_debiased"))
  out$n_t_confirm    <- length(grab(tx, "t_confirm"))
  out$n_t_gate       <- length(grab(tx, "t_gate"))

  # --- computational timing table (text form) ------------------------------
  i <- regexpr("Computational Timing", tx, fixed = TRUE)
  out$timing_block <- if (i > 0) substr(tx, i, i + 1500) else NA_character_

  # --- section (A) block, for manual reading --------------------------------
  j <- regexpr("(A) ForestSearch Using DINA for Screening", tx, fixed = TRUE)
  out$sectionA_block <- if (j > 0) substr(tx, j, j + 4000) else NA_character_

  out
}

if (!interactive() && !exists(".FS_SOURCED")) {
  args <- commandArgs(trailingOnly = TRUE)
  res <- lapply(args, extract_html)
  names(res) <- basename(args)
  saveRDS(res, "dev/replication-check/out/html_extracts.rds")
  for (nm in names(res)) {
    cat("\n==================================================================\n")
    cat(nm, "\n")
    cat("==================================================================\n")
    r <- res[[nm]]
    cat("subgroup_identified :", paste(r$subgroup_identified, collapse = "  |  "), "\n")
    cat("fs_overall          :", paste(r$fs_overall_sec, collapse = " | "), "\n")
    cat("fs_completed        :", paste(r$fs_completed, collapse = " | "), "\n")
    cat("boot_completed      :", paste(r$boot_completed, collapse = " | "), "\n")
    cat("consistency_used    :", paste(r$consistency_method_used, collapse = " | "), "\n")
    cat("mr_est              :", paste(r$mr_est, collapse = ", "), "\n")
    cat("sel_bias            :", paste(r$sel_bias, collapse = ", "), "\n")
    cat("harm_flag_lines     :\n"); cat(paste0("   ", r$harm_flag_lines, collapse = "\n"), "\n")
    cat("vocab: debias_gate=", r$n_debias_gate, " run_debias_gate=", r$n_run_debias,
        " gate_draws=", r$n_gate_draws, " | mr_inference=", r$n_mr_inference,
        " mr_draws=", r$n_mr_draws, " mr_harm_confirmed=", r$n_mr_harm_conf,
        " harm_flag_debiased=", r$n_harm_flag_db,
        " t_gate=", r$n_t_gate, " t_confirm=", r$n_t_confirm, "\n", sep = "")
  }
  cat("\nsaved: dev/replication-check/out/html_extracts.rds\n")
}
