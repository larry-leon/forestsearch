# =============================================================================
# Second-pass HTML extraction: the quantities the payload does not carry
# =============================================================================
# Anchors on console output rather than headings (headings also appear in the
# table of contents, which the first pass matched by mistake).

source("dev/replication-check/R/02_extract_html.R")

# window of text around the k-th fixed match of `anchor`
window <- function(tx, anchor, k, before = 2600, after = 1200) {
  p <- gregexpr(anchor, tx, fixed = TRUE)[[1]]
  if (p[1] == -1 || length(p) < k) return(NA_character_)
  substr(tx, max(1, p[k] - before), p[k] + after)
}

squish <- function(s) {
  s <- gsub("\r", "", s)
  s <- gsub("[ \t]+", " ", s)
  gsub("\n{3,}", "\n\n", s)
}

# parse a "CANDIDATE EVALUATION SUMMARY" console block into a data frame
parse_candidates <- function(blk) {
  if (is.na(blk)) return(NULL)
  ln <- strsplit(blk, "\n", fixed = TRUE)[[1]]
  ln <- trimws(sub("^#>", "", ln))
  # rows look like: 1 2.537 61 34 2 0.990 * * * {er <= 0} & {size <= 35}
  rx <- "^([0-9]+) +([0-9.]+) +([0-9]+) +([0-9]+) +([0-9]+) +([0-9.]+) +([*-]) +([*-]) +([*-]) +(.+)$"
  hit <- grep(rx, ln, perl = TRUE)
  if (!length(hit)) return(NULL)
  m <- regmatches(ln[hit], regexec(rx, ln[hit], perl = TRUE))
  do.call(rbind, lapply(m, function(x) data.frame(
    rank = as.integer(x[2]), effect = as.numeric(x[3]),
    N = as.integer(x[4]), E = as.integer(x[5]), K = as.integer(x[6]),
    Pcons = as.numeric(x[7]), on_frontier = x[8], in_band = x[9],
    selected = x[10], subgroup = trimws(x[11]),
    stringsAsFactors = FALSE)))
}

parse_eval_line <- function(blk) {
  if (is.na(blk)) return(NA_character_)
  m <- regmatches(blk, regexpr("Evaluated: [0-9]+ +Passed: [0-9]+ +On frontier: [0-9]+ +In band: [0-9]+ +Selected: m=[0-9]+", blk))
  if (!length(m)) NA_character_ else trimws(m)
}

extract2 <- function(path) {
  tx <- to_text(html_text(path))
  o <- list(file = basename(path))

  # ---- the two consistency-search fits: main FS (k=1), (A) screened (k=2) --
  for (k in 1:2) {
    tag <- c("fs", "fs_dina_screen")[k]
    w <- window(tx, "Subgroup identified:", k, before = 3200, after = 300)
    o[[paste0(tag, "_block")]] <- squish(w)
    o[[paste0(tag, "_subgroup")]] <- if (is.na(w)) NA_character_ else
      trimws(sub(".*Subgroup identified: (\\{[^\n]*).*", "\\1",
                 regmatches(w, regexpr("Subgroup identified: \\{[^\n]*", w))))
    o[[paste0(tag, "_candidates")]] <- parse_candidates(w)
    o[[paste0(tag, "_evalline")]]   <- parse_eval_line(w)
    o[[paste0(tag, "_secs")]] <- {
      s <- regmatches(w, regexpr("Seconds and minutes forestsearch overall = [0-9.]+", w))
      if (length(s)) as.numeric(sub(".*= ", "", s)) else NA_real_
    }
  }

  # ---- (A) reported identified subgroup + size ----------------------------
  o$A_identified <- trimws(grab(tx, "Identified subgroup \\(H\\):[^\n]{0,160}"))
  o$A_size       <- trimws(grab(tx, "Subgroup size[^\n]{0,160}"))
  o$n_pct_lines  <- trimws(grab(tx, "n *= *[0-9]+ *\\([0-9.]+%\\)[^\n]{0,80}", 20))

  # ---- harm confirmation, in document order -------------------------------
  gh <- grab(tx, "(De-biased gate|Gate harm flag|MR harm confirm[a-z]*|harm_flag_debiased|mr_harm_confirmed)[^\n]{0,220}")
  o$harm_lines <- trimws(gh[!grepl("%s|%\\.[0-9]f|%d", gh)])

  # ---- timings ------------------------------------------------------------
  o$fs_completed   <- trimws(grab(tx, "ForestSearch completed in [^\n]{0,60}"))
  o$boot_completed <- trimws(grab(tx, "Bootstrap completed in [^\n]{0,60}"))
  o$cores_line     <- trimws(grab(tx, "Using [0-9]+ of [0-9]+ total cores[^\n]{0,40}"))

  i <- regexpr("Computational Timing\n", tx)
  if (i < 0) i <- tail(gregexpr("Computational Timing", tx, fixed = TRUE)[[1]], 1)
  o$timing_block <- if (i > 0) squish(substr(tx, i, i + 2000)) else NA_character_

  # ---- CV / LOO must not have run -----------------------------------------
  o$cv_skipped  <- length(grab(tx, "K-Fold cross-validation skipped")) > 0
  o$loo_skipped <- length(grab(tx, "Leave-one-out .* skipped")) > 0

  # ---- vocabulary fingerprint ---------------------------------------------
  o$vocab <- c(
    debias_gate = length(grab(tx, "debias_gate")),
    run_debias_gate = length(grab(tx, "run_debias_gate")),
    gate_draws = length(grab(tx, "gate_draws")),
    harm_flag_debiased = length(grab(tx, "harm_flag_debiased")),
    t_gate = length(grab(tx, "t_gate")),
    mr_inference = length(grab(tx, "mr_inference")),
    mr_draws = length(grab(tx, "mr_draws")),
    mr_harm_confirmed = length(grab(tx, "mr_harm_confirmed")),
    t_confirm = length(grab(tx, "t_confirm")),
    mr_in_replicates = length(grab(tx, "mr_in_replicates")))

  o
}

if (!interactive() && !exists(".FS_SOURCED")) {
  a <- commandArgs(trailingOnly = TRUE)
  res <- lapply(a, extract2); names(res) <- basename(a)
  saveRDS(res, "dev/replication-check/out/html_extracts2.rds")
  for (nm in names(res)) {
    r <- res[[nm]]
    cat("\n############################################################\n")
    cat("#", nm, "\n############################################################\n")
    cat("cores          :", paste(r$cores_line, collapse=" | "), "\n")
    cat("CV skipped     :", r$cv_skipped, " LOO skipped:", r$loo_skipped, "\n")
    cat("FS  subgroup   :", r$fs_subgroup, "  (", r$fs_secs, "s )\n")
    cat("FS  eval       :", r$fs_evalline, "\n")
    cat("(A) subgroup   :", r$fs_dina_screen_subgroup, "  (", r$fs_dina_screen_secs, "s )\n")
    cat("(A) eval       :", r$fs_dina_screen_evalline, "\n")
    cat("(A) identified :", paste(r$A_identified, collapse=" | "), "\n")
    cat("fs_completed   :", paste(r$fs_completed, collapse=" | "), "\n")
    cat("boot_completed :", paste(r$boot_completed, collapse=" | "), "\n")
    cat("harm lines:\n"); cat(paste0("   - ", r$harm_lines, collapse="\n"), "\n")
    cat("vocab: "); print(r$vocab)
    cat("\nFS candidate table:\n"); print(r$fs_candidates)
    cat("\n(A) candidate table:\n"); print(r$fs_dina_screen_candidates)
  }
  cat("\nsaved: dev/replication-check/out/html_extracts2.rds\n")
}
