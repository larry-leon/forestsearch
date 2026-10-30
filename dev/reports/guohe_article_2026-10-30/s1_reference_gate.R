# s1_reference_gate.R -- reference-literal gate for the Guo & He article (task 2026-10-30, sec. 3).
# Run from the repo root: Rscript dev/reports/guohe_article_2026-10-30/s1_reference_gate.R
# Compares every literal in chunk `reference-data` of the placed article with the manuscript's
# S4 payload (read-only) and GH52_TABLE_7, at printed precision. Exits non-zero on any failure.

art_path <- "vignettes/articles/guohe-nested-continuum.qmd"
pay_path <- "quarto/GuoHe/_payloads/guohe_supp_section/guohe_supp_section_payload.rds"
gh_path  <- "quarto/GuoHe/guohe_sec52_sim.R"
out_dir  <- "dev/reports/guohe_article_2026-10-30/results"

src <- readLines(art_path, warn = FALSE)

# Text of an R assignment starting at the line matching `pat`: accumulate lines until it parses.
grab_expr <- function(lines, pat) {
  i <- grep(pat, lines)
  if (length(i) != 1L) stop("pattern ", pat, " matched ", length(i), " lines")
  for (j in i:length(lines)) {
    txt <- paste(lines[i:j], collapse = "\n")
    ok <- tryCatch({ parse(text = txt); TRUE }, error = function(e) FALSE)
    if (ok) return(txt)
  }
  stop("could not parse expression starting at ", pat)
}
# Body of a chunk by label: lines between its opening fence and the next closing fence, no #| lines.
chunk_body <- function(lines, label) {
  i <- grep(sprintf("^```\\{r %s[,}]", label), lines)
  if (length(i) != 1L) stop("chunk ", label, " found ", length(i), " times")
  ends <- which(lines == "```")
  j <- ends[ends > i][1]
  b <- lines[(i + 1L):(j - 1L)]
  b[!grepl("^#\\|", b)]
}

summ <- chunk_body(src, "summaries")
txt_methods <- c(grab_expr(summ, "^LOWER_METHODS <-"), grab_expr(summ, "^UPPER_METHODS <-"))
txt_ref <- paste(chunk_body(src, "reference-data"), collapse = "\n")
env <- new.env(parent = baseenv())
eval(parse(text = txt_methods), envir = env)
eval(parse(text = txt_ref), envir = env)

pay <- readRDS(pay_path)
need <- c("bound_t7", "ident_t7", "comp_t7", "joint_t7", "pre16", "meta")
if (!all(need %in% names(pay))) stop("payload lacks: ", paste(setdiff(need, names(pay)), collapse = ", "))
gh_src <- readLines(gh_path, warn = FALSE)
genv <- new.env(parent = baseenv())
eval(parse(text = grab_expr(gh_src, "^GH52_TABLE_7 <- ")), envir = genv)
GH52_TABLE_7 <- genv$GH52_TABLE_7

rows <- list()
cmp <- function(table, row, column, literal, payload, d) {
  stopifnot(length(literal) == length(payload), length(row) == length(literal))
  both_na <- is.na(literal) & is.na(payload)
  ad <- abs(literal - payload)
  pass <- both_na | (!is.na(ad) & ad <= 0.5 * 10^-d + 1e-9)
  rows[[length(rows) + 1L]] <<- data.frame(table = table, row = row, column = column,
    literal = literal, payload = payload, abs_diff = ad, pass = pass, stringsAsFactors = FALSE)
}
col_of <- function(df, nm) { if (!nm %in% names(df)) stop("payload column absent: ", nm); df[[nm]] }
beta_key <- function(b) sprintf("%.1f", b)

# ref_lower <-> bound_t7
rl <- env$ref_lower; bt <- pay$bound_t7
mmap <- c("Unadjusted" = "Naive", "Guo-He r = 1/12" = "G&H r=1/12",
          "MR (field)" = "MR (field)", "MR (IJ two-term)" = "MR (IJ)")
idx <- match(paste(beta_key(rl$beta2), mmap[rl$method]), paste(beta_key(bt$beta2), bt$method))
if (anyNA(idx)) stop("ref_lower rows unmatched in bound_t7")
rid <- paste0("beta2=", beta_key(rl$beta2), " | ", rl$method)
for (m in list(c("coverage", "cover", 3), c("wilson_lo", "cover_lo", 3), c("wilson_hi", "cover_hi", 3),
               c("margin", "mean_margin", 3), c("bias", "bias", 3), c("lower_hr", "mean_lower_hr", 2)))
  cmp("ref_lower", rid, m[1], rl[[m[1]]], col_of(bt, m[2])[idx], as.integer(m[3]))

# ref_select <-> ident_t7
rs <- env$ref_select; it <- pay$ident_t7
idx <- match(beta_key(rs$beta2), beta_key(it$beta2)); if (anyNA(idx)) stop("ref_select unmatched")
rid <- paste0("beta2=", beta_key(rs$beta2))
for (m in list(c("c_hat_median", "chat_med", 1), c("c_hat_q25", "chat_q25", 1), c("c_hat_q75", "chat_q75", 1),
               c("c_hat_min", "chat_min", 1), c("c_hat_max", "chat_max", 1), c("specificity", "spec_mean", 3),
               c("ppv", "ppv_mean", 3), c("n_sel_mean", "n_sel_mean", 1), c("n_sel_median", "n_sel_med", 0),
               c("gamma_mean", "gamma_mean", 3), c("dilution", "dilution", 3), c("p_hat_median", "phat_med", 3),
               c("p_hat_q25", "phat_q25", 3), c("p_hat_q75", "phat_q75", 3), c("share_at_30", "p_chat_30", 3),
               c("M_eff_median", "M_eff_med", 2), c("sd_marginal", "sd_marg", 3), c("sd_error", "sd_err", 3)))
  cmp("ref_select", rid, m[1], rs[[m[1]]], col_of(it, m[2])[idx], as.integer(m[3]))

# ref_location
rlo <- env$ref_location
idx <- match(beta_key(rlo$beta2), beta_key(it$beta2)); if (anyNA(idx)) stop("ref_location unmatched")
rid <- paste0("beta2=", beta_key(rlo$beta2))
cmp("ref_location", rid, "theta_H", rlo$theta_H, exp(it$beta2[idx]), 2L)
cmp("ref_location", rid, "hr_in_H_hat", rlo$hr_in_H_hat, exp(col_of(it, "gamma_mean")[idx]), 2L)

# ref_upper <-> comp_t7
ru <- env$ref_upper; ct <- pay$comp_t7
umap <- c("Unadjusted" = "naive", "MR (IJ two-term)" = "ij", "MR (field-s)" = "field_s")
idx <- match(paste(beta_key(ru$beta2), umap[ru$method]), paste(beta_key(ct$beta2), ct$method))
if (anyNA(idx)) stop("ref_upper rows unmatched in comp_t7")
rid <- paste0("beta2=", beta_key(ru$beta2), " | ", ru$method)
for (m in list(c("coverage", "cover", 3), c("wilson_lo", "cover_lo", 3), c("wilson_hi", "cover_hi", 3),
               c("upper_hr", "mean_upper_hr", 2)))
  cmp("ref_upper", rid, m[1], ru[[m[1]]], col_of(ct, m[2])[idx], as.integer(m[3]))

# ref_joint <-> joint_t7
rj <- env$ref_joint; jt <- pay$joint_t7
idx <- match(beta_key(rj$beta2), beta_key(jt$beta2)); if (anyNA(idx)) stop("ref_joint unmatched")
rid <- paste0("beta2=", beta_key(rj$beta2))
for (m in list(c("joint_coverage", "joint_s", 3), c("wilson_lo", "joint_s_lo", 3), c("wilson_hi", "joint_s_hi", 3),
               c("joint_unscaled", "joint", 3), c("unscaled_lo", "joint_lo", 3), c("unscaled_hi", "joint_hi", 3),
               c("corr_field_s", "corr_s_mean", 3), c("corr_unscaled", "corr_mean", 3),
               c("pair_lower_H_hr", "pair_lower_H_hr", 2), c("pair_upper_Hc_hr", "pair_upper_Hc_s_hr", 2)))
  cmp("ref_joint", rid, m[1], rj[[m[1]]], col_of(jt, m[2])[idx], as.integer(m[3]))

# ref_16 <-> pre16 (cells in payload order)
r16 <- env$ref_16; p16 <- pay$pre16
cell_order <- unique(p16$cell)
if (length(cell_order) != nrow(r16)) stop("ref_16 has ", nrow(r16), " rows; pre16 has ", length(cell_order), " cells")
rid <- paste0(cell_order, " | ", r16$design, " ", r16$cell)
for (pm in list(c("gh", "G&H r=1/12"), c("field", "MR (field)"))) {
  sub <- p16[p16$method == pm[2], ]
  idx <- match(cell_order, sub$cell); if (anyNA(idx)) stop("pre16 lacks ", pm[2], " for some cell")
  for (m in list(c("coverage", "cover", 3), c("margin", "margin", 3), c("bias", "bias", 3))) {
    ac <- paste0(pm[1], "_", m[1])
    cmp("ref_16", rid, ac, r16[[ac]], col_of(sub, m[2])[idx], as.integer(m[3]))
  }
}

# ref_gh_published <-> GH52_TABLE_7
rg <- env$ref_gh_published
idx <- match(beta_key(rg$beta2), beta_key(GH52_TABLE_7$beta2)); if (anyNA(idx)) stop("ref_gh_published unmatched")
rid <- paste0("beta2=", beta_key(rg$beta2))
cmp("ref_gh_published", rid, "gh_r_1_12", rg$gh_r_1_12, GH52_TABLE_7$r_1_12[idx], 3L)
cmp("ref_gh_published", rid, "unadjusted", rg$unadjusted, GH52_TABLE_7$naive[idx], 3L)

res <- do.call(rbind, rows)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(out_dir, "reference_gate.csv"), row.names = FALSE)

meta <- pay$meta
stated <- list(nrep = 2000L, draws = 5000L, field_R_out = 1000L, field_R_in = 500L)
agree <- vapply(names(stated), function(k) identical(as.integer(meta[[k]]), stated[[k]]), logical(1))
by_tab <- aggregate(pass ~ table, data = res, FUN = function(x) sprintf("%d/%d", sum(x), length(x)))
kv <- c(sprintf("reference_gate_comparisons = %d", nrow(res)),
        sprintf("reference_gate_passes = %d", sum(res$pass)),
        sprintf("reference_gate_failures = %d", sum(!res$pass)),
        sprintf("reference_gate_by_table: %s = %s", by_tab$table, by_tab$pass),
        sprintf("payload_meta_%s = %s  (article states %s; agree = %s)", names(stated),
                vapply(names(stated), function(k) as.character(meta[[k]]), ""),
                vapply(stated, as.character, ""), agree),
        sprintf("payload_meta_head_at_render = %s", meta$head_at_render),
        sprintf("payload_meta_forestsearch_version = %s", meta$forestsearch_version))
writeLines(kv, file.path(out_dir, "keyvalues.txt"))
cat(kv, sep = "\n")
if (any(!res$pass)) { print(res[!res$pass, ], row.names = FALSE); quit(status = 1L) }
if (!all(agree)) { cat("META DISAGREEMENT\n"); quit(status = 2L) }
cat("REFERENCE GATE: GREEN\n")
