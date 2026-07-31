# =============================================================================
# Comparison: new run vs baseline B (v2_2new) and baseline A (v2_2A_linux)
# =============================================================================
# Payload is the source of truth for FS / DINA / GRF rows and engine timings.
# HTML supplies analysis (A) (never written to the payload), the harm
# confirmation flags, and the search / bootstrap wall-clock lines.

.FS_SOURCED <- TRUE
source("dev/replication-check/R/03_extract2.R")

GB  <- "/home/larryleon/Documents/GitHub/forestsearch/quarto/applications/gbsg"
.a <- commandArgs(trailingOnly = TRUE)
NEW_DIR <- if (length(.a)) .a[1] else "dev/replication-check/render"
OUT_TAG <- if (length(.a) > 1) .a[2] else ""
.o <- function(f) file.path("dev/replication-check/out",
                            sub("\\.", paste0(OUT_TAG, "."), f))

paths <- list(
  new = list(payload = file.path(NEW_DIR, "gbsg_table2new_payload.rds"),
             html    = file.path(NEW_DIR, "analysis_gbsg_cox_multimethod_psi_v2_2new.html")),
  B   = list(payload = file.path(GB, "gbsg_table2new_payload.rds"),
             html    = file.path(GB, "analysis_gbsg_cox_multimethod_psi_v2_2new.html")),
  A   = list(payload = file.path(GB, "gbsg_table2Alinux_payload.rds"),
             html    = file.path(GB, "analysis_gbsg_cox_multimethod_psi_v2_2A_linux.html")))

P <- lapply(paths, function(x) readRDS(x$payload))
H <- lapply(paths, function(x) extract2(x$html))
saveRDS(list(P = P, H = H), .o("compare_inputs.rds"))

# ---------------------------------------------------------------------------
# helper: verdict label
# ---------------------------------------------------------------------------
vd <- function(a, b, tol = 0) {
  if (is.null(a) && is.null(b)) return("same")
  if (is.null(a) || is.null(b)) return("differs")
  if (is.numeric(a) && is.numeric(b)) {
    if (length(a) != length(b)) return("differs")
    return(if (isTRUE(all.equal(a, b, tolerance = tol))) "same" else "differs")
  }
  if (identical(as.character(a), as.character(b))) "same" else "differs"
}
fmt <- function(x, d = 3) if (is.numeric(x)) formatC(x, format = "f", digits = d) else as.character(x)

# ---------------------------------------------------------------------------
# 1. Manuscript table: 6 rows (method x region) x 13 columns
# ---------------------------------------------------------------------------
tab_rows <- function(p) {
  t <- p$table
  rownames(t) <- paste(t$method, t$region, sep = "/")
  t
}
TN <- tab_rows(P$new); TB <- tab_rows(P$B); TA <- tab_rows(P$A)
qcols <- c("n","pct","naive_est","naive_lo","naive_hi",
           "fb_est","fb_lo","fb_hi","mr_est","mr_lo","mr_hi")

cmp_table <- do.call(rbind, lapply(rownames(TB), function(r)
  do.call(rbind, lapply(qcols, function(q) data.frame(
    row = r, quantity = q,
    new = TN[r, q], B = TB[r, q], A = TA[r, q],
    vs_B = vd(TN[r, q], TB[r, q]), vs_A = vd(TN[r, q], TA[r, q]),
    stringsAsFactors = FALSE)))))

# ---------------------------------------------------------------------------
# 2. Selected subgroup labels (the primary criterion)
# ---------------------------------------------------------------------------
lab <- function(p, h) c(
  `FS main (fs)`               = p$labels$fs,
  `(A) DINA-screened`          = h$fs_dina_screen_subgroup,
  `(C) DINA-selected (fs_dina)`= p$labels$dina,
  `(D) GRF-selected (fs_grf)`  = p$labels$grf,
  `(B) standalone DINA`        = paste(p$labels$dina_standalone$label, collapse = ", "))
LN <- lab(P$new, H$new); LB <- lab(P$B, H$B); LA <- lab(P$A, H$A)

sizes <- function(p, h) {
  t <- p$table; H_ <- t[t$region == "H", ]
  c(`FS main (fs)` = sprintf("%d (%.1f%%)", H_$n[H_$method=="FS"], H_$pct[H_$method=="FS"]),
    `(A) DINA-screened` = {
      cd <- h$fs_dina_screen_candidates
      if (is.null(cd)) NA_character_ else {
        s <- cd[cd$selected == "*", ][1, ]
        sprintf("%d (%.1f%%)", s$N, 100 * s$N / p$meta$n_total) }},
    `(C) DINA-selected (fs_dina)` = sprintf("%d (%.1f%%)", H_$n[H_$method=="DINA"], H_$pct[H_$method=="DINA"]),
    `(D) GRF-selected (fs_grf)`   = sprintf("%d (%.1f%%)", H_$n[H_$method=="GRF"],  H_$pct[H_$method=="GRF"]),
    `(B) standalone DINA` = as.character(p$labels$dina_standalone$n))
}
SN <- sizes(P$new, H$new); SB <- sizes(P$B, H$B); SA <- sizes(P$A, H$A)

cmp_sg <- data.frame(
  analysis = names(LB),
  new = unname(LN), B = unname(LB), A = unname(LA),
  vs_B = mapply(vd, LN, LB), vs_A = mapply(vd, LN, LA),
  size_new = unname(SN), size_B = unname(SB), size_A = unname(SA),
  size_vs_B = mapply(vd, SN, SB), size_vs_A = mapply(vd, SN, SA),
  stringsAsFactors = FALSE, row.names = NULL)

# ---------------------------------------------------------------------------
# 3. MR harm confirmation
# ---------------------------------------------------------------------------
harm_flags <- function(h) {
  x <- h$harm_lines
  x <- x[grepl("de-biased HR|De-biased gate", x)]
  data.frame(line = substr(x, 1, 110), stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# 4. Timings
# ---------------------------------------------------------------------------
tm <- function(p, h, which) {
  t <- p$meta$timings[[which]]
  # the rename moved this field: pre-rename payloads carry `gate_sec`,
  # post-rename payloads carry `mr_sec`.  Same quantity.
  ms <- if (!is.null(t$mr_sec)) t$mr_sec else t$gate_sec
  data.frame(engine = which, boot_sec = t$boot_sec, mr_sec = ms,
             speedup = t$speedup, n_family = t$n_family,
             selection_bias = t$selection_bias, stringsAsFactors = FALSE)
}
cmp_time <- do.call(rbind, lapply(c("fs","dina","grf"), function(w) {
  a <- tm(P$new,,w); b <- tm(P$B,,w); cc <- tm(P$A,,w)
  data.frame(engine = w,
             quantity = c("boot_sec","mr_sec","speedup","n_family","selection_bias"),
             new = c(a$boot_sec,a$mr_sec,a$speedup,a$n_family,a$selection_bias),
             B   = c(b$boot_sec,b$mr_sec,b$speedup,b$n_family,b$selection_bias),
             A   = c(cc$boot_sec,cc$mr_sec,cc$speedup,cc$n_family,cc$selection_bias),
             stringsAsFactors = FALSE)
}))
cmp_time$vs_B <- ifelse(mapply(function(x,y) isTRUE(all.equal(x,y)), cmp_time$new, cmp_time$B), "same","differs")
cmp_time$vs_A <- ifelse(mapply(function(x,y) isTRUE(all.equal(x,y)), cmp_time$new, cmp_time$A), "same","differs")

wall <- data.frame(
  quantity = c("FS search (console)", "Main bootstrap (console)", "cores"),
  new = c(paste(H$new$fs_completed, collapse="|"), paste(H$new$boot_completed, collapse="|"), paste(H$new$cores_line, collapse="|")),
  B   = c(paste(H$B$fs_completed,   collapse="|"), paste(H$B$boot_completed,   collapse="|"), paste(H$B$cores_line,   collapse="|")),
  A   = c(paste(H$A$fs_completed,   collapse="|"), paste(H$A$boot_completed,   collapse="|"), paste(H$A$cores_line,   collapse="|")),
  stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# 5. Emit
# ---------------------------------------------------------------------------
out <- list(cmp_sg = cmp_sg, cmp_table = cmp_table, cmp_time = cmp_time,
            wall = wall,
            vocab = rbind(new = H$new$vocab, B = H$B$vocab, A = H$A$vocab),
            evallines = data.frame(
              which = c("FS main","(A) screened"),
              new = c(H$new$fs_evalline, H$new$fs_dina_screen_evalline),
              B   = c(H$B$fs_evalline,   H$B$fs_dina_screen_evalline),
              A   = c(H$A$fs_evalline,   H$A$fs_dina_screen_evalline),
              stringsAsFactors = FALSE),
            candidates = list(
              fs  = list(new = H$new$fs_candidates, B = H$B$fs_candidates, A = H$A$fs_candidates),
              scr = list(new = H$new$fs_dina_screen_candidates,
                         B = H$B$fs_dina_screen_candidates, A = H$A$fs_dina_screen_candidates)),
            harm = list(new = harm_flags(H$new), B = harm_flags(H$B), A = harm_flags(H$A)),
            cv_loo = data.frame(run = c("new","B","A"),
                                cv_skipped  = c(H$new$cv_skipped,  H$B$cv_skipped,  H$A$cv_skipped),
                                loo_skipped = c(H$new$loo_skipped, H$B$loo_skipped, H$A$loo_skipped)),
            meta = data.frame(run = c("new","B","A"),
                              built_at = sapply(P, function(p) format(p$built_at)),
                              version  = sapply(P, function(p) p$forestsearch_version),
                              B_boots  = sapply(P, function(p) p$meta$B),
                              n_total  = sapply(P, function(p) p$meta$n_total),
                              cores    = sapply(P, function(p) p$meta$cores),
                              draws    = sapply(P, function(p) {
                                d <- p$meta$mr_draws; if (is.null(d)) p$meta$gate_draws else d }),
                              stringsAsFactors = FALSE))
saveRDS(out, .o("comparison.rds"))

# ---- console report --------------------------------------------------------
cat("\n=============== META ===============\n"); print(out$meta, row.names = FALSE)
cat("\n=============== CV / LOO GATED OFF ===============\n"); print(out$cv_loo, row.names = FALSE)
cat("\n=============== SELECTED SUBGROUPS ===============\n")
print(out$cmp_sg[, c("analysis","new","B","A","vs_B","vs_A")], row.names = FALSE)
cat("\n--- sizes ---\n")
print(out$cmp_sg[, c("analysis","size_new","size_B","size_A","size_vs_B","size_vs_A")], row.names = FALSE)
cat("\n=============== EVAL LINES ===============\n"); print(out$evallines, row.names = FALSE)
cat("\n=============== MANUSCRIPT TABLE ===============\n")
print(out$cmp_table, row.names = FALSE)
cat("\nrows differing vs B:", sum(out$cmp_table$vs_B != "same"),
    " vs A:", sum(out$cmp_table$vs_A != "same"), "\n")
cat("\n=============== TIMINGS (payload) ===============\n"); print(out$cmp_time, row.names = FALSE)
cat("\n=============== WALL CLOCK (console) ===============\n"); print(out$wall, row.names = FALSE)
cat("\n=============== VOCABULARY ===============\n"); print(out$vocab)
cat("\n=============== HARM CONFIRMATION ===============\n")
for (k in c("new","B","A")) { cat("--", k, "--\n"); print(out$harm[[k]], row.names = FALSE) }
cat("\nsaved:", .o("comparison.rds"), "\n")
