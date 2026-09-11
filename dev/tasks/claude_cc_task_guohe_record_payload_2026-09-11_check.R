#!/usr/bin/env Rscript
# dev/tasks/claude_cc_task_guohe_record_payload_2026-09-11_check.R
# Run from the repository root AFTER the edited record has rendered:
#   Rscript dev/tasks/claude_cc_task_guohe_record_payload_2026-09-11_check.R
# Part 1 (GATE I): the new payload must reproduce the committed records wherever
# they overlap. Tolerance: 0.6 units in the last printed digit of the record's
# value (the 2026-09-05 record transcribed 3-4 dp displays and double-rounds in
# six entries; everything else matches exactly). Any failure -> exit status 1.
# Part 2: prints every payload element as markdown for the REPORT (computed,
# never typed).
p <- readRDS(file.path("quarto", "GuoHe", "_payloads", "guohe_supp_section",
                       "guohe_supp_section_payload.rds"))
tol_of <- function(s) { d <- nchar(sub("^[^.]*\\.?", "", sub("^[+-]", "", s))); 0.6 * 10^(-d) }
fail <- 0L; checked <- 0L
chk <- function(what, got, rec_str) {
  rec <- as.numeric(rec_str); tol <- tol_of(rec_str); ok <- abs(got - rec) <= tol
  checked <<- checked + 1L
  if (!isTRUE(ok)) { fail <<- fail + 1L
    cat(sprintf("  MISMATCH %-40s payload %.6f  record %s  (tol %.4g)\n", what, got, rec_str, tol)) }
}
# ---- REPORT_mr_field_vs_guohe_2026-09-05.md lines 12-27 (G&H at r = 1/30) ----
r0905 <- read.table(header = TRUE, colClasses = "character", text = "
cell ij_b fld_b gh_b ij_c fld_c gh_c ij_m fld_m gh_m
t35_beta2_00 +0.028 +0.010 +0.001 0.999 0.945 0.949 0.50 0.30 0.27
t35_beta2_01 +0.020 +0.004 -0.008 0.999 0.949 0.956 0.50 0.30 0.26
t35_beta2_02 +0.016 +0.003 -0.020 0.999 0.945 0.950 0.51 0.30 0.26
t35_beta2_03 +0.005 -0.003 -0.038 1.000 0.936 0.947 0.52 0.31 0.26
t35_beta2_04 -0.009 -0.012 -0.061 0.997 0.936 0.953 0.54 0.31 0.26
t35_beta2_05 -0.013 -0.011 -0.072 0.997 0.942 0.954 0.56 0.31 0.26
t6_k02 +0.026 +0.008 +0.000 0.998 0.944 0.950 0.50 0.30 0.27
t6_k06 +0.072 +0.030 +0.004 0.987 0.941 0.948 0.40 0.29 0.22
t6_k10 +0.088 +0.038 +0.004 0.982 0.944 0.953 0.37 0.29 0.20
t6_k12 +0.093 +0.041 +0.002 0.972 0.942 0.949 0.36 0.29 0.20
t7_beta2_00 +0.045 +0.030 +0.011 0.998 0.940 0.954 0.58 0.31 0.31
t7_beta2_01 +0.040 +0.026 +0.005 0.998 0.933 0.954 0.58 0.31 0.31
t7_beta2_02 +0.039 +0.025 -0.000 0.999 0.938 0.958 0.59 0.31 0.31
t7_beta2_03 +0.033 +0.021 -0.010 0.999 0.937 0.962 0.60 0.31 0.30
t7_beta2_04 +0.028 +0.017 -0.021 0.999 0.941 0.967 0.60 0.32 0.30
t7_beta2_05 +0.015 +0.007 -0.041 0.999 0.951 0.973 0.61 0.32 0.31")
x <- p$pre16
g <- function(cell, m, col) x[x$cell == cell & x$method == m, col]
for (i in seq_len(nrow(r0905))) { cl <- r0905$cell[i]
  for (m in list(c("MR (IJ)", "ij"), c("MR (field)", "fld"), c("G&H r=1/30", "gh"))) {
    chk(paste(cl, m[1], "bias"),   g(cl, m[1], "bias"),   r0905[[paste0(m[2], "_b")]][i])
    chk(paste(cl, m[1], "cover"),  g(cl, m[1], "cover"),  r0905[[paste0(m[2], "_c")]][i])
    chk(paste(cl, m[1], "margin"), g(cl, m[1], "margin"), r0905[[paste0(m[2], "_m")]][i]) } }
# ---- REPORT_mr_complement_vs_guohe_2026-09-09.md lines 167-172 and 183-188 ----
rc <- read.table(header = TRUE, colClasses = "character", text = "
cell fs_c fs_loc f_c f_loc js jn corr
t7_beta2_00 0.928 0.2922 0.926 0.2868 0.933 0.931 +0.018
t7_beta2_01 0.940 0.2915 0.935 0.2897 0.933 0.932 +0.021
t7_beta2_02 0.938 0.2823 0.939 0.2859 0.942 0.941 +0.021
t7_beta2_03 0.941 0.2713 0.940 0.2775 0.946 0.944 +0.018
t7_beta2_04 0.948 0.2827 0.950 0.2892 0.946 0.947 +0.015
t7_beta2_05 0.936 0.2753 0.939 0.2828 0.942 0.943 +0.009")
cm <- p$comp_t7; jt <- p$joint_t7
for (i in seq_len(nrow(rc))) { cl <- rc$cell[i]
  chk(paste(cl, "field-s cover"), cm$cover[cm$cell == cl & cm$method == "field_s"], rc$fs_c[i])
  chk(paste(cl, "field-s location"), cm$mean_upper[cm$cell == cl & cm$method == "field_s"], rc$fs_loc[i])
  chk(paste(cl, "field cover"), cm$cover[cm$cell == cl & cm$method == "field"], rc$f_c[i])
  chk(paste(cl, "field location"), cm$mean_upper[cm$cell == cl & cm$method == "field"], rc$f_loc[i])
  chk(paste(cl, "joint_s"), jt$joint_s[jt$cell == cl], rc$js[i])
  chk(paste(cl, "joint"), jt$joint[jt$cell == cl], rc$jn[i])
  chk(paste(cl, "corr"), jt$corr_mean[jt$cell == cl], rc$corr[i]) }
# ---- bound_t7 must equal the record's own B2 summary (SUM) by construction ----
b <- p$bound_t7
stopifnot(nrow(b) == 42L, all(is.finite(b$mean_lower)), all(b$n == 2000L))
cat(sprintf("GATE I: %d comparisons, %d outside tolerance\n", checked, fail))
if (fail > 0L) quit(status = 1L)

# ---- Part 2: markdown tables for the REPORT ----------------------------------
md <- function(d, digits = 4) {
  d[] <- lapply(d, function(v) {
    if (!is.numeric(v)) return(v)
    if (all(is.na(v) | v == round(v))) return(ifelse(is.na(v), "NA", format(v, trim = TRUE)))
    ifelse(is.na(v), "NA", formatC(v, format = "f", digits = digits)) })
  cat("| ", paste(names(d), collapse = " | "), " |\n|", paste(rep("---", ncol(d)), collapse = "|"), "|\n", sep = "")
  for (i in seq_len(nrow(d))) cat("| ", paste(unlist(d[i, ]), collapse = " | "), " |\n", sep = "")
  cat("\n")
}
for (nm in c("bound_t7", "ident_t7", "comp_t7", "joint_t7", "naive_convention", "pre16")) {
  cat("### Payload element `", nm, "`\n\n", sep = ""); md(p[[nm]]) }
cat("### Payload meta\n\n")
m <- p$meta; m$inputs_md5 <- NULL
for (k in names(m)) cat("- `", k, "`: ", paste(format(m[[k]]), collapse = " "), "\n", sep = "")
cat("\nInput bundle MD5s:\n\n"); md5 <- p$meta$inputs_md5
for (k in names(md5)) cat("- `", basename(k), "` ", md5[[k]], "\n", sep = "")
