# Section 5 re-run: do the quoted FS coverage ranges change?  Per harm cell, the four products
# of REPORT_fs_products_reconciliation_2026-09-12.md (field lower, field-s upper, joint, joint_s;
# detected replicates; definitions as scripts_dinamr/grfmr_numbers.R / blockC_numbers.R), before
# (committed bundle) and after (s5rerun).  Two cell sets: (a) the Section 5 set (p12x20 at 12.4%,
# cert20/e1stud at 31%); (b) the designated-comparator set of the reconciliation report (p12ext /
# tier2 at 12.4%, not re-run here, unchanged; cert20/e1stud at 31%, before -> after).
# usage (from gbsg_020/): Rscript s5rerun/ranges.R
cv <- function(r) { d <- r[r$detected %in% 1L, ]
  c(field_lo = mean(d$fld_H_lo1s <= d$betaHhat_H), fields_up = mean(d$fld_Hc_up1s_s >= d$betaHhat_Hc),
    joint = mean(d$betaHhat_H >= d$fld_joint_bonf_loH & d$betaHhat_Hc <= d$fld_joint_bonf_upHc),
    joint_s = mean(d$betaHhat_H >= d$fld_joint_s_bonf_loH & d$betaHhat_Hc <= d$fld_joint_s_bonf_upHc)) }
cells <- read.table("mrs5sweep/cells.txt", comment.char = "#",
  col.names = c("id","src","hr","n","z1q","cw","cwk","stem"), stringsAsFactors = FALSE)
cells <- cells[cells$hr > 1, ]
af <- function(c1) sprintf("results/fs_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_s5rerun_combined_1_2000.rds",
  round(100 * c1$hr), c1$n, if (c1$z1q != "-") "_z1q60" else "")
comp124 <- function(hr, n) sprintf("results/fs_maxeffCons_fb_mr_field_m1_h%03d_knoise0_n%d_%s_combined_1_2000.rds",
  round(100 * hr), n, if (abs(hr - 1.75) < 1e-9) "tier2" else "p12ext")
R <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) { c1 <- cells[i, ]
  b <- cv(readRDS(c1$stem)$results); a <- cv(readRDS(af(c1))$results)
  data.frame(cell = c1$id, src = c1$src, hr = c1$hr, n = c1$n, product = names(b), before = b, after = a) }))
rownames(R) <- NULL
C124 <- do.call(rbind, lapply(which(cells$src == "p12x20"), function(i) { c1 <- cells[i, ]
  v <- cv(readRDS(comp124(c1$hr, c1$n))$results)
  data.frame(cell = paste0(c1$id, "(comparator)"), src = "p12ext/tier2", hr = c1$hr, n = c1$n, product = names(v), before = v, after = v) }))
rownames(C124) <- NULL
rng <- function(D, col) tapply(D[[col]], D$product, function(v) sprintf("%.4f-%.4f", min(v), max(v)))
out <- c("## Per harm cell (4 dp), before -> after", "",
  "| cell | source | HR | n | field lower | field-s upper | joint | joint_s |", "|---|---|---|---|---|---|---|---|")
for (k in unique(R$cell)) { s <- R[R$cell == k, ]; f <- function(p) { x <- s[s$product == p, ]
    if (x$before == x$after) sprintf("%.4f", x$before) else sprintf("**%.4f -> %.4f**", x$before, x$after) }
  out <- c(out, sprintf("| %s | %s | %.2f | %d | %s | %s | %s | %s |", k, s$src[1], s$hr[1], s$n[1],
    f("field_lo"), f("fields_up"), f("joint"), f("joint_s"))) }
S5 <- R; D2 <- rbind(R[R$src != "p12x20", ], C124)
out <- c(out, "", "## Ranges over the 12 harm cells, before / after", "",
  "| cell set | field lower | field-s upper | joint | joint_s |", "|---|---|---|---|---|",
  sprintf("| (a) Section 5 set, before | %s | %s | %s | %s |", rng(S5,"before")["field_lo"], rng(S5,"before")["fields_up"], rng(S5,"before")["joint"], rng(S5,"before")["joint_s"]),
  sprintf("| (a) Section 5 set, after | %s | %s | %s | %s |", rng(S5,"after")["field_lo"], rng(S5,"after")["fields_up"], rng(S5,"after")["joint"], rng(S5,"after")["joint_s"]),
  sprintf("| (b) designated-comparator set, before | %s | %s | %s | %s |", rng(D2,"before")["field_lo"], rng(D2,"before")["fields_up"], rng(D2,"before")["joint"], rng(D2,"before")["joint_s"]),
  sprintf("| (b) designated-comparator set, after | %s | %s | %s | %s |", rng(D2,"after")["field_lo"], rng(D2,"after")["fields_up"], rng(D2,"after")["joint"], rng(D2,"after")["joint_s"]))
writeLines(out, "s5rerun/ranges.md"); cat(out, sep = "\n")

# ---- HR 1.00 cells: the bound-location figures of REPORT_grfmr_completion_2026-09-12.md (FS rows) ----
nul <- read.table("mrs5sweep/cells.txt", comment.char = "#",
  col.names = c("id","src","hr","n","z1q","cw","cwk","stem"), stringsAsFactors = FALSE)
nul <- nul[nul$hr == 1, ]
loc <- function(r) { d <- r[r$detected %in% 1L & is.finite(r$fld_H_lo1s) & is.finite(r$betaHhat_H), ]
  sprintf("%d | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f", nrow(d), median(d$fld_H_lo1s), median(d$betaHhat_H),
    median(d$fld_H_lo1s) / median(d$betaHhat_H), median(d$fld_H_lo1s / d$betaHhat_H),
    mean(d$fld_H_lo1s >= 1.00), mean(d$fld_H_lo1s >= 1.25)) }
o2 <- c("", "## HR 1.00 cells: field lower bound location (detected replicates), before / after", "",
  "| cell | source | n | arm | n_eval | median lower | median theta(Hhat) | bound/theta | paired ratio | share >= 1.00 | share >= 1.25 |",
  "|---|---|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(nul))) { c1 <- nul[i, ]
  o2 <- c(o2, sprintf("| %s | %s | %d | before | %s |", c1$id, c1$src, c1$n, loc(readRDS(c1$stem)$results)),
              sprintf("| %s | %s | %d | after | %s |", c1$id, c1$src, c1$n, loc(readRDS(af(c1))$results))) }
write(o2, "s5rerun/ranges.md", append = TRUE); cat(o2, sep = "\n")
