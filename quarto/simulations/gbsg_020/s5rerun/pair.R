# usage: Rscript pair.R <cell> <after.rds> <committed.rds> <out.rds>
a <- commandArgs(TRUE); cell <- a[1]
A <- readRDS(a[2]); B <- readRDS(a[3])
x <- A$results; b <- B$results
stopifnot(identical(A$meta$seed_base, B$meta$seed_base))
b <- b[match(x$sim_id, b$sim_id), ]; stopifnot(identical(x$sim_id, b$sim_id))
q <- c(est = "mr_H_est", lo = "fld_H_lo1s", up = "fld_Hc_up1s_s",
       bonf_lo = "fld_joint_s_bonf_loH", bonf_up = "fld_joint_s_bonf_upHc")
det_same <- identical(as.integer(x$detected), as.integer(b$detected))
d <- x$detected == 1L
same_lab <- (x$label == b$label) | (is.na(x$label) & is.na(b$label))
reg_chg <- d & !(same_lab & x$n_sel == b$n_sel)
adm_chg <- d & (x$admitted_n != b$admitted_n)
top1 <- function(z) trimws(sub("\\|.*$", "", z))
mr_top1_chg <- sum(d & top1(x$p_hat_top_labels) != top1(b$p_hat_top_labels))
errs <- sum(!is.na(x$err_msg) & nzchar(x$err_msg)) ; mrbad <- sum(d & !x$mr_ok)
st <- list()
for (k in names(q)) {
  dd <- x[[q[k]]][d] - b[[q[k]]][d]
  i <- which.max(abs(dd))
  st[[k]] <- data.frame(cell = cell, qty = k, col = q[k], n_det = sum(d),
    before = mean(b[[q[k]]][d]), after = mean(x[[q[k]]][d]),
    mean = mean(dd), median = median(dd), p05 = unname(quantile(dd, .05)), p95 = unname(quantile(dd, .95)),
    maxabs = dd[i], max_sim = x$sim_id[d][i], max_label = x$label[d][i], max_nsel = x$n_sel[d][i],
    nonfinite = sum(!is.finite(dd)), share_nonzero = mean(abs(dd) > 0),
    flips075 = sum((x[[q[k]]][d] >= 0.75) != (b[[q[k]]][d] >= 0.75)),
    flips125 = sum((x[[q[k]]][d] >= 1.25) != (b[[q[k]]][d] >= 1.25)),
    med_dist_thr = median(pmin(abs(b[[q[k]]][d] - 0.75), abs(b[[q[k]]][d] - 1.25))))
}
st <- do.call(rbind, st); rownames(st) <- NULL
res <- list(cell = cell, n = nrow(x), stats = st, det_same = det_same,
  det_rate_after = mean(x$detected == 1L), det_rate_before = mean(b$detected == 1L),
  region_changed = sum(reg_chg), mr_top1_changed = mr_top1_chg, admitted_changed = sum(adm_chg), n_det = sum(d),
  errors = errs, mr_not_ok = mrbad, meta_after = A$meta, meta_before = B$meta[c("forestsearch_version","n_workers","campaign_tag","harm_z1_quantile","harm_prevalence_super","target_hr_harm","n_sample","seed_base")],
  err_ids = x$sim_id[!is.na(x$err_msg) & nzchar(x$err_msg)])
saveRDS(res, a[4])
cat(sprintf("PAIR %s n=%d det_same=%s det_after=%.4f det_before=%.4f n_det=%d region_chg=%d mr_top1_chg=%d errors=%d mr_not_ok=%d\n",
  cell, nrow(x), det_same, res$det_rate_after, res$det_rate_before, sum(d), sum(reg_chg), mr_top1_chg, errs, mrbad))
print(st[, c("qty","before","after","mean","median","p05","p95","maxabs","max_sim","share_nonzero","nonfinite","flips075","flips125","med_dist_thr")], digits = 4)
if (!det_same) quit(status = 2)
