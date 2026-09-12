# Part A wall accounting (TASK_dinamr_blockC_grfprobe_2026-09-11).
#
# The authoritative wall is the campaign driver's own output, not a
# reconstruction: render.sh prints WALL_SECONDS per render and campaign.sh
# prints "CELL DONE: <cell>  wall=<n>s" per cell.  This script reads both out of
# the driver log and reconciles them against the bundles' timing columns:
#
#   driver_wall = batch renders + combine render          (measured)
#   compute     = sum(fit_mr_secs) / n_workers            (timing columns)
#   overhead    = batch renders - compute                 (the residual)
#
# The reconciliation is exact.  It also corrects the SHAPE of the overhead model
# used at Gate 1: projectC.R charged 3 x a uniform per-render overhead, but the
# combine render is a flat 9 s on every cell and the whole of the overhead sits
# in the two BATCH renders (the n_super = 1e5 DGM build plus compute_dgm_cde).
# The total was right; the shape was not.  Point DINAMR_PARTA_LOG at the log.
LOG <- Sys.getenv("DINAMR_PARTA_LOG", unset = "partA.log")
setwd(Sys.getenv("DINAMR_QMD_DIR", unset = ".."))
L <- readLines(LOG)
rw <- grep("^WALL_SECONDS=", L, value = TRUE)
R <- data.frame(secs = as.numeric(sub("^WALL_SECONDS=([0-9]+).*", "\\1", rw)),
                out  = sub(".*OUT=(.*)$", "\\1", rw), stringsAsFactors = FALSE)
R$cell <- sub("^dinamr_(.*)_(batch|combine)_[0-9]+$", "\\1", R$out)
R$kind <- sub("^dinamr_.*_(batch|combine)_[0-9]+$", "\\1", R$out)
dn <- grep("^CELL DONE", L, value = TRUE)
D <- data.frame(cell = sub("^CELL DONE: ([A-Za-z0-9_]+).*", "\\1", dn),
                driver_wall = as.numeric(sub(".*wall=([0-9]+)s.*", "\\1", dn)), stringsAsFactors = FALSE)
stem <- c(B_h175_n1500 = "h175_knoise0_n1500_z1q60",
          C124_h100_n500 = "h100_knoise0_n500",  C124_h100_n1000 = "h100_knoise0_n1000",
          C124_h100_n1500 = "h100_knoise0_n1500",
          C31_h100_n500  = "h100_knoise0_n500_z1q60", C31_h100_n1000 = "h100_knoise0_n1000_z1q60",
          C31_h100_n1500 = "h100_knoise0_n1500_z1q60")
proj <- c(B_h175_n1500 = 1.978, C124_h100_n500 = 0.3462, C124_h100_n1000 = 0.2690,
          C124_h100_n1500 = 0.1920, C31_h100_n500 = 0.7084, C31_h100_n1000 = 0.6716,
          C31_h100_n1500 = 0.6347)
out <- do.call(rbind, lapply(D$cell, function(k) {
  fs <- Sys.glob(sprintf("results/dina_effMaxSG_fb_mr_field_m1_%s_nb20_dinamr_res_*.rds", stem[[k]]))
  w <- sum(vapply(fs, function(f) sum(readRDS(f)$results$fit_mr_secs, na.rm=TRUE), 0))
  comp <- w/12
  batch <- sum(R$secs[R$cell==k & R$kind=="batch"]); comb <- sum(R$secs[R$cell==k & R$kind=="combine"])
  dw <- D$driver_wall[D$cell==k]
  data.frame(cell=k, driver_wall_s=dw, driver_wall_h=dw/3600,
             compute_s=round(comp,1), batch_renders_s=batch, combine_s=comb,
             batch_overhead_s=round(batch-comp,1),
             ovh_per_batch_s=round((batch-comp)/2,1),
             proj_h=proj[[k]], realized_over_proj=round((dw/3600)/proj[[k]],3)) }))
print(out, row.names=FALSE, digits=5)
cat(sprintf("\ncells done %d of 7 | wall so far %.3f h | projected for those %.3f h | ratio %.3f\n",
    nrow(out), sum(out$driver_wall_s)/3600, sum(out$proj_h),
    (sum(out$driver_wall_s)/3600)/sum(out$proj_h)))
cat(sprintf("Part A elapsed budget: %.3f h of a 9 h ceiling\n", sum(out$driver_wall_s)/3600))
