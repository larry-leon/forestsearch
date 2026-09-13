# Realized walls against the Gate 1 projection, grfmr completion
# (TASK_grfmr_completion_2026-09-12).  Reads the driver's own `CELL DONE ...
# wall=` lines (logs/grfmrC.driver.log) and projectionGC.rds.
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
p <- readRDS(file.path(SCRATCH, "projectionGC.rds"))$cells
L <- readLines(file.path(SCRATCH, "logs/grfmrC.driver.log"))
d <- grep("^CELL DONE", L, value = TRUE)
w <- data.frame(cell = sub("^CELL DONE: ([^ ]+).*$", "\\1", d),
                wall_s = as.numeric(sub("^.*wall=([0-9]+)s$", "\\1", d)))
x <- merge(p[, c("order", "cell", "h_uncal", "h_med")], w, by = "cell")
x <- x[order(x$order), ]
x$real_h <- x$wall_s / 3600
x$ratio_cal <- x$real_h / x$h_med
x$ratio_uncal <- x$real_h / x$h_uncal
print(x, digits = 4, row.names = FALSE)
cat(sprintf("\nTOTAL realized %d s = %.4f h ; projected %.4f h calibrated / %.4f h uncalibrated ; ratio %.4f / %.4f\n",
    sum(x$wall_s), sum(x$real_h), sum(x$h_med), sum(x$h_uncal),
    sum(x$real_h) / sum(x$h_med), sum(x$real_h) / sum(x$h_uncal)))
h2 <- x$order <= 2
cat(sprintf("two harm cells: realized %.4f h, projected %.4f h (ratio %.4f)\n",
    sum(x$real_h[h2]), sum(x$h_med[h2]), sum(x$real_h[h2]) / sum(x$h_med[h2])))
cat(sprintf("six null cells: realized %.4f h, projected %.4f h (ratio %.4f)\n",
    sum(x$real_h[!h2]), sum(x$h_med[!h2]), sum(x$real_h[!h2]) / sum(x$h_med[!h2])))
b <- grep("^WALL_SECONDS", L, value = TRUE)
r <- data.frame(render = sub("^.*OUT=grfmr_", "", b),
                secs = as.numeric(sub("^WALL_SECONDS=([0-9]+).*$", "\\1", b)))
cat(sprintf("\nper-render walls (%d renders):\n", nrow(r)))
print(r, row.names = FALSE)
st <- sub("^=== GRFMR C RUN START (.*) ; hard.*$", "\\1", grep("RUN START", L, value = TRUE))
en <- sub("^GRFMR C RUN COMPLETE ", "", grep("RUN COMPLETE", L, value = TRUE))
cat(sprintf("\nrun span: %s -> %s\n", st, en))
saveRDS(x, file.path(SCRATCH, "wallsGC.rds"))
