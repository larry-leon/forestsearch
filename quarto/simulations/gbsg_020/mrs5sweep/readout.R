# Section 5 sweep read-out: one row per completed cell, from mrs5sweep/pair/<cell>.rds.
# usage (from gbsg_020/): Rscript mrs5sweep/readout.R  -> markdown on stdout, mrs5sweep/readout.rds
cells <- read.table("mrs5sweep/cells.txt", comment.char = "#",
  col.names = c("id","src","hr","n","z1q","cw","cwk","stem"), stringsAsFactors = FALSE)
f5 <- function(v) sprintf("%+.5f", v)
rows <- list(); all <- list()
for (i in seq_len(nrow(cells))) {
  p <- file.path("mrs5sweep/pair", paste0(cells$id[i], ".rds")); if (!file.exists(p)) next
  r <- readRDS(p); s <- r$stats; all[[cells$id[i]]] <- r
  g <- function(q, k) s[s$qty == q, k]
  big <- s[which.max(abs(s$maxabs)), ]
  prev <- if (cells$z1q[i] == "-") "12.4%" else "31%"
  rows[[length(rows)+1]] <- sprintf("| %s | %s | %.2f | %d | %d/%d | %.4f / %.4f | %s | %s | %s | %s (sim %d, %s) | %d | %d | %d |",
    cells$id[i], prev, cells$hr[i], cells$n[i], r$n_det, r$n, r$det_rate_before, r$det_rate_after,
    paste(f5(c(g("est","mean"), g("est","median"), g("est","p05"), g("est","p95"))), collapse = " / "),
    paste(f5(c(g("lo","mean"), g("lo","median"), g("lo","p05"), g("lo","p95"))), collapse = " / "),
    paste(f5(c(g("up","mean"), g("up","median"), g("up","p05"), g("up","p95"))), collapse = " / "),
    f5(big$maxabs), big$max_sim, big$qty, r$region_changed, r$mr_top1_changed,
    sum(s$flips075[s$qty %in% c("lo","up","bonf_lo","bonf_up")]) + sum(s$flips125[s$qty %in% c("lo","up","bonf_lo","bonf_up")]))
}
cat("| cell | prev | HR | n | declaring | decl. rate before / after | est shift mean / med / p05 / p95 | field lower shift mean / med / p05 / p95 | field-s upper shift mean / med / p05 / p95 | largest single shift (sim, qty) | region chg | MR top-label chg | bound flips at 0.75/1.25 |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|---|---|\n"); cat(unlist(rows), sep = "\n"); cat("\n")
saveRDS(all, "mrs5sweep/readout.rds")
