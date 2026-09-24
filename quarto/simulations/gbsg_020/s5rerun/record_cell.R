# usage: Rscript record_cell.R <cell> <combined.rds> <wall_s> <workers> <out.tsv>
a <- commandArgs(TRUE); B <- readRDS(a[2]); x <- B$results
`%||%` <- function(u, v) if (is.null(u) || !length(u)) v else u
err <- !is.na(x$err_msg) & nzchar(x$err_msg)
ok <- nrow(x) == 2000L && identical(sort(x$sim_id), 1:2000) && !anyDuplicated(x$sim_id)
row <- data.frame(cell = a[1], wall_s = as.integer(a[3]), workers = as.integer(a[4]), meta_workers = B$meta$n_workers %||% NA,
  reps = nrow(x), simid_1_2000 = ok, errors = sum(err), err_ids = paste(x$sim_id[err], collapse = ";"),
  decl_rate = mean(x$detected == 1L), n_decl = sum(x$detected == 1L), campaign = B$meta$campaign_tag %||% NA,
  fs_version = B$meta$forestsearch_version %||% NA, seed_base = B$meta$seed_base %||% NA)
write.table(row, a[5], sep = "\t", row.names = FALSE, quote = FALSE, append = file.exists(a[5]), col.names = !file.exists(a[5]))
print(row)
