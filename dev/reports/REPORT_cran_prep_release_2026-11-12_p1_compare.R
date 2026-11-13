a <- readRDS(commandArgs(TRUE)[1]); b <- readRDS(commandArgs(TRUE)[2])
stopifnot(identical(names(a), names(b)))
ok <- vapply(names(a), function(n) identical(a[[n]], b[[n]]), logical(1))
for (n in names(a)) cat(sprintf("%-28s identical=%s\n", n, ok[[n]]))
cat(sprintf("TOTAL %d/%d identical\n", sum(ok), length(ok)))
# columns dropped before comparison: fs_minutes (wall time) only
