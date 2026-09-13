src <- readLines(commandArgs(trailingOnly = TRUE)[1], warn = FALSE)
st <- grep("^```[{]r", src); ok <- 0L
for (s in st) {
  e <- s + which(src[(s+1):length(src)] == "```")[1]
  r <- try(parse(text = src[(s+1):(e-1)]), silent = TRUE)
  if (inherits(r, "try-error")) cat("PARSE FAIL at line", s, ":", conditionMessage(attr(r, "condition")), "\n")
  else ok <- ok + 1L
}
cat("chunks parsed:", ok, "of", length(st), "\n")
