#!/usr/bin/env Rscript
# qmd_chunk_index.R — read-only static index of the R code chunks in a .qmd
#
# Usage:  Rscript qmd_chunk_index.R <file.qmd> <out.md>
#
# Base R only; no package is loaded and nothing in the document is evaluated.
# For every code chunk the script records: label, header options and `#|`
# options, line range, top-level objects assigned (or modified in place),
# free symbols (used but not assigned earlier in the same chunk), functions
# called, namespace-qualified calls (pkg::fun), numeric literals, RNG /
# draw / root-finding / integration calls, and file I/O calls. It then
# derives the chunk dependency graph (which earlier chunk last assigned each
# free symbol) and lists the objects that the prose reads through inline
# `r ...` code. Anything the parser cannot handle is reported, never guessed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: Rscript qmd_chunk_index.R <file.qmd> <out.md>")
qmd_path <- args[[1L]]
out_path <- args[[2L]]
stopifnot(file.exists(qmd_path))
lines <- readLines(qmd_path, warn = FALSE)

## ---- 1. locate chunks ------------------------------------------------------
open_re  <- "^\\s*(`{3,})\\{([A-Za-z=][A-Za-z0-9_.]*)(.*)\\}\\s*$"
close_re <- "^\\s*(`{3,})\\s*$"
chunks <- list()
i <- 1L
n <- length(lines)
while (i <= n) {
  if (grepl(open_re, lines[i])) {
    fence <- sub(open_re, "\\1", lines[i])
    lang  <- sub(open_re, "\\2", lines[i])
    hdr   <- trimws(sub(open_re, "\\3", lines[i]))
    j <- i + 1L
    while (j <= n && !(grepl(close_re, lines[j]) &&
                       nchar(sub(close_re, "\\1", lines[j])) >= nchar(fence))) j <- j + 1L
    body <- if (j > i + 1L) lines[(i + 1L):(j - 1L)] else character(0)
    chunks[[length(chunks) + 1L]] <- list(lang = lang, header = hdr, start = i,
                                          end = min(j, n), body = body,
                                          unterminated = j > n)
    i <- j + 1L
  } else i <- i + 1L
}

## ---- 2. per-chunk options and label --------------------------------------
parse_header <- function(hdr) {
  hdr <- sub("^,\\s*", "", hdr)
  if (!nzchar(hdr)) return(list(label = NA_character_, opts = character(0)))
  # split on commas outside quotes/parens (good enough for chunk headers)
  parts <- character(0); depth <- 0L; cur <- ""; inq <- ""
  for (ch in strsplit(hdr, "")[[1L]]) {
    if (nzchar(inq)) { cur <- paste0(cur, ch); if (ch == inq) inq <- ""; next }
    if (ch %in% c("'", '"')) { inq <- ch; cur <- paste0(cur, ch); next }
    if (ch %in% c("(", "[", "{")) depth <- depth + 1L
    if (ch %in% c(")", "]", "}")) depth <- depth - 1L
    if (ch == "," && depth == 0L) { parts <- c(parts, trimws(cur)); cur <- ""; next }
    cur <- paste0(cur, ch)
  }
  parts <- c(parts, trimws(cur))
  parts <- parts[nzchar(parts)]
  label <- NA_character_
  if (length(parts) && !grepl("=", parts[1L], fixed = TRUE)) {
    label <- parts[1L]; parts <- parts[-1L]
  }
  list(label = label, opts = parts)
}

for (k in seq_along(chunks)) {
  ch <- chunks[[k]]
  ph <- parse_header(ch$header)
  hash_opts <- grep("^\\s*#\\|", ch$body, value = TRUE)
  hash_opts <- trimws(sub("^\\s*#\\|\\s?", "", hash_opts))
  lab <- ph$label
  hl <- grep("^label:\\s*", hash_opts, value = TRUE)
  if (is.na(lab) && length(hl)) lab <- trimws(sub("^label:\\s*", "", hl[1L]))
  if (is.na(lab)) lab <- sprintf("unlabelled-%02d", k)
  chunks[[k]]$label <- lab
  chunks[[k]]$opts  <- c(ph$opts, hash_opts)
  chunks[[k]]$code  <- ch$body[!grepl("^\\s*#\\|", ch$body)]
}

## ---- 3. AST walkers ----------------------------------------------------------
ASSIGN_OPS <- c("<-", "<<-", "=", "->", "->>")
MODIFIER_FUNS <- c("$", "[[", "[", "@", "attr", "names", "dimnames", "levels",
                   "class", "dim", "colnames", "rownames", "body", "formals",
                   "environment", "attributes", "diag", "substr", "is.na")
RNG_FUNS <- c("set.seed", "rnorm", "runif", "rbinom", "rpois", "rexp", "rgamma",
              "rbeta", "rt", "rchisq", "sample", "sample.int", "rmvnorm", "mvrnorm",
              "RNGkind", "nextRNGStream", "nextRNGSubStream")
ANALYTIC_FUNS <- c("pnorm", "qnorm", "dnorm", "pmvnorm", "qmvnorm", "integrate",
                   "uniroot", "optimize", "optimise", "optim", "nlm", "nlminb",
                   "pt", "qt", "pchisq", "qchisq", "pbinom", "qbinom")
IO_FUNS <- c("readRDS", "saveRDS", "source", "load", "save", "read.csv",
             "read.table", "readr::read_csv", "fread", "file.path", "sys.source",
             "write.csv", "writeLines", "readLines", "here", "file.exists")

lhs_root <- function(e) {
  # returns list(root = symbol name or NA, text = deparsed LHS, modify = TRUE/FALSE)
  if (is.symbol(e)) return(list(root = as.character(e), text = as.character(e), modify = FALSE))
  if (is.character(e) && length(e) == 1L) return(list(root = e, text = e, modify = FALSE))
  if (is.call(e)) {
    f <- e[[1L]]
    fn <- if (is.symbol(f)) as.character(f) else paste(deparse(f), collapse = "")
    if (fn %in% MODIFIER_FUNS || grepl("<-$", fn)) {
      inner <- lhs_root(e[[2L]])
      return(list(root = inner$root, text = paste(deparse(e), collapse = ""), modify = TRUE))
    }
    if (fn %in% c("(", "{")) return(lhs_root(e[[2L]]))
  }
  list(root = NA_character_, text = paste(deparse(e), collapse = ""), modify = TRUE)
}

new_acc <- function() {
  a <- new.env()
  for (nm in c("assigned", "assigned_text", "uses", "funs", "ns", "nums", "strs",
               "rng", "analytic", "io", "io_text", "pkgs")) assign(nm, character(0), envir = a)
  a
}

is_empty <- function(x) identical(x, quote(expr = ))
walk_elems <- function(lst, acc, locals, from = 1L) {
  # iterate by index so the empty symbol (e.g. x[1, ]) is never bound to a variable
  if (length(lst) < from) return(invisible())
  for (i in from:length(lst)) if (!is_empty(lst[[i]]) && !is.null(lst[[i]])) walk(lst[[i]], acc, locals)
  invisible()
}

walk <- function(e, acc, locals = character(0)) {
  # acc: environment with character vectors: assigned, assigned_text, uses,
  #      funs, ns, nums, strs, rng, analytic, io, io_text, local_defs
  if (is.symbol(e)) {
    nm <- as.character(e)
    if (nzchar(nm) && !(nm %in% locals)) acc$uses <- c(acc$uses, nm)
    return(invisible())
  }
  if (is.numeric(e) && length(e) == 1L) { acc$nums <- c(acc$nums, format(e, digits = 15)); return(invisible()) }
  if (is.character(e) && length(e) == 1L) { acc$strs <- c(acc$strs, e); return(invisible()) }
  if (!is.call(e)) return(invisible())
  f <- e[[1L]]
  if (is.symbol(f)) {
    fn <- as.character(f)
    if (fn %in% ASSIGN_OPS && length(e) == 3L) {
      if (fn %in% c("->", "->>")) { lhs <- e[[3L]]; rhs <- e[[2L]] } else { lhs <- e[[2L]]; rhs <- e[[3L]] }
      # RHS first (uses), then LHS (assignment)
      walk(rhs, acc, locals)
      lr <- lhs_root(lhs)
      if (!is.na(lr$root)) {
        if (lr$modify) {
          # modifying an existing object reads it first
          if (!(lr$root %in% locals)) acc$uses <- c(acc$uses, lr$root)
          # index/attr arguments may use symbols too
          if (is.call(lhs) && !(as.character(lhs[[1L]])[1L] %in% c("$", "@")))
            walk_elems(as.list(lhs), acc, locals, from = 3L)
        }
        acc$assigned <- c(acc$assigned, lr$root)
        acc$assigned_text <- c(acc$assigned_text, lr$text)
      } else {
        walk(lhs, acc, locals)
      }
      return(invisible())
    }
    if (fn == "function") {
      fml <- e[[2L]]
      fnames <- if (is.null(fml)) character(0) else names(fml)
      # defaults may reference globals
      if (!is.null(fml)) walk_elems(as.list(fml), acc, c(locals, fnames))
      body_acc <- new_acc()
      walk(e[[3L]], body_acc, c(locals, fnames))
      inner_locals <- unique(c(fnames, body_acc$assigned))
      acc$uses <- c(acc$uses, setdiff(body_acc$uses, inner_locals))
      acc$funs <- c(acc$funs, setdiff(body_acc$funs, inner_locals))
      acc$ns <- c(acc$ns, body_acc$ns); acc$nums <- c(acc$nums, body_acc$nums)
      acc$strs <- c(acc$strs, body_acc$strs); acc$rng <- c(acc$rng, body_acc$rng)
      acc$analytic <- c(acc$analytic, body_acc$analytic); acc$io <- c(acc$io, body_acc$io)
      acc$io_text <- c(acc$io_text, body_acc$io_text); acc$pkgs <- c(acc$pkgs, body_acc$pkgs)
      return(invisible())
    }
    if (fn == "local" && length(e) >= 2L) {
      # local({ ... }) scopes its assignments; only its free symbols reach the chunk
      body_acc <- new_acc()
      walk(e[[2L]], body_acc, locals)
      acc$uses <- c(acc$uses, setdiff(body_acc$uses, body_acc$assigned))
      acc$funs <- c(acc$funs, "local", setdiff(body_acc$funs, body_acc$assigned))
      acc$ns <- c(acc$ns, body_acc$ns); acc$nums <- c(acc$nums, body_acc$nums)
      acc$strs <- c(acc$strs, body_acc$strs); acc$rng <- c(acc$rng, body_acc$rng)
      acc$analytic <- c(acc$analytic, body_acc$analytic); acc$io <- c(acc$io, body_acc$io)
      acc$io_text <- c(acc$io_text, body_acc$io_text); acc$pkgs <- c(acc$pkgs, body_acc$pkgs)
      return(invisible())
    }
    if (fn == "for" && length(e) == 4L) {
      acc$assigned <- c(acc$assigned, as.character(e[[2L]]))
      acc$assigned_text <- c(acc$assigned_text, paste0("for-var ", as.character(e[[2L]])))
      walk(e[[3L]], acc, locals); walk(e[[4L]], acc, c(locals, as.character(e[[2L]])))
      return(invisible())
    }
    if (fn == "assign" && length(e) >= 3L && is.character(e[[2L]])) {
      acc$assigned <- c(acc$assigned, e[[2L]])
      acc$assigned_text <- c(acc$assigned_text, paste0("assign(\"", e[[2L]], "\")"))
      walk_elems(as.list(e), acc, locals, from = 3L)
      return(invisible())
    }
    if (fn %in% c("library", "require", "requireNamespace", "loadNamespace") && length(e) >= 2L) {
      acc$funs <- c(acc$funs, fn)
      acc$pkgs <- c(acc$pkgs, paste(deparse(e[[2L]]), collapse = ""))
      return(invisible())
    }
    if (fn %in% c("::", ":::") && length(e) == 3L) {
      acc$ns <- c(acc$ns, paste0(deparse(e[[2L]]), fn, deparse(e[[3L]])))
      return(invisible())
    }
    if (fn %in% c("$", "@") && length(e) == 3L) {
      walk(e[[2L]], acc, locals)   # the name after $ is not a free symbol
      return(invisible())
    }
    if (!(fn %in% locals)) acc$funs <- c(acc$funs, fn)
    if (fn %in% RNG_FUNS) acc$rng <- c(acc$rng, fn)
    if (fn %in% ANALYTIC_FUNS) acc$analytic <- c(acc$analytic, fn)
    if (fn %in% IO_FUNS) {
      acc$io <- c(acc$io, fn)
      acc$io_text <- c(acc$io_text, paste(deparse(e, width.cutoff = 500L), collapse = " "))
    }
    walk_elems(as.list(e), acc, locals, from = 2L)
    return(invisible())
  }
  # call whose function is itself an expression, e.g. pkg::f(x) or obj$f(x)
  if (is.call(f) && is.symbol(f[[1L]]) && as.character(f[[1L]]) %in% c("::", ":::")) {
    q <- paste0(deparse(f[[2L]]), as.character(f[[1L]]), deparse(f[[3L]]))
    acc$ns <- c(acc$ns, q); acc$funs <- c(acc$funs, q)
    bare <- deparse(f[[3L]])
    if (bare %in% RNG_FUNS) acc$rng <- c(acc$rng, q)
    if (bare %in% ANALYTIC_FUNS) acc$analytic <- c(acc$analytic, q)
    if (bare %in% IO_FUNS) { acc$io <- c(acc$io, q); acc$io_text <- c(acc$io_text, paste(deparse(e, width.cutoff = 500L), collapse = " ")) }
    walk_elems(as.list(e), acc, locals, from = 2L)
    return(invisible())
  }
  walk_elems(as.list(e), acc, locals, from = 1L)
  invisible()
}


analyse_chunk <- function(code) {
  res <- list(parse_error = NA_character_, free = character(0), assigned = character(0),
              assigned_text = character(0), funs = character(0), ns = character(0),
              nums = character(0), rng = character(0), analytic = character(0),
              io_text = character(0), pkgs = character(0), n_expr = 0L)
  if (!length(code) || !any(nzchar(trimws(code)))) return(res)
  exprs <- tryCatch(parse(text = code, keep.source = FALSE),
                    error = function(e) { res$parse_error <<- conditionMessage(e); NULL })
  if (is.null(exprs)) return(res)
  res$n_expr <- length(exprs)
  defined <- character(0)
  acc_all <- new_acc()
  free <- character(0)
  for (ex in as.list(exprs)) {
    a <- new_acc()
    walk(ex, a)
    # free = used in this expression and not yet defined by an earlier expression
    # (a symbol both used and assigned inside ONE expression, e.g. x <- x + 1, counts as free)
    free <- c(free, setdiff(a$uses, defined))
    defined <- unique(c(defined, a$assigned))
    for (nm in ls(acc_all)) assign(nm, c(get(nm, envir = acc_all), get(nm, envir = a)), envir = acc_all)
  }
  res$free <- unique(free)
  res$assigned <- unique(acc_all$assigned)
  res$assigned_text <- unique(acc_all$assigned_text)
  fl <- unique(acc_all$funs)
  res$funs <- fl[grepl("^[A-Za-z.][A-Za-z0-9._]*$", fl) | grepl(":::?", fl)]
  res$funs <- setdiff(res$funs, c("if", "for", "while", "repeat", "break", "next", "function"))
  res$pkgs <- unique(acc_all$pkgs)
  res$ns <- unique(acc_all$ns)
  res$nums <- acc_all$nums
  res$rng <- unique(acc_all$rng)
  res$analytic <- unique(acc_all$analytic)
  res$io_text <- unique(acc_all$io_text)
  res
}

for (k in seq_along(chunks)) {
  chunks[[k]]$an <- if (chunks[[k]]$lang == "r") analyse_chunk(chunks[[k]]$code) else NULL
}

## ---- 4. dependency graph ---------------------------------------------------
# a symbol used freely in chunk k is attributed to the LAST earlier R chunk that assigned it
edges <- data.frame(chunk = character(0), depends_on = character(0), symbols = character(0),
                    stringsAsFactors = FALSE)
external <- list()
for (k in seq_along(chunks)) {
  an <- chunks[[k]]$an
  if (is.null(an)) next
  # function names that are not defined anywhere in the document are base/package calls;
  # only symbols used as *values* or as functions defined by an earlier chunk make edges
  cand <- unique(c(an$free, an$funs))
  ext <- character(0)
  by_src <- list()
  for (s in cand) {
    src <- NA_integer_
    if (k > 1L) for (j in (k - 1L):1L) {
      if (!is.null(chunks[[j]]$an) && s %in% chunks[[j]]$an$assigned) { src <- j; break }
    }
    if (is.na(src)) { if (s %in% an$free) ext <- c(ext, s) } else {
      key <- chunks[[src]]$label
      by_src[[key]] <- c(by_src[[key]], s)
    }
  }
  for (key in names(by_src))
    edges <- rbind(edges, data.frame(chunk = chunks[[k]]$label, depends_on = key,
                                     symbols = paste(sort(unique(by_src[[key]])), collapse = ", "),
                                     stringsAsFactors = FALSE))
  external[[chunks[[k]]$label]] <- sort(unique(ext))
}

## ---- 5. inline code in prose -------------------------------------------------
in_chunk <- rep(FALSE, n)
for (ch in chunks) in_chunk[ch$start:ch$end] <- TRUE
prose <- lines[!in_chunk]
inline <- unlist(regmatches(prose, gregexpr("`r [^`]+`", prose)))
inline_code <- sub("^`r ", "", sub("`$", "", inline))
inline_syms <- character(0)
inline_bad <- character(0)
for (ic in inline_code) {
  ex <- tryCatch(parse(text = ic, keep.source = FALSE), error = function(e) NULL)
  if (is.null(ex)) { inline_bad <- c(inline_bad, ic); next }
  a <- new_acc(); for (e1 in as.list(ex)) walk(e1, a)
  inline_syms <- c(inline_syms, a$uses)
}
inline_tab <- if (length(inline_syms)) sort(table(inline_syms), decreasing = TRUE) else NULL

## ---- 6. write the report -----------------------------------------------------
md <- character(0)
p <- function(...) md <<- c(md, paste0(...))
esc <- function(x) gsub("|", "\\|", x, fixed = TRUE)
p("# Static chunk index — `", basename(qmd_path), "`")
p("")
p("Generated by `qmd_chunk_index.R` (base R, static parse only; nothing evaluated).  ")
p("File: `", qmd_path, "` — ", n, " lines, ", length(chunks), " fenced code chunks (",
  sum(vapply(chunks, function(c) c$lang == "r", logical(1))), " R).  ")
p("Line numbers are valid only at the commit recorded in the provenance header of the report that embeds this index.")
p("")
p("## 1. Chunk index (document order)")
p("")
p("| # | label | lang | lines | n expr | options | assigns / modifies | free symbols (inputs) | pkg::calls | RNG / draw | analytic | parse |")
p("|---|---|---|---|---|---|---|---|---|---|---|---|")
for (k in seq_along(chunks)) {
  ch <- chunks[[k]]; an <- ch$an
  p("| ", k, " | `", esc(ch$label), "` | ", ch$lang, " | ", ch$start, "–", ch$end, " | ",
    if (is.null(an)) "" else an$n_expr, " | ", esc(paste(ch$opts, collapse = "; ")), " | ",
    if (is.null(an)) "" else esc(paste(sort(an$assigned), collapse = ", ")), " | ",
    if (is.null(an)) "" else esc(paste(sort(an$free), collapse = ", ")), " | ",
    if (is.null(an)) "" else esc(paste(sort(an$ns), collapse = ", ")), " | ",
    if (is.null(an)) "" else esc(paste(sort(an$rng), collapse = ", ")), " | ",
    if (is.null(an)) "" else esc(paste(sort(an$analytic), collapse = ", ")), " | ",
    if (is.null(an)) "n/a" else if (is.na(an$parse_error)) "ok" else paste0("ERROR: ", esc(gsub("\\s+", " ", an$parse_error))),
    if (isTRUE(ch$unterminated)) " (UNTERMINATED FENCE)" else "", " |")
}
p("")
p("## 2. Chunk dependency edges (free symbol → last earlier chunk that assigned it)")
p("")
if (nrow(edges)) {
  p("| chunk | depends on | via symbols |"); p("|---|---|---|")
  for (r in seq_len(nrow(edges))) p("| `", esc(edges$chunk[r]), "` | `", esc(edges$depends_on[r]), "` | ", esc(edges$symbols[r]), " |")
} else p("(no edges)")
p("")
p("```mermaid"); p("flowchart TD")
for (k in seq_along(chunks)) if (chunks[[k]]$lang == "r") p("  c", k, "[\"", gsub("\"", "'", chunks[[k]]$label), "\"]")
if (nrow(edges)) for (r in seq_len(nrow(edges))) {
  from <- which(vapply(chunks, function(c) c$label == edges$depends_on[r], logical(1)))[1L]
  to   <- which(vapply(chunks, function(c) c$label == edges$chunk[r], logical(1)))[1L]
  p("  c", from, " --> c", to)
}
p("```")
p("")
p("## 3. Symbols used but never assigned by any earlier chunk (external: package, base, or defined in prose/YAML)")
p("")
p("| chunk | external symbols |"); p("|---|---|")
for (key in names(external)) if (length(external[[key]])) p("| `", esc(key), "` | ", esc(paste(external[[key]], collapse = ", ")), " |")
p("")
p("## 4. Assignment targets in full (including in-place modifications such as `x$y <- …`)")
p("")
for (k in seq_along(chunks)) {
  an <- chunks[[k]]$an
  if (is.null(an) || !length(an$assigned_text)) next
  p("- `", esc(chunks[[k]]$label), "`: ", esc(paste(an$assigned_text, collapse = " ; ")))
}
p("")
p("## 5. Numeric literals by chunk (raw; classify their roles by hand)")
p("")
p("| chunk | literal | count |"); p("|---|---|---|")
for (k in seq_along(chunks)) {
  an <- chunks[[k]]$an
  if (is.null(an) || !length(an$nums)) next
  tb <- table(an$nums)
  for (v in names(tb)) p("| `", esc(chunks[[k]]$label), "` | ", v, " | ", tb[[v]], " |")
}
p("")
p("## 6. Namespace-qualified calls (all chunks)")
p("")
allns <- sort(unique(unlist(lapply(chunks, function(c) if (is.null(c$an)) NULL else c$an$ns))))
if (length(allns)) for (q in allns) {
  where <- vapply(chunks, function(c) !is.null(c$an) && q %in% c$an$ns, logical(1))
  p("- `", q, "` — chunks: ", paste(vapply(chunks[where], function(c) c$label, character(1)), collapse = ", "))
} else p("(none)")
p("")
p("## 6b. Packages attached / loaded (library, require, requireNamespace, loadNamespace)")
p("")
for (k in seq_along(chunks)) {
  an <- chunks[[k]]$an
  if (is.null(an) || !length(an$pkgs)) next
  p("- `", esc(chunks[[k]]$label), "`: ", esc(paste(an$pkgs, collapse = ", ")))
}
p("")
p("## 7. File I/O calls (verbatim)")
p("")
for (k in seq_along(chunks)) {
  an <- chunks[[k]]$an
  if (is.null(an) || !length(an$io_text)) next
  for (t in an$io_text) p("- `", esc(chunks[[k]]$label), "`: `", esc(t), "`")
}
p("")
p("## 8. Functions called, by chunk (bare names; base/package/document-defined not distinguished here)")
p("")
for (k in seq_along(chunks)) {
  an <- chunks[[k]]$an
  if (is.null(an) || !length(an$funs)) next
  p("- `", esc(chunks[[k]]$label), "`: ", esc(paste(sort(an$funs), collapse = ", ")))
}
p("")
p("## 9. Objects read by inline `r …` code in prose")
p("")
p(length(inline_code), " inline code spans; ", length(inline_bad), " failed to parse.")
p("")
if (!is.null(inline_tab)) {
  p("| symbol | inline uses |"); p("|---|---|")
  for (s in names(inline_tab)) p("| `", esc(s), "` | ", inline_tab[[s]], " |")
}
if (length(inline_bad)) { p(""); p("Unparseable inline spans:"); for (b in inline_bad) p("- `", esc(b), "`") }
p("")
writeLines(md, out_path)
cat("wrote", out_path, "-", length(chunks), "chunks;",
    sum(vapply(chunks, function(c) !is.null(c$an) && !is.na(c$an$parse_error), logical(1))),
    "parse errors\n")
