# compare.R -- legacy vs lean-fit re-render comparison (task cc_task_rerender_compare_2026-10-26, §4)
# Run from dev/reports/rerender_2026-10-26/:  Rscript compare.R
# Inputs:  legacy_<id>.html, new_<id>.html   (id in V1, V2, C, B)
# Outputs: outputs_<id>_{legacy,new}.txt, diff_<id>.txt, keyvalues_<id>.csv,
#          diff_counts.csv, deprecation_counts.csv, membership_B.txt
source("extract_lib.R")

ids <- c("V1", "V2", "C", "B")

# ---- prose lines (asis output rendered as markdown: B's gate comparisons) ----
extract_prose <- function(path) {
  h <- read_html(path)
  xml_remove(xml_find_all(h, "//pre|//script|//style|//div[contains(@class,'cell-output')]|//details|//nav"))
  t <- xml_text(xml_find_all(h, "//body//p|//body//li|//body//tr"))
  t <- trimws(gsub("\\s+", " ", t))
  t[nzchar(t)]
}

# ---- key-quantity specs: quantity, regex (one capture), source, occurrence ----
S <- function(q, re, k = 1L, src = "out") data.frame(quantity = q, re = re, k = k, src = src)
specs <- list(
  V1 = rbind(
    S("selected subgroup", "^Definition: (.*)$"),
    S("N", "^N: (\\S+)$"), S("HR", "^HR: (\\S+)$"), S("Pcons", "^Pcons: (\\S+)$"),
    S("candidates evaluated", "^Candidates evaluated: (\\d+)"),
    S("candidates passed", "^Candidates passed: (\\d+)"),
    S("candidate family (family_status)", "^Candidate family: (.*)$"),
    S("admission set", "^Admission set: (.*)$"),
    S("max subgroup HR", "^Maximum subgroup HR estimate: (\\S+)"),
    S("naive HR", "^naive HR on Hhat +(\\S+)"),
    S("de-biased HR", "^de-biased HR on Hhat +(\\S+)"),
    S("field one-sided lower (H)", "one-sided 95% lower bound on HR +(\\S+)"),
    S("field-s one-sided upper (Hc)", "one-sided 95% upper bound on HR +(\\S+)"),
    S("Bonferroni H lower", "Joint \\(Bonferroni\\): +H lower ([0-9.]+)"),
    S("Bonferroni Hc upper", "Joint \\(Bonferroni\\).*Hc upper ([0-9.]+)"),
    S("IJ two-sided (H)", "Two-sided \\(IJ, secondary\\): H (\\(.*\\))"),
    S("se_field", "se_field +([0-9.]+)"), S("se_field_s", "se_field_s +([0-9.]+)"),
    S("se_ij", "se_ij +([0-9.]+)"),
    S("p-hat(H)", "p-hat\\(H\\) = ([0-9.]+)"),
    S("top-3 re-selection mass", "top-3 re-selection mass: +(.*)$"),
    S("n_family", "over a family of (\\d+) candidates"),
    S("sens_H", "^sens_H += ([0-9.]+)"), S("ppv_H", "^ppv_H += ([0-9.]+)"),
    S("containment q10/q50/q90", "^containment q10 / q50 / q90 = (.*)$"),
    S("share containing all of Hhat", "share of draws containing all of Hhat = ([0-9.]+)"),
    S("recovery draws", "containing all of Hhat = .*over (\\d+) draws")),
  V2 = rbind(
    S("GRF subgroup", "^GRF subgroup: (.*)$"),
    S("Cox-LASSO selected", "^Cox-LASSO selected: (.*)$"),
    S("LASSO omitted", "^Omitted: (.*)$"),
    S("candidate factors", "^Candidate factors: (\\d+)"),
    S("combinations evaluated", "^Combinations evaluated: (\\d+)"),
    S("passed effect threshold", "^Passed effect threshold .*: (\\d+)"),
    S("candidates found", "^Found (\\d+) subgroup"),
    S("consistency evaluated", "^Evaluated: (\\d+)"),
    S("consistency passed", "^Evaluated: \\d+ +Passed: (\\d+)"),
    S("Pcons (selected)", "^% consistency criteria met= (\\S+)"),
    S("selected subgroup", "^Subgroup identified: (.*)$"),
    S("subgroup size", "^Subgroup size: (.*)$"),
    S("max HR estimate", "Maximum HR estimate = (\\S+)"),
    S("ITT row", "^ITT \\| (.*)$"),
    S("Questionable (H) row", "^Questionable1 \\| (.*)$"),
    S("Recommend (Hc) row", "^Recommend \\| (.*)$"),
    S("bootstrap success", "^Successful subgroup ID: +(.*)$"),
    S("SG_CIs H row (bias-corrected)", "^Qstnbl3 \\| (.*)$"),
    S("SG_CIs Hc row (bias-corrected)", "^Recmnd \\| (.*)$"),
    S("H bias-corrected HR", "Bias-corrected HR \\| (\\S+)", 1L),
    S("Hc bias-corrected HR", "Bias-corrected HR \\| (\\S+)", 2L),
    S("H bootstrap CV", "Bootstrap CV \\(%\\) \\| (\\S+)", 1L),
    S("Hc bootstrap CV", "Bootstrap CV \\(%\\) \\| (\\S+)", 2L),
    S("subgroup identified in bootstrap", "Subgroup identified in ([0-9.]+)% of bootstrap"),
    S("CV Sensitivity (H)", "^Sensitivity \\(H\\) \\| .* \\| (\\S+)$"),
    S("CV Sensitivity (Hc)", "^Sensitivity \\(Hc\\) \\| .* \\| (\\S+)$"),
    S("CV PPV (H)", "^PPV \\(H\\) \\| .* \\| (\\S+)$"),
    S("CV PPV (Hc)", "^PPV \\(Hc\\) \\| .* \\| (\\S+)$"),
    S("CV Any Found", "^Any Found \\| .* \\| (\\S+)$"),
    S("CV Exact Match", "^Exact Match \\| .* \\| (\\S+)$"),
    S("CV At Least 1", "^At Least 1 \\| .* \\| (\\S+)$")),
  C = rbind(
    S("N", "^Analysis N = (\\d+)"),
    S("ITT MD", "^unadjusted ITT mean difference = (.*)$"),
    S("selected subgroup", "^H-hat = (.*?);"),
    S("n(H-hat)", "n\\(H-hat\\) = (\\d+)"),
    S("T_obs (MD, oriented)", "T_obs = ([0-9.]+)"),
    S("Pcons", "p\\.consistency = ([0-9.]+)"),
    S("band line", "^max oriented MD (.*)$"),
    S("H Naive", "^Naive \\| (.*)$", 1L),
    S("H MR (IJ two-term)", "^MR \\(IJ two-term\\) \\| (.*)$", 1L),
    S("H MR (field)", "^MR \\(field\\) \\| (.*)$"),
    S("Hc Naive", "^Naive \\| (.*)$", 2L),
    S("Hc MR (IJ two-term)", "^MR \\(IJ two-term\\) \\| (.*)$", 2L),
    S("Hc MR (field, unstudentized)", "^MR \\(field, unstudentized\\) \\| (.*)$"),
    S("Hc MR (field-s)", "^MR \\(field-s\\) \\| (.*)$"),
    S("regime diagnostic", "Regime diagnostic: (.*)$"),
    S("joint subtitle", "^unstudentized: γ = (.*)$"),
    S("separate one-sided pair", "^Separate one-sided 95% field bounds \\(unstudentized\\) \\| (.*)$"),
    S("Bonferroni pair (unstudentized)", "^Bonferroni \\(γ = 0.025 each\\) \\(unstudentized\\) \\| (.*)$"),
    S("Bonferroni pair (field-s)", "^Bonferroni \\(γ = 0.025 each\\), field-s \\| (.*)$"),
    S("p-hat(H)", "^p̂\\(Ĥ\\) \\[selected: .*\\] \\| (.*)$"),
    S("top re-selected", "^Top re-selected: (.*)$"),
    S("2nd re-selected", "^2nd: (.*)$"), S("3rd re-selected", "^3rd: (.*)$"),
    S("selection rate; n_family", "^Selection rate .* \\| (.*)$"),
    S("field Lambda*-mean", "^Field Λ\\*-mean .* \\| (.*)$"),
    S("field lambda-SD ratios", "^Field λ-SD vs IJ SE.* \\| (.*)$"),
    S("field draw usage", "^Field draw usage .* \\| (.*)$"),
    S("complement fits; lambda-SDc ratio", "^Complement fits per gate call.* \\| (.*)$")),
  B = rbind(
    S("ITT OR line", "^OR: +(.*)$"),
    S("FS: combinations evaluated", "^Combinations evaluated: (\\d+)", 1L),
    S("FS: passed effect screen", "^Passed effect threshold .*: (\\d+)", 1L),
    S("FS: near-duplicates removed", "^Removed (\\d+) near-duplicate", 1L),
    S("FS: consistency pool (Evaluated)", "^Evaluated: (\\d+)", 1L),
    S("FS: consistency passed", "^Evaluated: \\d+ +Passed: (\\d+)", 1L),
    S("FS: on frontier", "On frontier: (\\d+) +In band", 1L),
    S("FS: selected subgroup (search print)", "^Subgroup identified: (.*)$", 1L),
    S("FS: MR harm confirmation", "^MR harm confirmation: (.*)$", 1L),
    S("FS: candidate subgroups found (table)", "Candidate subgroups found = (\\d+)", 1L),
    S("FS: max OR estimate", "Maximum OR estimate = (\\S+)", 1L),
    S("FS: ITT row", "^ITT \\| (.*)$", 1L),
    S("FS: H row (naive)", "^Questionable1 \\| (.*)$", 1L),
    S("FS: Hc row (naive)", "^Recommend \\| (.*)$", 1L),
    S("FS: identified subgroup (H)", "^Identified subgroup \\(H\\): (.*)$", 1L),
    S("FS: subgroup size", "^Subgroup size: (.*)$", 1L),
    S("FS: Pareto selected row", "^★ \\| (.*)$", 1L),
    S("FS vs GRF agreement", "^Agreement: (.*?) Cohen"),
    S("FS vs GRF kappa", "Cohen.s κ: ([0-9.]+)"),
    S("FS: bootstrap success", "^Successful subgroup ID: +(.*)$", 1L),
    S("FS: SG_CIs H (bootstrap)", "^Qstnbl3 \\| (.*)$", 1L),
    S("FS: SG_CIs Hc (bootstrap)", "^Recmnd \\| (.*)$", 1L),
    S("FS: SG_CIs H (MR gate)", "^Qstnbl3 \\| (.*)$", 2L),
    S("FS: SG_CIs Hc (MR gate)", "^Recmnd \\| (.*)$", 2L),
    S("FS: H bias-corrected OR", "Bias-corrected OR \\| (\\S+)", 1L),
    S("FS: Hc bias-corrected OR", "Bias-corrected OR \\| (\\S+)", 2L),
    S("FS: H bootstrap CV", "Bootstrap CV \\(%\\) \\| (\\S+)", 1L),
    S("FS: boot exact match with original", "Exact match with original +(.*)$", 1L),
    S("FS: boot partial match", "Partial match \\(shared factor\\) +(.*)$", 1L),
    S("FS: gate de-biased OR (T2 | T1)", "^de-biased OR \\(95% CI\\) (.*)$", 1L, "prose"),
    S("FS: gate complement de-biased OR (T2 | T1)", "^complement de-biased OR \\(95% CI\\) (.*)$", 1L, "prose"),
    S("FS: candidate family (family_status)", "^Candidate family: ([0-9].*)$", 1L, "prose"),
    S("FS: scale cross-check", "^Scale cross-check.*: (.*)$", 1L, "prose"),
    S("FS: gate harm confirmed", "^Gate mr_harm_confirmed = ([A-Z]+)", 1L, "prose"),
    S("FS: gate selection bias (log-OR)", "selection bias ([-0-9.]+) \\(log-OR\\)", 1L, "prose"),
    S("FS: IJ SE vs robust SE", "(IJ SE [0-9.]+ vs robust SE [0-9.]+)", 1L, "prose"),
    S("FS: gate complement", "Complement \\(n = (.*)\\.?$", 1L, "prose"),
    S("DINA-screen: GRF cut", "^GRF subgroup: (.*)$", 1L),
    S("DINA-screen: combinations evaluated", "^Combinations evaluated: (\\d+)", 2L),
    S("DINA-screen: passed effect screen", "^Passed effect threshold .*: (\\d+)", 2L),
    S("DINA-screen: consistency pool", "^Evaluated: (\\d+)", 2L),
    S("DINA-screen: consistency passed", "^Evaluated: \\d+ +Passed: (\\d+)", 2L),
    S("DINA-screen: selected subgroup", "^Subgroup identified: (.*)$", 2L),
    S("DINA-screen: subgroup size", "^Subgroup size: (.*)$", 2L),
    S("standalone DINA: selection 1", "^DINA-selected subgroup: ([^{].*)$", 1L),
    S("standalone DINA: selection 2", "^DINA-selected subgroup: ([^{].*)$", 2L),
    S("standalone DINA: direction/N 1", "^Direction: (.*)$", 1L),
    S("standalone DINA: direction/N 2", "^Direction: (.*)$", 2L),
    S("FS-dina: candidates searched", "^Candidates searched: +(\\d+)"),
    S("FS-dina: candidates qualifying", "^Candidates qualifying .*: (\\d+)"),
    S("FS-dina: SELECTED", "^SELECTED: (.*)$"),
    S("FS-dina: selected subgroup (H)", "^DINA-selected subgroup \\(H\\): (.*)$"),
    S("FS-dina: subgroup size", "^Subgroup size: (.*)$", 3L),
    S("FS-dina: bootstrap success", "^Successful subgroup ID: +(.*)$", 2L),
    S("FS-dina: SG_CIs H (bootstrap)", "^Qstnbl3 \\| (.*)$", 3L),
    S("FS-dina: SG_CIs Hc (bootstrap)", "^Recmnd \\| (.*)$", 3L),
    S("FS-dina: gate de-biased OR (T2 | T1)", "^de-biased OR \\(95% CI\\) (.*)$", 2L, "prose"),
    S("FS-dina: gate complement de-biased OR (T2 | T1)", "^complement de-biased OR \\(95% CI\\) (.*)$", 2L, "prose"),
    S("FS-dina: gate candidate family", "candidate family = (\\d{2,})", 1L, "prose"),
    S("FS-dina: gate selection bias (log-OR)", "selection bias ([-0-9.]+) \\(log-OR\\)", 2L, "prose"),
    S("FS-grf: SELECTED", "^SELECTED \\(depth 2\\): (.*)$"),
    S("FS-grf: selected subgroup (H)", "^GRF-selected subgroup \\(H\\): (.*)$"),
    S("FS-grf: subgroup size", "^Subgroup size: (.*)$", 4L),
    S("FS-grf: bootstrap success", "^Successful subgroup ID: +(.*)$", 3L),
    S("FS-grf: SG_CIs H (bootstrap)", "^Qstnbl3 \\| (.*)$", 5L),
    S("FS-grf: SG_CIs Hc (bootstrap)", "^Recmnd \\| (.*)$", 5L),
    S("FS-grf: gate de-biased OR (T2 | T1)", "^de-biased OR \\(95% CI\\) (.*)$", 3L, "prose"),
    S("FS-grf: gate complement de-biased OR (T2 | T1)", "^complement de-biased OR \\(95% CI\\) (.*)$", 3L, "prose"),
    S("FS-grf: gate candidate family", "candidate family = (\\d{2,})", 2L, "prose"),
    S("FS-grf: gate selection bias (log-OR)", "selection bias ([-0-9.]+) \\(log-OR\\)", 3L, "prose"),
    S("summary: definition", "^Definition: (.*)$", 1L),
    S("summary: H size", "^Size: (.*)$", 1L),
    S("summary: Hc size", "^Size: (.*)$", 2L))
)

grab <- function(lines, re, k) {
  m <- regmatches(lines, regexec(re, lines, perl = TRUE))
  hits <- vapply(Filter(length, m), `[`, "", 2)
  if (length(hits) >= k) trimws(hits[k]) else NA_character_
}
clauses <- function(x) sort(trimws(strsplit(gsub("^\\s+|\\s+$", "", x), " & ", fixed = TRUE)[[1]]))
decimals <- function(x) {
  n <- regmatches(x, gregexpr("-?[0-9]+\\.[0-9]+", x))[[1]]
  if (!length(n)) "integer/text" else paste(unique(nchar(sub(".*\\.", "", n))), collapse = "/")
}
verdict <- function(a, b) {
  if (is.na(a) && is.na(b)) return("absent in both")
  if (is.na(a)) return("new only")
  if (is.na(b)) return("legacy only")
  if (identical(a, b)) return("identical")
  if (grepl(" & ", a) && grepl(" & ", b) && identical(clauses(a), clauses(b))) return("same clauses, reordered")
  "changed"
}

dcount <- list(); depr <- list()
for (id in ids) {
  fL <- sprintf("legacy_%s.html", id); fN <- sprintf("new_%s.html", id)
  if (!file.exists(fN)) { message(id, ": not rendered"); next }
  oL <- extract_outputs(fL); oN <- extract_outputs(fN)
  writeLines(oL, sprintf("outputs_%s_legacy.txt", id)); writeLines(oN, sprintf("outputs_%s_new.txt", id))
  d <- suppressWarnings(system2("diff", c("-u", sprintf("outputs_%s_legacy.txt", id), sprintf("outputs_%s_new.txt", id)), stdout = TRUE))
  writeLines(d, sprintf("diff_%s.txt", id))
  body <- d[-(1:2)]
  dcount[[id]] <- data.frame(id = id, lines_legacy = length(oL), lines_new = length(oN),
                             lines_removed = sum(grepl("^-", body)), lines_added = sum(grepl("^\\+", body)),
                             identical_lines = length(oL) - sum(grepl("^-", body)))
  pL <- extract_prose(fL); pN <- extract_prose(fN)
  sp <- specs[[id]]
  kv <- do.call(rbind, lapply(seq_len(nrow(sp)), function(i) {
    srcL <- if (sp$src[i] == "prose") pL else oL; srcN <- if (sp$src[i] == "prose") pN else oN
    a <- grab(srcL, sp$re[i], sp$k[i]); b <- grab(srcN, sp$re[i], sp$k[i])
    data.frame(quantity = sp$quantity[i], legacy = a, new = b, result = verdict(a, b),
               printed_precision = if (is.na(a)) NA else decimals(a), source = sp$src[i])
  }))
  write.csv(kv, sprintf("keyvalues_%s.csv", id), row.names = FALSE)
  cnt <- function(f, re) { x <- readLines(f, warn = FALSE); c(lines = sum(grepl(re, x)), occurrences = sum(lengths(regmatches(x, gregexpr(re, x))))) }
  depr[[id]] <- data.frame(id = id, rbind(c(legacy_grep_c = cnt(fL, "deprecated")[["lines"]], legacy_occ = cnt(fL, "deprecated")[["occurrences"]],
                                            new_grep_c = cnt(fN, "deprecated")[["lines"]], new_occ = cnt(fN, "deprecated")[["occurrences"]],
                                            new_vi_grf_min_warn_occ = cnt(fN, "vi\\.grf\\.min.{0,40}deprecated|deprecated.{0,80}vi\\.grf\\.min")[["occurrences"]])))
  cat(sprintf("%s: legacy %d / new %d lines; -%d +%d; key quantities %d (identical %d)\n", id, length(oL), length(oN),
              dcount[[id]]$lines_removed, dcount[[id]]$lines_added, nrow(kv), sum(kv$result == "identical")))
}
write.csv(do.call(rbind, dcount), "diff_counts.csv", row.names = FALSE)
write.csv(do.call(rbind, depr), "deprecation_counts.csv", row.names = FALSE)
print(do.call(rbind, depr))
