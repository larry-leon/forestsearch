#!/usr/bin/env bash
# s2_text_gates.sh -- v2 sec. 5 text gates on the existing render, output only.
# Run from the repo root. Method: delete every <pre ...>...</pre> element from the HTML
# (26 echoed code blocks, class "sourceCode r", and 2 bare <pre> printed-output blocks),
# then convert with pandoc --wrap=none -t plain. article.txt is the whole page, for the record.
set -u
H=vignettes/articles/guohe-nested-continuum.html
Q=vignettes/articles/guohe-nested-continuum.qmd
R=dev/reports/guohe_article_2026-10-30/results
pandoc --wrap=none -t plain "$H" > "$R/article.txt"
Rscript -e 'a <- commandArgs(TRUE); h <- paste(readLines(a[1], warn = FALSE), collapse = "\n");
  n <- lengths(regmatches(h, gregexpr("<pre[ >]", h)));
  h2 <- gsub("(?s)<pre[ >].*?</pre>", "", h, perl = TRUE);
  stopifnot(!grepl("<pre[ >]", h2)); writeLines(h2, a[2]); cat("pre_elements_removed =", n, "\n")' \
  "$H" "$R/article_nopre.html"
pandoc --wrap=none -f html -t plain "$R/article_nopre.html" > "$R/article_output.txt"
rm "$R/article_nopre.html"

gate() { # id, regex (ERE), expectation (awk condition on n)
  local o w
  o=$(grep -oE "$2" "$R/article_output.txt" | wc -l); w=$(grep -oE "$2" "$R/article.txt" | wc -l)
  local ok; ok=$(awk -v n="$o" "BEGIN{print ($3) ? \"PASS\" : \"FAIL\"}")
  echo "$1 output_only = $o ; whole_page_v1 = $w ; diff = $((w - o)) ; expect: $3 ; $ok"
}
gate tg1_walkthrough_manuscript_settings "Replicate 1, beta_2 = 0.3, at the manuscript.s settings" "n==1"
gate tg2a_own_lower "This article.s run: one-sided 95% lower bounds" "n==1"
gate tg2b_own_select "This article.s run: what the argmax selects" "n==1"
gate tg2c_own_pair "This article.s run: the Bonferroni pair" "n==1"
gate tg3_wallclock_2_workers "Wall-clock of this article.s run on 2 workers" "n==1"
gate tg4a_0.952 "0\.952 \(0\.941, 0\.960\)" "n>=2"
gate tg4b_0.928 "0\.928 \(0\.916, 0\.939\)" "n>=1"
gate tg4c_S23_row1 "0\.945 / 0\.302 / \+0\.010" "n==1"
gate tg5a_NaN "NaN" "n==0"
o=$(grep -cP '\bNA \(' "$R/article_output.txt"); w=$(grep -cP '\bNA \(' "$R/article.txt")
echo "tg5b_NA_paren output_only = $o ; whole_page_v1 = $w ; diff = $((w - o)) ; expect: n==0 ; $([ "$o" -eq 0 ] && echo PASS || echo FAIL)"
gate tg5c_Inf "(^|[^[:alpha:]])-?Inf([^[:alpha:]]|$)" "n==0"
o=$(grep -c "Guo and He.s argmax and the multiplier resampling both carry the same candidate" "$R/article_output.txt")
w=$(grep -c "Guo and He.s argmax and the multiplier resampling both carry the same candidate" "$R/article.txt")
echo "tg6_selection_agreement output_only = $o ; whole_page_v1 = $w ; diff = $((w - o)) ; expect: n==1 ; $([ "$o" -eq 1 ] && echo PASS || echo FAIL)"
L=$(grep -n '^## Limits' "$Q" | cut -d: -f1); M=$(grep -n "^## Running the manuscript" "$Q" | cut -d: -f1)
grep -nE '^\s*[0-9]+[.)]( |$)' "$Q" > "$R/numbered_lines.txt"
awk -F: -v l="$L" -v m="$M" '{ if ($1>l && $1<m) i++; else x++ } END { printf "tg7_numbered_list (source) inside_Limits = %d ; outside = %d ; expect: outside==0 ; %s\n", i, x, (x+0==0) ? "PASS" : "FAIL" }' "$R/numbered_lines.txt"
o=$(grep -c 'Bias of the point estimate against beta(H-hat)' "$R/article_output.txt"); w=$(grep -c 'Bias of the point estimate against beta(H-hat)' "$R/article.txt")
echo "tg8a_fig_lower_caption output_only = $o ; whole_page_v1 = $w ; diff = $((w - o)) ; expect: n>=1 ; $([ "$o" -ge 1 ] && echo PASS || echo FAIL)"
o=$(grep -c 'The complement and the joint claim, by cell' "$R/article_output.txt"); w=$(grep -c 'The complement and the joint claim, by cell' "$R/article.txt")
echo "tg8b_fig_complement_caption output_only = $o ; whole_page_v1 = $w ; diff = $((w - o)) ; expect: n>=1 ; $([ "$o" -ge 1 ] && echo PASS || echo FAIL)"
