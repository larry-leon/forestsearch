#!/usr/bin/env bash
# REPORT_cran_prep_code_fixes_2026-11-12_final_gates.sh -- "After the fixes" step 8.
# Usage (from anywhere): bash REPORT_cran_prep_code_fixes_2026-11-12_final_gates.sh [staged]
# Recomputes the final gates from git alone (no R):
#   G1 git status --porcelain is empty (with "staged": only the report files may be staged)
#   G2 both sibling repositories are at their start SHAs, on main, with empty status, and no
#      file in them (outside .git) is newer than the start commit
#   G3 git diff <start> --stat lists only files the report names, nothing under vignettes/
#      or quarto/
#   G4 each fix commit touches only its own files (code, its man page, its test file)
#   G5 the branch head is not on any remote-tracking ref (nothing pushed); no merge commit
set -u
REPO=$HOME/Documents/GitHub/forestsearch
START=76a62a71
REPORT=$REPO/dev/reports/REPORT_cran_prep_code_fixes_2026-11-12.md
DETAIL=$REPO/dev/reports/REPORT_cran_prep_code_fixes_2026-11-12_detail.md
MODE=${1:-final}
cd "$REPO" || exit 2
ok=1
say() { printf '%s\n' "$*"; }

say "== G1 working tree ($MODE)"
st=$(git status --porcelain)
if [ "$MODE" = staged ]; then
  bad=$(printf '%s\n' "$st" | grep -v '^A  dev/reports/REPORT_cran_prep_code_fixes_2026-11-12' | grep -v '^$' || true)
  [ -z "$bad" ] && say "PASS: only report files staged ($(printf '%s\n' "$st" | grep -c '^A '))" || { say "FAIL:"; say "$bad"; ok=0; }
else
  [ -z "$st" ] && say "PASS: git status --porcelain is empty" || { say "FAIL:"; say "$st"; ok=0; }
fi

say "== G2 sibling repositories"
start_time=$(git log -1 --format=%cI $START)
for r in fs-post-selection:5d165acce104a336e11093f1881cfdb572348462 fs-glms-interpretable:dc94a17e0e93653d5d4a9fb83bacca93bb8f87ea; do
  name=${r%%:*}; sha=${r#*:}; d=$HOME/Documents/GitHub/$name
  b=$(git -C "$d" branch --show-current); h=$(git -C "$d" rev-parse HEAD); s=$(git -C "$d" status --porcelain | wc -l)
  newer=$(find "$d" -path "$d/.git" -prune -o -type f -newermt "$start_time" -print | wc -l)
  if [ "$b" = main ] && [ "$h" = "$sha" ] && [ "$s" = 0 ] && [ "$newer" = 0 ]; then say "PASS: $name main $h, status empty, 0 files newer than $start_time"
  else say "FAIL: $name branch=$b head=$h status_lines=$s newer=$newer"; ok=0; fi
done

say "== G3 git diff $START --stat (committed + staged)"
files=$(git diff --cached --name-only $START)
printf '%s\n' "$files" | sed 's/^/  /'
n=$(printf '%s\n' "$files" | grep -c .)
say "  $n files"
vq=$(printf '%s\n' "$files" | grep -E '^(vignettes|quarto)/' || true)
[ -z "$vq" ] && say "PASS: nothing under vignettes/ or quarto/" || { say "FAIL: $vq"; ok=0; }
if [ -f "$REPORT" ]; then
  miss=""
  for f in $files; do grep -qF "$f" "$REPORT" "$DETAIL" 2>/dev/null || grep -qF "$(basename "$f")" "$REPORT" || miss="$miss $f"; done
  [ -z "$miss" ] && say "PASS: every path is named in the report" || { say "FAIL: not named:$miss"; ok=0; }
fi

say "== G4 per-commit scope"
declare -A ALLOW=(
  [0dc99947]="R/summarize_bootstrap_subgroups.R tests/testthat/test-bootstrap-summary-bins.R"
  [1c2565b3]="R/forestsearch_cross_validation.R tests/testthat/test-cv-fold-na-assignments.R"
  [3612d0ba]="R/forestsearch_cross_validation.R tests/testthat/test-cv-exact-match.R"
  [4aaab2f0]="R/test_hte_crump.R tests/testthat/test-hte-crump-confounders.R"
  [f5b65516]="R/forestsearch_main.R R/get_FSdata_helpers.R man/lasso_selection.Rd tests/testthat/test-lasso-intercept.R"
  [65831b54]="R/fpr_calibration.R man/fpr_calibration.Rd tests/testthat/test-fpr-calibration-fallback.R"
  [ea5f82d9]="R/get_fsdata.R tests/testthat/test-duplicate-cut-names.R"
  [a30d6e7c]="R/mrct_case_study.R R/mrct_simulation.R man/mrct_case_study.Rd man/mrct_region_sims.Rd tests/testthat/test-mrct-default-confounders.R"
  [f2e8de47]="R/grf_subg_harm_glm.R tests/testthat/test-grf-count-orientation.R"
  [25ecaa0f]="R/mr_estimates_table.R R/plot_sg_glm_outcomes.R R/plot_subgroup_results_forestplot.R tests/testthat/test-display-orientation.R"
  [89ad9c19]="R/glm_effect_estimators.R tests/testthat/test-rd-standard-error.R"
  [685e0c5d]="R/forestsearch_main.R tests/testthat/test-dina-gaussian-floor.R"
  [a4cce90a]="R/forestsearch_cross_validation.R R/forestsearch_helpers.R R/run_simulation_analysis.R tests/testthat/test-missing-membership-messages.R"
  [61c15088]="R/forestsearch_main.R tests/testthat/test-dina-screening-tau-sign.R"
  [88302212]="DESCRIPTION NEWS.md"
  [a577b59d]="NEWS.md"
  [76a62a71]="dev/tasks/TASK_cran_prep_code_fixes_2026-11-12.md"
)
for c in $(git rev-list --reverse $START^..HEAD); do
  s=$(git rev-parse --short=8 $c); got=$(git show --name-only --format= $c | sort | tr '\n' ' ' | sed 's/ $//')
  if [ -n "${ALLOW[$s]+x}" ]; then
    want=$(printf '%s\n' ${ALLOW[$s]} | sort | tr '\n' ' ' | sed 's/ $//')
    [ "$got" = "$want" ] && say "PASS: $s $got" || { say "FAIL: $s got [$got] want [$want]"; ok=0; }
  else
    extra=$(printf '%s\n' $got | grep -v '^dev/reports/REPORT_cran_prep_code_fixes_2026-11-12' || true)
    [ -z "$extra" ] && say "PASS: $s (report commit) $got" || { say "FAIL: $s unexpected files: $extra"; ok=0; }
  fi
done

say "== G5 nothing pushed or merged"
on_remote=$(git branch -r --contains HEAD 2>/dev/null | wc -l)
merges=$(git rev-list --merges $START..HEAD | wc -l)
[ "$on_remote" = 0 ] && [ "$merges" = 0 ] && say "PASS: HEAD on no remote-tracking ref; 0 merge commits since $START" || { say "FAIL: remote refs containing HEAD=$on_remote merges=$merges"; ok=0; }

say ""
[ $ok = 1 ] && say "FINAL GATES: PASS" || say "FINAL GATES: FAIL"
