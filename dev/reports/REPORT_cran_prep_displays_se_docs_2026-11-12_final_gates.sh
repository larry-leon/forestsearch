#!/usr/bin/env bash
# REPORT_cran_prep_displays_se_docs_2026-11-12_final_gates.sh -- "After the parts" step 8 (part 5).
# Usage (from anywhere): bash REPORT_cran_prep_displays_se_docs_2026-11-12_final_gates.sh [staged]
# Recomputes the final gates from git alone (no R):
#   G1 git status --porcelain is empty (with "staged": only the report files may be staged)
#   G2 both sibling repositories are at their start SHAs, on main, with empty status, and no
#      file in them (outside .git) is newer than the start commit
#   G3 git diff <start> --stat lists only files the report names; under vignettes/ only the two
#      C7 files; nothing under quarto/
#   G4 each commit since the start touches only its own files (lists taken from the gated patches)
#   G5 the branch head is not on any remote-tracking ref (nothing pushed); no merge commit
#   G6 the documentation gates, recomputed per commit against its parent: comment-only for every
#      R/ and tests/ file of C4, C5, C6; C7 changes only prose lines outside code chunks
set -u
REPO=$HOME/Documents/GitHub/forestsearch
START=b06de4e8
PFX=dev/reports/REPORT_cran_prep_displays_se_docs_2026-11-12
REPORT=$REPO/$PFX.md
DETAIL=${REPORT%.md}_detail.md
MODE=${1:-final}
cd "$REPO" || exit 2
ok=1
say() { printf '%s\n' "$*"; }

say "== G1 working tree ($MODE)"
st=$(git status --porcelain)
if [ "$MODE" = staged ]; then
  bad=$(printf '%s\n' "$st" | grep -v "^A  $PFX" | grep -v '^$' || true)
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
q=$(printf '%s\n' "$files" | grep -E '^quarto/' || true)
v=$(printf '%s\n' "$files" | grep -E '^vignettes/' | sort | tr '\n' ' ' | sed 's/ $//')
[ -z "$q" ] && say "PASS: nothing under quarto/" || { say "FAIL: $q"; ok=0; }
[ "$v" = "vignettes/articles/biomarker_effects.qmd vignettes/articles/paper_simulations.qmd" ] && say "PASS: under vignettes/ only the two C7 files" || { say "FAIL: vignettes/ files: $v"; ok=0; }
if [ -f "$REPORT" ]; then
  miss=""
  for f in $files; do grep -qF "$f" "$REPORT" "$DETAIL" 2>/dev/null || grep -qF "$(basename "$f")" "$REPORT" || miss="$miss $f"; done
  [ -z "$miss" ] && say "PASS: every path is named in the report" || { say "FAIL: not named:$miss"; ok=0; }
fi

say "== G4 per-commit scope"
declare -A ALLOW=(
  [2984dbee]="R/forestsearch_main.R tests/testthat/test-threshold-remap-messages.R"
  [f4c7dd2d]="R/summary_utility_functions.R tests/testthat/test-sg-tables-orientation.R"
  [a2cd36db]="R/forestsearch_cross_validation.R tests/testthat/test-cv-tables-orientation.R"
  [fb21396f]="R/frontier_cis.R tests/testthat/test-frontier-cis-orientation.R"
  [8e87cf30]="R/glm_effect_estimators.R tests/testthat/test-rd-adjusted-se.R"
  [9ec67b7e]="R/format_subgroup_summary_tables.R tests/testthat/test-consistency-highlight-labels.R"
  [c3323a59]="R/summary_utility_functions.R tests/testthat/test-grf-footnote-scale.R"
  [ba7c4ea2]="man/calibrate_k_treat.Rd man/compute_dgm_cde.Rd man/CV_sgs.Rd man/cv_summary_tables.Rd man/dot-glm_resample_supported.Rd man/find_k_inter_for_target_hr.Rd man/forestsearch_KfoldOut.Rd man/forestsearch.Rd man/generate_bootstrap_with_noise.Rd man/make_effect_estimator.Rd man/mr_estimates_table.Rd man/plot_sg_glm_outcomes.Rd man/plot_subgroup_results_forestplot.Rd man/setup_gbsg_dgm.Rd man/sg_tables.Rd man/summarize_bootstrap_subgroups.Rd man/test_hte_from_forestsearch.Rd R/calibrate_k_treat.R R/cv_summary_tables.R R/dina_subgroup.R R/find_k_inter_main.R R/forestsearch_cross_validation.R R/forestsearch_main.R R/get_fsdata.R R/glm_effect_estimators.R R/mr_estimates_table.R R/oc_analyses.R R/plot_sg_glm_outcomes.R R/plot_subgroup_results_forestplot.R R/setup_gbsg_dgm.R R/subgroup_consistency_helpers.R R/summarize_bootstrap_subgroups.R R/summary_utility_functions.R R/synthetic_data_perturbation.R R/test_hte_crump.R"
  [5073fd67]="man/evaluate_consistency_twostage.Rd man/evaluate_subgroup_consistency.Rd man/format_bootstrap_table.Rd man/fs_fdr_report.Rd man/pareto_frontier_table.Rd man/plot.fs_forestplot.Rd man/plot.fs_sg_plot.Rd man/plot_sg_distribution.Rd man/print.fs_forestplot.Rd man/print.fs_kfold.Rd man/print.fs_sg_plot.Rd man/print.fs_tenfold.Rd man/run_subgroup_sims.Rd man/sg_consistency_out.Rd man/sg_tables.Rd man/summarize_simulation_results.Rd R/bootstrap_summaries_helpers.R R/forestsearch_cross_validation.R R/fs_fdr_report.R R/oc_analyses.R R/pareto_frontier_table.R R/plot_sg_distribution.R R/plot_sg_results.R R/render_forestplot.R R/run_subgroup_sims.R R/subgroup_consistency_helpers.R R/summary_utility_functions.R"
  [ffc55486]="tests/testthat/test-bootstrap-summary-bins.R tests/testthat/test-dina-gaussian-floor.R tests/testthat/test-display-orientation.R tests/testthat/test-duplicate-cut-names.R tests/testthat/test-grf-count-orientation.R tests/testthat/test-hte-crump-confounders.R tests/testthat/test-lasso-intercept.R tests/testthat/test-mrct-default-confounders.R tests/testthat/test-rd-standard-error.R"
  [835fac1b]="vignettes/articles/biomarker_effects.qmd vignettes/articles/paper_simulations.qmd"
  [fb80dd8f]="DESCRIPTION NEWS.md"
  [2a8be3bc]="NEWS.md"
  [b06de4e8]="dev/tasks/TASK_cran_prep_displays_se_docs_2026-11-12.md"
)
for c in $(git rev-list --reverse $START^..HEAD); do
  s=$(git rev-parse --short=8 $c); got=$(git show --name-only --format= $c | sort | tr '\n' ' ' | sed 's/ $//')
  if [ -n "${ALLOW[$s]+x}" ]; then
    want=$(printf '%s\n' ${ALLOW[$s]} | sort | tr '\n' ' ' | sed 's/ $//')
    [ "$got" = "$want" ] && say "PASS: $s ($(printf '%s\n' $got | wc -l) files)" || { say "FAIL: $s got [$got] want [$want]"; ok=0; }
  else
    extra=$(printf '%s\n' $got | grep -v "^$PFX" || true)
    [ -z "$extra" ] && say "PASS: $s (report commit) $got" || { say "FAIL: $s unexpected files: $extra"; ok=0; }
  fi
done

say "== G5 nothing pushed or merged"
on_remote=$(git branch -r --contains HEAD 2>/dev/null | wc -l)
merges=$(git rev-list --merges $START..HEAD | wc -l)
[ "$on_remote" = 0 ] && [ "$merges" = 0 ] && say "PASS: HEAD on no remote-tracking ref; 0 merge commits since $START" || { say "FAIL: remote refs containing HEAD=$on_remote merges=$merges"; ok=0; }

say "== G6 documentation gates (each commit against its parent)"
for c in ba7c4ea2 5073fd67 ffc55486; do
  bad=""
  for f in $(git diff --name-only $c^ $c -- R tests); do
    out=$(git diff -U0 $c^ $c -- "$f" | grep -E '^[+-]' | grep -vE '^(\+\+\+|---)' | grep -vE "^[+-][[:space:]]*#")
    [ -n "$out" ] && bad="$bad $f"
  done
  nf=$(git diff --name-only $c^ $c -- R tests | wc -l)
  [ -z "$bad" ] && say "PASS: $c comment-only in $nf R/ and tests/ files" || { say "FAIL: $c non-comment lines in:$bad"; ok=0; }
done
c7=$(python3 - "$REPO" 835fac1b <<'PY'
import re, subprocess, sys
repo, c = sys.argv[1], sys.argv[2]
def git(*a): return subprocess.run(['git', *a], capture_output=True, text=True, cwd=repo).stdout
def inchunk(lines):
    st, inside = [], False
    for l in lines:
        if re.match(r'^\s*```', l): st.append(True); inside = not inside; continue
        st.append(inside)
    return st
bad = 0; n = 0
for f in git('diff', '--name-only', c + '^', c).split():
    new = git('show', f'{c}:{f}').splitlines(); old = git('show', f'{c}^:{f}').splitlines()
    ni, oi = inchunk(new), inchunk(old)
    for h in re.finditer(r'^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@', git('diff', '-U0', c + '^', c, '--', f), re.M):
        os_, oc, ns, nc = int(h[1]), int(h[2] or 1), int(h[3]), int(h[4] or 1)
        for i in range(ns, ns + nc): n += 1; bad += ni[i - 1]
        for i in range(os_, os_ + oc): n += 1; bad += oi[i - 1]
print(f'{n} {bad}')
PY
)
set -- $c7
[ "$2" = 0 ] && say "PASS: 835fac1b $1 changed lines (old + new sides), none inside a code chunk or on a fence" || { say "FAIL: 835fac1b $2 of $1 changed lines inside a code chunk"; ok=0; }

say ""
[ $ok = 1 ] && say "FINAL GATES: PASS" || say "FINAL GATES: FAIL"
