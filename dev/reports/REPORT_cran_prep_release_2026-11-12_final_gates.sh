#!/usr/bin/env bash
# Final gates -- TASK_cran_prep_release_2026-11-12, "Report and close" step 3.
# Usage: final_gates.sh <start sha>
# 1. git status --porcelain is empty in forestsearch.
# 2. Both sibling repositories are unchanged: same branch, same HEAD as at
#    Step 0, git status --porcelain empty.
# 3. git diff <start> --stat lists only the files the report names, and
#    nothing under vignettes/ or quarto/.
# 4. Nothing pushed, merged or tagged: no new tag, origin/cran-prep and
#    origin/master unchanged from Step 0 (2beb79c4), no merge commit since the
#    start SHA.
set -u
START=$1
REPO=$HOME/Documents/GitHub/forestsearch
ok=1
say() { printf '%s\n' "$*"; }
cd "$REPO" || exit 2
say "head: $(git rev-parse HEAD) ($(git branch --show-current))"
st=$(git status --porcelain)
say "G1 git status --porcelain empty: $([ -z "$st" ] && echo PASS || { ok=0; echo "FAIL"; printf '%s\n' "$st"; })"
for r in fs-post-selection:main:5d165acce104a336e11093f1881cfdb572348462 \
         fs-glms-interpretable:main:dc94a17e0e93653d5d4a9fb83bacca93bb8f87ea; do
  name=${r%%:*}; rest=${r#*:}; br=${rest%%:*}; sha=${rest#*:}
  d=$HOME/Documents/GitHub/$name
  cb=$(git -C "$d" branch --show-current); ch=$(git -C "$d" rev-parse HEAD); cs=$(git -C "$d" status --porcelain | wc -l)
  if [ "$cb" = "$br" ] && [ "$ch" = "$sha" ] && [ "$cs" -eq 0 ]; then res=PASS; else res=FAIL; ok=0; fi
  say "G2 $name: branch $cb, HEAD $ch, porcelain lines $cs -> $res"
done
say "G3 git diff $START --stat:"
git diff "$START" --stat | sed 's/^/    /'
files=$(git diff "$START" --name-only)
bad=$(printf '%s\n' "$files" | command grep -E '^(vignettes|quarto)/' || true)
say "G3 nothing under vignettes/ or quarto/: $([ -z "$bad" ] && echo PASS || { ok=0; echo "FAIL: $bad"; })"
allowed=$(cat <<'EOF'
dev/tasks/TASK_cran_prep_release_2026-11-12.md
R/run_simulation_analysis.R
man/run_simulation_analysis.Rd
tests/testthat/test-sim-effect-scales.R
NEWS.md
DESCRIPTION
cran-comments.md
EOF
)
extra=$(printf '%s\n' "$files" | command grep -v -x -F -f <(printf '%s\n' "$allowed") | command grep -v -E '^dev/reports/REPORT_cran_prep_release_2026-11-12' || true)
say "G3 every changed file is one the report names: $([ -z "$extra" ] && echo PASS || { ok=0; echo "FAIL: $extra"; })"
tags_new=$(git tag --contains "$START" 2>/dev/null)
merges=$(git rev-list --merges "$START"..HEAD | wc -l)
ocp=$(git rev-parse origin/cran-prep); om=$(git rev-parse origin/master)
say "G4 tags containing the start SHA: ${tags_new:-none}; merge commits since start: $merges; origin/cran-prep $ocp; origin/master $om"
if [ -n "$tags_new" ] || [ "$merges" -ne 0 ] || [ "$ocp" != "2beb79c45a50620932cd0f53f2c4e7e87df25579" ] || [ "$om" != "2beb79c45a50620932cd0f53f2c4e7e87df25579" ]; then ok=0; say "G4 FAIL"; else say "G4 PASS"; fi
say "FINAL GATES: $([ $ok -eq 1 ] && echo PASS || echo FAIL)"
