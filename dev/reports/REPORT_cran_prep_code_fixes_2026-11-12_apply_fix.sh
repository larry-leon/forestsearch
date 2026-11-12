#!/usr/bin/env bash
# apply_fix.sh <FIX> <filter> -- main-session red/green on the repository tree.
# 1. apply only the patch's tests/ hunks, run the fix's test (RED, expect failures)
# 2. apply the rest of the patch, run the test again (GREEN, expect 0 failures)
# Output: /tmp/cran_prep_code_fixes_2026-11-12/main/<FIX>/{red,green}.txt
set -u
F=$1; FILT=$2
S=/tmp/cran_prep_code_fixes_2026-11-12; P=$S/patches/$F.patch; O=$S/main/$F
REPO=$HOME/Documents/GitHub/forestsearch
mkdir -p $O; cd $REPO
[ -z "$(git status --porcelain)" ] || { echo "tree not clean"; exit 2; }
if git apply --check $P; then A="git apply"; else echo "plain apply fails; using -3"; A="git apply -3"; $A --check $P || { echo "patch does not apply"; exit 3; }; fi
$A --include='tests/*' $P || exit 4
t0=$(date +%s.%N)
NOT_CRAN=true timeout 300 Rscript --vanilla -e "devtools::test(filter = '$FILT', reporter = 'summary', stop_on_failure = FALSE)" > $O/red.txt 2>&1
t1=$(date +%s.%N)
$A --exclude='tests/*' $P || exit 5
NOT_CRAN=true timeout 300 Rscript --vanilla -e "devtools::test(filter = '$FILT', reporter = 'summary', stop_on_failure = FALSE)" > $O/green.txt 2>&1
t2=$(date +%s.%N)
printf "red wall %.1f s; green wall %.1f s\n" $(echo "$t1 - $t0" | bc) $(echo "$t2 - $t1" | bc) | tee $O/wall.txt
echo "--- RED tail"; tail -15 $O/red.txt; echo "--- GREEN tail"; tail -8 $O/green.txt
git status --porcelain
