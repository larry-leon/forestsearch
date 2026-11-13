#!/usr/bin/env bash
# commit_fix.sh <FIX> <message-file> -- stage the patch's paths, check the staged diff equals the
# gated patch (ignoring 'index' and hunk-header lines), print the stat, commit (part 5).
set -eu
F=$1; MSG=$2; S=/tmp/cran_prep_displays_se_docs_2026-11-12; P=$S/patches/$F.patch
mkdir -p $S/main/$F
cd $HOME/Documents/GitHub/forestsearch
paths=$(grep '^diff --git' $P | sed 's#^diff --git a/\(.*\) b/.*#\1#')
git add -- $paths
if diff <(git diff --cached --binary | grep -v '^index \|^@@') <(grep -v '^index \|^@@' $P) > $S/main/$F/patch_equal.txt; then echo "staged diff == gated patch"; else echo "STAGED DIFF DIFFERS FROM PATCH"; head -20 $S/main/$F/patch_equal.txt; exit 9; fi
git diff --cached --stat | tee $S/main/$F/diffstat.txt
[ -z "$(git status --porcelain | grep -v '^[MA] ')" ] || { echo "unstaged/untracked leftovers:"; git status --porcelain; exit 8; }
git commit -q -F $MSG
git log --oneline -1 | tee $S/main/$F/commit.txt
