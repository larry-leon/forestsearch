#!/usr/bin/env bash
# comment_gate.sh <pre-ref> [repo-dir] -- the task's comment-only gate on every changed file under
# R/ and tests/ between <pre-ref> and the working tree (staged + unstaged). Empty output = PASS.
set -u
PRE=$1; D=${2:-$HOME/Documents/GitHub/forestsearch}; cd $D
fail=0
for f in $(git diff --name-only $PRE -- R tests); do
  out=$(git diff -U0 $PRE -- "$f" | grep -E '^[+-]' | grep -vE '^(\+\+\+|---)' | grep -vE "^[+-][[:space:]]*#")
  if [ -n "$out" ]; then echo "NON-COMMENT LINES in $f:"; echo "$out"; fail=1; else echo "comment-only: $f"; fi
done
[ $fail = 0 ] && echo "COMMENT GATE: PASS" || echo "COMMENT GATE: FAIL"
