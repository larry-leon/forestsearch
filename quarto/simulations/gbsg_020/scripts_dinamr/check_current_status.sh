#!/bin/zsh
# Post-condition for the gbsg_020 closeout step:
#   the pin stated in current_status.md == HEAD at the time it was committed.
#
# current_status.md is committed ALONE, as the child of the commit it describes,
# so "HEAD at commit time" is that file's parent commit.  Two modes:
#
#   check_current_status.sh          before committing: stated pin == HEAD now
#   check_current_status.sh --commit after committing : stated pin == HEAD~1,
#                                    and the closeout commit touched only
#                                    current_status.md
#
# Exit 0 = PASS, 1 = FAIL.  Run from anywhere inside the repo.
set -e
cd "$(git rev-parse --show-toplevel)"
F=quarto/simulations/gbsg_020/current_status.md
[[ -f $F ]] || { echo "FAIL: $F does not exist"; exit 1; }
STATED=$(grep -m1 -o '`[0-9a-f]\{7,40\}`' $F | tr -d '`')
[[ -n $STATED ]] || { echo "FAIL: no pin found in $F"; exit 1; }
if [[ "$1" == "--commit" ]]; then
  WANT=$(git rev-parse --short HEAD~1)
  FILES=$(git show --pretty=format: --name-only HEAD | grep -v '^$' | sort -u)
  if [[ "$FILES" != "$F" ]]; then
    echo "FAIL: the closeout commit touched more than $F:"; echo "$FILES"; exit 1
  fi
else
  WANT=$(git rev-parse --short HEAD)
fi
if [[ "$STATED" == "$WANT" ]]; then
  echo "PASS: stated pin $STATED == $([[ $1 == --commit ]] && echo 'HEAD~1' || echo 'HEAD') $WANT"
else
  echo "FAIL: stated pin $STATED != $WANT"; exit 1
fi
