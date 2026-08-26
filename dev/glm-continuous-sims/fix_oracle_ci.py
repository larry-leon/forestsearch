#!/usr/bin/env python3
"""
fix_oracle_ci_v2.py -- propagate the oracle confidence-interval fix by
TRANSPLANT from a reference document that already carries it.

WHY v2
    v1 carried hand-authored replacement text.  The dry run then found that
    quarto/simulations/actg175/continuous/sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_50.qmd
    ALREADY contains this fix, in its own textual form, with a comment recording
    the same measurement independently ("mean CI length -3.92 x SE ... coverage
    identically 0").  Applying v1 would have left five documents carrying two
    different implementations of one fix.  v2 therefore reads the corrected
    blocks out of that file at run time and transplants them verbatim, so every
    document ends byte-identical in these two blocks.

THE DEFECT (unchanged from v1)
    .oracle_md_on() built c(est, lo, hi, se) with orientation already applied to
    est, so lo < hi held by construction.  It then reordered positionally with
    [c(1,3,2,4)] -> (est, hi, lo, se) and repaired the values BY NAME via sort().
    The consumer assigned POSITIONALLY, so or_*_lo received hi and or_*_hi
    received lo.  Measured on the committed n=500 and n=700 payloads:
    or_H_lo > or_H_hi in 1000/1000 rows; mean stored width -77.67, |width| =
    2 * 1.96 * mean(se) exactly; oracle coverage 0.000 as stored, 0.942 swapped
    back, against a normal-theory check of 0.942.

SAFETY
    - The reference's two blocks are extracted at run time and PRINTED, so the
      transplanted text is visible in the log rather than trusted.
    - A target is edited only if its .oracle_md_on() agrees with the reference's
      byte for byte UP TO the divergence point.  That gate means the transplant
      cannot mask an unrelated difference in the function body.
    - Dry run by default; refuses to write while any eligible file is blocked;
      the write is the last operation.

USAGE
    python3 fix_oracle_ci_v2.py --root . --reference <path/to/fixed.qmd>
    python3 fix_oracle_ci_v2.py --root . --reference <path/to/fixed.qmd> --apply
"""

import argparse
import os
import re
import subprocess
import sys

MARKER = ".oracle_md_on"

# The defective code, byte for byte.  Prose that merely quotes "[c(1,3,2,4)]"
# in a comment does not match, because of the "-> v".
ANCHOR_FN = (
    '  c(est = s*est, lo = s*est - 1.96*se, hi = s*est + 1.96*se, se = se)[c(1,3,2,4)] -> v\n'
    '  if (s < 0) v[c("lo","hi")] <- sort(unname(v[c("lo","hi")]))\n'
    '  v\n'
)
ANCHOR_USE = (
    '  rec[c("or_H_est","or_H_lo","or_H_hi","or_H_se")] <- as.list(.oracle_md(df))\n'
    '  rec[c("or_Hc_est","or_Hc_lo","or_Hc_hi","or_Hc_se")] <- as.list(.oracle_md_c(df))\n'
)
DEFECT_CODE = "[c(1,3,2,4)] -> v"

FN_START = ".oracle_md_on <- function"
OR_COLS = ["or_H_est", "or_H_lo", "or_H_hi", "or_H_se",
           "or_Hc_est", "or_Hc_lo", "or_Hc_hi", "or_Hc_se"]

EXCLUDE = re.compile(
    r"(^|/)\.claude/|(^|/)\.git/|_legacy(/|$)|archive|gbsg_redux|\.Rcheck(/|$)")


def tracked(root):
    try:
        out = subprocess.run(["git", "-C", root, "ls-files"],
                             capture_output=True, text=True, check=True).stdout
        return set(out.splitlines())
    except Exception as exc:                                   # noqa: BLE001
        print(f"!! could not list tracked files: {exc}", file=sys.stderr)
        return set()


def extract_fn_block(text):
    """The whole `.oracle_md_on <- function ... \\n}\\n` block."""
    i = text.find(FN_START)
    if i < 0:
        return None
    j = text.find("\n}\n", i)
    if j < 0:
        return None
    return text[i:j + 3]


def extract_use_block(text):
    """The recorder block: the contiguous run that writes the or_* columns.

    Seeded on an actual `.oracle_md(df)` CALL (not the one-line function
    definitions, which read `.oracle_md    <- function(df)`), then expanded
    outward only across adjacent recorder lines.  Two earlier attempts failed
    here and both are worth naming: taking the whole surrounding non-blank run
    reached back through simulate_from_glm_dgm() and would have transplanted
    unrelated code; taking min..max over every line mentioning `.oracle_md` or
    an or_* column swept in the `rec` initialiser and the function definitions.
    """
    lines = text.split("\n")

    def is_recorder(l):
        return (".oracle_md" in l or ".oH" in l or ".oHc" in l
                or "or_H_" in l or "or_Hc_" in l)

    seeds = [k for k, l in enumerate(lines)
             if ".oracle_md(df)" in l and "<- function" not in l]
    for h in seeds:
        a = h
        while a > 0 and is_recorder(lines[a - 1]) and lines[a - 1].strip():
            a -= 1
        b = h
        while b < len(lines) - 1 and is_recorder(lines[b + 1]) and lines[b + 1].strip():
            b += 1
        block = "\n".join(lines[a:b + 1]) + "\n"
        if all(c in block for c in OR_COLS) and "=" not in block.split("\n")[0].split("<-")[0]:
            return block
        if all(c in block for c in OR_COLS):
            return block
    return None


def common_prefix(a, b):
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return a[:i]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    ref_text = open(args.reference, encoding="utf-8").read()

    ref_fn = extract_fn_block(ref_text)
    ref_use = extract_use_block(ref_text)
    if ref_fn is None or ref_use is None:
        print("!! could not extract both blocks from the reference.")
        return 1
    if DEFECT_CODE in ref_fn or ANCHOR_USE in ref_text:
        print("!! the reference still carries the defect; wrong file.")
        return 1

    print("=" * 72)
    print(f"REFERENCE: {args.reference}")
    print("=" * 72)
    print("--- .oracle_md_on() to be transplanted ---")
    print(ref_fn)
    print("--- recorder block to be transplanted ---")
    print(ref_use)
    print("=" * 72 + "\n")

    tracked_set = tracked(root)
    found = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in (".git", "node_modules")]
        for fn in filenames:
            if not fn.endswith((".qmd", ".Rmd", ".R")):
                continue
            full = os.path.join(dirpath, fn)
            try:
                text = open(full, encoding="utf-8", errors="replace").read()
            except OSError:
                continue
            if MARKER in text:
                found.append((os.path.relpath(full, root), full, text))

    rows, staged, blocked = [], [], []
    for rel, full, text in sorted(found):
        if os.path.abspath(full) == os.path.abspath(args.reference):
            rows.append((rel, "reference", "-", "reference (unchanged)"))
            continue
        if EXCLUDE.search(rel):
            rows.append((rel, "excluded", "-", "excluded"))
            continue
        if rel not in tracked_set:
            rows.append((rel, "untracked", "-", "untracked"))
            continue
        if DEFECT_CODE not in text and ANCHOR_USE not in text:
            rows.append((rel, "eligible", "-", "already-fixed"))
            continue
        if text.count(ANCHOR_FN) != 1 or text.count(ANCHOR_USE) != 1:
            rows.append((rel, "eligible", "-",
                         f"ANCHOR-MISMATCH fn={text.count(ANCHOR_FN)} "
                         f"use={text.count(ANCHOR_USE)}"))
            blocked.append(rel)
            continue

        tgt_fn = extract_fn_block(text)
        if tgt_fn is None:
            rows.append((rel, "eligible", "-", "NO-FN-BLOCK"))
            blocked.append(rel)
            continue

        pre = common_prefix(ref_fn, tgt_fn)
        # The divergence must begin at the defective tail, not earlier.
        diverges_at_defect = tgt_fn[len(pre):].lstrip().startswith("c(est = s*est")
        if not diverges_at_defect:
            rows.append((rel, "eligible", f"{len(pre)}",
                         "PREFIX-DIVERGES-EARLY"))
            blocked.append(rel)
            continue

        rows.append((rel, "eligible", f"{len(pre)}", "READY"))
        staged.append((full, text, tgt_fn))

    w = max([len(r[0]) for r in rows], default=10)
    print(f"{'file':<{w}}  {'class':<10} {'pfx':>6}  status")
    print("-" * (w + 30))
    for rel, cls, pfx, status in rows:
        print(f"{rel:<{w}}  {cls:<10} {pfx:>6}  {status}")
    print(f"\nREADY: {len(staged)}   blocked: {len(blocked)}")
    if blocked:
        print("\nBLOCKED -- not edited:")
        for rel in blocked:
            print(f"  {rel}")

    if not args.apply:
        print("\n(dry run -- nothing written; re-run with --apply)")
        return 0
    if blocked:
        print("\n!! refusing to write while any eligible file is blocked.")
        return 1

    for full, text, tgt_fn in staged:
        new = text.replace(tgt_fn, ref_fn).replace(ANCHOR_USE, ref_use)
        assert DEFECT_CODE not in new, f"defect survived: {full}"
        assert ANCHOR_USE not in new, f"recorder anchor survived: {full}"
        assert extract_fn_block(new) == ref_fn, f"fn block not identical: {full}"
        assert extract_use_block(new) == ref_use, f"use block not identical: {full}"
        with open(full, "w", encoding="utf-8") as fh:
            fh.write(new)
        print(f"written: {full}")

    print(f"\n{len(staged)} file(s) written, all byte-identical to the reference "
          f"in both blocks.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
