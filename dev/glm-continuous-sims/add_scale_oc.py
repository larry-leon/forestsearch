#!/usr/bin/env python3
"""
add_scale_oc.py -- store the computed scale and the OC summary in every payload.

WHAT IT ADDS, to each batch document:
  1. `dgm_scale <- fs_dgm_scale(dgm)` right after the scoring frame is built.
  2. `scale` and `oc` blocks in BOTH payload writes -- the batch branch and the
     combine branch.  Missing the combine branch would leave pooled bundles
     silently without them.
  3. A self-check printing the scale-predicted oracle SD against the measured
     one.  That check, had it existed, would have shown the 2.4% anchor drift
     in the payload in August rather than after a two-day investigation.

WHAT IT DOES NOT TOUCH: the estimation-table rendering.  Doc `batch_1_50` uses
a different, older layout (relative bias) that fs_mr_oc_summary() cannot
reproduce, so table rendering is a separate task.

SAFETY: dry run by default; every file gated on all four anchors matching byte
for byte exactly once; refuses to write while any eligible file is blocked; the
write is the last operation.

USAGE
    python3 add_scale_oc.py --root . --only <substring> [--apply]
    python3 add_scale_oc.py --root . [--apply]          # every eligible document

The five batch documents differ only in n and in run size, so one is revised
first as the prototype.  This script is idempotent and reports an already
patched file as `already-done`, so a later run without --only picks up whichever
documents remain.
"""
import argparse, os, re, subprocess, sys

MARKER = "fs_build_eval_frame"

A_SCALE = '  eval_df <- fs_build_eval_frame(dgm, outcome_type = "continuous")\n'
A_SCALE_ALT = 'eval_df <- fs_build_eval_frame(dgm, outcome_type = "continuous")\n'

R_SCALE = '''eval_df <- fs_build_eval_frame(dgm, outcome_type = "continuous")

# Computed sampling scale, replacing the measured anchor.  fs_dgm_scale()
# enumerates sigma^2 + V_g[mu_0] + C_g[mu_0, tau] + V_g[tau]/2 on df_super, so
# the scale is a constant of the fixture rather than a number read off one run.
# The anchor it replaces (13.6786 at n = 1000) implied a bracket of 16,119
# against a residual variance of 16,256 -- below the noise floor, which on the
# true region is impossible since the bracket there is sigma^2 + V_Q[mu_0].
dgm_scale <- fs_dgm_scale(dgm)
'''

A_BATCH_OPEN = '\nsaveRDS(list(results = results, truth = truth,\n'
A_BATCH_CLOSE = 'built_at = Sys.time())), out)\n'
R_BATCH_OPEN = '\n.payload <- list(results = results, truth = truth,\n'
R_BATCH_CLOSE = '''built_at = Sys.time()))

# Scale and OC summary travel WITH the payload.  Computed here, where the DGM
# and the results are both in hand, so a document's tables and its payload
# cannot disagree and the analytic documents read numbers instead of
# transcribing literals.
.payload$scale <- dgm_scale
.payload$oc    <- fs_mr_oc_summary(.payload)
saveRDS(.payload[c("results", "truth", "scale", "oc", "meta")], out)

# Self-check.  The oracle refits on the TRUE region, so its replicate spread
# estimates sigma[betahat_n(Q)] directly.  This is the comparison that was
# missing while the anchor sat 2.4% low across two designs.
if (sum(results$status %in% "DETECTED") > 1L) {
  .pred <- fs_scale_se(dgm_scale, n = n_sample, region = "Q", jensen = TRUE)
  .orow <- .payload$oc$estimation
  .orow <- .orow[.orow$block == "H" & .orow$estimator == "oracle", ]
  cat(sprintf(paste0("SCALE CHECK  predicted SD %.4f | measured SD %.4f (%+.2f%%)",
                     " | mean fitted SE %.4f (%+.2f%%)\\n"),
              .pred, .orow$sd_emp, 100 * (.orow$sd_emp / .pred - 1),
              .orow$se_hat, 100 * (.orow$se_hat / .pred - 1)))
}
'''

A_COMB_OPEN = '    saveRDS(list(results = results, truth = truth,\n'
A_COMB_CLOSE = '                             built_at = Sys.time())), cpath)\n'
R_COMB_OPEN = '    .cpayload <- list(results = results, truth = truth,\n'
R_COMB_CLOSE = '''                             built_at = Sys.time()))
    # Pooled bundles carry scale and oc too; omitting them here would make a
    # combined payload silently poorer than the batches it merges.
    .cpayload$scale <- dgm_scale
    .cpayload$oc    <- fs_mr_oc_summary(.cpayload)
    saveRDS(.cpayload[c("results", "truth", "scale", "oc", "meta")], cpath)
'''

EXCLUDE = re.compile(r"(^|/)\.claude/|(^|/)\.git/|_legacy(/|$)|archive|gbsg_redux|\.Rcheck(/|$)")
PAIRS = [(A_SCALE_ALT, R_SCALE), (A_BATCH_OPEN, R_BATCH_OPEN),
         (A_BATCH_CLOSE, R_BATCH_CLOSE), (A_COMB_OPEN, R_COMB_OPEN),
         (A_COMB_CLOSE, R_COMB_CLOSE)]

def tracked(root):
    try:
        return set(subprocess.run(["git","-C",root,"ls-files"],capture_output=True,
                                  text=True,check=True).stdout.splitlines())
    except Exception:
        return set()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="."); ap.add_argument("--apply", action="store_true")
    ap.add_argument("--only", default=None,
                    help="substring; restrict to matching paths. Omit for all.")
    a = ap.parse_args()
    root = os.path.abspath(a.root); tr = tracked(root)

    found=[]
    for dp,dn,fn in os.walk(root):
        dn[:] = [d for d in dn if d not in (".git","node_modules")]
        for f in fn:
            if not f.endswith((".qmd",".Rmd")): continue
            full=os.path.join(dp,f)
            try: t=open(full,encoding="utf-8",errors="replace").read()
            except OSError: continue
            if MARKER in t and "oracle_md_on" in t: found.append((os.path.relpath(full,root),full,t))

    rows=[];staged=[];blocked=[]
    for rel,full,t in sorted(found):
        if a.only and a.only not in rel:
            rows.append((rel,"skipped","not selected by --only")); continue
        if EXCLUDE.search(rel): rows.append((rel,"excluded","-")); continue
        if rel not in tr:       rows.append((rel,"untracked","-")); continue
        if "dgm_scale <- fs_dgm_scale(dgm)" in t:
            rows.append((rel,"eligible","already-done")); continue
        counts=[t.count(x) for x,_ in PAIRS]
        if counts==[1,1,1,1,1]:
            rows.append((rel,"eligible","READY")); staged.append((full,t))
        else:
            rows.append((rel,"eligible",f"ANCHOR-MISMATCH {counts}")); blocked.append(rel)

    w=max([len(r[0]) for r in rows],default=10)
    print(f"{'file':<{w}}  {'class':<10} status"); print("-"*(w+28))
    for rel,cls,st in rows: print(f"{rel:<{w}}  {cls:<10} {st}")
    print(f"\nREADY: {len(staged)}   blocked: {len(blocked)}")
    for b in blocked: print(f"  BLOCKED {b}")

    if not a.apply:
        print("\n(dry run -- nothing written; re-run with --apply)"); return 0
    if blocked:
        print("\n!! refusing to write while any eligible file is blocked."); return 1

    for full,t in staged:
        new=t
        for anc,rep in PAIRS:
            assert new.count(anc)==1, f"{anc[:40]} in {full}"
            new=new.replace(anc,rep)
        assert "dgm_scale <- fs_dgm_scale(dgm)" in new
        assert new.count("fs_mr_oc_summary(")==2, f"expected 2 oc calls: {full}"
        assert "saveRDS(.payload[" in new and "saveRDS(.cpayload[" in new
        open(full,"w",encoding="utf-8").write(new); print(f"written: {full}")
    print(f"\n{len(staged)} file(s) written.")
    return 0

if __name__=="__main__": sys.exit(main())
