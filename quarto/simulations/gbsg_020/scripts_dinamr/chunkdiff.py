import io, re, sys, difflib
def chunks(p):
    s = io.open(p, encoding="utf-8").read().split("\n")
    out, cur, lab = {}, None, None
    for ln in s:
        m = re.match(r'^```\{r ([A-Za-z0-9_-]*)', ln)
        if m:
            lab = m.group(1) or "(unlabelled)"; cur = []; continue
        if cur is not None and ln.strip() == "```":
            out[lab] = cur; cur = None; continue
        if cur is not None: cur.append(ln)
    return out
A = chunks("summary_cert20.qmd"); B = chunks("summary_dinamr.qmd")
ret = [k for k in A if k in B]
for k in ret:
    a, b = A[k], B[k]
    d = [l for l in difflib.unified_diff(a, b, lineterm="", n=0) if l[:1] in "+-" and l[:3] not in ("---","+++")]
    print("="*78); print("CHUNK  %s   (cert20 %d lines -> dinamr %d lines, %d changed lines)" % (k, len(a), len(b), len(d)))
    for l in d: print("   " + l[:160])
