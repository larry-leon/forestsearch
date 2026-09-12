#!/usr/bin/env python3
"""Apply the exact-string edit spec to quarto/GuoHe/guohe_supp_section.qmd and
check the post-conditions. Exact matching only: every 'replace' edit's old text
must occur exactly once in the current text; every 'span' edit's anchors must
each occur exactly once and the text between them must match its SHA-256.
Any failure -> non-zero exit before anything is written.

Usage (from the repo root):
  python3 dev/tasks/claude_cc_task_guohe_record_payload_2026-09-11_apply.py \
      dev/tasks/claude_cc_task_guohe_record_payload_2026-09-11_edits.json [--check-only]
"""
import hashlib, json, re, sys

def sha(s): return hashlib.sha256(s.encode("utf-8")).hexdigest()

FIXED_1 = ("Table 8 of Guo and He (2021) reports an upper confidence bound on the "
           "log-hazard-ratio scale; this is the identical one-sided lower-bound "
           "construction under their Section 4 convention βᵢ = −log HR, applied to "
           "the argmax-selected subgroup, and does not constitute an upper-bound "
           "capability for a second, complementary subgroup.")
FIXED_2 = ("The Guo–He correction is defined through the maximum functional over the "
           "supplied candidate family. The complement Ĥᶜ is not a member of that "
           "family and is not the argmax of any functional, so no analogous "
           "correction exists within their framework; the absence of a complement "
           "bound, and hence of a joint two-subgroup claim, is structural rather "
           "than an implementation gap.")
FORBIDDEN = ["stable-pick", "data-built", "CERT_", "cert_upper", "cert_joint",
             "The wider figures", '"dev", "notes"', "Classification is not the right lens"]
FORBIDDEN_CI = ["certif"]
LABELS = ["setup", "b2-build", "b2-table", "b2-adaptive-note", "b2-consolidated",
          "b3-table", "b5", "b5-chat-plot", "b6", "export-payload", "provenance",
          "provenance-platform"]

def blockquote_text(txt):
    """Join each run of '> ' lines into one whitespace-normalised string."""
    out, cur = [], []
    for line in txt.split("\n"):
        if line.startswith("> "):
            cur.append(line[2:].strip())
        elif cur:
            out.append(" ".join(cur)); cur = []
    if cur: out.append(" ".join(cur))
    return [re.sub(r"\s+", " ", s).strip() for s in out]

UNTOUCHED = ["b2-build", "b2-table", "b2-adaptive-note", "b2-consolidated",
             "b3-table", "b5-chat-plot", "provenance-platform"]

def chunks(txt):
    """label -> full chunk text (header line through closing fence)."""
    out = {}
    for m in re.finditer(r"^```\{r ([A-Za-z0-9_-]+)\}\s*\n.*?^```\s*$", txt, flags=re.M | re.S):
        out[m.group(1)] = m.group(0)
    return out

def postconditions(txt, before=None):
    res = {}
    if before is not None:
        cb, ca = chunks(before), chunks(txt)
        for lab in UNTOUCHED:
            res[f"untouched_identical[{lab}]"] = int(lab in cb and cb.get(lab) == ca.get(lab))
    bq = blockquote_text(txt)
    res["fixed_sentence_1_count"] = sum(s == FIXED_1 for s in bq)
    res["fixed_sentence_2_count"] = sum(s == FIXED_2 for s in bq)
    res["sec_b3_label_count"] = txt.count("{#sec-b3}")
    for f in FORBIDDEN: res[f"forbidden[{f}]"] = txt.count(f)
    for f in FORBIDDEN_CI: res[f"forbidden_ci[{f}]"] = txt.lower().count(f)
    for lab in LABELS:
        res[f"chunk[{lab}]"] = len(re.findall(r"^```\{r " + re.escape(lab) + r"\}\s*$", txt, flags=re.M))
    ok = (res["fixed_sentence_1_count"] == 1 and res["fixed_sentence_2_count"] == 1
          and res["sec_b3_label_count"] == 1
          and all(v == 0 for k, v in res.items() if k.startswith("forbidden"))
          and all(v == 1 for k, v in res.items() if k.startswith("chunk["))
          and all(v == 1 for k, v in res.items() if k.startswith("untouched_identical[")))
    return ok, res

def main():
    spec = json.load(open(sys.argv[1], encoding="utf-8"))
    check_only = "--check-only" in sys.argv
    path = spec["target"]
    txt = open(path, encoding="utf-8").read()
    original = txt
    if sha(txt) != spec["target_sha256"]:
        sys.exit(f"STOP: {path} sha256 {sha(txt)} != spec {spec['target_sha256']} "
                 "(the edits were written against a different text)")
    for e in spec["edits"]:
        if e["kind"] == "replace":
            n = txt.count(e["old"])
            if n != 1: sys.exit(f"STOP: {e['id']} old text occurs {n} times (need exactly 1)")
            txt = txt.replace(e["old"], e["new"], 1)
        elif e["kind"] == "span":
            ns, ne = txt.count(e["start"]), txt.count(e["end"])
            if ns != 1 or ne != 1: sys.exit(f"STOP: {e['id']} anchors occur {ns}/{ne} times (need 1/1)")
            i = txt.index(e["start"]); j = txt.index(e["end"], i)
            if sha(txt[i:j]) != e["span_sha256"]:
                sys.exit(f"STOP: {e['id']} span sha256 mismatch (span changed since the spec was written)")
            txt = txt[:i] + e["new"] + txt[j:]
        else:
            sys.exit(f"STOP: unknown edit kind {e['kind']}")
        print(f"applied {e['id']:>4}  {e['purpose']}")
    ok, res = postconditions(txt, before=original)
    for k, v in res.items(): print(f"  {k:<48} {v}")
    if not ok: sys.exit("STOP: post-conditions failed; file NOT written")
    if sha(txt) != spec["result_sha256"]:
        sys.exit(f"STOP: edited text sha256 {sha(txt)} != expected {spec['result_sha256']}")
    print(f"  result sha256 matches the spec: {spec['result_sha256']}")
    if check_only:
        print("check-only: all edits apply and all post-conditions hold; file not written"); return
    open(path, "w", encoding="utf-8").write(txt)
    print(f"WROTE {path}  sha256 {sha(txt)}")

if __name__ == "__main__":
    main()
