#!/usr/bin/env python3
"""Apply the article-alignment edit set to vignettes/articles/ (2026-10-30).

    python3 -I apply_alignment_edits.py alignment_edits_2026-10-30.json vignettes/articles [--check]

For every file in the edit set:
  * if its md5 already equals `md5_after`  -> "already applied", nothing written (idempotent);
  * if its md5 is not `md5_before`         -> FAIL (the file is not the master c2967362 text);
  * otherwise every anchor (`old`) must occur exactly once and every replacement
    (`new`) zero times before the edit is applied; after all edits the md5 must
    equal `md5_after`, or nothing is written.
With --check, nothing is written in any case.  Exit status 0 only when every file
is applied or already applied; 1 otherwise.  Writes are UTF-8, newline-preserving.
"""
import hashlib, json, os, sys

def md5(s: str) -> str:
    return hashlib.md5(s.encode("utf-8")).hexdigest()

def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__); return 2
    edits_path, art_dir = sys.argv[1], sys.argv[2]
    check_only = "--check" in sys.argv[3:]
    with open(edits_path, encoding="utf-8") as fh:
        spec = json.load(fh)
    status_ok = True
    for fname, meta in spec["files"].items():
        path = os.path.join(art_dir, fname)
        with open(path, encoding="utf-8", newline="") as fh:
            text = fh.read()
        h = md5(text)
        edits = [e for e in spec["edits"] if e["file"] == fname]
        if h == meta["md5_after"]:
            print(f"{fname}: already applied ({len(edits)} edits; md5 {h[:8]})")
            continue
        if h != meta["md5_before"]:
            print(f"{fname}: FAIL md5 {h[:8]} is neither master ({meta['md5_before'][:8]}) "
                  f"nor applied ({meta['md5_after'][:8]})")
            status_ok = False
            continue
        bad = []
        for e in edits:
            n_old, n_new = text.count(e["old"]), text.count(e["new"])
            if n_old != 1 or n_new != 0:
                bad.append((e["id"], n_old, n_new))
        if bad:
            for eid, n_old, n_new in bad:
                print(f"{fname}: FAIL anchor {eid}: old occurs {n_old} (want 1), new occurs {n_new} (want 0)")
            status_ok = False
            continue
        out = text
        for e in edits:
            out = out.replace(e["old"], e["new"], 1)
        for e in edits:
            if out.count(e["new"]) != 1:
                print(f"{fname}: FAIL post-check {e['id']}: new occurs {out.count(e['new'])} after apply")
                status_ok = False
        h2 = md5(out)
        if h2 != meta["md5_after"]:
            print(f"{fname}: FAIL result md5 {h2[:8]} != expected {meta['md5_after'][:8]}")
            status_ok = False
            continue
        if check_only:
            print(f"{fname}: OK (check only) {len(edits)} edits; {h[:8]} -> {h2[:8]}")
            continue
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write(out)
        print(f"{fname}: applied {len(edits)} edits; md5 {h[:8]} -> {h2[:8]}")
    print("RESULT:", "OK" if status_ok else "FAIL")
    return 0 if status_ok else 1

if __name__ == "__main__":
    sys.exit(main())
