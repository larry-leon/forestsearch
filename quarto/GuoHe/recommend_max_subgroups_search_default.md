# Recommendation — `max_subgroups_search`: default value and truncation visibility

**Context.** On GBSG, `sg_focus = "maxSG"` at `max_subgroups_search = 50`
returned **no subgroup**. At 2000 it returned `{er <= 0}` (n = 82, naive HR
1.9514) after evaluating only 90 candidates. Nothing else changed. The `NULL` was
candidate-pool truncation, not absence of a harmful subgroup — and it was
**silent**.

---

## 1. The mechanism, and why it is easy to miss

`subgroup_consistency_main.R:565-567`:

```r
maxsgs    <- min(nrow(found.hrs), stop_Kgroups)   # stop_Kgroups = max_subgroups_search
found.hrs <- found.hrs[seq_len(maxsgs), ]
```

Unconditional. It is **not** tied to `stop_threshold`, so disabling early
stopping does not disable it.

The truncation is applied to the `sort_subgroups_preview()` ordering, which is
`sg_focus`-dependent:

| `sg_focus` | preview order | consequence |
|---|---|---|
| `hr` | `-HR` | most harmful first — the pool's head *is* what a harm search wants |
| `maxSG` | `-n` | largest first — harmful subgroups, which tend to be **small**, sit in the tail |
| `hrMaxSG`/`hrMinSG` | effect band | early stopping separately reset to `NULL`, full pool evaluated |

So the exposure is **inversely related to the focus's alignment with the
search's objective**. `maxSG` needs the deepest pool precisely because its own
ordering reaches the relevant candidates last. On GBSG the harm subgroup sits at
n = 61, one above the `n.min = 60` floor — about as deep in a `-n` ordering as it
is possible to be.

---

## 2. Recommendation 1 (higher priority): make truncation announce itself

Raising the default reduces the chance of this failure. It does not remove it,
and it does not make it visible when it does occur. **A silent truncation is the
actual defect.**

The existing guardrail (`forestsearch_main.R:1610-1620`) does not cover this:

```r
if (details &&
    .normalize_sg_focus(sg_focus) %in% c("hrMaxSG", "hrMinSG") &&
    max_subgroups_search < 30) { ... }
```

Three reasons it stayed silent here: it fires only for `hrMaxSG`/`hrMinSG` (not
`maxSG`); only below a **fixed** threshold of 30 (the run used 50); and only when
`details = TRUE`. Meanwhile truncation removed the answer.

**Proposed:** warn whenever truncation actually binds, regardless of value,
focus, or `details`:

```r
  if (nrow(found.hrs) > stop_Kgroups) {
    warning(sprintf(
      paste0("max_subgroups_search = %s truncated the candidate pool from %d to %d. ",
             "Candidates ranked below %d under the sg_focus = '%s' preview ordering ",
             "were not evaluated. Size-focused foci ('maxSG'/'minSG') order by ",
             "subgroup size, so small high-effect subgroups rank last and are the ",
             "first to be truncated."),
      format(stop_Kgroups), nrow(found.hrs), stop_Kgroups, stop_Kgroups, sg_focus),
      call. = FALSE)
  }
```

This is condition-based rather than threshold-based, so it cannot be outgrown by
a future default. Had it existed, the `maxSG` `NULL` would have been
self-explaining.

Keep the existing advisory note as well — it carries focus-specific guidance the
generic message does not.

---

## 3. Recommendation 2: raise the defaults

Current: `forestsearch(max_subgroups_search = 30)`,
`subgroup.consistency(stop_Kgroups = 10)`.

Both are small relative to realistic pools. With `maxk = 2` and $L$ cut columns
the enumerated space is $L + \binom{L}{2}$; GBSG's 66 dummy-expanded cuts give
2211, of which 1744 clear `n.min`. Continuous covariates cut at several quantiles
each make $L$ in the dozens routine, so pools in the hundreds-to-thousands are
the norm, not the exception. A cap of 30 evaluates under 2% of that.

**Recommended: `max_subgroups_search = 200`, `stop_Kgroups = 200`.**

Your suggestion of 100-200 is right; 200 is the safer end, for two reasons.

**The cost is usually bounded by early stopping, not by the cap.** With
`stop_threshold` active the search halts at the first passer — the depth-2000 run
evaluated **90** candidates, not 2000. So a larger cap mostly buys *availability*,
not work. Cost scales with the cap only when no candidate passes, or when early
stopping is off (`hrMaxSG`/`hrMinSG`, or `stop_threshold = NULL`).

**Where it does scale, `consistency_method` dominates.** Under `"resample"` each
candidate is one fit plus the closed form $2\Phi((\hat\beta-c)/\sigma_D)-1$ —
200 candidates is seconds. Under `"split"` it is `fs.splits` (default 1000) × 2
refits per candidate, so 200 candidates is ~400k fits. **The default should be
raised in tandem with a note that deep pools assume `"resample"`.**

Going beyond ~200 by default is not obviously worthwhile: the marginal
protection falls off once the cap exceeds the depth at which passers occur (90
here), while worst-case cost keeps rising.

---

## 4. Compatibility

Changing a default **changes results** for anyone relying on it — subgroups may
differ where truncation was previously binding. This is a behaviour change, not
a bug fix, so:

- **NEWS entry**, stating explicitly that fits using the old default may select
  different subgroups, with the GBSG `maxSG` case as the worked example;
- **version bump** at least at the patch level;
- consider a **transitional release** where the truncation warning ships first
  and the default changes one release later — so users see *why* their results
  move before they move.

---

## 5. Interaction with the `maxeff` patch

None. `sg_focus = "maxeff"` overrides `max_subgroups_search` to `Inf`
unconditionally, so it is unaffected by whatever default is chosen. The two
changes are independent and can land in either order.

---

## 6. Summary

| | change | priority |
|---|---|---|
| 1 | warn whenever `nrow(found.hrs) > stop_Kgroups` | **high** — converts a silent failure into a visible one |
| 2 | `max_subgroups_search`, `stop_Kgroups` -> 200 | medium — reduces frequency; requires NEWS + version bump |
| 3 | document that deep pools assume `consistency_method = "resample"` | medium |
| 4 | keep the existing `hrMaxSG`/`hrMinSG` advisory | low |

Recommendation 1 is worth doing even if the default never changes.
