# NOTE — `guohe_adaptive_r()`'s r̂ is substantially seed-determined on the t7 design (2026-09-09)

Raised by the Gate 1b pilot and confirmed by a seed-offset test; both recorded in
`quarto/GuoHe/REPORT_guohe_adaptive_gate1b_2026-09-09.md`. **Record only — no repair, and the
seed derivation in the committed driver is untouched.**

## The observation

- Gate 1b pilot, cell `t7_beta2_03`, 10 replicates, `r_grid = c(1/3, 1/12)`, `v = 5`, `B = 2000`,
  `orient = +1`: **r̂ = 1/12 on 5/5 even replicate indices and 0/5 odd** — perfect alternation.
- Under a fair coin the probability of a perfect parity split in either direction is
  **2 × (1/2)¹⁰ = 1/512 ≈ 0.002 (0.2%)**.
- Cell `t7_beta2_00` leaned the same way without being perfect: 1/12 on 4/5 even, 2/5 odd.
- Every selection was a near-tie. Objective differences (1/12 − 1/3) had median |diff| 2.18e-03
  and max 8.75e-03, against an objective spread across replicates of ≈ 0.28 — i.e. **the gaps that
  decide r̂ are under 2% of the objective's own spread.**

## The seed derivation

- `quarto/GuoHe/guohe_sec52_adaptive_run.R`: `GHA_SEED_OFFSET <- 600000L`, and per replicate
  `seed_ad = base + m + GHA_SEED_OFFSET`, with `base = 1000000L + as.integer(sum(utf8ToInt(id)) * 100003L)`.
- **Consecutive replicates therefore differ by exactly 1 in the adaptive seed.**
- The offset was chosen to be distinct from `gh52_one_rep()`'s `boot_seed = seed + 500000L` so the
  streams could not collide. That objective is met; the parity effect is a separate matter.

## The confirming test

Same 10 replicates of `t7_beta2_03`, everything identical except the offset `+600000L` →
`+700000L`. The committed driver was not modified; the test replicated its per-replicate logic in
a scratch script. Bundle: `quarto/GuoHe/guohe_adaptive_t7_beta2_03_grid2_seed700k.rds`.

| | r̂ = 1/12 on even m | r̂ = 1/12 on odd m |
|---|---|---|
| offset `+600000L` | **5/5** | **0/5** |
| offset `+700000L` | 4/5 | 4/5 |

- **The parity alignment moved with the offset** — perfect at 600k, absent at 700k.
- **r̂ changed on 5 of 10 replicates from the seed change alone**, with the data, candidate family
  and selection untouched (selection keys reproduced 10/10).
- Per-replicate cost essentially unchanged: 367.9 s vs 357.3 s.

## Mechanism

- `guohe_adaptive_r()` draws its v-fold assignment from the supplied seed, and draws
  **independently across r** — no common random numbers on the candidate grid
  (`quarto/GuoHe/guohe_reproduction_RUN.md:131-136` flags this as a known risk).
- With no CRN, the between-candidate comparison carries the full Monte Carlo noise of two
  independent sets of draws. On t7 that noise is comparable to the true objective gap between
  1/3 and 1/12, so the comparison is a near-tie in essentially every replicate.
- A near-tie is decided by whatever the fold assignment happens to be, and the fold assignment is
  a deterministic function of `base + m + offset`. Consecutive Mersenne-Twister seeds are a known
  source of correlated early draws, which is a credible route from "consecutive m" to
  "parity-aligned folds". **The route is plausible and consistent with the evidence; it was not
  isolated further.**

## What it implies for a full Adaptive column

- An Adaptive column at 2000 replicates on this design would report an r̂ distribution that is
  **substantially an artefact of the seed grid, not a measurement of adaptive selection.**
- Its coverage would be a blend of the fixed-r columns in a mixing ratio set by the seed
  derivation. Changing the offset would change the column without changing any data.
- The Gate 1b cost measurement stands independently: the adaptive path costs **×10.28** the
  fixed-r bound (356.5 vs 34.7 core-s per replicate). That premium buys a selection that, on this
  design and grid, is close to arbitrary.
- This is consistent with, and gives a concrete mechanism for, **Guo & He's own Table 6 caution**
  on the adaptive procedure.

## Scope of the claim

- Established on **t7 only**, two cells, 10 replicates each, with a **two-element** grid
  `c(1/3, 1/12)`. Nothing here is claimed about the published four-element grid, about Tables 3–6,
  or about designs where the objective separates the candidates.
- The finding is about **r̂'s stability**, not about the correctness of `guohe_algorithm3()` at any
  fixed r. The fixed-r bounds are unaffected.

## Disposition

- **No repair in this task.** `R/guohe_adaptive_r.R` is not modified; the seed derivation in
  `guohe_sec52_adaptive_run.R` is not modified.
- **No committed result is invalidated.** The Adaptive column was never produced, and the
  reproduction's Tables 3–6 Adaptive columns are out of scope here (see also
  `dev/notes/NOTE_adaptive_B_inert_2026-09-09.md`).
- Options, if it is ever pursued: common random numbers across the candidate grid; a seed stride
  larger than 1 across replicates; or reporting a fixed r and saying why. **Larry's call.**
