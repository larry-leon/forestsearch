# REPORT — Guo & He record: factual corrections and the manuscript payload

Task: `dev/tasks/claude_cc_task_guohe_record_payload_2026-09-11.md` (issued 2026-09-11).
Executed 2026-09-11, unattended, by Claude Code. No `R/` change; no compute.

## Provenance

| item | value |
|---|---|
| machine | `pop-os` (Linux 7.0.11-76070011-generic) |
| repo / branch | `larry-leon/forestsearch` / `feature/glm-extension` |
| BASE (HEAD before the §1 commit) | `1d9401cb7765a6fdcf8239382a78598ec0c0e1db` |
| final HEAD (record + payload commit) | `bd1ce4bf` (this report is committed after it) |
| HEAD contains `1d9401cb` | yes (BASE *is* `1d9401cb`) |
| R | R version 4.6.1 (2026-06-24) |
| forestsearch | 0.3.5 |
| Quarto binary | `/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto` |
| Quarto version | 1.9.38 |
| locale | `LANG=en_US.UTF-8`, `LC_ALL=en_US.UTF-8` |

Companion SHA-256s, verified on copy from `~/Downloads` into `dev/tasks/` (§1):

| file | SHA-256 | matches table |
|---|---|---|
| `..._edits.json` | `8c0c09b8f15c7b2f1a137d22f9374893eafb53598e0c09c64239f89c1c850b82` | yes |
| `..._apply.py` | `3913c2017da21ab9e405d1a22e8d9b08d502c3832a7f4b58689fd2c2ad05edf2` | yes |
| `..._check.R` | `0789fccb2046d2280235287f5aa6dc372b92728a037d8f0c83a21cf8716957a2` | yes |

Renders (both with the binary above, UTF-8 locale):

| render | exit | wall | `WARNING` lines |
|---|---|---|---|
| §3 baseline (unedited record) | 0 | 9 s | 0 |
| §5.1 edited record | 0 | 11 s | 0 |

Working-tree dirt in this task's paths at kickoff (`quarto/GuoHe/`, `dev/tasks/`): none.
Six pre-existing untracked files elsewhere (`quarto/extreme_subgroups/`, `quarto/simulations/gbsg_020/`) left alone.

## Stage 0 — preconditions and source quotes

**2.1 Input hash.** `quarto/GuoHe/guohe_supp_section.qmd`
= `0058800765baa846721e562b59756be2cf875fd9b2e23d319cd65e0b00f12a4a`, git blob `4afa7e42` — matches the spec. PASS.

**2.2 Inputs present.** All 30 required bundles present under `quarto/GuoHe/`: the 16
`mr_field_vs_guohe_*` (t35_beta2_00..05, t6_k02/k06/k10/k12, t7_beta2_00..05), the 6
`guohe_repro_t7_beta2_00..05`, the 6 `mr_field_complement_vs_guohe_t7_beta2_00..05`, and the 2
`guohe_adaptive_t7_beta2_0{0,3}_fixedr00833`. None absent. PASS.

**2.3 Quotes by content** (current tree, with line numbers):

`R/fs_mr_inference.R`:

| line | text |
|---|---|
| 695 | `  se_wald    <- sdv[sel]` |
| 1019 | `                    lower = to_eff(beta_naive - z975 * se_wald),` |
| 810 | `                        lower = to_eff(bnc - z975 * sec),` |
| 816 | `                        se = sec_used, se_ij = se_ijc_rep$se, se_wald = sec,` |

`quarto/GuoHe/mr_field_complement_vs_guohe_run.R`:

| line | text |
|---|---|
| 228 | `    cgate_naive_est = .lg(gc$naive$est),` |
| 229 | `    cgate_est = .lg(gc$debiased$est),` |
| 230 | `    cgate_se_ij = .nz(gc$debiased$se_ij),` |
| 231 | `    cgate_se_wald = .nz(gc$debiased$se_wald),` |

All eight quotes present. The `comp_t7` rows `naive` and `ij` are therefore **confirmed**, not unconfirmed.

**2.4 Quarto.** The record's setup comment (line 40) says: "RENDERING: render with RStudio's
bundled Quarto binary on the executing machine and record the binary path used in the render
record." Path and version recorded in Provenance above.

## Stage 1 — baseline render (GATE 1)

Exit 0, wall 9 s, zero `WARNING` lines, `Output created: guohe_supp_section.html`.
The HTML is untracked (`quarto/.gitignore:62` → `GuoHe/*.html`). Nothing committed in this stage. PASS.

## Stage 2 — apply the edits (GATE 2)

### `--check-only`, verbatim (exit 0)

```
applied  E01  YAML title/subtitle: analysis record; no 'data-built' for the t7 family; Adaptive not a reported column
applied  E02  setup comment: this is the analysis record; the manuscript section is drafted from its payload
applied  E03  setup: remove the E1 range-citation block and the dev/notes read (span verified by SHA-256)
applied  E04  opening: 'nested, data-built' -> 'nested'
applied  E05  opening: say what the order-statistic family is
applied  E06  opening: D1 rationale without 'data-built'
applied  E07  B2 callout: replace the unsupported adjacency mechanism with Guo & He's own Section 2.5 reading
applied  E08  B3 fixed sentence 1, D1: 'hazard-ratio scale' -> 'log-hazard-ratio scale' (Guo & He Table 8 caption)
applied  E09  B3: render the complement-truth explanation (was only in the never-taken branch) and scope the structural argument
applied  E10  B5 table: specificity/PPV shown as '—' at beta2 = 0 (no harm region); add mean gamma_c-hat
applied  E11  B5 caption: say why beta2 = 0 shows '—' and what mean gamma_c-hat is
applied  E12  B5 callout: identification has one dimension here; the null cell as reference; overshoot reaches the estimand
applied  E13  B6: limit 1 re-grounded on this design's records; E1 ranges and dev/notes read removed; limit 3 acknowledges G&H Section 6; limit 5 added; cert-note chunk removed
applied  E14  insert the 'Manuscript payload' section (export-payload chunk) before Provenance
applied  E15  provenance: item label
applied  E16  provenance: item value
  untouched_identical[b2-build]                    1
  untouched_identical[b2-table]                    1
  untouched_identical[b2-adaptive-note]            1
  untouched_identical[b2-consolidated]             1
  untouched_identical[b3-table]                    1
  untouched_identical[b5-chat-plot]                1
  untouched_identical[provenance-platform]         1
  fixed_sentence_1_count                           1
  fixed_sentence_2_count                           1
  sec_b3_label_count                               1
  forbidden[stable-pick]                           0
  forbidden[data-built]                            0
  forbidden[CERT_]                                 0
  forbidden[cert_upper]                            0
  forbidden[cert_joint]                            0
  forbidden[The wider figures]                     0
  forbidden["dev", "notes"]                        0
  forbidden[Classification is not the right lens]  0
  forbidden_ci[certif]                             0
  chunk[setup]                                     1
  chunk[b2-build]                                  1
  chunk[b2-table]                                  1
  chunk[b2-adaptive-note]                          1
  chunk[b2-consolidated]                           1
  chunk[b3-table]                                  1
  chunk[b5]                                        1
  chunk[b5-chat-plot]                              1
  chunk[b6]                                        1
  chunk[export-payload]                            1
  chunk[provenance]                                1
  chunk[provenance-platform]                       1
  result sha256 matches the spec: 15abd1209a4ea42432389a6e8d8c7ece1972c4a39f7f4632ac0d65a7b33c6866
check-only: all edits apply and all post-conditions hold; file not written
```

### write run, verbatim (exit 0)

```
applied  E01  YAML title/subtitle: analysis record; no 'data-built' for the t7 family; Adaptive not a reported column
applied  E02  setup comment: this is the analysis record; the manuscript section is drafted from its payload
applied  E03  setup: remove the E1 range-citation block and the dev/notes read (span verified by SHA-256)
applied  E04  opening: 'nested, data-built' -> 'nested'
applied  E05  opening: say what the order-statistic family is
applied  E06  opening: D1 rationale without 'data-built'
applied  E07  B2 callout: replace the unsupported adjacency mechanism with Guo & He's own Section 2.5 reading
applied  E08  B3 fixed sentence 1, D1: 'hazard-ratio scale' -> 'log-hazard-ratio scale' (Guo & He Table 8 caption)
applied  E09  B3: render the complement-truth explanation (was only in the never-taken branch) and scope the structural argument
applied  E10  B5 table: specificity/PPV shown as '—' at beta2 = 0 (no harm region); add mean gamma_c-hat
applied  E11  B5 caption: say why beta2 = 0 shows '—' and what mean gamma_c-hat is
applied  E12  B5 callout: identification has one dimension here; the null cell as reference; overshoot reaches the estimand
applied  E13  B6: limit 1 re-grounded on this design's records; E1 ranges and dev/notes read removed; limit 3 acknowledges G&H Section 6; limit 5 added; cert-note chunk removed
applied  E14  insert the 'Manuscript payload' section (export-payload chunk) before Provenance
applied  E15  provenance: item label
applied  E16  provenance: item value
  untouched_identical[b2-build]                    1
  untouched_identical[b2-table]                    1
  untouched_identical[b2-adaptive-note]            1
  untouched_identical[b2-consolidated]             1
  untouched_identical[b3-table]                    1
  untouched_identical[b5-chat-plot]                1
  untouched_identical[provenance-platform]         1
  fixed_sentence_1_count                           1
  fixed_sentence_2_count                           1
  sec_b3_label_count                               1
  forbidden[stable-pick]                           0
  forbidden[data-built]                            0
  forbidden[CERT_]                                 0
  forbidden[cert_upper]                            0
  forbidden[cert_joint]                            0
  forbidden[The wider figures]                     0
  forbidden["dev", "notes"]                        0
  forbidden[Classification is not the right lens]  0
  forbidden_ci[certif]                             0
  chunk[setup]                                     1
  chunk[b2-build]                                  1
  chunk[b2-table]                                  1
  chunk[b2-adaptive-note]                          1
  chunk[b2-consolidated]                           1
  chunk[b3-table]                                  1
  chunk[b5]                                        1
  chunk[b5-chat-plot]                              1
  chunk[b6]                                        1
  chunk[export-payload]                            1
  chunk[provenance]                                1
  chunk[provenance-platform]                       1
  result sha256 matches the spec: 15abd1209a4ea42432389a6e8d8c7ece1972c4a39f7f4632ac0d65a7b33c6866
WROTE quarto/GuoHe/guohe_supp_section.qmd  sha256 15abd1209a4ea42432389a6e8d8c7ece1972c4a39f7f4632ac0d65a7b33c6866
```

Both exit 0. The applier verified the input hash, all 16 edits by exact single-string matching
(E03 and E13 as SHA-256-verified spans), both fixed B3 sentences exactly once (sentence 1 in its
D1 "log-hazard-ratio scale" form), `{#sec-b3}` once, zero occurrences of every retired term and
the `dev/notes` read, every chunk label once, the seven untouched chunks byte-identical, and the
result SHA-256
`15abd1209a4ea42432389a6e8d8c7ece1972c4a39f7f4632ac0d65a7b33c6866` — equal to the spec value.
Independently re-hashed after the write: same value. GATE 2 PASS.

## Stage 3 — render and identity

**5.1 Edited render.** Exit 0, wall 11 s, zero `WARNING` lines. Neither `stopifnot()` guard in B6
fired. The payload exists at
`quarto/GuoHe/_payloads/guohe_supp_section/guohe_supp_section_payload.rds`.
The rendered HTML contains "Five limits belong on the record" (3 occurrences). Clean render. PASS.

**5.3 Payload SHA-256:**
`b5f4cc2e820d080ae54db9acb85c9a4ce284ab23e30715fa8dd904e5604ba19a`

**5.2 GATE 3b.** `Rscript dev/tasks/claude_cc_task_guohe_record_payload_2026-09-11_check.R`
exited 0 with empty stderr. First line:

```
GATE I: 186 comparisons, 0 outside tolerance
```

PASS. The remainder of `/tmp/guohe_payload_check.md`, verbatim:

### Payload element `bound_t7`

| cell | method | mean_est | bias | sd_emp | mean_margin | cover | cover_lo | cover_hi | n | mean_lower | beta2 | mean_lower_hr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| t7_beta2_00 | Naive | 0.1219 | 0.1219 | 0.1795 | 0.3049 | 0.8660 | 0.8504 | 0.8802 | 2000 | -0.1830 | 0.0000 | 0.8328 |
| t7_beta2_00 | G&H r=1/3 | 0.0404 | 0.0404 | 0.1788 | 0.3065 | 0.9355 | 0.9239 | 0.9455 | 2000 | -0.2661 | 0.0000 | 0.7664 |
| t7_beta2_00 | G&H r=1/12 | 0.0138 | 0.0138 | 0.1766 | 0.3060 | 0.9515 | 0.9412 | 0.9601 | 2000 | -0.2922 | 0.0000 | 0.7466 |
| t7_beta2_00 | G&H r=1/21 | 0.0120 | 0.0120 | 0.1765 | 0.3060 | 0.9525 | 0.9423 | 0.9610 | 2000 | -0.2940 | 0.0000 | 0.7453 |
| t7_beta2_00 | G&H r=1/30 | 0.0114 | 0.0114 | 0.1764 | 0.3060 | 0.9535 | 0.9434 | 0.9619 | 2000 | -0.2946 | 0.0000 | 0.7449 |
| t7_beta2_00 | MR (field) | 0.0298 | 0.0298 | 0.1856 | 0.3095 | 0.9395 | 0.9282 | 0.9491 | 2000 | -0.2797 | 0.0000 | 0.7560 |
| t7_beta2_00 | MR (IJ) | 0.0449 | 0.0449 | 0.1833 | 0.5772 | 0.9975 | 0.9942 | 0.9989 | 2000 | -0.5323 | 0.0000 | 0.5873 |
| t7_beta2_01 | Naive | 0.1981 | 0.1163 | 0.1864 | 0.3103 | 0.8700 | 0.8545 | 0.8840 | 2000 | -0.1121 | 0.1000 | 0.8939 |
| t7_beta2_01 | G&H r=1/3 | 0.1162 | 0.0344 | 0.1870 | 0.3081 | 0.9395 | 0.9282 | 0.9491 | 2000 | -0.1919 | 0.1000 | 0.8254 |
| t7_beta2_01 | G&H r=1/12 | 0.0894 | 0.0076 | 0.1839 | 0.3053 | 0.9515 | 0.9412 | 0.9601 | 2000 | -0.2159 | 0.1000 | 0.8058 |
| t7_beta2_01 | G&H r=1/21 | 0.0876 | 0.0058 | 0.1836 | 0.3051 | 0.9540 | 0.9439 | 0.9623 | 2000 | -0.2175 | 0.1000 | 0.8045 |
| t7_beta2_01 | G&H r=1/30 | 0.0870 | 0.0052 | 0.1836 | 0.3051 | 0.9540 | 0.9439 | 0.9623 | 2000 | -0.2180 | 0.1000 | 0.8041 |
| t7_beta2_01 | MR (field) | 0.1074 | 0.0256 | 0.1959 | 0.3099 | 0.9335 | 0.9217 | 0.9436 | 2000 | -0.2025 | 0.1000 | 0.8167 |
| t7_beta2_01 | MR (IJ) | 0.1222 | 0.0404 | 0.1920 | 0.5824 | 0.9980 | 0.9949 | 0.9992 | 2000 | -0.4602 | 0.1000 | 0.6312 |
| t7_beta2_02 | Naive | 0.2805 | 0.1115 | 0.1850 | 0.3174 | 0.8825 | 0.8676 | 0.8959 | 2000 | -0.0369 | 0.2000 | 0.9638 |
| t7_beta2_02 | G&H r=1/3 | 0.1993 | 0.0302 | 0.1860 | 0.3110 | 0.9445 | 0.9336 | 0.9537 | 2000 | -0.1117 | 0.2000 | 0.8943 |
| t7_beta2_02 | G&H r=1/12 | 0.1711 | 0.0021 | 0.1821 | 0.3053 | 0.9575 | 0.9477 | 0.9655 | 2000 | -0.1342 | 0.2000 | 0.8744 |
| t7_beta2_02 | G&H r=1/21 | 0.1693 | 0.0002 | 0.1818 | 0.3050 | 0.9580 | 0.9483 | 0.9659 | 2000 | -0.1357 | 0.2000 | 0.8731 |
| t7_beta2_02 | G&H r=1/30 | 0.1686 | -0.0004 | 0.1817 | 0.3048 | 0.9580 | 0.9483 | 0.9659 | 2000 | -0.1362 | 0.2000 | 0.8726 |
| t7_beta2_02 | MR (field) | 0.1941 | 0.0250 | 0.1951 | 0.3127 | 0.9375 | 0.9260 | 0.9473 | 2000 | -0.1186 | 0.2000 | 0.8882 |
| t7_beta2_02 | MR (IJ) | 0.2078 | 0.0387 | 0.1917 | 0.5905 | 0.9985 | 0.9956 | 0.9995 | 2000 | -0.3828 | 0.2000 | 0.6820 |
| t7_beta2_03 | Naive | 0.3638 | 0.1029 | 0.1886 | 0.3226 | 0.8950 | 0.8808 | 0.9077 | 2000 | 0.0412 | 0.3000 | 1.0421 |
| t7_beta2_03 | G&H r=1/3 | 0.2833 | 0.0224 | 0.1908 | 0.3131 | 0.9510 | 0.9406 | 0.9596 | 2000 | -0.0299 | 0.3000 | 0.9706 |
| t7_beta2_03 | G&H r=1/12 | 0.2535 | -0.0073 | 0.1854 | 0.3048 | 0.9610 | 0.9516 | 0.9686 | 2000 | -0.0513 | 0.3000 | 0.9500 |
| t7_beta2_03 | G&H r=1/21 | 0.2515 | -0.0093 | 0.1850 | 0.3043 | 0.9610 | 0.9516 | 0.9686 | 2000 | -0.0528 | 0.3000 | 0.9486 |
| t7_beta2_03 | G&H r=1/30 | 0.2508 | -0.0100 | 0.1848 | 0.3041 | 0.9615 | 0.9521 | 0.9691 | 2000 | -0.0533 | 0.3000 | 0.9481 |
| t7_beta2_03 | MR (field) | 0.2823 | 0.0214 | 0.2025 | 0.3139 | 0.9370 | 0.9255 | 0.9468 | 2000 | -0.0316 | 0.3000 | 0.9689 |
| t7_beta2_03 | MR (IJ) | 0.2942 | 0.0333 | 0.1979 | 0.5978 | 0.9985 | 0.9956 | 0.9995 | 2000 | -0.3036 | 0.3000 | 0.7381 |
| t7_beta2_04 | Naive | 0.4512 | 0.0937 | 0.1949 | 0.3255 | 0.9010 | 0.8871 | 0.9133 | 2000 | 0.1258 | 0.4000 | 1.1340 |
| t7_beta2_04 | G&H r=1/3 | 0.3713 | 0.0138 | 0.1972 | 0.3151 | 0.9535 | 0.9434 | 0.9619 | 2000 | 0.0562 | 0.4000 | 1.0578 |
| t7_beta2_04 | G&H r=1/12 | 0.3391 | -0.0184 | 0.1912 | 0.3049 | 0.9650 | 0.9560 | 0.9722 | 2000 | 0.0343 | 0.4000 | 1.0349 |
| t7_beta2_04 | G&H r=1/21 | 0.3370 | -0.0206 | 0.1907 | 0.3042 | 0.9665 | 0.9577 | 0.9735 | 2000 | 0.0327 | 0.4000 | 1.0333 |
| t7_beta2_04 | G&H r=1/30 | 0.3362 | -0.0213 | 0.1906 | 0.3040 | 0.9670 | 0.9582 | 0.9740 | 2000 | 0.0322 | 0.4000 | 1.0327 |
| t7_beta2_04 | MR (field) | 0.3747 | 0.0172 | 0.2090 | 0.3159 | 0.9405 | 0.9293 | 0.9500 | 2000 | 0.0587 | 0.4000 | 1.0605 |
| t7_beta2_04 | MR (IJ) | 0.3851 | 0.0276 | 0.2041 | 0.6013 | 0.9990 | 0.9964 | 0.9997 | 2000 | -0.2162 | 0.4000 | 0.8056 |
| t7_beta2_05 | Naive | 0.5357 | 0.0764 | 0.1965 | 0.3289 | 0.9180 | 0.9052 | 0.9292 | 2000 | 0.2068 | 0.5000 | 1.2297 |
| t7_beta2_05 | G&H r=1/3 | 0.4572 | -0.0021 | 0.1989 | 0.3187 | 0.9610 | 0.9516 | 0.9686 | 2000 | 0.1386 | 0.5000 | 1.1486 |
| t7_beta2_05 | G&H r=1/12 | 0.4220 | -0.0373 | 0.1923 | 0.3061 | 0.9715 | 0.9633 | 0.9779 | 2000 | 0.1159 | 0.5000 | 1.1229 |
| t7_beta2_05 | G&H r=1/21 | 0.4196 | -0.0397 | 0.1918 | 0.3053 | 0.9725 | 0.9644 | 0.9788 | 2000 | 0.1143 | 0.5000 | 1.1210 |
| t7_beta2_05 | G&H r=1/30 | 0.4188 | -0.0406 | 0.1916 | 0.3051 | 0.9730 | 0.9649 | 0.9792 | 2000 | 0.1137 | 0.5000 | 1.1204 |
| t7_beta2_05 | MR (field) | 0.4661 | 0.0068 | 0.2119 | 0.3196 | 0.9510 | 0.9406 | 0.9596 | 2000 | 0.1465 | 0.5000 | 1.1577 |
| t7_beta2_05 | MR (IJ) | 0.4746 | 0.0153 | 0.2068 | 0.6081 | 0.9990 | 0.9964 | 0.9997 | 2000 | -0.1335 | 0.5000 | 0.8751 |

### Payload element `ident_t7`

| cell | beta2 | chat_min | chat_q25 | chat_med | chat_q75 | chat_max | p_chat_30 | spec_mean | ppv_mean | n_sel_mean | n_sel_med | gamma_mean | dilution | phat_med | phat_q25 | phat_q75 | phat_mean | M_eff_med | sd_marg | sd_err |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| t7_beta2_00 | 0.0000 | 30 | 32.3553 | 40.0314 | 52.3669 | 59.9990 | 0.0625 | NA | NA | 212.3965 | 202 | 0.0000 | 0.0000 | 0.1956 | 0.1350 | 0.2757 | 0.2211 | 2.1205 | 0.1856 | 0.1856 |
| t7_beta2_01 | 0.1000 | 30 | 31.1927 | 36.3286 | 48.7361 | 59.9979 | 0.0905 | 0.7936 | 0.7881 | 201.8010 | 183 | 0.0818 | 0.0182 | 0.2092 | 0.1421 | 0.2995 | 0.2340 | 2.1338 | 0.1959 | 0.1905 |
| t7_beta2_02 | 0.2000 | 30 | 30.6701 | 33.5437 | 42.0559 | 59.9967 | 0.1140 | 0.8478 | 0.8357 | 188.3420 | 171 | 0.1691 | 0.0309 | 0.2290 | 0.1586 | 0.3243 | 0.2555 | 2.1404 | 0.1951 | 0.1879 |
| t7_beta2_03 | 0.3000 | 30 | 30.3824 | 32.2793 | 37.5390 | 59.9774 | 0.1455 | 0.8899 | 0.8743 | 177.8440 | 165 | 0.2608 | 0.0392 | 0.2442 | 0.1659 | 0.3475 | 0.2726 | 2.1637 | 0.2025 | 0.1930 |
| t7_beta2_04 | 0.4000 | 30 | 30.2156 | 31.4403 | 35.2615 | 59.9758 | 0.1675 | 0.9124 | 0.8976 | 171.8925 | 161 | 0.3575 | 0.0425 | 0.2718 | 0.1828 | 0.3926 | 0.3010 | 2.1968 | 0.2090 | 0.1969 |
| t7_beta2_05 | 0.5000 | 30 | 30.0767 | 30.9360 | 33.4944 | 59.9169 | 0.2130 | 0.9388 | 0.9237 | 165.2535 | 158 | 0.4593 | 0.0407 | 0.3088 | 0.2039 | 0.4417 | 0.3346 | 2.2348 | 0.2119 | 0.2021 |

### Payload element `comp_t7`

| cell | beta2 | method | cover | cover_lo | cover_hi | mean_upper | mean_upper_hr |
|---|---|---|---|---|---|---|---|
| t7_beta2_00 | 0.0000 | naive | 0.8625 | 0.8467 | 0.8769 | 0.2120 | 1.2362 |
| t7_beta2_00 | 0.0000 | ij | 1.0000 | 0.9981 | 1.0000 | 0.5705 | 1.7691 |
| t7_beta2_00 | 0.0000 | field | 0.9260 | 0.9137 | 0.9367 | 0.2868 | 1.3322 |
| t7_beta2_00 | 0.0000 | field_s | 0.9280 | 0.9158 | 0.9385 | 0.2922 | 1.3394 |
| t7_beta2_01 | 0.1000 | naive | 0.8810 | 0.8661 | 0.8945 | 0.2170 | 1.2423 |
| t7_beta2_01 | 0.1000 | ij | 0.9990 | 0.9964 | 0.9997 | 0.5634 | 1.7567 |
| t7_beta2_01 | 0.1000 | field | 0.9350 | 0.9233 | 0.9450 | 0.2897 | 1.3360 |
| t7_beta2_01 | 0.1000 | field_s | 0.9395 | 0.9282 | 0.9491 | 0.2915 | 1.3384 |
| t7_beta2_02 | 0.2000 | naive | 0.8950 | 0.8808 | 0.9077 | 0.2162 | 1.2414 |
| t7_beta2_02 | 0.2000 | ij | 0.9990 | 0.9964 | 0.9997 | 0.5496 | 1.7326 |
| t7_beta2_02 | 0.2000 | field | 0.9385 | 0.9271 | 0.9482 | 0.2859 | 1.3310 |
| t7_beta2_02 | 0.2000 | field_s | 0.9380 | 0.9266 | 0.9478 | 0.2823 | 1.3262 |
| t7_beta2_03 | 0.3000 | naive | 0.8960 | 0.8819 | 0.9086 | 0.2120 | 1.2362 |
| t7_beta2_03 | 0.3000 | ij | 0.9985 | 0.9956 | 0.9995 | 0.5342 | 1.7060 |
| t7_beta2_03 | 0.3000 | field | 0.9395 | 0.9282 | 0.9491 | 0.2775 | 1.3198 |
| t7_beta2_03 | 0.3000 | field_s | 0.9410 | 0.9298 | 0.9505 | 0.2713 | 1.3116 |
| t7_beta2_04 | 0.4000 | naive | 0.9135 | 0.9004 | 0.9250 | 0.2301 | 1.2588 |
| t7_beta2_04 | 0.4000 | ij | 0.9995 | 0.9972 | 0.9999 | 0.5444 | 1.7235 |
| t7_beta2_04 | 0.4000 | field | 0.9500 | 0.9396 | 0.9587 | 0.2892 | 1.3354 |
| t7_beta2_04 | 0.4000 | field_s | 0.9480 | 0.9374 | 0.9569 | 0.2827 | 1.3266 |
| t7_beta2_05 | 0.5000 | naive | 0.9170 | 0.9041 | 0.9283 | 0.2306 | 1.2594 |
| t7_beta2_05 | 0.5000 | ij | 0.9985 | 0.9956 | 0.9995 | 0.5371 | 1.7111 |
| t7_beta2_05 | 0.5000 | field | 0.9385 | 0.9271 | 0.9482 | 0.2828 | 1.3269 |
| t7_beta2_05 | 0.5000 | field_s | 0.9365 | 0.9250 | 0.9464 | 0.2753 | 1.3169 |

### Payload element `joint_t7`

| cell | beta2 | joint_s | joint_s_lo | joint_s_hi | joint | joint_lo | joint_hi | corr_mean | corr_s_mean | pair_lower_H | pair_lower_H_hr | pair_upper_Hc_s | pair_upper_Hc_s_hr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| t7_beta2_00 | 0.0000 | 0.9330 | 0.9212 | 0.9431 | 0.9310 | 0.9190 | 0.9413 | 0.0183 | 0.0153 | -0.3431 | 0.7095 | 0.3536 | 1.4241 |
| t7_beta2_01 | 0.1000 | 0.9330 | 0.9212 | 0.9431 | 0.9320 | 0.9201 | 0.9422 | 0.0212 | 0.0200 | -0.2654 | 0.7669 | 0.3512 | 1.4208 |
| t7_beta2_02 | 0.2000 | 0.9420 | 0.9309 | 0.9514 | 0.9410 | 0.9298 | 0.9505 | 0.0212 | 0.0215 | -0.1817 | 0.8338 | 0.3394 | 1.4041 |
| t7_beta2_03 | 0.3000 | 0.9460 | 0.9352 | 0.9551 | 0.9435 | 0.9325 | 0.9528 | 0.0177 | 0.0194 | -0.0945 | 0.9098 | 0.3271 | 1.3869 |
| t7_beta2_04 | 0.4000 | 0.9460 | 0.9352 | 0.9551 | 0.9465 | 0.9358 | 0.9555 | 0.0150 | 0.0170 | -0.0035 | 0.9965 | 0.3372 | 1.4010 |
| t7_beta2_05 | 0.5000 | 0.9415 | 0.9303 | 0.9510 | 0.9425 | 0.9314 | 0.9519 | 0.0091 | 0.0116 | 0.0842 | 1.0879 | 0.3288 | 1.3893 |

### Payload element `naive_convention`

| cell | beta2 | cover_replication_se | cover_engine_se | mean_se_replication | mean_se_engine |
|---|---|---|---|---|---|
| t7_beta2_00 | 0.0000 | 0.8660 | 0.8620 | 0.1854 | 0.1838 |
| t7_beta2_01 | 0.1000 | 0.8700 | 0.8655 | 0.1886 | 0.1870 |
| t7_beta2_02 | 0.2000 | 0.8825 | 0.8800 | 0.1930 | 0.1913 |
| t7_beta2_03 | 0.3000 | 0.8950 | 0.8920 | 0.1961 | 0.1944 |
| t7_beta2_04 | 0.4000 | 0.9010 | 0.8995 | 0.1979 | 0.1959 |
| t7_beta2_05 | 0.5000 | 0.9180 | 0.9140 | 0.2000 | 0.1981 |

### Payload element `pre16`

| cell | method | cover | cover_lo | cover_hi | margin | bias | n | design |
|---|---|---|---|---|---|---|---|---|
| t35_beta2_00 | Naive | 0.9015 | 0.8877 | 0.9138 | 0.3067 | 0.1017 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_00 | G&H r=1/12 | 0.9460 | 0.9352 | 0.9551 | 0.2651 | 0.0029 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_00 | G&H r=1/30 | 0.9490 | 0.9385 | 0.9578 | 0.2649 | 0.0008 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_00 | MR (field) | 0.9445 | 0.9336 | 0.9537 | 0.3024 | 0.0103 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_00 | MR (IJ) | 0.9985 | 0.9956 | 0.9995 | 0.5015 | 0.0280 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_01 | Naive | 0.9120 | 0.8988 | 0.9236 | 0.3049 | 0.0906 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_01 | G&H r=1/12 | 0.9525 | 0.9423 | 0.9610 | 0.2630 | -0.0063 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_01 | G&H r=1/30 | 0.9555 | 0.9456 | 0.9637 | 0.2627 | -0.0084 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_01 | MR (field) | 0.9485 | 0.9379 | 0.9574 | 0.3014 | 0.0036 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_01 | MR (IJ) | 0.9985 | 0.9956 | 0.9995 | 0.5016 | 0.0200 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_02 | Naive | 0.9150 | 0.9020 | 0.9264 | 0.3051 | 0.0768 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_02 | G&H r=1/12 | 0.9480 | 0.9374 | 0.9569 | 0.2616 | -0.0171 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_02 | G&H r=1/30 | 0.9495 | 0.9390 | 0.9583 | 0.2613 | -0.0198 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_02 | MR (field) | 0.9445 | 0.9336 | 0.9537 | 0.3039 | 0.0034 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_02 | MR (IJ) | 0.9990 | 0.9964 | 0.9997 | 0.5116 | 0.0155 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_03 | Naive | 0.9105 | 0.8972 | 0.9222 | 0.3055 | 0.0553 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_03 | G&H r=1/12 | 0.9450 | 0.9341 | 0.9542 | 0.2605 | -0.0352 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_03 | G&H r=1/30 | 0.9465 | 0.9358 | 0.9555 | 0.2600 | -0.0384 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_03 | MR (field) | 0.9365 | 0.9250 | 0.9464 | 0.3068 | -0.0030 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_03 | MR (IJ) | 0.9995 | 0.9972 | 0.9999 | 0.5240 | 0.0045 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_04 | Naive | 0.9205 | 0.9078 | 0.9316 | 0.3055 | 0.0295 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_04 | G&H r=1/12 | 0.9500 | 0.9396 | 0.9587 | 0.2599 | -0.0570 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_04 | G&H r=1/30 | 0.9525 | 0.9423 | 0.9610 | 0.2591 | -0.0609 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_04 | MR (field) | 0.9360 | 0.9244 | 0.9459 | 0.3091 | -0.0117 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_04 | MR (IJ) | 0.9970 | 0.9935 | 0.9986 | 0.5390 | -0.0094 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_05 | Naive | 0.9335 | 0.9217 | 0.9436 | 0.3059 | 0.0148 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_05 | G&H r=1/12 | 0.9530 | 0.9428 | 0.9614 | 0.2598 | -0.0676 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_05 | G&H r=1/30 | 0.9540 | 0.9439 | 0.9623 | 0.2588 | -0.0722 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_05 | MR (field) | 0.9415 | 0.9303 | 0.9510 | 0.3109 | -0.0113 | 2000 | Tables 3-5: two complementary subgroups |
| t35_beta2_05 | MR (IJ) | 0.9970 | 0.9935 | 0.9986 | 0.5558 | -0.0129 | 2000 | Tables 3-5: two complementary subgroups |
| t6_k02 | Naive | 0.8970 | 0.8829 | 0.9096 | 0.3069 | 0.1012 | 2000 | Table 6: k disjoint null subgroups |
| t6_k02 | G&H r=1/12 | 0.9480 | 0.9374 | 0.9569 | 0.2651 | 0.0021 | 2000 | Table 6: k disjoint null subgroups |
| t6_k02 | G&H r=1/30 | 0.9500 | 0.9396 | 0.9587 | 0.2650 | 0.0001 | 2000 | Table 6: k disjoint null subgroups |
| t6_k02 | MR (field) | 0.9440 | 0.9330 | 0.9533 | 0.3020 | 0.0081 | 2000 | Table 6: k disjoint null subgroups |
| t6_k02 | MR (IJ) | 0.9980 | 0.9949 | 0.9992 | 0.5009 | 0.0264 | 2000 | Table 6: k disjoint null subgroups |
| t6_k06 | Naive | 0.7465 | 0.7270 | 0.7651 | 0.3086 | 0.2370 | 2000 | Table 6: k disjoint null subgroups |
| t6_k06 | G&H r=1/12 | 0.9445 | 0.9336 | 0.9537 | 0.2180 | 0.0075 | 2000 | Table 6: k disjoint null subgroups |
| t6_k06 | G&H r=1/30 | 0.9475 | 0.9368 | 0.9564 | 0.2178 | 0.0040 | 2000 | Table 6: k disjoint null subgroups |
| t6_k06 | MR (field) | 0.9410 | 0.9298 | 0.9505 | 0.2945 | 0.0302 | 2000 | Table 6: k disjoint null subgroups |
| t6_k06 | MR (IJ) | 0.9870 | 0.9810 | 0.9911 | 0.3984 | 0.0715 | 2000 | Table 6: k disjoint null subgroups |
| t6_k10 | Naive | 0.5990 | 0.5774 | 0.6203 | 0.3092 | 0.2894 | 2000 | Table 6: k disjoint null subgroups |
| t6_k10 | G&H r=1/12 | 0.9485 | 0.9379 | 0.9574 | 0.2025 | 0.0074 | 2000 | Table 6: k disjoint null subgroups |
| t6_k10 | G&H r=1/30 | 0.9525 | 0.9423 | 0.9610 | 0.2023 | 0.0037 | 2000 | Table 6: k disjoint null subgroups |
| t6_k10 | MR (field) | 0.9440 | 0.9330 | 0.9533 | 0.2907 | 0.0382 | 2000 | Table 6: k disjoint null subgroups |
| t6_k10 | MR (IJ) | 0.9820 | 0.9752 | 0.9870 | 0.3666 | 0.0884 | 2000 | Table 6: k disjoint null subgroups |
| t6_k12 | Naive | 0.5555 | 0.5336 | 0.5772 | 0.3095 | 0.3053 | 2000 | Table 6: k disjoint null subgroups |
| t6_k12 | G&H r=1/12 | 0.9465 | 0.9358 | 0.9555 | 0.1980 | 0.0056 | 2000 | Table 6: k disjoint null subgroups |
| t6_k12 | G&H r=1/30 | 0.9490 | 0.9385 | 0.9578 | 0.1978 | 0.0020 | 2000 | Table 6: k disjoint null subgroups |
| t6_k12 | MR (field) | 0.9420 | 0.9309 | 0.9514 | 0.2898 | 0.0410 | 2000 | Table 6: k disjoint null subgroups |
| t6_k12 | MR (IJ) | 0.9715 | 0.9633 | 0.9779 | 0.3584 | 0.0931 | 2000 | Table 6: k disjoint null subgroups |
| t7_beta2_00 | Naive | 0.8660 | 0.8504 | 0.8802 | 0.3049 | 0.1219 | 2000 | Table 7: nested continuum |
| t7_beta2_00 | G&H r=1/12 | 0.9515 | 0.9412 | 0.9601 | 0.3060 | 0.0138 | 2000 | Table 7: nested continuum |
| t7_beta2_00 | G&H r=1/30 | 0.9535 | 0.9434 | 0.9619 | 0.3060 | 0.0114 | 2000 | Table 7: nested continuum |
| t7_beta2_00 | MR (field) | 0.9395 | 0.9282 | 0.9491 | 0.3095 | 0.0298 | 2000 | Table 7: nested continuum |
| t7_beta2_00 | MR (IJ) | 0.9975 | 0.9942 | 0.9989 | 0.5772 | 0.0449 | 2000 | Table 7: nested continuum |
| t7_beta2_01 | Naive | 0.8700 | 0.8545 | 0.8840 | 0.3103 | 0.1163 | 2000 | Table 7: nested continuum |
| t7_beta2_01 | G&H r=1/12 | 0.9515 | 0.9412 | 0.9601 | 0.3053 | 0.0076 | 2000 | Table 7: nested continuum |
| t7_beta2_01 | G&H r=1/30 | 0.9540 | 0.9439 | 0.9623 | 0.3051 | 0.0052 | 2000 | Table 7: nested continuum |
| t7_beta2_01 | MR (field) | 0.9335 | 0.9217 | 0.9436 | 0.3099 | 0.0256 | 2000 | Table 7: nested continuum |
| t7_beta2_01 | MR (IJ) | 0.9980 | 0.9949 | 0.9992 | 0.5824 | 0.0404 | 2000 | Table 7: nested continuum |
| t7_beta2_02 | Naive | 0.8825 | 0.8676 | 0.8959 | 0.3174 | 0.1115 | 2000 | Table 7: nested continuum |
| t7_beta2_02 | G&H r=1/12 | 0.9575 | 0.9477 | 0.9655 | 0.3053 | 0.0021 | 2000 | Table 7: nested continuum |
| t7_beta2_02 | G&H r=1/30 | 0.9580 | 0.9483 | 0.9659 | 0.3048 | -0.0004 | 2000 | Table 7: nested continuum |
| t7_beta2_02 | MR (field) | 0.9375 | 0.9260 | 0.9473 | 0.3127 | 0.0250 | 2000 | Table 7: nested continuum |
| t7_beta2_02 | MR (IJ) | 0.9985 | 0.9956 | 0.9995 | 0.5905 | 0.0387 | 2000 | Table 7: nested continuum |
| t7_beta2_03 | Naive | 0.8950 | 0.8808 | 0.9077 | 0.3226 | 0.1029 | 2000 | Table 7: nested continuum |
| t7_beta2_03 | G&H r=1/12 | 0.9610 | 0.9516 | 0.9686 | 0.3048 | -0.0073 | 2000 | Table 7: nested continuum |
| t7_beta2_03 | G&H r=1/30 | 0.9615 | 0.9521 | 0.9691 | 0.3041 | -0.0100 | 2000 | Table 7: nested continuum |
| t7_beta2_03 | MR (field) | 0.9370 | 0.9255 | 0.9468 | 0.3139 | 0.0214 | 2000 | Table 7: nested continuum |
| t7_beta2_03 | MR (IJ) | 0.9985 | 0.9956 | 0.9995 | 0.5978 | 0.0333 | 2000 | Table 7: nested continuum |
| t7_beta2_04 | Naive | 0.9010 | 0.8871 | 0.9133 | 0.3255 | 0.0937 | 2000 | Table 7: nested continuum |
| t7_beta2_04 | G&H r=1/12 | 0.9650 | 0.9560 | 0.9722 | 0.3049 | -0.0184 | 2000 | Table 7: nested continuum |
| t7_beta2_04 | G&H r=1/30 | 0.9670 | 0.9582 | 0.9740 | 0.3040 | -0.0213 | 2000 | Table 7: nested continuum |
| t7_beta2_04 | MR (field) | 0.9405 | 0.9293 | 0.9500 | 0.3159 | 0.0172 | 2000 | Table 7: nested continuum |
| t7_beta2_04 | MR (IJ) | 0.9990 | 0.9964 | 0.9997 | 0.6013 | 0.0276 | 2000 | Table 7: nested continuum |
| t7_beta2_05 | Naive | 0.9180 | 0.9052 | 0.9292 | 0.3289 | 0.0764 | 2000 | Table 7: nested continuum |
| t7_beta2_05 | G&H r=1/12 | 0.9715 | 0.9633 | 0.9779 | 0.3061 | -0.0373 | 2000 | Table 7: nested continuum |
| t7_beta2_05 | G&H r=1/30 | 0.9730 | 0.9649 | 0.9792 | 0.3051 | -0.0406 | 2000 | Table 7: nested continuum |
| t7_beta2_05 | MR (field) | 0.9510 | 0.9406 | 0.9596 | 0.3196 | 0.0068 | 2000 | Table 7: nested continuum |
| t7_beta2_05 | MR (IJ) | 0.9990 | 0.9964 | 0.9997 | 0.6081 | 0.0153 | 2000 | Table 7: nested continuum |

### Payload meta

- `nrep`: 2000
- `mcse`: 0.004873397
- `draws`: 5000
- `field_R_out`: 1000
- `field_R_in`: 500
- `r_fixed`: 1/12 (G&H r=1/12: r2 / gh_r2 columns)
- `scale`: oriented log-HR; *_hr = exp(mean log-scale bound)
- `complement_truth`: 0
- `head_at_render`: d975a51eda0420775f6b6e4f55b424380e5657e8
- `built_at`: 2026-09-11 10:45:35
- `R`: R version 4.6.1 (2026-06-24)
- `forestsearch_version`: 0.3.5

Input bundle MD5s:

- `mr_field_vs_guohe_t35_beta2_00.rds` edd5ed5f1d9f621020115bb6c48ad6ff
- `mr_field_vs_guohe_t35_beta2_01.rds` c31b304b97a62d2a06e446d9f3b48cdf
- `mr_field_vs_guohe_t35_beta2_02.rds` d8b0bec4f6761cf832893a49d47c79ef
- `mr_field_vs_guohe_t35_beta2_03.rds` df144ae66e093279b4f258cf618e344e
- `mr_field_vs_guohe_t35_beta2_04.rds` 52bb369059723f1142e728f85d976027
- `mr_field_vs_guohe_t35_beta2_05.rds` 8640fa09a5c1d7ebaa32430179a11cab
- `mr_field_vs_guohe_t6_k02.rds` 4b8152f1d77a47faeca271852f031dd1
- `mr_field_vs_guohe_t6_k06.rds` dc9b76e3b230376f62217aaedfe8a50b
- `mr_field_vs_guohe_t6_k10.rds` 8e02fcfab4181803d48c01eb52f83506
- `mr_field_vs_guohe_t6_k12.rds` 95ff0ba61dc9c2b927a856d9947806c5
- `mr_field_vs_guohe_t7_beta2_00.rds` d60527e8f769f9302ed25c9c7651de52
- `mr_field_vs_guohe_t7_beta2_01.rds` 8263d0d69dc064974c1505105bca8317
- `mr_field_vs_guohe_t7_beta2_02.rds` 5afe483e7f83bd6528980f1c9eb5f079
- `mr_field_vs_guohe_t7_beta2_03.rds` 145cfbcdbcf6391edd59fa8c54a01783
- `mr_field_vs_guohe_t7_beta2_04.rds` c4c80365c7cd9fcc042d98b055145c51
- `mr_field_vs_guohe_t7_beta2_05.rds` 699290395e1752823e4327221d58ec43
- `guohe_repro_t7_beta2_00.rds` 339f8b85fa1f5234f8b052d791e7f6dc
- `guohe_repro_t7_beta2_01.rds` 8f22d714b6a3b00502eb690dbf24e77d
- `guohe_repro_t7_beta2_02.rds` 67f2e08227e7f481650962e0f65ab560
- `guohe_repro_t7_beta2_03.rds` 6df7dbfbb0f1bc2f67c6c4f8f8acfcec
- `guohe_repro_t7_beta2_04.rds` 9a288e513de74ddebc9b23d586367e92
- `guohe_repro_t7_beta2_05.rds` 112e8d6ad25096af904f3fe2a2e801ff
- `mr_field_complement_vs_guohe_t7_beta2_00.rds` 77cd1c1b26ec8ebf41eedf5c351284b7
- `mr_field_complement_vs_guohe_t7_beta2_01.rds` f6d6132afde44e145ed22b83bbe79797
- `mr_field_complement_vs_guohe_t7_beta2_02.rds` f1aa2cf4edf00470f5f87fd475bf2a7a
- `mr_field_complement_vs_guohe_t7_beta2_03.rds` 8fe9bbc3a3cf5d6bbbdc6bc2e287bcb8
- `mr_field_complement_vs_guohe_t7_beta2_04.rds` 2fa534fc55a3313743e33f77930cbfb2
- `mr_field_complement_vs_guohe_t7_beta2_05.rds` 4b534b4b8743f06ddbfef1296d8f8ccc

## OPEN ITEMS

None. Every gate passed on the first attempt:

- §1 provenance and first commit — PASS (branch, HEAD, no dirt, three companion hashes matched).
- GATE 0 — PASS (input hash, 30 bundles present, all eight §2.3 quotes found).
- GATE 1 baseline render — PASS (exit 0, no warnings).
- GATE 2 applier — PASS (both runs exit 0, output hash equals the spec).
- GATE 3 edited render — PASS (exit 0, no warnings, payload written, B6 sentinel present, no guard fired).
- GATE 3b payload check — PASS (exit 0, 186 comparisons, 0 outside tolerance).

Two facts recorded without interpretation, neither a gap:

1. The payload's `head_at_render` is `d975a51eda0420775f6b6e4f55b424380e5657e8` — the §1
   task-files commit, which was HEAD when the edited record rendered. The record and payload were
   committed afterwards, at `bd1ce4bf`.
2. BASE is `1d9401cb` itself: HEAD was exactly the required commit at kickoff, so
   `git merge-base --is-ancestor 1d9401cb HEAD` held trivially.

No task is proposed; nothing is blocked.

## Commits, BASE..HEAD

```
bd1ce4bf guohe record: 16 exact-string corrections and the manuscript payload
d975a51e guohe record/payload: task document and companions as received
```

(The commit of this report follows, and is not listed above.)

Files touched, all by explicit path: the four task files under `dev/tasks/`,
`quarto/GuoHe/guohe_supp_section.qmd`,
`quarto/GuoHe/_payloads/guohe_supp_section/guohe_supp_section_payload.rds`, and this report.
No `R/` change, no push, no touch of `feature/glm-extension-mac`, no manuscript text.
