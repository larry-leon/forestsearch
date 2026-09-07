# scripts_mdf1 — session scripts behind the Mac continuous/field records (2026-09-07)

Committed as records of what ran, not as a toolkit: each script hardcodes the session's scratchpad path (`SP=/private/tmp/claude-501/.../scratchpad`) and the repository path, and was run from `quarto/simulations/actg175/continuous/`. Cited by `REPORT_continuous_field_stage0_2026-09-07.md`, `REPORT_continuous_field_stage1_2026-09-07.md`, `REPORT_continuous_field_gate2_2026-09-07.md` and `REPORT_continuous_field_2026-09-07.md`.

Stage 2 campaign (`mdf1`):
- `run_cell.sh` — per-cell driver: batch 1 (sim 1–1000) → batch 2 (1001–2000) → combine render of `sim_fs_maxeffCons_mr_field_md_template.qmd`, each under the hard timeout, with the 1.5×-projection and cumulative-ceiling checks and the wall log (`MD N FBJOIN THRESH_SEC SAMPLE_MEM`).
- `tmo.sh` — portable hard timeout (Perl alarm; exit 124), used because macOS has no `timeout`.
- `mem_sampler.sh` — 5-second sampler of the summed RSS of all R/Quarto processes while the first render ran (peak memory in the Gate 2 record).
- `gate2_check.R` — the Gate 2 checker: completeness, save-guard paths, pairing proof by sim_id against the committed bundle (≤ 1e-8 relative, rules as term sets, enumerated rows classified as label tie / target move / MR-numerics / selection flip; written to `gate2_flips.txt` beside each cell's bundles), new columns finite, interval invariants, bound identities, γ range, p̂ range.
- `stage3_aggregate.R` — cross-cell aggregation from the combined bundles (Table-2 layout, bound location, joint pair, regime diagnostics, identity-scale display; the markdown tables and figures in the Stage 3 report). The summary document `summary_continuous_field_mdf1.qmd` carries the same code inline.

Stage 0 / Stage 1 (cited by those records):
- `time_one_rep.R` — 0d single-replicate timing of the twin and the replicate-1 identity peek.
- `fixture_baseline_save.R`, `identity_postchange.R` — the 28 pre-/post-change `fs_sim_bias_coverage()` tables and their `identical()` comparison for the `scale` argument.
- `fixture_check.R` — the 14-point bias-coverage fixture check on the committed s7/map1 bundles.
- `id_renders.sh`, `identity_check.R` — the eight 5-replicate identity renders and their comparison with the committed bundles (field invariants included).
- `rep1_field.R`, `reselection_check.R` — replicate 1 under `ci_method = "field"` and the re-selection-map reproduction of the observed Ĥ on the unperturbed effects.
- `cal_renders.sh` — the worker calibration renders (13 and 10 workers at n = 500; 13 at n = 700).
