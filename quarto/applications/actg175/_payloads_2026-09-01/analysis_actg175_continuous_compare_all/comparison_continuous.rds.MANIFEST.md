# Manifest — comparison_continuous.rds (0.3.5 render, parked out of tree)

The payload itself exceeds GitHub's 100 MB hard limit and is therefore not tracked
(same exclusion the 08-31 reference set carries, commit `8ad2c3f6`). The original
file, byte-identical to the fingerprint below, is parked at
`~/fs_parked_2026-09-01/_payloads_2026-09-01/analysis_actg175_continuous_compare_all/comparison_continuous.rds`,
beside the 08-31 parked copy (`~/fs_parked_2026-08-31/…`).

| field | value |
|---|---|
| filename | `comparison_continuous.rds` |
| byte size | 111,266,019 |
| sha256 (computed in-place before the move) | `341fa8d9d575a5e70c05a4eab8b1191e3af83c4d5c151a8c1a5082fa7e210a25` |
| forestsearch_version | 0.3.5 — the object has no self-describing version field (top-level names `combos,fs,ci_tab,plots,plot_grid,plot_combined,plot_combined_subsets,combined_skip_reason,console,diagnostics,errors`; none inside `fs[[i]]`/`args_call_all` either, matching the 08-31 report §7's "no version field expected"); vintage attributed from the render: written by the `analysis_actg175_continuous_compare_all.qmd` render of 2026-09-01 under installed forestsearch 0.3.5 (`REPORT_applications_render_0.3.5_2026-09-01.md` §1.1, §2) |
| built_at | no such field in the object; file mtime `2026-09-01 14:23:23 PDT`, i.e. during that render's 14:18:16–14:23:30 window (render log, report §2) |

Comparison A verdict for this file (report §3.1): CLEAN — every COMPARED leaf
`identical()` against the parked 08-31 (0.3.1) original.
