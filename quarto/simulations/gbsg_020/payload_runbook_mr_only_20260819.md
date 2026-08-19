# Standalone MR-only payload campaign — runbook rev 3 (2026-08-19)

REV 3: campaign_tag = "standalone1" added to all six drivers. The campaign writes ONLY
to *_seedtab_s500_standalone1 directories: nothing pre-existing is touched, read, or
reused (the legacy knoise0 dir with FB contents stays exactly as-is for its other
summaries), resumability cannot pick up stale cells, and the portable set is one glob.

## Campaign: 6 drivers, 10 cells, n_sims = 500 -> 5,000 fresh replicates (all FS-only, MR-only)
h10_knoise0 / h10_knoise3_fsonly / h10_knoise6: n = 500,1000 (2 cells each)
h15_knoise0_fsonly: n = 500,1000 (2 cells) | h15/h20_knoise3_fsonly_n1500: 1 cell each
Cross-checked: the 13 artifacts the four fragments require are all inside the 16 the
campaign produces (cells + grids), path-for-path. Fragments are retargeted to the
*_standalone1 tags (rev 3 copies delivered — use these, they supersede all earlier).

## Order of operations
1. Kill any running legacy-chain renders (byte-duplicate outputs in wrong namespace).
2. Place all six rev-3 drivers in quarto/simulations/gbsg_020/ (overwrite rev 2).
3. Run the 6-render loop (terminal or CC), quarto binary:
   /usr/lib/rstudio/resources/app/bin/quarto/bin/quarto
4. Port: cp -r mr_sweep/*_standalone1 <manuscript>/_payloads/   (one glob, six dirs)
5. Drop rev-3 fragments into the manuscript tree; render. Captions ("500 replicates
   per size") remain correct. _supp_sim_synthesis stays on its FB t1_t2 bundles — not
   part of this campaign.

## Notes
- The knoise3_fsonly driver now writes its OWN standalone dir (no longer shares the
  committed template's dir — that earlier caveat is dissolved).
- campaign_tag = "" in any driver reproduces the legacy tag exactly (documented inline).
- Repo hygiene: the new dirs are untracked until you choose a results checkpoint
  commit; the drivers themselves can be committed scoped whenever convenient.
