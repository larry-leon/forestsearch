# REPORT — declcal re-run with the clinically specified null level c0 (`declcalc0`)

Task: `dev/tasks/TASK_declcal_c0_campaign_2026-09-22.md`. Branch `feature/glm-extension`, pop-os, 64 workers. No `R/` change. No push.

Full tables with Wilson 95% intervals for every entry: `quarto/simulations/gbsg_020/scripts_dinamr/logs/declcalc0_findings.txt`, written by `scripts_dinamr/declcalc0_findings.R`, with its object in `scripts_dinamr/declcalc0_findings.rds`. The tables below are copied from that file.

**Implied p\*.** Every κ̂ row in this report carries its implied p\* = 2·pnorm(κ̂) − 1. This footnote applies to all of them.<sup>†</sup>

<sup>†</sup> **A literal re-run of the fixed-p\* screen at an implied p\* needs `pconsistency.digits >= 5`.** The implied p\* values are quoted to 5 decimals. At the default `pconsistency.digits = 2`, the consistency proportion is rounded to 2 decimals before it is compared with p\*. For example, p\* = 0.99586 then acts as "rounded proportion = 1.00", so the executed cutoff is not κ̂.

## 0. Pins, build, commits

- **§0 assertion.** The installed package passes it, re-checked at report time:
  - `declaration_c0` is a formal of `forestsearch:::fs_mr_inference` with default `NULL`;
  - `forestsearch:::fs_declaration_calibration` accepts `c0`.
- **Installed build.** `forestsearch` 0.3.5.9000, R 4.6.1, built 2026-09-23 02:32:08 UTC. Every one of the 10 payload metas records this same build, host pop-os and 64 workers.
- **HEAD at each stage:**
  - smoke: `a46bf9b7`;
  - first Stage 2 run: `dbc14acb`;
  - resumed Stage 2 run: `ee64f11d`.
  - No commit in between touches `R/`.
- **Comparators.** The committed `results/declcal_{inull,power}_<cell>_res_1_2000.rds` at `81752681`, and `scripts_dinamr/logs/declcal_fixedk_practical.txt`.

**Commits (task §8):**

| # | commit | what |
|---|---|---|
| 1 | `1f7dc5b3` | task document |
| 2 | `a46bf9b7` | scripts: `declcalc0_run.R`, `declcalc0.sh`, cells files. The named-line changes are listed in each file's header |
| 3 | `dbc14acb` | Stage 1 smoke record (B5, 5 reps, identity gate green) |
| 4a | `ee64f11d` | B1 and B2 payloads and logs, from the stopped Stage 2 run. Committed by explicit paths; the partial B3 was not committed |
| 4b | `02852113` | Block B remainder: B3–B6 payloads and logs. Committed by the driver |
| 4c | `8a44bee3` | Block C: C1–C4 payloads and logs. Committed by the driver |
| 5 | `7ac9c277` | report (first version), findings script and tables, driver logs, resume cells files, closeout-script registration |
| 6 | `55125a6d` | closeout: `current_status.md` at pin `7ac9c277` |
| 5′ | the commit that revises this file | Wilson intervals in the primary tables, the three free checks (§5.5), tuned fixed-p\* labels; findings script extended |
| 6′ | the closeout commit after it | `current_status.md` regenerated at the new pin |

## 1. Run record

**Stopped Stage 2 run** (`logs/declcalc0_stage2_driver.log`):

- Started 19:35:17 at `dbc14acb`, with a hard cap of 6,300 s.
- B1 completed in 322 s and B2 in 394 s. Both passed the identity gate on 2000 of 2000 replicates.
- B3 had checkpointed 500 replicates (chunk 1 at 19:49) when the run stopped. It stopped at about 900 s elapsed.
- The driver log ends mid-cell, with no halt record and no gate line.

**Resumed Stage 2 run** (`logs/declcalc0_stage2_resume_driver.log`):

- **Driver.** The same `declcalc0.sh`, unchanged.
- **Cells.** They came from `scripts_dinamr/resume/declcalc0_{inull,power}.cells`:
  - the files keep the block basenames, so the campaign tags, payload paths and log paths are unchanged;
  - the inull file lists B3–B6, and the power file is a copy of the C1–C4 block file.
- **B3.** It restarted from replicate 1. The runner does not resume from checkpoints, so B3's partial payload was overwritten.
- **Caps:**
  - per replicate: 116 s;
  - per cell: 1% abort/error gate;
  - campaign hard cap: 5,400 s. That is the task's 6,300 s cap less the ~900 s the stopped run had already spent, so the cap still binds the campaign as a whole.
- **Result.** Started 21:15:08. **Completed 8 of 8 cells in 3,454 s**, with 0 aborts and 0 errors. The longest single replicate took 18.3 s, against the 116 s cap.

**Campaign compute.** The smoke plus about 900 s plus 3,454 s is **≈ 4,350 s ≈ 1.2 h**. The projection was 4,200 s, the hard cap 6,300 s, and Larry's ceiling 3 h.

**Per-cell elapsed (s):**

| B1 | B2 | B3 | B4 | B5 | B6 | C1 | C2 | C3 | C4 |
|---|---|---|---|---|---|---|---|---|---|
| 288 | 364 | 441 | 305 | 375 | 445 | 393 | 461 | 415 | 485 |

## 2. Gates, with measured values

| cell | ok reps | identity: reps disagreeing on any of 13 columns | Mstar_c0 monotone in c0 (all draws) | κ̂₀.₀₅(c0) ≤ κ̂₀.₀₅ | fidelity (declared_conv = search indicator) | abort / error |
|---|---|---|---|---|---|---|
| B1 | 2000 | 0 | 2000 / 2000 | 2000 / 2000 | 2000 / 2000 | 0 |
| B2 | 2000 | 0 | 2000 / 2000 | 2000 / 2000 | 2000 / 2000 | 0 |
| B3 | 2000 | 0 | 2000 / 2000 | 2000 / 2000 | 2000 / 2000 | 0 |
| B4 | 2000 | 0 | 2000 / 2000 | 2000 / 2000 | 2000 / 2000 | 0 |
| B5 | 2000 | 0 | 2000 / 2000 | 2000 / 2000 | 2000 / 2000 | 0 |
| B6 | 2000 | 0 | 2000 / 2000 | 2000 / 2000 | 2000 / 2000 | 0 |
| C1 | 2000 | 0 | 2000 / 2000 | 2000 / 2000 | 2000 / 2000 | 0 |
| C2 | 2000 | 0 | 2000 / 2000 | 2000 / 2000 | 2000 / 2000 | 0 |
| C3 | 2000 | 0 | 2000 / 2000 | 2000 / 2000 | 2000 / 2000 | 0 |
| C4 | 2000 | 0 | 2000 / 2000 | 2000 / 2000 | 2000 / 2000 | 0 |

- **Identity gate (§4): PASS on 20,000 of 20,000 replicates.** Thirteen columns match the committed `declcal` payloads exactly: `max_T_pre`, `max_T_post`, `G_pre`, `G_post`, `declared_conv`, `declared_conv_exact`, `Mstar_q90/q95/q99`, `kappa_hat_05/10` and `declared_cal05/10`. `Mstar_c0` is therefore the only new information.
- **Other gates:**
  - smoke (B5, 5 reps): identity green;
  - per-replicate cap: 0 trips;
  - per-cell 1% gate: 0 trips;
  - campaign hard cap: not reached.

## 3. Primary tables (§6.1)

**How to read these tables:**

- **Entries:** proportion x / 2000 with its Wilson 95% interval, for exact and approx rows alike.
- **"exact" rows** are this campaign: the per-replicate rule `max_T_pre >= κ̂_α(c0)`, with B_cal 500.
- **"approx" rows** sit directly beneath the exact rows. They come from `REPORT_declcal_c0_approx_2026-09-22` (`e9de87c7`): a fixed plug-in cutoff at the median κ̂(c0) over 40 captures at that n (B 2000), applied to `max_T_pre`.
- **"p\* 0.90 exec."** is `declared_conv`, the rounded rule on the post-reduction family.
- **"c0 = c2"** is the committed `declared_cal05` or `declared_cal10`.

### α 0.05

| cell | DGM | n | p\* 0.90 as executed | c0 0.70 | c0 0.75 | c0 0.80 | c0 0.85 | c0 = c2 (committed) |
|---|---|---|---|---|---|---|---|---|
| B1 | uniform benefit HR 0.657 | 500 | 0.0950 [0.0829, 0.1086] | 0.0080 [0.0049, 0.0130] | 0.0035 [0.0017, 0.0072] | 0.0010 [0.0003, 0.0036] | 0.0005 [0.0001, 0.0028] | 0.0005 [0.0001, 0.0028] |
| B1 approx | plug-in, fixed cutoff | 500 | — | 0.0065 [0.0038, 0.0111] | 0.0025 [0.0011, 0.0058] | 0.0015 [0.0005, 0.0044] | 0.0005 [0.0001, 0.0028] | — |
| B2 | uniform benefit HR 0.657 | 1000 | 0.0325 [0.0256, 0.0412] | 0.0065 [0.0038, 0.0111] | 0.0020 [0.0008, 0.0051] | 0.0005 [0.0001, 0.0028] | 0.0005 [0.0001, 0.0028] | 0.0000 [0.0000, 0.0019] |
| B2 approx | plug-in, fixed cutoff | 1000 | — | 0.0065 [0.0038, 0.0111] | 0.0030 [0.0014, 0.0065] | 0.0010 [0.0003, 0.0036] | 0.0005 [0.0001, 0.0028] | — |
| B3 | uniform benefit HR 0.657 | 1500 | 0.0070 [0.0042, 0.0117] | 0.0045 [0.0024, 0.0085] | 0.0030 [0.0014, 0.0065] | 0.0015 [0.0005, 0.0044] | 0.0000 [0.0000, 0.0019] | 0.0000 [0.0000, 0.0019] |
| B3 approx | plug-in, fixed cutoff | 1500 | — | 0.0040 [0.0020, 0.0079] | 0.0020 [0.0008, 0.0051] | 0.0015 [0.0005, 0.0044] | 0.0000 [0.0000, 0.0019] | — |
| B4 | uniform benefit HR 0.721 | 500 | 0.2325 [0.2145, 0.2515] | 0.0315 [0.0247, 0.0401] | 0.0130 [0.0089, 0.0190] | 0.0040 [0.0020, 0.0079] | 0.0010 [0.0003, 0.0036] | 0.0005 [0.0001, 0.0028] |
| B4 approx | plug-in, fixed cutoff | 500 | — | 0.0270 [0.0208, 0.0351] | 0.0080 [0.0049, 0.0130] | 0.0030 [0.0014, 0.0065] | 0.0010 [0.0003, 0.0036] | — |
| B5 | uniform benefit HR 0.721 | 1000 | 0.1110 [0.0980, 0.1255] | 0.0285 [0.0221, 0.0367] | 0.0100 [0.0065, 0.0154] | 0.0040 [0.0020, 0.0079] | 0.0010 [0.0003, 0.0036] | 0.0000 [0.0000, 0.0019] |
| B5 approx | plug-in, fixed cutoff | 1000 | — | 0.0255 [0.0194, 0.0334] | 0.0100 [0.0065, 0.0154] | 0.0035 [0.0017, 0.0072] | 0.0005 [0.0001, 0.0028] | — |
| B6 | uniform benefit HR 0.721 | 1500 | 0.0325 [0.0256, 0.0412] | 0.0210 [0.0156, 0.0283] | 0.0065 [0.0038, 0.0111] | 0.0020 [0.0008, 0.0051] | 0.0010 [0.0003, 0.0036] | 0.0000 [0.0000, 0.0019] |
| B6 approx | plug-in, fixed cutoff | 1500 | — | 0.0185 [0.0135, 0.0254] | 0.0040 [0.0020, 0.0079] | 0.0020 [0.0008, 0.0051] | 0.0010 [0.0003, 0.0036] | — |
| C1 | planted harm HR 1.5 | 1000 | 0.6985 [0.6780, 0.7182] | 0.5025 [0.4806, 0.5244] | 0.3730 [0.3521, 0.3944] | 0.2705 [0.2515, 0.2904] | 0.1880 [0.1715, 0.2057] | 0.0495 [0.0408, 0.0599] |
| C1 approx | plug-in, fixed cutoff | 1000 | — | 0.5025 [0.4806, 0.5244] | 0.3650 [0.3442, 0.3863] | 0.2685 [0.2495, 0.2884] | 0.1925 [0.1758, 0.2104] | — |
| C2 | planted harm HR 1.5 | 1500 | 0.7905 [0.7721, 0.8078] | 0.7470 [0.7275, 0.7656] | 0.6200 [0.5985, 0.6410] | 0.4910 [0.4691, 0.5129] | 0.3625 [0.3417, 0.3838] | 0.0935 [0.0815, 0.1071] |
| C2 approx | plug-in, fixed cutoff | 1500 | — | 0.7505 [0.7311, 0.7690] | 0.6190 [0.5975, 0.6400] | 0.4785 [0.4567, 0.5004] | 0.3515 [0.3309, 0.3727] | — |
| C3 | planted harm HR 2.0 | 1000 | 0.9635 [0.9544, 0.9709] | 0.8960 [0.8819, 0.9086] | 0.8340 [0.8171, 0.8497] | 0.7530 [0.7336, 0.7714] | 0.6760 [0.6552, 0.6962] | 0.3550 [0.3343, 0.3762] |
| C3 approx | plug-in, fixed cutoff | 1000 | — | 0.8945 [0.8803, 0.9072] | 0.8325 [0.8155, 0.8482] | 0.7570 [0.7377, 0.7753] | 0.6740 [0.6531, 0.6942] | — |
| C4 | planted harm HR 2.0 | 1500 | 0.9880 [0.9822, 0.9919] | 0.9830 [0.9763, 0.9878] | 0.9660 [0.9571, 0.9731] | 0.9275 [0.9153, 0.9381] | 0.8820 [0.8671, 0.8954] | 0.6410 [0.6197, 0.6617] |
| C4 approx | plug-in, fixed cutoff | 1500 | — | 0.9850 [0.9787, 0.9895] | 0.9650 [0.9560, 0.9722] | 0.9295 [0.9174, 0.9399] | 0.8780 [0.8629, 0.8916] | — |

- **Size:** every uniform-benefit rate at every c0 is below 0.05. The largest is 0.0315 (B4, c0 0.70, Wilson 0.0247–0.0401).
- **Mechanism:** within each B design, the rate falls as c0 rises and falls as n grows.
- **HR 0.721 vs 0.657:** B4–B6 run higher than B1–B3 at every c0. Their benefit is weaker, which puts it nearer each c0.
- **Power:** every c0 recovers most of the power that the committed c0 = c2 calibration gives up. At HR 1.5, n 1000, power is 0.19–0.50 against 0.05.
- **Exact vs approx:** the exact per-replicate rule agrees with the plug-in approximation to within 0.0125 in every entry. The largest gaps are C2 at c0 0.80 (0.4910 vs 0.4785) and C2 at 0.85 (0.3625 vs 0.3515). The approximation reads slightly low in the B cells: 0.0315 vs 0.0270 at B4, c0 0.70.

### α 0.10

| cell | DGM | n | p\* 0.90 as executed | c0 0.70 | c0 0.75 | c0 0.80 | c0 0.85 | c0 = c2 (committed) |
|---|---|---|---|---|---|---|---|---|
| B1 | uniform benefit HR 0.657 | 500 | 0.0950 [0.0829, 0.1086] | 0.0240 [0.0181, 0.0317] | 0.0095 [0.0061, 0.0148] | 0.0035 [0.0017, 0.0072] | 0.0015 [0.0005, 0.0044] | 0.0005 [0.0001, 0.0028] |
| B1 approx | plug-in, fixed cutoff | 500 | — | 0.0220 [0.0164, 0.0294] | 0.0085 [0.0053, 0.0136] | 0.0040 [0.0020, 0.0079] | 0.0015 [0.0005, 0.0044] | — |
| B2 | uniform benefit HR 0.657 | 1000 | 0.0325 [0.0256, 0.0412] | 0.0170 [0.0122, 0.0237] | 0.0075 [0.0046, 0.0123] | 0.0025 [0.0011, 0.0058] | 0.0005 [0.0001, 0.0028] | 0.0000 [0.0000, 0.0019] |
| B2 approx | plug-in, fixed cutoff | 1000 | — | 0.0175 [0.0126, 0.0242] | 0.0060 [0.0034, 0.0105] | 0.0030 [0.0014, 0.0065] | 0.0010 [0.0003, 0.0036] | — |
| B3 | uniform benefit HR 0.657 | 1500 | 0.0070 [0.0042, 0.0117] | 0.0120 [0.0081, 0.0178] | 0.0035 [0.0017, 0.0072] | 0.0025 [0.0011, 0.0058] | 0.0005 [0.0001, 0.0028] | 0.0000 [0.0000, 0.0019] |
| B3 approx | plug-in, fixed cutoff | 1500 | — | 0.0100 [0.0065, 0.0154] | 0.0035 [0.0017, 0.0072] | 0.0020 [0.0008, 0.0051] | 0.0005 [0.0001, 0.0028] | — |
| B4 | uniform benefit HR 0.721 | 500 | 0.2325 [0.2145, 0.2515] | 0.0705 [0.0601, 0.0826] | 0.0350 [0.0278, 0.0440] | 0.0165 [0.0118, 0.0231] | 0.0055 [0.0031, 0.0098] | 0.0005 [0.0001, 0.0028] |
| B4 approx | plug-in, fixed cutoff | 500 | — | 0.0570 [0.0477, 0.0680] | 0.0300 [0.0234, 0.0384] | 0.0120 [0.0081, 0.0178] | 0.0030 [0.0014, 0.0065] | — |
| B5 | uniform benefit HR 0.721 | 1000 | 0.1110 [0.0980, 0.1255] | 0.0645 [0.0545, 0.0761] | 0.0255 [0.0194, 0.0334] | 0.0115 [0.0077, 0.0172] | 0.0035 [0.0017, 0.0072] | 0.0000 [0.0000, 0.0019] |
| B5 approx | plug-in, fixed cutoff | 1000 | — | 0.0565 [0.0472, 0.0675] | 0.0240 [0.0181, 0.0317] | 0.0100 [0.0065, 0.0154] | 0.0030 [0.0014, 0.0065] | — |
| B6 | uniform benefit HR 0.721 | 1500 | 0.0325 [0.0256, 0.0412] | 0.0525 [0.0436, 0.0632] | 0.0155 [0.0109, 0.0219] | 0.0035 [0.0017, 0.0072] | 0.0020 [0.0008, 0.0051] | 0.0000 [0.0000, 0.0019] |
| B6 approx | plug-in, fixed cutoff | 1500 | — | 0.0470 [0.0386, 0.0572] | 0.0145 [0.0101, 0.0207] | 0.0040 [0.0020, 0.0079] | 0.0015 [0.0005, 0.0044] | — |
| C1 | planted harm HR 1.5 | 1000 | 0.6985 [0.6780, 0.7182] | 0.6165 [0.5950, 0.6376] | 0.5010 [0.4791, 0.5229] | 0.3695 [0.3486, 0.3909] | 0.2650 [0.2461, 0.2848] | 0.0825 [0.0712, 0.0954] |
| C1 approx | plug-in, fixed cutoff | 1000 | — | 0.6055 [0.5839, 0.6267] | 0.4980 [0.4761, 0.5199] | 0.3635 [0.3427, 0.3848] | 0.2655 [0.2466, 0.2853] | — |
| C2 | planted harm HR 1.5 | 1500 | 0.7905 [0.7721, 0.8078] | 0.8235 [0.8062, 0.8396] | 0.7225 [0.7025, 0.7417] | 0.5985 [0.5768, 0.6198] | 0.4585 [0.4368, 0.4804] | 0.1445 [0.1298, 0.1606] |
| C2 approx | plug-in, fixed cutoff | 1500 | — | 0.8225 [0.8051, 0.8386] | 0.7185 [0.6984, 0.7378] | 0.6005 [0.5789, 0.6218] | 0.4615 [0.4397, 0.4834] | — |
| C3 | planted harm HR 2.0 | 1000 | 0.9635 [0.9544, 0.9709] | 0.9425 [0.9314, 0.9519] | 0.8960 [0.8819, 0.9086] | 0.8330 [0.8160, 0.8487] | 0.7555 [0.7362, 0.7738] | 0.4475 [0.4258, 0.4694] |
| C3 approx | plug-in, fixed cutoff | 1000 | — | 0.9380 [0.9266, 0.9478] | 0.8915 [0.8771, 0.9044] | 0.8310 [0.8139, 0.8468] | 0.7560 [0.7367, 0.7743] | — |
| C4 | planted harm HR 2.0 | 1500 | 0.9880 [0.9822, 0.9919] | 0.9915 [0.9864, 0.9947] | 0.9820 [0.9752, 0.9870] | 0.9635 [0.9544, 0.9709] | 0.9210 [0.9084, 0.9320] | 0.7255 [0.7055, 0.7446] |
| C4 approx | plug-in, fixed cutoff | 1500 | — | 0.9915 [0.9864, 0.9947] | 0.9810 [0.9740, 0.9861] | 0.9600 [0.9505, 0.9677] | 0.9210 [0.9084, 0.9320] | — |

- **Size:** every uniform-benefit rate is below 0.10. The largest is 0.0705 (B4, c0 0.70, Wilson 0.0601–0.0826).
- **The calibrated screen is looser than p\* 0.90 at n 1500 and c0 0.70.** There the median κ̂₀.₁₀ is 1.506 (p\* 0.86786), below the executed cutoff of about 1.62. It declares more often than p\* 0.90 as executed in:
  - B3: 0.0120 vs 0.0070;
  - B6: 0.0525 vs 0.0325;
  - C2: 0.8235 vs 0.7905;
  - C4: 0.9915 vs 0.9880.
- **Exact vs approx:** the approximation understates the exact B rates by up to 0.0135 (B4, c0 0.70). The power entries agree to within 0.011.

## 4. The calibration's own quantities per c0 (§6.2)

- **κ̂₀.₀₅ and implied p\*:** medians with IQRs.
- **fw_1645, fw_1621:** means of `mean(Mstar_c0 > cutoff)`.
- **n_admitted:** the median over all replicates. The last column is the share of declaring replicates (cal05) with exactly one admitted subgroup.

| cell | c0 | κ̂₀.₀₅ median (IQR) | implied p\*<sup>†</sup> median (IQR) | mean fw_1645 | mean fw_1621 | n_admitted median | declaring reps | share n_admitted = 1 |
|---|---|---|---|---|---|---|---|---|
| B1 | 0.70 | 2.476 (2.416–2.546) | 0.98670 (0.98429–0.98911) | 0.3329 | 0.3465 | 0 | 16 | 0.5625 (9/16) |
| B1 | 0.75 | 2.676 (2.617–2.739) | 0.99255 (0.99112–0.99383) | 0.4547 | 0.4701 | 0 | 7 | 0.5714 (4/7) |
| B1 | 0.80 | 2.867 (2.811–2.927) | 0.99586 (0.99507–0.99658) | 0.5808 | 0.5961 | 0 | 2 | 0.5000 (1/2) |
| B1 | 0.85 | 3.058 (3.005–3.116) | 0.99777 (0.99734–0.99816) | 0.6974 | 0.7115 | 0 | 1 | 0/1 |
| B1 | c2 (committed) | 3.606 (3.556–3.658) | 0.99969 (0.99962–0.99975) | 0.9242 | 0.9302 | 0 | 1 | 0/1 |
| B2 | 0.70 | 2.125 (2.068–2.180) | 0.96640 (0.96132–0.97073) | 0.1745 | 0.1841 | 0 | 13 | 0.6923 (9/13) |
| B2 | 0.75 | 2.379 (2.326–2.432) | 0.98266 (0.98000–0.98497) | 0.2967 | 0.3099 | 0 | 4 | 4/4 |
| B2 | 0.80 | 2.628 (2.580–2.679) | 0.99140 (0.99013–0.99261) | 0.4514 | 0.4670 | 0 | 1 | 1/1 |
| B2 | 0.85 | 2.872 (2.827–2.920) | 0.99592 (0.99530–0.99650) | 0.6160 | 0.6318 | 0 | 1 | 1/1 |
| B2 | c2 (committed) | 3.604 (3.563–3.654) | 0.99969 (0.99963–0.99974) | 0.9405 | 0.9457 | 0 | 0 | — |
| B3 | 0.70 | 1.808 (1.754–1.865) | 0.92937 (0.92051–0.93775) | 0.0814 | 0.0866 | 0 | 9 | 0.6667 (6/9) |
| B3 | 0.75 | 2.112 (2.062–2.165) | 0.96535 (0.96081–0.96965) | 0.1730 | 0.1825 | 0 | 6 | 6/6 |
| B3 | 0.80 | 2.409 (2.362–2.458) | 0.98399 (0.98181–0.98604) | 0.3201 | 0.3343 | 0 | 3 | 3/3 |
| B3 | 0.85 | 2.702 (2.654–2.747) | 0.99310 (0.99205–0.99399) | 0.5098 | 0.5260 | 0 | 0 | — |
| B3 | c2 (committed) | 3.592 (3.545–3.635) | 0.99967 (0.99961–0.99972) | 0.9395 | 0.9448 | 0 | 0 | — |
| B4 | 0.70 | 2.426 (2.370–2.489) | 0.98475 (0.98220–0.98720) | 0.3113 | 0.3248 | 0 | 63 | 0.5556 (35/63) |
| B4 | 0.75 | 2.630 (2.576–2.692) | 0.99145 (0.99002–0.99291) | 0.4348 | 0.4500 | 0 | 26 | 0.5769 (15/26) |
| B4 | 0.80 | 2.830 (2.777–2.889) | 0.99534 (0.99452–0.99613) | 0.5652 | 0.5808 | 0 | 8 | 0.5000 (4/8) |
| B4 | 0.85 | 3.023 (2.975–3.079) | 0.99750 (0.99707–0.99792) | 0.6874 | 0.7017 | 0 | 2 | 1/2 |
| B4 | c2 (committed) | 3.589 (3.542–3.638) | 0.99967 (0.99960–0.99972) | 0.9247 | 0.9307 | 0 | 1 | 0/1 |
| B5 | 0.70 | 2.072 (2.019–2.127) | 0.96173 (0.95655–0.96657) | 0.1571 | 0.1659 | 0 | 57 | 0.5088 (29/57) |
| B5 | 0.75 | 2.335 (2.284–2.387) | 0.98045 (0.97760–0.98303) | 0.2766 | 0.2895 | 0 | 20 | 0.7000 (14/20) |
| B5 | 0.80 | 2.591 (2.543–2.641) | 0.99044 (0.98901–0.99173) | 0.4323 | 0.4481 | 0 | 8 | 0.7500 (6/8) |
| B5 | 0.85 | 2.840 (2.798–2.889) | 0.99548 (0.99485–0.99614) | 0.6023 | 0.6184 | 0 | 2 | 2/2 |
| B5 | c2 (committed) | 3.591 (3.549–3.639) | 0.99967 (0.99961–0.99973) | 0.9406 | 0.9457 | 0 | 0 | — |
| B6 | 0.70 | 1.754 (1.699–1.807) | 0.92051 (0.91075–0.92932) | 0.0703 | 0.0751 | 0 | 42 | 0.4762 (20/42) |
| B6 | 0.75 | 2.066 (2.014–2.119) | 0.96117 (0.95603–0.96593) | 0.1564 | 0.1653 | 0 | 13 | 0.6154 (8/13) |
| B6 | 0.80 | 2.368 (2.322–2.420) | 0.98212 (0.97977–0.98450) | 0.3015 | 0.3153 | 0 | 4 | 4/4 |
| B6 | 0.85 | 2.669 (2.623–2.719) | 0.99238 (0.99128–0.99346) | 0.4936 | 0.5102 | 0 | 2 | 2/2 |
| B6 | c2 (committed) | 3.580 (3.535–3.625) | 0.99966 (0.99959–0.99971) | 0.9397 | 0.9450 | 0 | 0 | — |
| C1 | 0.70 | 2.072 (2.018–2.132) | 0.96169 (0.95640–0.96697) | 0.1565 | 0.1653 | 1 | 1005 | 0.0945 (95/1005) |
| C1 | 0.75 | 2.332 (2.282–2.386) | 0.98030 (0.97748–0.98295) | 0.2750 | 0.2879 | 0 | 746 | 0.0925 (69/746) |
| C1 | 0.80 | 2.587 (2.539–2.638) | 0.99031 (0.98887–0.99167) | 0.4310 | 0.4467 | 0 | 541 | 0.1035 (56/541) |
| C1 | 0.85 | 2.837 (2.792–2.888) | 0.99545 (0.99476–0.99613) | 0.6014 | 0.6174 | 0 | 376 | 0.1090 (41/376) |
| C1 | c2 (committed) | 3.592 (3.548–3.641) | 0.99967 (0.99961–0.99973) | 0.9416 | 0.9468 | 0 | 99 | 0.0808 (8/99) |
| C2 | 0.70 | 1.755 (1.704–1.809) | 0.92081 (0.91155–0.92948) | 0.0704 | 0.0751 | 3 | 1494 | 0.0382 (57/1494) |
| C2 | 0.75 | 2.065 (2.016–2.117) | 0.96107 (0.95621–0.96573) | 0.1559 | 0.1647 | 2 | 1240 | 0.0508 (63/1240) |
| C2 | 0.80 | 2.366 (2.321–2.416) | 0.98203 (0.97972–0.98429) | 0.2992 | 0.3130 | 0 | 982 | 0.0479 (47/982) |
| C2 | 0.85 | 2.666 (2.621–2.712) | 0.99232 (0.99122–0.99332) | 0.4922 | 0.5087 | 0 | 725 | 0.0441 (32/725) |
| C2 | c2 (committed) | 3.580 (3.536–3.625) | 0.99966 (0.99959–0.99971) | 0.9412 | 0.9463 | 0 | 187 | 0.0535 (10/187) |
| C3 | 0.70 | 2.061 (2.008–2.119) | 0.96071 (0.95535–0.96587) | 0.1533 | 0.1620 | 7 | 1792 | 0.0195 (35/1792) |
| C3 | 0.75 | 2.324 (2.271–2.378) | 0.97985 (0.97684–0.98259) | 0.2710 | 0.2838 | 4 | 1668 | 0.0294 (49/1668) |
| C3 | 0.80 | 2.577 (2.532–2.630) | 0.99004 (0.98865–0.99146) | 0.4271 | 0.4430 | 3 | 1506 | 0.0325 (49/1506) |
| C3 | 0.85 | 2.830 (2.787–2.879) | 0.99535 (0.99467–0.99601) | 0.5989 | 0.6150 | 2 | 1352 | 0.0392 (53/1352) |
| C3 | c2 (committed) | 3.590 (3.545–3.635) | 0.99967 (0.99961–0.99972) | 0.9422 | 0.9473 | 0 | 710 | 0.0451 (32/710) |
| C4 | 0.70 | 1.744 (1.692–1.801) | 0.91877 (0.90940–0.92824) | 0.0686 | 0.0732 | 13 | 1966 | 0.0061 (12/1966) |
| C4 | 0.75 | 2.056 (2.006–2.108) | 0.96023 (0.95518–0.96500) | 0.1530 | 0.1618 | 9 | 1932 | 0.0057 (11/1932) |
| C4 | 0.80 | 2.360 (2.312–2.410) | 0.98173 (0.97921–0.98403) | 0.2957 | 0.3093 | 6 | 1855 | 0.0097 (18/1855) |
| C4 | 0.85 | 2.660 (2.613–2.706) | 0.99218 (0.99101–0.99320) | 0.4892 | 0.5057 | 4 | 1764 | 0.0085 (15/1764) |
| C4 | c2 (committed) | 3.575 (3.533–3.621) | 0.99965 (0.99959–0.99971) | 0.9417 | 0.9469 | 2 | 1282 | 0.0187 (24/1282) |

<sup>†</sup> See the note at the top of this report. A literal re-run at the implied p\* needs `pconsistency.digits >= 5`.

## 5. Readings (§6.3)

### 5.1 The trade, with c0 as rows

- **Worst B:** the worst uniform-benefit false-declaration rate over B1–B6, and the cell where it occurs.
- **C columns:** power against the planted harm region.
- **"approx" rows:** these come from `REPORT_declcal_c0_approx_2026-09-22` and sit beneath the exact rows.

**α 0.05**

| screen | worst B (cell) | HR 1.5 n 1000 | HR 1.5 n 1500 | HR 2.0 n 1000 | HR 2.0 n 1500 |
|---|---|---|---|---|---|
| calibrated, c0 0.70 | 0.0315 (B4) | 0.5025 | 0.7470 | 0.8960 | 0.9830 |
| — approx, c0 0.70 | 0.0270 (B4) | 0.5025 | 0.7505 | 0.8945 | 0.9850 |
| calibrated, c0 0.75 | 0.0130 (B4) | 0.3730 | 0.6200 | 0.8340 | 0.9660 |
| — approx, c0 0.75 | 0.0100 (B5) | 0.3650 | 0.6190 | 0.8325 | 0.9650 |
| calibrated, c0 0.80 | 0.0040 (B4) | 0.2705 | 0.4910 | 0.7530 | 0.9275 |
| — approx, c0 0.80 | 0.0035 (B5) | 0.2685 | 0.4785 | 0.7570 | 0.9295 |
| calibrated, c0 0.85 | 0.0010 (B4) | 0.1880 | 0.3625 | 0.6760 | 0.8820 |
| — approx, c0 0.85 | 0.0010 (B4) | 0.1925 | 0.3515 | 0.6740 | 0.8780 |
| calibrated, c0 = c2 = 1.00 (committed) | 0.0005 (B1) | 0.0495 | 0.0935 | 0.3550 | 0.6410 |
| tuned fixed p\* 0.9545 (k 2.0), post-reduction (committed practical) | 0.0725 (B4) | 0.5255 | 0.6475 | 0.9095 | 0.9705 |
| tuned fixed p\* 0.9545 (k 2.0), pre-reduction | 0.0955 (B4) | 0.5305 | 0.6500 | 0.9120 | 0.9710 |
| p\* 0.90 as executed | 0.2325 (B4) | 0.6985 | 0.7905 | 0.9635 | 0.9880 |

**α 0.10**

| screen | worst B (cell) | HR 1.5 n 1000 | HR 1.5 n 1500 | HR 2.0 n 1000 | HR 2.0 n 1500 |
|---|---|---|---|---|---|
| calibrated, c0 0.70 | 0.0705 (B4) | 0.6165 | 0.8235 | 0.9425 | 0.9915 |
| — approx, c0 0.70 | 0.0570 (B4) | 0.6055 | 0.8225 | 0.9380 | 0.9915 |
| calibrated, c0 0.75 | 0.0350 (B4) | 0.5010 | 0.7225 | 0.8960 | 0.9820 |
| — approx, c0 0.75 | 0.0300 (B4) | 0.4980 | 0.7185 | 0.8915 | 0.9810 |
| calibrated, c0 0.80 | 0.0165 (B4) | 0.3695 | 0.5985 | 0.8330 | 0.9635 |
| — approx, c0 0.80 | 0.0120 (B4) | 0.3635 | 0.6005 | 0.8310 | 0.9600 |
| calibrated, c0 0.85 | 0.0055 (B4) | 0.2650 | 0.4585 | 0.7555 | 0.9210 |
| — approx, c0 0.85 | 0.0030 (B4) | 0.2655 | 0.4615 | 0.7560 | 0.9210 |
| calibrated, c0 = c2 = 1.00 (committed) | 0.0005 (B1) | 0.0825 | 0.1445 | 0.4475 | 0.7255 |
| tuned fixed p\* 0.9545 (k 2.0), post-reduction (committed practical) | 0.0725 (B4) | 0.5255 | 0.6475 | 0.9095 | 0.9705 |
| tuned fixed p\* 0.9545 (k 2.0), pre-reduction | 0.0955 (B4) | 0.5305 | 0.6500 | 0.9120 | 0.9710 |
| p\* 0.90 as executed | 0.2325 (B4) | 0.6985 | 0.7905 | 0.9635 | 0.9880 |

- **Source of the fixed-k rows.** They are parsed from `logs/declcal_fixedk_practical.txt` (k = 2.0000). The post-reduction column is reproduced from the committed payloads in all 10 of 10 cells. The fixed-k rows do not depend on α.
- **Where the worst B rate falls.** It is B4 (HR 0.721, n 500) in every calibrated exact row. That is the weakest benefit at the smallest n.
- **c0 0.70 at α 0.05 against fixed k 2.0 (post-reduction):**
  - worst B rate: 0.0315 against 0.0725;
  - power at n 1000: about 0.02 lower (C1 0.5025 vs 0.5255; C3 0.8960 vs 0.9095);
  - power at n 1500: higher (C2 0.7470 vs 0.6475; C4 0.9830 vs 0.9705).
  - The reason is that κ̂(c0) falls with n while k is fixed at 2.0.
- **Stepping c0 from 0.70 to 0.85 at α 0.05.** The worst B rate falls from 0.0315 to 0.0010. HR 1.5 power at n 1000 falls from 0.50 to 0.19. The steepest part of the trade is HR 1.5 at n 1000.
- **Exact vs approx.** The exact rule reproduces the approximate table's ordering and magnitudes. The approximation runs slightly optimistic on size: at α 0.10, c0 0.70, the exact worst B rate is 0.0705 against an approximate 0.0570.

### 5.2 Is κ̂(c0) configuration-invariant?

- **Rows:** the median of the per-cell medians at each n, over the B and C cells at that n.
- **Spread:** the min–max of the per-cell medians, in parentheses.
- **Comparison:** the approximate table's median sits beside each exact value, each with its implied p\*<sup>†</sup>.

**κ̂₀.₀₅**

| n | c0 | cells | exact median (min–max) | exact implied p\* | approx median | approx implied p\* |
|---|---|---|---|---|---|---|
| 500 | 0.70 | B1,B4 | 2.451 (2.426–2.476) | 0.98576 | 2.448 | 0.98563 |
| 500 | 0.75 | B1,B4 | 2.653 (2.630–2.676) | 0.99202 | 2.657 | 0.99211 |
| 500 | 0.80 | B1,B4 | 2.848 (2.830–2.867) | 0.99561 | 2.848 | 0.99560 |
| 500 | 0.85 | B1,B4 | 3.040 (3.023–3.058) | 0.99764 | 3.046 | 0.99768 |
| 500 | c2 (committed) | B1,B4 | 3.598 (3.589–3.606) | 0.99968 | 3.613 | 0.99970 |
| 1000 | 0.70 | B2,B5,C1,C3 | 2.072 (2.061–2.125) | 0.96171 | 2.080 | 0.96249 |
| 1000 | 0.75 | B2,B5,C1,C3 | 2.333 (2.324–2.379) | 0.98038 | 2.343 | 0.98086 |
| 1000 | 0.80 | B2,B5,C1,C3 | 2.589 (2.577–2.628) | 0.99038 | 2.586 | 0.99030 |
| 1000 | 0.85 | B2,B5,C1,C3 | 2.838 (2.830–2.872) | 0.99547 | 2.835 | 0.99542 |
| 1000 | c2 (committed) | B2,B5,C1,C3 | 3.591 (3.590–3.604) | 0.99967 | 3.586 | 0.99966 |
| 1500 | 0.70 | B3,B6,C2,C4 | 1.755 (1.744–1.808) | 0.92066 | 1.758 | 0.92124 |
| 1500 | 0.75 | B3,B6,C2,C4 | 2.065 (2.056–2.112) | 0.96112 | 2.077 | 0.96221 |
| 1500 | 0.80 | B3,B6,C2,C4 | 2.367 (2.360–2.409) | 0.98208 | 2.391 | 0.98318 |
| 1500 | 0.85 | B3,B6,C2,C4 | 2.667 (2.660–2.702) | 0.99235 | 2.690 | 0.99285 |
| 1500 | c2 (committed) | B3,B6,C2,C4 | 3.580 (3.575–3.592) | 0.99966 | 3.588 | 0.99967 |

**κ̂₀.₁₀**

| n | c0 | cells | exact median (min–max) | exact implied p\* | approx median | approx implied p\* |
|---|---|---|---|---|---|---|
| 500 | 0.70 | B1,B4 | 2.192 (2.170–2.214) | 0.97163 | 2.200 | 0.97220 |
| 500 | 0.75 | B1,B4 | 2.394 (2.374–2.414) | 0.98333 | 2.402 | 0.98370 |
| 500 | 0.80 | B1,B4 | 2.590 (2.572–2.608) | 0.99040 | 2.600 | 0.99069 |
| 500 | 0.85 | B1,B4 | 2.782 (2.768–2.797) | 0.99461 | 2.793 | 0.99477 |
| 500 | c2 (committed) | B1,B4 | 3.338 (3.331–3.345) | 0.99916 | 3.341 | 0.99917 |
| 1000 | 0.70 | B2,B5,C1,C3 | 1.823 (1.811–1.871) | 0.93165 | 1.841 | 0.93434 |
| 1000 | 0.75 | B2,B5,C1,C3 | 2.086 (2.075–2.129) | 0.96298 | 2.095 | 0.96379 |
| 1000 | 0.80 | B2,B5,C1,C3 | 2.342 (2.331–2.378) | 0.98081 | 2.348 | 0.98113 |
| 1000 | 0.85 | B2,B5,C1,C3 | 2.592 (2.586–2.623) | 0.99045 | 2.595 | 0.99053 |
| 1000 | c2 (committed) | B2,B5,C1,C3 | 3.343 (3.339–3.354) | 0.99917 | 3.333 | 0.99914 |
| 1500 | 0.70 | B3,B6,C2,C4 | 1.506 (1.497–1.560) | 0.86786 | 1.514 | 0.86991 |
| 1500 | 0.75 | B3,B6,C2,C4 | 1.818 (1.811–1.865) | 0.93091 | 1.827 | 0.93225 |
| 1500 | 0.80 | B3,B6,C2,C4 | 2.123 (2.116–2.162) | 0.96623 | 2.129 | 0.96675 |
| 1500 | 0.85 | B3,B6,C2,C4 | 2.424 (2.416–2.455) | 0.98465 | 2.434 | 0.98505 |
| 1500 | c2 (committed) | B3,B6,C2,C4 | 3.332 (3.329–3.340) | 0.99914 | 3.337 | 0.99914 |

- **Invariance is approximate.** κ̂(c0) is invariant to about ±0.03 at fixed n and c0. The widest per-cell spread is 0.064 (n 1000, c0 0.70). The unshifted κ̂ was invariant to within 0.017.
- **Where the spread comes from.** It is systematic, not noise. The HR 0.657 cells (B2, B3) sit highest at n 1000 and 1500: B2 is 2.125 against 2.061–2.072 for the rest at c0 0.70. B5, C1 and C3 agree to within 0.011. The shift (c2 − c0)/σ_D(g) depends on σ_D, which the design moves slightly.
- **Unlike the unshifted κ̂, κ̂(c0) depends strongly on n.** At c0 0.70, κ̂₀.₀₅ falls from 2.45 to 1.76 between n 500 and n 1500. The equivalent fixed p\* therefore changes with n: from 0.986 to 0.921 at α 0.05.
- **Exact vs approx medians.** The exact and approximate κ̂ medians agree to within 0.024 (n 1500, c0 0.80).

### 5.3 Is the re-selection question live?

This is the share of declaring replicates (cal05) with `n_admitted_cal05_c0 == 1`, pooled.

| c0 | B: declaring reps | B: share = 1 | C: declaring reps | C: share = 1 | C: median n_admitted among declaring |
|---|---|---|---|---|---|
| 0.70 | 200 | 0.5400 | 6257 | 0.0318 | 7 |
| 0.75 | 76 | 0.6711 | 5586 | 0.0344 | 6 |
| 0.80 | 26 | 0.7308 | 4884 | 0.0348 | 5 |
| 0.85 | 8 | 0.7500 | 4217 | 0.0334 | 4 |

- **The re-selection question is live at every c0 in the grid.** Under planted harm, about 97% of declaring replicates admit more than one subgroup at the calibrated cutoff. The median is 4–7 admitted subgroups, so which subgroup the calibrated screen would select is not determined by the declaration.
- **The committed c0 = c2 calibration had the same property.** Its C-cell shares of 1 are 0.02–0.08.
- **Under uniform benefit,** about half of the few declarations admit a single subgroup.
- **Why this matters.** The executed search selects among subgroups that pass the rounded p\* 0.90 rule, not the calibrated rule. `sg_focus` on the calibrated set is out of scope (task §9). The column above prices that question.

### 5.4 Does fw_c0 now track a realized rate?

| cell | uniform marginal HR | c0 | mean fw_1645_c0 | realized pre-family exact rate (max_T_pre ≥ 1.6449) | diff (paired MC SE) | mean fw_1621_c0 | executed rate | diff (paired MC SE) |
|---|---|---|---|---|---|---|---|---|
| B4 | 0.721 | 0.70 | 0.3113 | 0.2500 | +0.0613 (0.0098) | 0.3248 | 0.2325 | +0.0923 (0.0096) |
| B5 | 0.721 | 0.70 | 0.1571 | 0.1050 | +0.0521 (0.0069) | 0.1659 | 0.1110 | +0.0549 (0.0071) |
| B6 | 0.721 | 0.70 | 0.0703 | 0.0300 | +0.0403 (0.0038) | 0.0751 | 0.0325 | +0.0426 (0.0040) |
| B1 | 0.657 | 0.70 | 0.3329 | 0.1220 | +0.2109 (0.0075) | 0.3465 | 0.0950 | +0.2515 (0.0067) |
| B2 | 0.657 | 0.70 | 0.1745 | 0.0325 | +0.1420 (0.0040) | 0.1841 | 0.0325 | +0.1516 (0.0040) |
| B3 | 0.657 | 0.70 | 0.0814 | 0.0070 | +0.0744 (0.0019) | 0.0866 | 0.0070 | +0.0796 (0.0019) |

All 24 (cell, c0) rows are in the findings log.

- **Nearer, but not tracking.** In B4–B6, whose marginal Cox HR of 0.721 is nearest c0 0.70, `fw_1645_0.70` overstates the realized pre-family rate by 0.04–0.06. That is 6 to 11 paired SEs.
- **Compared with the unshifted diagnostic,** this is much closer: the unshifted fw is 0.92–0.94 against realized rates of 0.007–0.25. It still overstates in every cell.
- **The ordering is right.** The diagnostic ranks the B4–B6 cells correctly (0.31 > 0.16 > 0.07 against 0.25 > 0.105 > 0.03).
- **The overstatement grows as c0 moves away from the true effect.** At c0 0.85 it is 0.44–0.50.
- **Vocabulary.** The null designs are fixed by marginal Cox HR. Their uniform patient-level HR is 0.6566 at the 0.721 point and 0.5829 at the 0.657 point (`current_status.md` §1). c0 0.70 therefore lies above the patient-level HR of the B4–B6 designs. The diagnostic's residual overstatement is consistent with the true null sitting further from c2 than c0 does. That reading is offered as an explanation, not tested here.

### 5.5 The three free checks (declcal task §7.3, per c0)

These are findings, not gates. Block A is not re-run, so check (a) uses the B cells.

#### (a) Eq. 8 at the shifted field as an estimator

- **Pairings:** mean `fw_1621_c0` against the as-executed rate, and mean `fw_1645_c0` against the exact-z pre-family rate (`max_T_pre ≥ 1.6449`).
- **Standard error:** paired Monte Carlo SE = sd / √2000.

| cell | uniform HR | c0 | mean fw_1645_c0 | realized pre-family exact rate (max_T_pre >= 1.6449) | diff | paired MC SE | mean fw_1621_c0 | executed rate (post, rounded) | diff | paired MC SE |
|---|---|---|---|---|---|---|---|---|---|---|
| B1 | 0.657 | 0.70 | 0.3329 | 0.1220 | +0.2109 | 0.0075 | 0.3465 | 0.0950 | +0.2515 | 0.0067 |
| B1 | 0.657 | 0.75 | 0.4547 | 0.1220 | +0.3327 | 0.0074 | 0.4701 | 0.0950 | +0.3751 | 0.0067 |
| B1 | 0.657 | 0.80 | 0.5808 | 0.1220 | +0.4588 | 0.0074 | 0.5961 | 0.0950 | +0.5011 | 0.0067 |
| B1 | 0.657 | 0.85 | 0.6974 | 0.1220 | +0.5754 | 0.0074 | 0.7115 | 0.0950 | +0.6165 | 0.0066 |
| B2 | 0.657 | 0.70 | 0.1745 | 0.0325 | +0.1420 | 0.0040 | 0.1841 | 0.0325 | +0.1516 | 0.0040 |
| B2 | 0.657 | 0.75 | 0.2967 | 0.0325 | +0.2642 | 0.0040 | 0.3099 | 0.0325 | +0.2774 | 0.0040 |
| B2 | 0.657 | 0.80 | 0.4514 | 0.0325 | +0.4189 | 0.0040 | 0.4670 | 0.0325 | +0.4345 | 0.0040 |
| B2 | 0.657 | 0.85 | 0.6160 | 0.0325 | +0.5835 | 0.0040 | 0.6318 | 0.0325 | +0.5993 | 0.0040 |
| B3 | 0.657 | 0.70 | 0.0814 | 0.0070 | +0.0744 | 0.0019 | 0.0866 | 0.0070 | +0.0796 | 0.0019 |
| B3 | 0.657 | 0.75 | 0.1730 | 0.0070 | +0.1660 | 0.0019 | 0.1825 | 0.0070 | +0.1755 | 0.0019 |
| B3 | 0.657 | 0.80 | 0.3201 | 0.0070 | +0.3131 | 0.0020 | 0.3343 | 0.0070 | +0.3273 | 0.0020 |
| B3 | 0.657 | 0.85 | 0.5098 | 0.0070 | +0.5028 | 0.0020 | 0.5260 | 0.0070 | +0.5190 | 0.0020 |
| B4 | 0.721 | 0.70 | 0.3113 | 0.2500 | +0.0613 | 0.0098 | 0.3248 | 0.2325 | +0.0923 | 0.0096 |
| B4 | 0.721 | 0.75 | 0.4348 | 0.2500 | +0.1848 | 0.0098 | 0.4500 | 0.2325 | +0.2175 | 0.0096 |
| B4 | 0.721 | 0.80 | 0.5652 | 0.2500 | +0.3152 | 0.0098 | 0.5808 | 0.2325 | +0.3483 | 0.0095 |
| B4 | 0.721 | 0.85 | 0.6874 | 0.2500 | +0.4374 | 0.0097 | 0.7017 | 0.2325 | +0.4692 | 0.0095 |
| B5 | 0.721 | 0.70 | 0.1571 | 0.1050 | +0.0521 | 0.0069 | 0.1659 | 0.1110 | +0.0549 | 0.0071 |
| B5 | 0.721 | 0.75 | 0.2766 | 0.1050 | +0.1716 | 0.0069 | 0.2895 | 0.1110 | +0.1785 | 0.0071 |
| B5 | 0.721 | 0.80 | 0.4323 | 0.1050 | +0.3273 | 0.0069 | 0.4481 | 0.1110 | +0.3371 | 0.0071 |
| B5 | 0.721 | 0.85 | 0.6023 | 0.1050 | +0.4973 | 0.0069 | 0.6184 | 0.1110 | +0.5074 | 0.0070 |
| B6 | 0.721 | 0.70 | 0.0703 | 0.0300 | +0.0403 | 0.0038 | 0.0751 | 0.0325 | +0.0426 | 0.0040 |
| B6 | 0.721 | 0.75 | 0.1564 | 0.0300 | +0.1264 | 0.0039 | 0.1653 | 0.0325 | +0.1328 | 0.0040 |
| B6 | 0.721 | 0.80 | 0.3015 | 0.0300 | +0.2715 | 0.0039 | 0.3153 | 0.0325 | +0.2828 | 0.0040 |
| B6 | 0.721 | 0.85 | 0.4936 | 0.0300 | +0.4636 | 0.0039 | 0.5102 | 0.0325 | +0.4777 | 0.0040 |

- **Material gap, recorded as a finding.** `fw_c0` overstates the realized rate in all 24 rows. The smallest gaps are at c0 0.70 in B4–B6 (+0.04 to +0.09), which are the cells whose marginal HR is nearest c0.
- **Why.** Eq. 8 at the shifted field is the family-wise size *if the null sat at c0 for every candidate*. Under uniform benefit the true patient-level effect is further from c2 than c0, so the realized rate is lower (§5.4).

#### (b) How strict the calibration is: `pstar_implied_05_c0`<sup>†</sup>

The last column is the share of replicates whose κ̂₀.₀₅(c0) is below the executed cutoff of 1.621. In those replicates the calibrated rule is looser than p\* 0.90 as executed.

| cell | c0 | min | 5% | 25% | 50% | 75% | 95% | max | fraction > 0.90 | fraction < p*(executed cutoff 1.621) |
|---|---|---|---|---|---|---|---|---|---|---|
| B1 | 0.70 | 0.97096 | 0.98026 | 0.98429 | 0.98670 | 0.98911 | 0.99204 | 0.99663 | 1.0000 | 0.0000 |
| B1 | 0.75 | 0.98230 | 0.98893 | 0.99112 | 0.99255 | 0.99383 | 0.99541 | 0.99794 | 1.0000 | 0.0000 |
| B1 | 0.80 | 0.99002 | 0.99378 | 0.99507 | 0.99586 | 0.99658 | 0.99747 | 0.99884 | 1.0000 | 0.0000 |
| B1 | 0.85 | 0.99468 | 0.99665 | 0.99734 | 0.99777 | 0.99816 | 0.99860 | 0.99925 | 1.0000 | 0.0000 |
| B2 | 0.70 | 0.93889 | 0.95339 | 0.96132 | 0.96640 | 0.97073 | 0.97663 | 0.98806 | 1.0000 | 0.0000 |
| B2 | 0.75 | 0.96748 | 0.97573 | 0.98000 | 0.98266 | 0.98497 | 0.98812 | 0.99324 | 1.0000 | 0.0000 |
| B2 | 0.80 | 0.98320 | 0.98791 | 0.99013 | 0.99140 | 0.99261 | 0.99417 | 0.99611 | 1.0000 | 0.0000 |
| B2 | 0.85 | 0.99176 | 0.99422 | 0.99530 | 0.99592 | 0.99650 | 0.99725 | 0.99837 | 1.0000 | 0.0000 |
| B3 | 0.70 | 0.87164 | 0.90822 | 0.92051 | 0.92937 | 0.93775 | 0.94930 | 0.96554 | 0.9840 | 0.0055 |
| B3 | 0.75 | 0.93283 | 0.95391 | 0.96081 | 0.96535 | 0.96965 | 0.97575 | 0.98208 | 1.0000 | 0.0000 |
| B3 | 0.80 | 0.96973 | 0.97824 | 0.98181 | 0.98399 | 0.98604 | 0.98896 | 0.99167 | 1.0000 | 0.0000 |
| B3 | 0.85 | 0.98551 | 0.99045 | 0.99205 | 0.99310 | 0.99399 | 0.99524 | 0.99638 | 1.0000 | 0.0000 |
| B4 | 0.70 | 0.96941 | 0.97773 | 0.98220 | 0.98475 | 0.98720 | 0.99038 | 0.99429 | 1.0000 | 0.0000 |
| B4 | 0.75 | 0.98315 | 0.98756 | 0.99002 | 0.99145 | 0.99291 | 0.99456 | 0.99691 | 1.0000 | 0.0000 |
| B4 | 0.80 | 0.99070 | 0.99320 | 0.99452 | 0.99534 | 0.99613 | 0.99703 | 0.99846 | 1.0000 | 0.0000 |
| B4 | 0.85 | 0.99465 | 0.99632 | 0.99707 | 0.99750 | 0.99792 | 0.99841 | 0.99910 | 1.0000 | 0.0000 |
| B5 | 0.70 | 0.93067 | 0.94855 | 0.95655 | 0.96173 | 0.96657 | 0.97260 | 0.98241 | 1.0000 | 0.0000 |
| B5 | 0.75 | 0.96431 | 0.97318 | 0.97760 | 0.98045 | 0.98303 | 0.98634 | 0.99042 | 1.0000 | 0.0000 |
| B5 | 0.80 | 0.98155 | 0.98671 | 0.98901 | 0.99044 | 0.99173 | 0.99343 | 0.99558 | 1.0000 | 0.0000 |
| B5 | 0.85 | 0.99140 | 0.99363 | 0.99485 | 0.99548 | 0.99614 | 0.99697 | 0.99807 | 1.0000 | 0.0000 |
| B6 | 0.70 | 0.87348 | 0.89657 | 0.91075 | 0.92051 | 0.92932 | 0.94191 | 0.95775 | 0.9165 | 0.0410 |
| B6 | 0.75 | 0.92945 | 0.94832 | 0.95603 | 0.96117 | 0.96593 | 0.97204 | 0.98062 | 1.0000 | 0.0000 |
| B6 | 0.80 | 0.96559 | 0.97573 | 0.97977 | 0.98212 | 0.98450 | 0.98762 | 0.99231 | 1.0000 | 0.0000 |
| B6 | 0.85 | 0.98499 | 0.98941 | 0.99128 | 0.99238 | 0.99346 | 0.99474 | 0.99664 | 1.0000 | 0.0000 |
| C1 | 0.70 | 0.93558 | 0.94801 | 0.95640 | 0.96169 | 0.96697 | 0.97369 | 0.98345 | 1.0000 | 0.0000 |
| C1 | 0.75 | 0.96452 | 0.97293 | 0.97748 | 0.98030 | 0.98295 | 0.98655 | 0.99057 | 1.0000 | 0.0000 |
| C1 | 0.80 | 0.98206 | 0.98666 | 0.98887 | 0.99031 | 0.99167 | 0.99348 | 0.99558 | 1.0000 | 0.0000 |
| C1 | 0.85 | 0.99142 | 0.99363 | 0.99476 | 0.99545 | 0.99613 | 0.99697 | 0.99799 | 1.0000 | 0.0000 |
| C2 | 0.70 | 0.87316 | 0.89637 | 0.91155 | 0.92081 | 0.92948 | 0.94293 | 0.95719 | 0.9220 | 0.0425 |
| C2 | 0.75 | 0.93435 | 0.94792 | 0.95621 | 0.96107 | 0.96573 | 0.97228 | 0.98086 | 1.0000 | 0.0000 |
| C2 | 0.80 | 0.96941 | 0.97570 | 0.97972 | 0.98203 | 0.98429 | 0.98742 | 0.99181 | 1.0000 | 0.0000 |
| C2 | 0.85 | 0.98464 | 0.98941 | 0.99122 | 0.99232 | 0.99332 | 0.99470 | 0.99648 | 1.0000 | 0.0000 |
| C3 | 0.70 | 0.93408 | 0.94699 | 0.95535 | 0.96071 | 0.96587 | 0.97292 | 0.98276 | 1.0000 | 0.0000 |
| C3 | 0.75 | 0.96635 | 0.97260 | 0.97684 | 0.97985 | 0.98259 | 0.98620 | 0.99072 | 1.0000 | 0.0000 |
| C3 | 0.80 | 0.98187 | 0.98639 | 0.98865 | 0.99004 | 0.99146 | 0.99321 | 0.99496 | 1.0000 | 0.0000 |
| C3 | 0.85 | 0.99032 | 0.99358 | 0.99467 | 0.99535 | 0.99601 | 0.99689 | 0.99795 | 1.0000 | 0.0000 |
| C4 | 0.70 | 0.86335 | 0.89482 | 0.90940 | 0.91877 | 0.92824 | 0.94147 | 0.95744 | 0.9100 | 0.0510 |
| C4 | 0.75 | 0.92832 | 0.94716 | 0.95518 | 0.96023 | 0.96500 | 0.97154 | 0.98086 | 1.0000 | 0.0000 |
| C4 | 0.80 | 0.96534 | 0.97536 | 0.97921 | 0.98173 | 0.98403 | 0.98716 | 0.99181 | 1.0000 | 0.0000 |
| C4 | 0.85 | 0.98587 | 0.98923 | 0.99101 | 0.99218 | 0.99320 | 0.99458 | 0.99675 | 1.0000 | 0.0000 |

- **Every calibrated replicate is stricter than p\* 0.90 except at c0 0.70 and n 1500.** There, the implied p\* is at most 0.90 in 1.6% (B3) to 9.0% (C4) of replicates. κ̂₀.₀₅(c0) is below 1.621 in 0.6% (B3), 4.1% (B6), 4.3% (C2) and 5.1% (C4) of them.
- **The median implied p\* ranges from 0.919 (C4, c0 0.70) to 0.998 (B1, c0 0.85).** Compare about 0.9997 for c0 = c2 in every cell.

#### (c) Can the calibrated rule declare where the conventional one did not?

The table shows only the (cell, c0) rows with a non-zero count; the other 23 rows are zero in every column (full table in the log).
- **Rule definitions:**
  - conv = `declared_conv`, the rounded post-reduction rule;
  - exact = `declared_conv_exact` (`max_T_post ≥ 1.6449`).
- **Mechanisms:**
  - *via family*: `max_T_post < κ̂`. The admitting subgroup is in the pre-reduction family, which the calibrated rule reads, but not in the post-reduction family that the executed screen evaluated.
  - *via cutoff*: `max_T_post ≥ κ̂`, yet the rounded rule did not admit it. This means κ̂ is below 1.621.

| cell | c0 | cal05 = 1 & conv = 0 | cal10 = 1 & conv = 0 | cal05 = 1 & exact = 0 | cal10 = 1 & exact = 0 | min kappa_hat_05_c0 | min kappa_hat_10_c0 | cal05 & !conv: via family / via cutoff | cal10 & !conv: via family / via cutoff |
|---|---|---|---|---|---|---|---|---|---|
| B1 | 0.70 | 1 | 10 | 1 | 10 | 2.183 | 1.942 | 1 / 0 | 10 / 0 |
| B1 | 0.75 | 0 | 2 | 0 | 2 | 2.372 | 2.150 | 0 / 0 | 2 / 0 |
| B2 | 0.70 | 0 | 1 | 0 | 1 | 1.873 | 1.665 | 0 / 0 | 1 / 0 |
| B3 | 0.70 | 0 | 11 | 0 | 12 | 1.521 | 1.339 | 0 / 0 | 2 / 9 |
| B4 | 0.70 | 3 | 13 | 3 | 13 | 2.162 | 1.917 | 3 / 0 | 13 / 0 |
| B4 | 0.75 | 0 | 4 | 0 | 4 | 2.390 | 2.158 | 0 / 0 | 4 / 0 |
| B5 | 0.70 | 2 | 2 | 2 | 2 | 1.816 | 1.624 | 2 / 0 | 2 / 0 |
| B5 | 0.75 | 1 | 2 | 1 | 2 | 2.100 | 1.896 | 1 / 0 | 2 / 0 |
| B6 | 0.70 | 1 | 40 | 1 | 45 | 1.528 | 1.293 | 0 / 1 | 3 / 37 |
| C1 | 0.70 | 3 | 4 | 4 | 5 | 1.849 | 1.619 | 3 / 0 | 4 / 0 |
| C1 | 0.75 | 2 | 3 | 2 | 4 | 2.103 | 1.866 | 2 / 0 | 3 / 0 |
| C1 | 0.80 | 0 | 2 | 0 | 2 | 2.367 | 2.141 | 0 / 0 | 2 / 0 |
| C2 | 0.70 | 4 | 66 | 5 | 76 | 1.527 | 1.292 | 3 / 1 | 2 / 64 |
| C2 | 0.75 | 0 | 2 | 0 | 2 | 1.841 | 1.618 | 0 / 0 | 2 / 0 |
| C3 | 0.70 | 1 | 2 | 2 | 3 | 1.839 | 1.624 | 1 / 0 | 2 / 0 |
| C3 | 0.75 | 0 | 1 | 0 | 2 | 2.124 | 1.901 | 0 / 0 | 1 / 0 |
| C4 | 0.70 | 0 | 7 | 0 | 7 | 1.488 | 1.289 | 0 / 0 | 0 / 7 |

- **Yes. Unlike the committed c0 = c2 calibration, which gave zero in every cell, the c0 screens do declare where p\* 0.90 did not.** Pooled over all 40 (cell, c0) rows: cal05 & !conv 18 (via family 16, via cutoff 2); cal10 & !conv 172 (via family 55, via cutoff 117).
- **At α 0.05 the mechanism is almost entirely the family, not the cutoff.** The calibrated rule reads `max_T_pre`, which includes candidates that the post-reduction floors (d0/d1) remove.
- **At α 0.10, c0 0.70 and n 1500, the cutoff dominates.** κ̂₀.₁₀ falls to about 1.29–1.34 there (B3 9, B6 37, C2 64, C4 7).
- **These counts are small against the declaring totals, except C2 at α 0.10, c0 0.70 (66 of 1,647).** They are nonetheless real: the calibrated screen at c0 is not a sub-event of the executed screen.

## 6. OPEN ITEMS

- **The first Stage 2 run stopped without a record.** It stopped mid-B3, at about 900 s elapsed, with no halt record, no gate line and no watchdog message. It was far below the 6,300 s cap. The cause is not recorded in the tree. The resumed run re-ran B3 from replicate 1, so no replicate from the stopped run is in any committed payload other than B1 and B2.
- **Resume cells files.** `scripts_dinamr/resume/declcalc0_{inull,power}.cells` are new files outside the task's §7 list. They exist only so the unchanged driver would run the remaining cells under the same tags. They are committed with this report.
- **`declcalc0_findings.R` has changes beyond a straight transplant of `declcal_findings.R`.** It reads the approximate-table payload `results/declcal_c0approx_res.rds` for the approx rows and adds the implied-p\* footnote. It was not in the task's §7 transplant list and is committed with this report.
- **The implied p\* are not runnable at the default rounding.** Any operational use of a fixed p\* taken from these κ̂ values must set `pconsistency.digits >= 5`. At the default of 2, every implied p\* above 0.995 collapses to "rounded proportion = 1.00".
- **Mixed B_cal across the comparison.** The exact rows use B_cal 500, and the approximate rows B 2000. The approximate report measured the B 2000 vs B 500 κ̂ difference at SD 0.057 per replicate. The median-level agreement in §5.2, within 0.024, is consistent with that.
- **Unrelated untracked files are left untouched:** `scripts_dinamr/logs/nullmr_findings.err` and the `actg175/` files. None is in any commit of this task.
