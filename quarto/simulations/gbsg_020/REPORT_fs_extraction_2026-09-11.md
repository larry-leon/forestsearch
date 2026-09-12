# REPORT — FS extraction for the manuscript

- **Date:** 2026-09-12 (task dated 2026-09-11)
- **Task:** `dev/tasks/TASK_grfmr_campaign_2026-09-11.md`, **Part C**
- **Producer:** `quarto/simulations/gbsg_020/scripts_dinamr/fs_extraction.R`
- **Inputs:** committed FS comparator bundles only — `tier2` / `p12ext` at 12.4%, `cert20` /
  `e1stud` at 31%. **No re-run, no new simulation, no recorder change.** Every number below comes
  from a column that is already in those bundles.

## Coverage and the criterion that is not shared

All **18** committed FS cells are covered: both prevalences × HR 1.50 / 1.75 / 1.00 × n 500 / 1000
/ 1500. The two prevalences were **not run under the same criterion**, and the comparator's
`sg_focus` and ε are carried beside every row for that reason:

| prevalence | campaigns | `sg_focus` | ε |
|---|---|---|---|
| 12.4% | `p12ext` (HR 1.50, and HR 1.00 at n 1000/1500); `tier2` (HR 1.75, and HR 1.00 at n 500) | `maxeffCons` | 0.10 |
| 31% | `e1stud` (n 500 at HR 1.50/1.75); `cert20` (all others) | `effMaxSG` | 0.20 |

Every comparison that crosses the two prevalences therefore carries a criterion confound in
addition to whatever it is meant to show.

## How the 2×2 was recovered, and what a Wilson interval is put on

`sens` / `spec` / `ppv` / `npv` are recorded **per replicate**
(`sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:811–822`), as `tp/(tp+fn)`, `tn/(tn+fp)`,
`tp/(tp+fp)`, `tn/(tn+fn)` over that replicate's own subjects. The 2×2 counts themselves are not
recorder columns, but they are **exactly recoverable**: `n_sel = tp + fp`, `n_true = tp + fn`, and
`tp = sens × n_true = ppv × n_sel`, so `tp`, `fp`, `fn` and `tn = n_sample − tp − fp − fn` all
follow. The two expressions for `tp` agree to **5.68e-14** over every row of every cell, so the
reconstruction is exact and not an approximation.

A **Wilson interval is put on the pooled subject-level counts** summed over replicates, which is
the quantity a Wilson interval is defined for. The **replicate-mean** rate — the form the DINA
tables use — is reported beside it without an interval, because a Wilson interval on a mean of
proportions would not be one. All rates and all location figures are computed on the **detected**
replicates, the set on which the quantities exist; `n_eval` is stated per cell.

---

## 1. Classification against the planted region

### 1a. Pooled rates, Wilson intervals

| cell | campaign | focus | ε | n_eval | detection | sensitivity | specificity |
|---|---|---|---|---|---|---|---|
| 12.4% HR 1.50 n 500 | p12ext | maxeffCons | 0.10 | 1822 | 0.9110 | 0.5770 [0.5741, 0.5799] | 0.9145 [0.9139, 0.9151] |
| 12.4% HR 1.50 n 1000 | p12ext | maxeffCons | 0.10 | 1948 | 0.9740 | 0.6969 [0.6950, 0.6987] | 0.9585 [0.9582, 0.9588] |
| 12.4% HR 1.50 n 1500 | p12ext | maxeffCons | 0.10 | 1976 | 0.9880 | 0.8194 [0.8182, 0.8206] | 0.9751 [0.9749, 0.9753] |
| 12.4% HR 1.75 n 500 | tier2 | maxeffCons | 0.10 | 1900 | 0.9500 | 0.6650 [0.6623, 0.6677] | 0.9297 [0.9292, 0.9303] |
| 12.4% HR 1.75 n 1000 | tier2 | maxeffCons | 0.10 | 1990 | 0.9950 | 0.7866 [0.7850, 0.7882] | 0.9714 [0.9712, 0.9717] |
| 12.4% HR 1.75 n 1500 | tier2 | maxeffCons | 0.10 | 1998 | 0.9990 | 0.8826 [0.8816, 0.8837] | 0.9831 [0.9829, 0.9832] |
| 12.4% HR 1.00 n 500 | tier2 | maxeffCons | 0.10 | 1361 | 0.6805 | 0.3526 [0.3494, 0.3558] | 0.8698 [0.8689, 0.8706] |
| 12.4% HR 1.00 n 1000 | p12ext | maxeffCons | 0.10 | 1319 | 0.6595 | 0.4127 [0.4103, 0.4151] | 0.9118 [0.9113, 0.9123] |
| 12.4% HR 1.00 n 1500 | p12ext | maxeffCons | 0.10 | 1248 | 0.6240 | 0.5284 [0.5264, 0.5304] | 0.9306 [0.9302, 0.9310] |
| 31% HR 1.50 n 500 | e1stud | effMaxSG | 0.20 | 1999 | 0.9995 | 0.5893 [0.5876, 0.5911] | 0.8934 [0.8927, 0.8941] |
| 31% HR 1.50 n 1000 | cert20 | effMaxSG | 0.20 | 2000 | 1.0000 | 0.7073 [0.7062, 0.7084] | 0.9276 [0.9272, 0.9281] |
| 31% HR 1.50 n 1500 | cert20 | effMaxSG | 0.20 | 2000 | 1.0000 | 0.8389 [0.8382, 0.8397] | 0.9256 [0.9252, 0.9259] |
| 31% HR 1.75 n 500 | e1stud | effMaxSG | 0.20 | 1999 | 0.9995 | 0.6590 [0.6573, 0.6606] | 0.9189 [0.9183, 0.9196] |
| 31% HR 1.75 n 1000 | cert20 | effMaxSG | 0.20 | 2000 | 1.0000 | 0.7683 [0.7673, 0.7694] | 0.9475 [0.9471, 0.9479] |
| 31% HR 1.75 n 1500 | cert20 | effMaxSG | 0.20 | 2000 | 1.0000 | 0.8738 [0.8731, 0.8745] | 0.9466 [0.9463, 0.9469] |
| 31% HR 1.00 n 500 | cert20 | effMaxSG | 0.20 | 1841 | 0.9205 | 0.3871 [0.3853, 0.3889] | 0.8278 [0.8269, 0.8288] |
| 31% HR 1.00 n 1000 | cert20 | effMaxSG | 0.20 | 1909 | 0.9545 | 0.4600 [0.4587, 0.4613] | 0.8527 [0.8521, 0.8533] |
| 31% HR 1.00 n 1500 | cert20 | effMaxSG | 0.20 | 1918 | 0.9590 | 0.5814 [0.5804, 0.5824] | 0.8430 [0.8425, 0.8435] |

| cell | PPV | NPV |
|---|---|---|
| 12.4% HR 1.50 n 500 | 0.4890 [0.4863, 0.4916] | 0.9385 [0.9379, 0.9390] |
| 12.4% HR 1.50 n 1000 | 0.7038 [0.7020, 0.7056] | 0.9572 [0.9568, 0.9575] |
| 12.4% HR 1.50 n 1500 | 0.8234 [0.8222, 0.8246] | 0.9744 [0.9742, 0.9746] |
| 12.4% HR 1.75 n 500 | 0.5727 [0.5701, 0.5754] | 0.9514 [0.9510, 0.9519] |
| 12.4% HR 1.75 n 1000 | 0.7957 [0.7941, 0.7973] | 0.9699 [0.9696, 0.9701] |
| 12.4% HR 1.75 n 1500 | 0.8807 [0.8796, 0.8817] | 0.9834 [0.9832, 0.9835] |
| 12.4% HR 1.00 n 500 | 0.2770 [0.2744, 0.2797] | 0.9047 [0.9039, 0.9055] |
| 12.4% HR 1.00 n 1000 | 0.3989 [0.3966, 0.4012] | 0.9163 [0.9158, 0.9168] |
| 12.4% HR 1.00 n 1500 | 0.5193 [0.5173, 0.5214] | 0.9329 [0.9325, 0.9333] |
| 31% HR 1.50 n 500 | 0.7090 [0.7073, 0.7108] | 0.8315 [0.8307, 0.8324] |
| 31% HR 1.50 n 1000 | 0.8119 [0.8108, 0.8129] | 0.8777 [0.8772, 0.8783] |
| 31% HR 1.50 n 1500 | 0.8328 [0.8320, 0.8335] | 0.9286 [0.9283, 0.9290] |
| 31% HR 1.75 n 500 | 0.7817 [0.7801, 0.7833] | 0.8594 [0.8586, 0.8602] |
| 31% HR 1.75 n 1000 | 0.8659 [0.8650, 0.8669] | 0.9026 [0.9021, 0.9031] |
| 31% HR 1.75 n 1500 | 0.8784 [0.8778, 0.8791] | 0.9444 [0.9441, 0.9447] |
| 31% HR 1.00 n 500 | 0.4980 [0.4959, 0.5001] | 0.7538 [0.7528, 0.7548] |
| 31% HR 1.00 n 1000 | 0.5795 [0.5781, 0.5809] | 0.7815 [0.7809, 0.7822] |
| 31% HR 1.00 n 1500 | 0.6207 [0.6197, 0.6218] | 0.8200 [0.8195, 0.8206] |

### 1b. Replicate-mean rates, and mean |Ĥ| against |H|

| cell | sens | spec | PPV | NPV | mean \|Ĥ\| | mean \|H\| | ratio of means | median paired ratio |
|---|---|---|---|---|---|---|---|---|
| 12.4% HR 1.50 n 500 | 0.5732 | 0.9148 | 0.5082 | 0.9382 | 73.25 | 62.07 | 1.1801 | 1.1241 |
| 12.4% HR 1.50 n 1000 | 0.6975 | 0.9586 | 0.7059 | 0.9573 | 122.78 | 124.00 | 0.9901 | 1.0000 |
| 12.4% HR 1.50 n 1500 | 0.8198 | 0.9752 | 0.8200 | 0.9747 | 185.20 | 186.10 | 0.9951 | 1.0058 |
| 12.4% HR 1.75 n 500 | 0.6594 | 0.9300 | 0.5902 | 0.9512 | 72.03 | 62.04 | 1.1611 | 1.1120 |
| 12.4% HR 1.75 n 1000 | 0.7867 | 0.9715 | 0.7943 | 0.9701 | 122.51 | 123.93 | 0.9885 | 1.0000 |
| 12.4% HR 1.75 n 1500 | 0.8829 | 0.9831 | 0.8767 | 0.9836 | 186.44 | 186.03 | 1.0022 | 1.0280 |
| 12.4% HR 1.00 n 500 | 0.3520 | 0.8699 | 0.2920 | 0.9044 | 78.90 | 61.99 | 1.2728 | 1.1935 |
| 12.4% HR 1.00 n 1000 | 0.4127 | 0.9118 | 0.4127 | 0.9162 | 128.52 | 124.23 | 1.0346 | 0.9909 |
| 12.4% HR 1.00 n 1500 | 0.5295 | 0.9306 | 0.5285 | 0.9330 | 189.61 | 186.36 | 1.0174 | 0.9834 |
| 31% HR 1.50 n 500 | 0.5903 | 0.8934 | 0.7063 | 0.8372 | 127.12 | 152.95 | 0.8312 | 0.7933 |
| 31% HR 1.50 n 1000 | 0.7078 | 0.9276 | 0.8040 | 0.8868 | 266.81 | 306.25 | 0.8712 | 0.9119 |
| 31% HR 1.50 n 1500 | 0.8391 | 0.9255 | 0.8374 | 0.9365 | 463.08 | 459.69 | 1.0074 | 1.0618 |
| 31% HR 1.75 n 500 | 0.6603 | 0.9189 | 0.7785 | 0.8659 | 128.93 | 152.95 | 0.8430 | 0.8207 |
| 31% HR 1.75 n 1000 | 0.7689 | 0.9475 | 0.8600 | 0.9114 | 271.74 | 306.25 | 0.8873 | 0.9495 |
| 31% HR 1.75 n 1500 | 0.8740 | 0.9465 | 0.8827 | 0.9509 | 457.27 | 459.69 | 0.9947 | 1.0546 |
| 31% HR 1.00 n 500 | 0.3872 | 0.8279 | 0.4937 | 0.7570 | 118.99 | 153.07 | 0.7774 | 0.7296 |
| 31% HR 1.00 n 1000 | 0.4599 | 0.8526 | 0.5646 | 0.7884 | 243.09 | 306.23 | 0.7938 | 0.7396 |
| 31% HR 1.00 n 1500 | 0.5810 | 0.8430 | 0.6177 | 0.8288 | 430.64 | 459.75 | 0.9367 | 0.9142 |

**What the classification block says, descriptively.** Sensitivity and PPV rise monotonically with
n at every (prevalence, HR) combination, and rise with HR at every (prevalence, n). The size of the
identified region tracks the planted one closely from n 1000 upward at 12.4% (ratio of means 0.988
–1.035) but runs **18% large at n 500** there, while at 31% it runs **17% small at n 500** and
converges from below. The null cells are not a different mechanism, only a weaker one: sensitivity
0.35–0.58 and PPV 0.28–0.62, still rising with n.

---

## 2. Bound location on the HR scale

Same columns as the DINA location tables (`scripts_dinamr/blockA_rest.R:63–76`): `fld_H_lo1s` (the
field lower bound), `betaHhat_H` (the realized θ(Ĥ) on the super-population), `fld_H_est2`,
`nv_H_est`, and `truth$marg_H`.

### 2a. Levels

| cell | focus | ε | n_eval | median lower bound | median θ(Ĥ) | median field est2 | median naive | planted marginal θ(H) |
|---|---|---|---|---|---|---|---|---|
| 12.4% HR 1.50 n 500 | maxeffCons | 0.10 | 1822 | 0.5019 | 0.9766 | 0.9382 | 1.745 | 1.5086 |
| 12.4% HR 1.50 n 1000 | maxeffCons | 0.10 | 1948 | 0.6729 | 1.2885 | 1.0840 | 1.601 | 1.5086 |
| 12.4% HR 1.50 n 1500 | maxeffCons | 0.10 | 1976 | 0.8059 | 1.3993 | 1.2096 | 1.522 | 1.5086 |
| 12.4% HR 1.75 n 500 | maxeffCons | 0.10 | 1900 | 0.5575 | 1.1878 | 1.0472 | 1.877 | 1.7691 |
| 12.4% HR 1.75 n 1000 | maxeffCons | 0.10 | 1990 | 0.8098 | 1.6120 | 1.3092 | 1.792 | 1.7691 |
| 12.4% HR 1.75 n 1500 | maxeffCons | 0.10 | 1998 | 1.0135 | 1.7390 | 1.5108 | 1.748 | 1.7691 |
| 12.4% HR 1.00 n 500 | maxeffCons | 0.10 | 1361 | 0.4276 | 0.7053 | 0.7928 | 1.580 | 1.0005 |
| 12.4% HR 1.00 n 1000 | maxeffCons | 0.10 | 1319 | 0.5116 | 0.7448 | 0.8138 | 1.356 | 1.0005 |
| 12.4% HR 1.00 n 1500 | maxeffCons | 0.10 | 1248 | 0.5571 | 0.8029 | 0.8113 | 1.214 | 1.0005 |
| 31% HR 1.50 n 500 | effMaxSG | 0.20 | 1999 | 0.6482 | 1.2733 | 1.1505 | 1.876 | 1.4990 |
| 31% HR 1.50 n 1000 | effMaxSG | 0.20 | 2000 | 0.8329 | 1.3417 | 1.2642 | 1.663 | 1.4990 |
| 31% HR 1.50 n 1500 | effMaxSG | 0.20 | 2000 | 0.9535 | 1.3778 | 1.3135 | 1.510 | 1.4990 |
| 31% HR 1.75 n 500 | effMaxSG | 0.20 | 1999 | 0.7594 | 1.5196 | 1.3588 | 2.102 | 1.7462 |
| 31% HR 1.75 n 1000 | effMaxSG | 0.20 | 2000 | 0.9997 | 1.6894 | 1.5144 | 1.895 | 1.7462 |
| 31% HR 1.75 n 1500 | effMaxSG | 0.20 | 2000 | 1.1595 | 1.6965 | 1.5719 | 1.753 | 1.7462 |
| 31% HR 1.00 n 500 | effMaxSG | 0.20 | 1841 | 0.4591 | 0.8241 | 0.8108 | 1.466 | 0.9999 |
| 31% HR 1.00 n 1000 | effMaxSG | 0.20 | 1909 | 0.5512 | 0.8494 | 0.8326 | 1.232 | 0.9999 |
| 31% HR 1.00 n 1500 | effMaxSG | 0.20 | 1918 | 0.6032 | 0.8709 | 0.8279 | 1.091 | 0.9999 |

### 2b. The gap, three ways, and the location shares with Wilson intervals

| cell | focus | ε | bound − θ | bound / θ | median paired ratio | **share ≥ 1.00** | **share ≥ 1.25** |
|---|---|---|---|---|---|---|---|
| 12.4% HR 1.50 n 500 | maxeffCons | 0.10 | −0.4747 | 0.5140 | 0.5206 | 0.0274 [0.0209, 0.0360] | 0.0115 [0.0076, 0.0176] |
| 12.4% HR 1.50 n 1000 | maxeffCons | 0.10 | −0.6156 | 0.5222 | 0.5732 | 0.1129 [0.0996, 0.1278] | 0.0385 [0.0308, 0.0480] |
| 12.4% HR 1.50 n 1500 | maxeffCons | 0.10 | −0.5934 | 0.5759 | 0.6168 | 0.2126 [0.1951, 0.2311] | 0.0557 [0.0464, 0.0667] |
| 12.4% HR 1.75 n 500 | maxeffCons | 0.10 | −0.6303 | 0.4694 | 0.4909 | 0.0784 [0.0672, 0.0914] | 0.0263 [0.0200, 0.0345] |
| 12.4% HR 1.75 n 1000 | maxeffCons | 0.10 | −0.8022 | 0.5024 | 0.5495 | 0.2603 [0.2415, 0.2800] | 0.1070 [0.0942, 0.1214] |
| 12.4% HR 1.75 n 1500 | maxeffCons | 0.10 | −0.7256 | 0.5828 | 0.6301 | 0.5220 [0.5001, 0.5439] | 0.2152 [0.1978, 0.2338] |
| 12.4% HR 1.00 n 500 | maxeffCons | 0.10 | −0.2777 | 0.6062 | 0.6036 | 0.0029 [0.0011, 0.0075] | 0.0000 [0.0000, 0.0028] |
| 12.4% HR 1.00 n 1000 | maxeffCons | 0.10 | −0.2332 | 0.6869 | 0.6933 | 0.0076 [0.0041, 0.0139] | 0.0023 [0.0008, 0.0067] |
| 12.4% HR 1.00 n 1500 | maxeffCons | 0.10 | −0.2459 | 0.6938 | 0.7026 | 0.0064 [0.0033, 0.0126] | 0.0008 [0.0001, 0.0045] |
| 31% HR 1.50 n 500 | effMaxSG | 0.20 | −0.6251 | 0.5091 | 0.5303 | 0.0875 [0.0759, 0.1007] | 0.0275 [0.0212, 0.0356] |
| 31% HR 1.50 n 1000 | effMaxSG | 0.20 | −0.5088 | 0.6208 | 0.6202 | 0.2280 [0.2101, 0.2469] | 0.0470 [0.0386, 0.0572] |
| 31% HR 1.50 n 1500 | effMaxSG | 0.20 | −0.4243 | 0.6921 | 0.6920 | 0.3990 [0.3778, 0.4206] | 0.0795 [0.0684, 0.0922] |
| 31% HR 1.75 n 500 | effMaxSG | 0.20 | −0.7602 | 0.4997 | 0.5158 | 0.2026 [0.1856, 0.2208] | 0.0695 [0.0592, 0.0815] |
| 31% HR 1.75 n 1000 | effMaxSG | 0.20 | −0.6897 | 0.5918 | 0.6211 | 0.4995 [0.4776, 0.5214] | 0.1830 [0.1667, 0.2005] |
| 31% HR 1.75 n 1500 | effMaxSG | 0.20 | −0.5369 | 0.6835 | 0.7049 | 0.7975 [0.7793, 0.8145] | 0.3325 [0.3122, 0.3535] |
| 31% HR 1.00 n 500 | effMaxSG | 0.20 | −0.3650 | 0.5571 | 0.5637 | 0.0060 [0.0033, 0.0107] | 0.0005 [0.0001, 0.0031] |
| 31% HR 1.00 n 1000 | effMaxSG | 0.20 | −0.2982 | 0.6489 | 0.6558 | 0.0105 [0.0068, 0.0161] | 0.0016 [0.0005, 0.0046] |
| 31% HR 1.00 n 1500 | effMaxSG | 0.20 | −0.2677 | 0.6926 | 0.7058 | 0.0094 [0.0059, 0.0148] | 0.0010 [0.0003, 0.0038] |

**What the location block says, descriptively.**

- The bound sits **below** the realized θ(Ĥ) in every one of the 18 cells; `bound − θ` is negative
  throughout and `bound / θ` runs 0.469–0.694. The ratio and the median paired ratio agree to
  within about 0.05 everywhere, so the ordering is not an artefact of taking medians separately.
- On the **harm** cells the share of lower bounds at or above 1.00 rises steeply with n and with
  HR: at 12.4%, 0.027 → 0.113 → 0.213 at HR 1.50 and 0.078 → 0.260 → 0.522 at HR 1.75; at 31%,
  0.088 → 0.228 → 0.399 and 0.203 → 0.500 → 0.798. The share at or above 1.25 follows the same
  ordering at roughly a third of the level.
- On the **null** cells the same shares stay at **0.003–0.011** (≥ 1.00) and **0.000–0.002**
  (≥ 1.25) at both prevalences and every n, and do **not** rise with n — at 12.4% they are flat
  within their Wilson intervals and at 31% they rise only from 0.0060 to 0.0094, an interval
  overlap. This is the contrast the location columns exist to carry: the bound moves with planted
  harm and does not move with n in its absence.
- The naive median runs above θ(Ĥ) in every cell and falls toward it with n, which is the optimism
  the de-biasing addresses; the field `est2` sits between the bound and θ(Ĥ) throughout.

---

## Quantities not derivable from the committed columns

Recorded here rather than obtained by changing a recorder, as the task directs:

- **The 2×2 counts are not stored**, only the four rates plus `n_sel` and `n_true`. They are
  exactly recoverable (checked at 5.68e-14), so nothing was lost — but a cell whose replicate had
  `n_true = 0` or `n_sel = 0` would have `NA` rates and no recoverable counts. No such replicate
  occurs in these bundles.
- **Classification on the non-detected replicates is undefined**, not missing: `sens`/`spec`/`ppv`/
  `npv` are written only when a subgroup exists. Every rate above is conditional on detection, and
  the detection rate is stated beside it in §1a so the conditioning is visible.
- **Per-subject membership is not stored**, so no classification metric finer than the four rates
  (for instance, agreement between two cells' Ĥ on the same draw) is derivable. No recorder change
  is proposed for it.

## Provenance

- Script: `scripts_dinamr/fs_extraction.R`; saved object `scripts_dinamr/fs_extraction.rds`.
- Bundles: 18 of 18 grid cells present on disk and read; none missing.
