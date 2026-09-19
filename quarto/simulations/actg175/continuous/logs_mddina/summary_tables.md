
=== TABLE 1 [section cells] Combined bundles read (2,000 replicates per cell where present); the commit is the one that added the bundle. (conditional on the proposed family)
|cell|status|
|---|---|
|md40 n500|2000 reps, 1996 declared (0.9980), pkg 0.3.5, host pop-os, identifier dina, dina_select_statistic effect, dina_args list(), commit 1856c107|
|md120 n500|2000 reps, 2000 declared (1.0000), pkg 0.3.5, host pop-os, identifier dina, dina_select_statistic effect, dina_args list(), commit c9f1b948|
|null n500 (no subgroup; homogeneous +26)|2000 reps, 1987 declared (0.9935), pkg 0.3.5, host pop-os, identifier dina, dina_select_statistic effect, dina_args list(), commit 9051b702|
|md40 n700|2000 reps, 1997 declared (0.9985), pkg 0.3.5, host pop-os, identifier dina, dina_select_statistic effect, dina_args list(), commit 21e48ee5|

=== TABLE 2 [section cells] Declaration per cell: replicates, declared, rate with its Monte Carlo SE, and FS’s (mdsgnb20) rate on the same seeds from the committed extract. The null cell has no subgroup; its declarations are read against a homogeneous +26. (conditional on the proposed family)
|cell|replicates|declared|rate|mc_se|fs_rate|
|---|---|---|---|---|---|
|md40 n500|2000|1996|0.9980|0.0010|0.9990|
|md120 n500|2000|2000|1.0000|0.0000|1.0000|
|null n500 (no subgroup; homogeneous +26)|2000|1987|0.9935|0.0018|0.9965|
|md40 n700|2000|1997|0.9985|0.0009|0.9995|

=== TABLE 3 [section ĥ-table] Hhat: declaration rate, bias on the MD scale and in SD units, SE/SD, one-sided lower and two-sided coverage. (conditional on the proposed family)
|cell|estimator|n|declaration rate|bias (MD)|bias (SD units)|SE/SD|one-sided LOWER coverage (Wilson)|two-sided coverage (Wilson)|
|---|---|---|---|---|---|---|---|---|
|md40 n500|naive|1996|0.998|61.464|3.591|1.500|0.111 (0.098, 0.125)|0.244 (0.226, 0.264)|
|md40 n500|oracle|1996|0.998|-0.662|-0.033|0.991|0.952 (0.942, 0.961)|0.943 (0.932, 0.952)|
|md40 n500|MR (IJ)|1996|0.998|19.566|1.084|1.919|0.976 (0.968, 0.982)|0.993 (0.988, 0.996)|
|md40 n500|MR (field)|1996|0.998|11.615|0.600|1.278|0.937 (0.926, 0.947)|0.971 (0.963, 0.978)|
|md120 n500|naive|2000|1.000|31.462|1.727|1.225|0.596 (0.574, 0.617)|0.705 (0.685, 0.725)|
|md120 n500|oracle|2000|1.000|-0.684|-0.034|0.991|0.953 (0.942, 0.961)|0.943 (0.932, 0.952)|
|md120 n500|MR (IJ)|2000|1.000|-2.764|-0.122|1.419|0.991 (0.985, 0.994)|0.984 (0.977, 0.988)|
|md120 n500|MR (field)|2000|1.000|-5.940|-0.238|1.010|0.965 (0.955, 0.972)|0.909 (0.896, 0.921)|
|null n500 (no subgroup; homogeneous +26)|naive|1987|0.994|62.973|3.602|1.493|0.115 (0.101, 0.130)|0.233 (0.215, 0.252)|
|null n500 (no subgroup; homogeneous +26)|oracle|0|0.994|NaN|NaN|NA|NaN (NA, NA)|NaN (NA, NA)|
|null n500 (no subgroup; homogeneous +26)|MR (IJ)|1987|0.994|21.827|1.209|1.963|0.977 (0.969, 0.983)|0.992 (0.987, 0.995)|
|null n500 (no subgroup; homogeneous +26)|MR (field)|1987|0.994|13.800|0.712|1.291|0.931 (0.919, 0.941)|0.967 (0.959, 0.974)|
|md40 n700|naive|1997|0.999|63.496|3.705|1.472|0.071 (0.061, 0.083)|0.168 (0.152, 0.185)|
|md40 n700|oracle|1997|0.999|0.100|0.006|0.968|0.941 (0.930, 0.950)|0.943 (0.932, 0.952)|
|md40 n700|MR (IJ)|1997|0.999|19.906|1.092|1.789|0.960 (0.951, 0.968)|0.985 (0.979, 0.989)|
|md40 n700|MR (field)|1997|0.999|11.646|0.596|1.243|0.930 (0.918, 0.941)|0.965 (0.957, 0.973)|

=== TABLE 4 [section ĥᶜ-table] Hhat^c: declaration rate, bias on the MD scale and in SD units, SE/SD, one-sided upper and two-sided coverage; MR (field-s) evaluated, MR (field) beside it. (conditional on the proposed family)
|cell|estimator|n|declaration rate|bias (MD)|bias (SD units)|SE/SD|one-sided UPPER coverage (Wilson)|two-sided coverage (Wilson)|
|---|---|---|---|---|---|---|---|---|
|md40 n500|naive|1996|0.998|-16.964|-1.425|1.105|0.659 (0.638, 0.679)|0.778 (0.759, 0.795)|
|md40 n500|oracle|1996|0.998|-0.234|-0.016|1.007|0.958 (0.948, 0.966)|0.953 (0.943, 0.962)|
|md40 n500|MR (IJ)|1996|0.998|-6.728|-0.519|1.844|0.991 (0.986, 0.994)|1.000 (0.998, 1.000)|
|md40 n500|MR (field)|1996|0.998|-4.630|-0.348|0.970|0.894 (0.880, 0.907)|0.930 (0.918, 0.940)|
|md40 n500|MR (field-s)|1996|0.998|-4.638|-0.350|0.975|0.897 (0.883, 0.909)|0.933 (0.922, 0.943)|
|md120 n500|naive|2000|1.000|-12.253|-0.853|0.978|0.788 (0.770, 0.806)|0.868 (0.852, 0.882)|
|md120 n500|oracle|2000|1.000|-0.301|-0.021|1.002|0.956 (0.946, 0.964)|0.953 (0.942, 0.961)|
|md120 n500|MR (IJ)|2000|1.000|-2.197|-0.145|1.654|0.995 (0.990, 0.997)|0.999 (0.996, 1.000)|
|md120 n500|MR (field)|2000|1.000|-0.674|-0.043|0.878|0.934 (0.923, 0.945)|0.927 (0.915, 0.938)|
|md120 n500|MR (field-s)|2000|1.000|-0.704|-0.045|0.902|0.937 (0.925, 0.947)|0.936 (0.925, 0.946)|
|null n500 (no subgroup; homogeneous +26)|naive|1987|0.994|-16.580|-1.420|1.122|0.671 (0.650, 0.692)|0.788 (0.770, 0.806)|
|null n500 (no subgroup; homogeneous +26)|oracle|1987|0.994|-0.262|-0.023|1.028|0.955 (0.945, 0.963)|0.958 (0.948, 0.966)|
|null n500 (no subgroup; homogeneous +26)|MR (IJ)|1987|0.994|-6.824|-0.528|1.853|0.992 (0.988, 0.995)|0.999 (0.997, 1.000)|
|null n500 (no subgroup; homogeneous +26)|MR (field)|1987|0.994|-4.809|-0.362|0.969|0.890 (0.875, 0.903)|0.931 (0.919, 0.941)|
|null n500 (no subgroup; homogeneous +26)|MR (field-s)|1987|0.994|-4.817|-0.364|0.972|0.893 (0.879, 0.906)|0.933 (0.921, 0.943)|
|md40 n700|naive|1997|0.999|-11.642|-1.127|1.041|0.720 (0.700, 0.739)|0.822 (0.805, 0.838)|
|md40 n700|oracle|1997|0.999|0.201|0.017|1.006|0.955 (0.945, 0.964)|0.950 (0.940, 0.959)|
|md40 n700|MR (IJ)|1997|0.999|-4.207|-0.382|1.816|0.995 (0.991, 0.997)|0.998 (0.996, 0.999)|
|md40 n700|MR (field)|1997|0.999|-2.659|-0.237|0.943|0.910 (0.897, 0.922)|0.926 (0.914, 0.937)|
|md40 n700|MR (field-s)|1997|0.999|-2.664|-0.237|0.946|0.911 (0.898, 0.923)|0.928 (0.916, 0.938)|

=== TABLE 5 [section bound-location-the-d3-ladder] Harm block: one-sided 95% LOWER bound (oriented MD) – location quantiles and the share at or above each tau (‘harm of at least tau supported’), with MC SEs; field | oracle, IJ and naive for orientation. (conditional on the proposed family)
|cell|estimator|n|mean|q05|q25|median|q75|q95|P(L>=0)|P(L>=10)|P(L>=20)|P(L>=30)|P(L>=40)|P(L>=60)|P(L>=80)|P(L>=100)|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|MR (field)|1996|1.2|-27.3|-12.2|0.4|12.8|33.3|0.508 (0.011)|0.301 (0.010)|0.157 (0.008)|0.070 (0.006)|0.029 (0.004)|0.004 (0.001)|0.001 (0.001)|0.000 (0.000)|
|md40 n500|oracle|1996|6.8|-27.0|-5.9|7.3|20.1|39.8|0.638 (0.011)|0.443 (0.011)|0.252 (0.010)|0.119 (0.007)|0.048 (0.005)|0.003 (0.001)|0.000 (0.000)|0.000 (0.000)|
|md40 n500|MR (IJ)|1996|-5.8|-34.4|-18.0|-6.2|5.8|23.6|0.370 (0.011)|0.184 (0.009)|0.076 (0.006)|0.030 (0.004)|0.011 (0.002)|0.001 (0.001)|0.000 (0.000)|0.000 (0.000)|
|md40 n500|naive|1996|50.8|24.5|39.8|50.7|61.3|77.5|0.996 (0.001)|0.990 (0.002)|0.971 (0.004)|0.905 (0.007)|0.747 (0.010)|0.279 (0.010)|0.039 (0.004)|0.004 (0.001)|
|md120 n500|MR (field)|2000|49.1|10.9|32.3|47.4|64.5|89.9|0.986 (0.003)|0.956 (0.005)|0.892 (0.007)|0.782 (0.009)|0.627 (0.011)|0.311 (0.010)|0.107 (0.007)|0.025 (0.003)|
|md120 n500|oracle|2000|86.7|53.0|74.1|87.3|100.1|119.8|1.000 (0.000)|1.000 (0.000)|0.999 (0.001)|0.997 (0.001)|0.986 (0.003)|0.907 (0.006)|0.638 (0.011)|0.252 (0.010)|
|md120 n500|MR (IJ)|2000|41.8|7.0|28.1|41.2|55.4|77.7|0.979 (0.003)|0.931 (0.006)|0.852 (0.008)|0.718 (0.010)|0.521 (0.011)|0.193 (0.009)|0.038 (0.004)|0.003 (0.001)|
|md120 n500|naive|2000|92.3|64.6|80.5|91.7|103.7|121.6|1.000 (0.000)|1.000 (0.000)|1.000 (0.000)|1.000 (0.000)|1.000 (0.000)|0.972 (0.004)|0.760 (0.010)|0.318 (0.010)|
|null n500 (no subgroup; homogeneous +26)|MR (field)|1987|-2.4|-32.1|-15.4|-3.2|9.7|29.8|0.431 (0.011)|0.243 (0.010)|0.120 (0.007)|0.050 (0.005)|0.017 (0.003)|0.003 (0.001)|0.001 (0.001)|0.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|oracle|0|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|
|null n500 (no subgroup; homogeneous +26)|MR (IJ)|1987|-10.2|-38.6|-22.8|-10.3|1.5|19.8|0.278 (0.010)|0.133 (0.008)|0.050 (0.005)|0.018 (0.003)|0.005 (0.002)|0.001 (0.001)|0.000 (0.000)|0.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|naive|1987|46.3|18.7|35.5|46.3|57.4|73.5|0.993 (0.002)|0.979 (0.003)|0.943 (0.005)|0.836 (0.008)|0.656 (0.011)|0.200 (0.009)|0.021 (0.003)|0.002 (0.001)|
|md40 n700|MR (field)|1997|1.6|-26.3|-11.0|-0.4|13.0|34.6|0.494 (0.011)|0.290 (0.010)|0.162 (0.008)|0.075 (0.006)|0.035 (0.004)|0.008 (0.002)|0.001 (0.001)|0.000 (0.000)|
|md40 n700|oracle|1997|12.5|-15.4|0.8|12.1|24.0|41.5|0.764 (0.009)|0.553 (0.011)|0.322 (0.010)|0.165 (0.008)|0.059 (0.005)|0.004 (0.001)|0.000 (0.000)|0.000 (0.000)|
|md40 n700|MR (IJ)|1997|-2.2|-28.8|-13.8|-3.3|8.1|27.8|0.421 (0.011)|0.222 (0.009)|0.107 (0.007)|0.044 (0.005)|0.017 (0.003)|0.002 (0.001)|0.000 (0.000)|0.000 (0.000)|
|md40 n700|naive|1997|53.6|29.3|43.7|53.1|62.7|79.9|0.996 (0.001)|0.995 (0.002)|0.985 (0.003)|0.945 (0.005)|0.825 (0.009)|0.316 (0.010)|0.050 (0.005)|0.007 (0.002)|

=== TABLE 6 [section bound-location-the-d3-ladder] Complement block: one-sided 95% UPPER bound (oriented MD) – location quantiles and the share at or below each tau (‘harm of at most tau supported’), with MC SEs; field-s | oracle, unstudentized field, IJ and naive for orientation. (conditional on the proposed family)
|cell|estimator|n|mean|q05|q25|median|q75|q95|P(U<=0)|P(U<=10)|P(U<=20)|P(U<=30)|P(U<=40)|P(U<=60)|P(U<=80)|P(U<=100)|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|MR (field-s)|1996|47.4|25.8|39.0|47.7|56.2|69.2|0.000 (0.000)|0.006 (0.002)|0.024 (0.003)|0.094 (0.007)|0.271 (0.010)|0.830 (0.008)|0.996 (0.001)|1.000 (0.000)|
|md40 n500|oracle|1996|49.7|27.1|39.9|50.0|59.5|72.5|0.001 (0.001)|0.006 (0.002)|0.020 (0.003)|0.079 (0.006)|0.252 (0.010)|0.758 (0.010)|0.984 (0.003)|1.000 (0.000)|
|md40 n500|MR (field)|1996|47.4|25.1|38.9|47.7|56.6|69.5|0.000 (0.000)|0.007 (0.002)|0.026 (0.004)|0.096 (0.007)|0.275 (0.010)|0.830 (0.008)|0.996 (0.001)|1.000 (0.000)|
|md40 n500|MR (IJ)|1996|63.4|42.8|55.1|63.7|71.5|84.6|0.000 (0.000)|0.000 (0.000)|0.000 (0.000)|0.007 (0.002)|0.033 (0.004)|0.392 (0.011)|0.900 (0.007)|0.999 (0.001)|
|md40 n500|naive|1996|35.5|16.1|27.8|35.6|43.3|55.8|0.002 (0.001)|0.020 (0.003)|0.092 (0.006)|0.318 (0.010)|0.658 (0.011)|0.984 (0.003)|1.000 (0.000)|1.000 (0.000)|
|md120 n500|MR (field-s)|2000|63.6|39.7|53.3|64.0|74.0|88.0|0.000 (0.000)|0.001 (0.000)|0.006 (0.002)|0.015 (0.003)|0.053 (0.005)|0.398 (0.011)|0.868 (0.008)|0.993 (0.002)|
|md120 n500|oracle|2000|49.7|26.8|39.8|49.9|59.5|72.5|0.001 (0.001)|0.006 (0.002)|0.021 (0.003)|0.081 (0.006)|0.253 (0.010)|0.758 (0.010)|0.985 (0.003)|1.000 (0.000)|
|md120 n500|MR (field)|2000|63.1|38.4|52.8|63.5|73.9|88.0|0.000 (0.000)|0.001 (0.000)|0.008 (0.002)|0.018 (0.003)|0.065 (0.005)|0.413 (0.011)|0.871 (0.007)|0.993 (0.002)|
|md120 n500|MR (IJ)|2000|80.3|57.0|70.5|80.5|90.1|104.1|0.000 (0.000)|0.000 (0.000)|0.000 (0.000)|0.001 (0.000)|0.008 (0.002)|0.082 (0.006)|0.486 (0.011)|0.918 (0.006)|
|md120 n500|naive|2000|52.2|29.6|43.0|52.4|61.7|74.3|0.000 (0.000)|0.005 (0.002)|0.015 (0.003)|0.054 (0.005)|0.186 (0.009)|0.711 (0.010)|0.981 (0.003)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|MR (field-s)|1987|42.6|21.2|34.1|42.6|51.3|64.7|0.001 (0.001)|0.010 (0.002)|0.045 (0.005)|0.170 (0.008)|0.423 (0.011)|0.903 (0.007)|0.999 (0.001)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|oracle|1987|45.2|26.8|37.5|45.2|52.6|64.0|0.000 (0.000)|0.002 (0.001)|0.014 (0.003)|0.088 (0.006)|0.330 (0.011)|0.896 (0.007)|0.999 (0.001)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|MR (field)|1987|42.6|21.1|34.0|42.6|51.6|64.9|0.001 (0.001)|0.011 (0.002)|0.048 (0.005)|0.170 (0.008)|0.426 (0.011)|0.900 (0.007)|0.999 (0.001)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|MR (IJ)|1987|58.8|38.5|50.4|58.8|67.0|80.5|0.000 (0.000)|0.000 (0.000)|0.002 (0.001)|0.014 (0.003)|0.063 (0.005)|0.543 (0.011)|0.947 (0.005)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|naive|1987|31.2|12.5|23.6|31.1|38.9|51.4|0.007 (0.002)|0.030 (0.004)|0.159 (0.008)|0.458 (0.011)|0.783 (0.009)|0.994 (0.002)|1.000 (0.000)|1.000 (0.000)|
|md40 n700|MR (field-s)|1997|45.6|27.0|38.0|45.8|53.3|64.2|0.000 (0.000)|0.002 (0.001)|0.013 (0.002)|0.079 (0.006)|0.303 (0.010)|0.904 (0.007)|0.998 (0.001)|1.000 (0.000)|
|md40 n700|oracle|1997|46.5|27.3|38.5|46.6|54.6|66.5|0.000 (0.000)|0.002 (0.001)|0.012 (0.002)|0.085 (0.006)|0.294 (0.010)|0.872 (0.007)|0.996 (0.001)|1.000 (0.000)|
|md40 n700|MR (field)|1997|45.6|27.0|37.9|45.7|53.3|64.1|0.000 (0.000)|0.002 (0.001)|0.013 (0.003)|0.079 (0.006)|0.307 (0.010)|0.904 (0.007)|0.998 (0.001)|1.000 (0.000)|
|md40 n700|MR (IJ)|1997|59.5|41.9|52.2|59.4|66.9|77.9|0.000 (0.000)|0.000 (0.000)|0.001 (0.001)|0.005 (0.001)|0.036 (0.004)|0.520 (0.011)|0.970 (0.004)|1.000 (0.000)|
|md40 n700|naive|1997|36.9|20.5|30.0|36.8|43.8|54.1|0.001 (0.001)|0.005 (0.002)|0.047 (0.005)|0.251 (0.010)|0.624 (0.011)|0.989 (0.002)|1.000 (0.000)|1.000 (0.000)|

=== TABLE 7 [section joint-pair] Joint (Hhat lower, Hhat^c upper) coverage of (beta(Hhat), beta(Hhat^c)) on declared replicates where both bounds exist: the field-s Bonferroni pair (evaluated), the unstudentized Bonferroni pair beside it, and the separate 95% pair; share_both = declared replicates carrying both bounds; margins in MD units from beta-tilde. (conditional on the proposed family)
|cell|pair|declared|both_bounds|share_both|joint|joint_mc_se|cov_H|cov_Hc|margin_H|margin_Hc|
|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|Bonferroni field-s (gamma = 0.025)|1996|1996|1|0.914 (0.901, 0.926)|0.006|0.972|0.941|59.103|27.392|
|md40 n500|Bonferroni unstudentized (gamma = 0.025)|1996|1996|1|0.911 (0.898, 0.923)|0.006|0.972|0.937|59.103|27.390|
|md40 n500|separate 95% bounds: field lower, field-s upper|1996|1996|1|0.840 (0.823, 0.856)|0.008|0.937|0.897|49.926|23.346|
|md120 n500|Bonferroni field-s (gamma = 0.025)|2000|2000|1|0.950 (0.940, 0.959)|0.005|0.983|0.967|54.754|28.896|
|md120 n500|Bonferroni unstudentized (gamma = 0.025)|2000|2000|1|0.945 (0.934, 0.954)|0.005|0.983|0.962|54.754|28.267|
|md120 n500|separate 95% bounds: field lower, field-s upper|2000|2000|1|0.903 (0.890, 0.916)|0.007|0.965|0.937|45.698|24.523|
|null n500 (no subgroup; homogeneous +26)|Bonferroni field-s (gamma = 0.025)|1987|1987|1|0.910 (0.897, 0.922)|0.006|0.969|0.941|59.719|27.218|
|null n500 (no subgroup; homogeneous +26)|Bonferroni unstudentized (gamma = 0.025)|1987|1987|1|0.909 (0.895, 0.921)|0.006|0.969|0.939|59.719|27.251|
|null n500 (no subgroup; homogeneous +26)|separate 95% bounds: field lower, field-s upper|1987|1987|1|0.832 (0.815, 0.848)|0.008|0.931|0.893|50.475|23.165|
|md40 n700|Bonferroni field-s (gamma = 0.025)|1997|1997|1|0.911 (0.898, 0.923)|0.006|0.966|0.944|59.135|22.295|
|md40 n700|Bonferroni unstudentized (gamma = 0.025)|1997|1997|1|0.909 (0.896, 0.921)|0.006|0.966|0.943|59.135|22.296|
|md40 n700|separate 95% bounds: field lower, field-s upper|1997|1997|1|0.846 (0.829, 0.861)|0.008|0.930|0.911|49.855|18.989|

=== TABLE 8 [section joint-pair] Calibrated split (diagnostics): mean gamma and corr(Lambda, Lambdac) for the studentized (_s) and unstudentized pairs. (conditional on the proposed family)
|cell|metric|value|mc_se|n|
|---|---|---|---|---|
|md40 n500|gamma_mean_s|0.0252|0.0000|1996|
|md40 n500|corr_s|0.0879|0.0013|1996|
|md40 n500|gamma_mean|0.0251|0.0000|1996|
|md40 n500|corr|0.0905|0.0014|1996|
|md120 n500|gamma_mean_s|0.0252|0.0000|2000|
|md120 n500|corr_s|0.0111|0.0017|2000|
|md120 n500|gamma_mean|0.0252|0.0000|2000|
|md120 n500|corr|0.0148|0.0017|2000|
|null n500 (no subgroup; homogeneous +26)|gamma_mean_s|0.0252|0.0000|1987|
|null n500 (no subgroup; homogeneous +26)|corr_s|0.0855|0.0013|1987|
|null n500 (no subgroup; homogeneous +26)|gamma_mean|0.0251|0.0000|1987|
|null n500 (no subgroup; homogeneous +26)|corr|0.0878|0.0013|1987|
|md40 n700|gamma_mean_s|0.0252|0.0000|1997|
|md40 n700|corr_s|0.0879|0.0012|1997|
|md40 n700|gamma_mean|0.0251|0.0000|1997|
|md40 n700|corr|0.0901|0.0013|1997|

=== TABLE 9 [section identification-dina-beside-fs-and-grf] DINA (this campaign) beside GRF (mdgrf) and FS (mdsgnb20), the same seeds: declaration, size of Hhat, true positives, sensitivity, PPV, proposed, family and admitted sizes, and the paired size comparison with FS. (conditional on the proposed family)
|cell|quantity|DINA|GRF_mdgrf|FS_mdsgnb20|
|---|---|---|---|---|
|md40 n500|declaration rate|0.9980|1.0000|0.9990|
|md40 n500|mean size of Hhat (n_sel)|108.53|109.46|111.62|
|md40 n500|mean true positives (sens x n_true)|42.39|-|-|
|md40 n500|sensitivity (mean over declared; NA where n_true = 0)|0.2462|0.2497|0.2679|
|md40 n500|PPV (mean over declared)|0.3888|0.3909|0.4122|
|md40 n500|mean proposed family (dina_proposed_n)|3341.7|-|-|
|md40 n500|mean family size K (n_family, MR’s kept family)|3341.7|1222.4|-|
|md40 n500|mean admitted set (admitted_n)|2631.0|624.6|-|
|md40 n500|size of Hhat DINA larger / equal / smaller than FS (paired by sim_id, both declared)|852 / 93 / 1050|-|-|
|md120 n500|declaration rate|1.0000|1.0000|1.0000|
|md120 n500|mean size of Hhat (n_sel)|146.59|143.61|152.74|
|md120 n500|mean true positives (sens x n_true)|113.50|-|-|
|md120 n500|sensitivity (mean over declared; NA where n_true = 0)|0.6596|0.5989|0.7116|
|md120 n500|PPV (mean over declared)|0.7666|0.7149|0.7980|
|md120 n500|mean proposed family (dina_proposed_n)|5940.0|-|-|
|md120 n500|mean family size K (n_family, MR’s kept family)|5940.0|1222.4|-|
|md120 n500|mean admitted set (admitted_n)|5599.1|1059.8|-|
|md120 n500|size of Hhat DINA larger / equal / smaller than FS (paired by sim_id, both declared)|757 / 78 / 1165|-|-|
|null n500 (no subgroup; homogeneous +26)|declaration rate|0.9935|1.0000|0.9965|
|null n500 (no subgroup; homogeneous +26)|mean size of Hhat (n_sel)|105.15|106.14|107.49|
|null n500 (no subgroup; homogeneous +26)|mean true positives (sens x n_true)|NA|-|-|
|null n500 (no subgroup; homogeneous +26)|sensitivity (mean over declared; NA where n_true = 0)|NA|NA|NA|
|null n500 (no subgroup; homogeneous +26)|PPV (mean over declared)|0.0000|0.0000|0.0000|
|null n500 (no subgroup; homogeneous +26)|mean proposed family (dina_proposed_n)|2644.1|-|-|
|null n500 (no subgroup; homogeneous +26)|mean family size K (n_family, MR’s kept family)|2644.1|1222.4|-|
|null n500 (no subgroup; homogeneous +26)|mean admitted set (admitted_n)|1959.0|508.5|-|
|null n500 (no subgroup; homogeneous +26)|size of Hhat DINA larger / equal / smaller than FS (paired by sim_id, both declared)|867 / 104 / 1010|-|-|
|md40 n700|declaration rate|0.9985|1.0000|0.9995|
|md40 n700|mean size of Hhat (n_sel)|113.84|116.33|117.99|
|md40 n700|mean true positives (sens x n_true)|44.92|-|-|
|md40 n700|sensitivity (mean over declared; NA where n_true = 0)|0.1862|0.1873|0.2041|
|md40 n700|PPV (mean over declared)|0.3895|0.3863|0.4101|
|md40 n700|mean proposed family (dina_proposed_n)|4007.0|-|-|
|md40 n700|mean family size K (n_family, MR’s kept family)|4007.0|1354.5|-|
|md40 n700|mean admitted set (admitted_n)|3202.8|711.4|-|
|md40 n700|size of Hhat DINA larger / equal / smaller than FS (paired by sim_id, both declared)|873 / 80 / 1043|-|-|

=== TABLE 10 [section the-three-identifiers] FS, GRF and DINA on the same four cells and seeds: declaration, one-sided coverage of the field lower bound on Hhat and the field-s upper bound on Hhat^c, field-s Bonferroni joint coverage, mean size of Hhat, sensitivity and PPV; FS and GRF copied from md_field_metrics.csv and md_grf_metrics.csv. The null cell has no subgroup (homogeneous +26); its sensitivity is undefined. (conditional on the proposed family)
|cell|identifier|declaration rate|field lower coverage, Hhat|field-s upper coverage, Hhat^c|Bonferroni joint (field-s)|mean size of Hhat|sensitivity|PPV|
|---|---|---|---|---|---|---|---|---|
|md40 n500|FS (mdsgnb20)|0.9990|0.9479|0.9139|0.9279|111.6241|0.2679|0.4122|
|md40 n500|GRF (mdgrf)|1.0000|0.9390|0.9200|0.9195|109.4590|0.2497|0.3909|
|md40 n500|DINA (mddina)|0.9980|0.9374|0.8968|0.9143|108.5326|0.2462|0.3888|
|md120 n500|FS (mdsgnb20)|1.0000|0.9670|0.9405|0.9490|152.7385|0.7116|0.7980|
|md120 n500|GRF (mdgrf)|1.0000|0.9560|0.9315|0.9425|143.6125|0.5989|0.7149|
|md120 n500|DINA (mddina)|1.0000|0.9645|0.9370|0.9500|146.5870|0.6596|0.7666|
|null n500 (no subgroup; homogeneous +26)|FS (mdsgnb20)|0.9965|0.9473|0.9172|0.9318|107.4922|NA|0.0000|
|null n500 (no subgroup; homogeneous +26)|GRF (mdgrf)|1.0000|0.9405|0.9200|0.9300|106.1385|NA|0.0000|
|null n500 (no subgroup; homogeneous +26)|DINA (mddina)|0.9935|0.9311|0.8933|0.9104|105.1500|NA|0.0000|
|md40 n700|FS (mdsgnb20)|0.9995|0.9395|0.9280|0.9275|117.9905|0.2041|0.4101|
|md40 n700|GRF (mdgrf)|1.0000|0.9430|0.9310|0.9345|116.3330|0.1873|0.3863|
|md40 n700|DINA (mddina)|0.9985|0.9304|0.9114|0.9109|113.8353|0.1862|0.3895|

=== TABLE 11 [section regime-diagnostics] Per cell: declared replicates, p-hat(Hhat) mean and share < 0.5 (tie regime), SD(beta-tilde^c)/mean naive SE^c, lambda-SD^c/naive SE^c for the field and for field-s, IJ SE / empirical SD per block, seconds per replicate (fit + MR + field; field; complement). (conditional on the proposed family)
|cell|n_det|p_hat_mean|p_hat_lt05|sd_btc_naive|lamc_naive|lamc_s_naive|ij_sd_H|ij_sd_Hc|fit_secs|field_secs|comp_secs|
|---|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|1996|0.080|0.996|0.985|0.979|0.982|1.919|1.844|134.235|69.798|6.999|
|md120 n500|2000|0.083|0.999|1.076|0.972|0.998|1.419|1.654|304.939|149.201|17.650|
|null n500 (no subgroup; homogeneous +26)|1987|0.092|0.990|0.986|0.981|0.983|1.963|1.853|100.276|54.362|4.951|
|md40 n700|1997|0.076|0.997|1.024|0.985|0.987|1.789|1.816|199.655|88.036|12.555|

=== TABLE 12 [section the-display-identity-scale] fs_sim_bias_coverage(scale = ‘identity’): b = retained bias / SD, r = mean SE / SD; refs = Phi(1.645 r -/+ b) and Phi(1.96 r - b) - Phi(-1.96 r - b). Hhat^c ‘fld_s’ rows come from a copy of the results with fld_Hc_s renamed to fld_Hc (field-s). (conditional on the proposed family)
|cell|block|estimator|n|bias (MD)|SD|mean SE|b|r|1-sided cov|1-sided ref|2-sided cov|2-sided ref|
|---|---|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|H|naive|1996|61.464|17.117|25.674|3.591|1.500|0.111|0.131|0.244|0.257|
|md40 n500|H|mr|1996|19.566|18.048|34.641|1.084|1.919|0.976|0.981|0.993|0.996|
|md40 n500|H|fld|1996|11.615|19.362|24.750|0.600|1.278|0.937|0.934|0.971|0.971|
|md40 n500|Hc|naive|1996|-16.964|11.907|13.160|-1.425|1.105|0.659|0.653|0.778|0.771|
|md40 n500|Hc|mr|1996|-6.728|12.960|23.904|-0.519|1.844|0.991|0.994|1.000|0.999|
|md40 n500|Hc|fld|1996|-4.630|13.289|12.889|-0.348|0.970|0.894|0.894|0.930|0.927|
|md40 n500|Hc|fld_s|1996|-4.638|13.263|12.927|-0.350|0.975|0.897|0.895|0.933|0.929|
|md120 n500|H|naive|2000|31.462|18.216|22.308|1.727|1.225|0.596|0.613|0.705|0.750|
|md120 n500|H|mr|2000|-2.764|22.680|32.192|-0.122|1.419|0.991|0.993|0.984|0.994|
|md120 n500|H|fld|2000|-5.940|24.997|25.239|-0.238|1.010|0.965|0.971|0.909|0.946|
|md120 n500|Hc|naive|2000|-12.253|14.367|14.051|-0.853|0.978|0.788|0.775|0.868|0.854|
|md120 n500|Hc|mr|2000|-2.197|15.122|25.016|-0.145|1.654|0.995|0.995|0.999|0.999|
|md120 n500|Hc|fld|2000|-0.674|15.555|13.657|-0.043|0.878|0.934|0.919|0.927|0.914|
|md120 n500|Hc|fld_s|2000|-0.704|15.549|14.019|-0.045|0.902|0.937|0.925|0.936|0.922|
|null n500 (no subgroup; homogeneous +26)|H|naive|1987|62.973|17.485|26.102|3.602|1.493|0.115|0.126|0.233|0.250|
|null n500 (no subgroup; homogeneous +26)|H|mr|1987|21.827|18.048|35.421|1.209|1.963|0.977|0.978|0.992|0.996|
|null n500 (no subgroup; homogeneous +26)|H|fld|1987|13.800|19.374|25.012|0.712|1.291|0.931|0.921|0.967|0.965|
|null n500 (no subgroup; homogeneous +26)|Hc|naive|1987|-16.580|11.679|13.101|-1.420|1.122|0.671|0.665|0.788|0.782|
|null n500 (no subgroup; homogeneous +26)|Hc|mr|1987|-6.824|12.914|23.930|-0.528|1.853|0.992|0.994|0.999|0.999|
|null n500 (no subgroup; homogeneous +26)|Hc|fld|1987|-4.809|13.267|12.858|-0.362|0.969|0.890|0.891|0.931|0.926|
|null n500 (no subgroup; homogeneous +26)|Hc|fld_s|1987|-4.817|13.244|12.873|-0.364|0.972|0.893|0.892|0.933|0.927|
|md40 n700|H|naive|1997|63.496|17.137|25.225|3.705|1.472|0.071|0.100|0.168|0.206|
|md40 n700|H|mr|1997|19.906|18.235|32.616|1.092|1.789|0.960|0.968|0.985|0.992|
|md40 n700|H|fld|1997|11.646|19.546|24.288|0.596|1.243|0.930|0.926|0.965|0.966|
|md40 n700|Hc|naive|1997|-11.642|10.330|10.753|-1.127|1.041|0.720|0.721|0.822|0.819|
|md40 n700|Hc|mr|1997|-4.207|11.010|19.997|-0.382|1.816|0.995|0.995|0.998|0.999|
|md40 n700|Hc|fld|1997|-2.659|11.229|10.586|-0.237|0.943|0.910|0.906|0.926|0.928|
|md40 n700|Hc|fld_s|1997|-2.664|11.219|10.617|-0.237|0.946|0.911|0.906|0.928|0.929|

EXTRACT LINE: X <- do.call(rbind, EX)
gi <- X$identifier == "dina"
X$value[gi] <- signif(X$value[gi], 8); X$mc_se[gi] <- signif(X$mc_se[gi], 6); X$wilson_lo[gi] <- signif(X$wilson_lo[gi], 6); X$wilson_hi[gi] <- signif(X$wilson_hi[gi], 6)
write.csv(X, file.path(out_dir, "md_dina_metrics.csv"), row.names = FALSE, na = "")
cols <- c(
"# COLUMNS — md_dina_metrics.csv (campaign mddina; TASK_md_dina_campaign_2026-09-17 §3.2)", "",
"Written by `summary_continuous_field_mddina.qmd` from the same objects its tables print. The schema of `md_grf_metrics.csv` (see `COLUMNS_md_grf.md`, itself `md_field_metrics.csv`'s plus `identifier`). One row per cell × block × estimator × metric (× tau for bound-location rows).", "",
"## Columns", "",
"- `campaign`: `mddina`; `mdsgnb20` for the FS comparator rows; `mdgrf` for the GRF comparator rows.",
"- `identifier`: `dina` (this campaign); `fs` for the FS comparator rows, copied from the committed `md_field_metrics.csv`; `grf` for the GRF comparator rows, copied from the committed `md_grf_metrics.csv`. Comparator rows keep their own `commit` and are not recomputed.",
"- `cell`: `md40 n500`, `md120 n500`, `null n500 (no subgroup; homogeneous +26)`, `md40 n700`.",
"- `block`: `H` (the selected subgroup Ĥ), `Hc` (its complement Ĥᶜ), `joint` (the pair), `all` (timing).",
"- `estimator`, `metric`, `tau`, `value`, `mc_se`, `wilson_lo`, `wilson_hi`, `n`: as in `COLUMNS_md_field.md`.",
"- `commit`: the commit that added the cell's combined bundle (for `fs` and `grf` rows, the value in the committed extract).", "",
"## Estimator codes", "",
"`naive`, `oracle`, `mr`, `fld`, `fld_s`, `bonf_s`, `bonf`, `separate_s`, `calibrated`, `all`: as in `COLUMNS_md_field.md`, applied to the DINA-selected Ĥ.", "",
"## Metrics", "",
"As in `COLUMNS_md_grf.md` (itself `COLUMNS_md_field.md`), with these identification rows:", "",
"- `mean_size_hhat`: mean |Ĥ| over declared replicates (the recorder's `n_sel`; its `n_harm` is the same count, the size of Ĥ, not a true-positive count).",
"- `mean_true_positives`: mean number of truly harmed patients in Ĥ, sensitivity × `n_true`, over declared replicates with `n_true` > 0 (undefined in the null cell).",
"- `sensitivity_mean`, `ppv_mean`: the recorder's `sens`, `ppv` (`.classify`), mean over declared replicates.",
"- `mean_proposed_n`: mean of DINA's proposed family (`dina_proposed_n`: candidates with oriented tau-hat at or above the floor 30 and at least `n.min` members); `mean_n_family`: mean of MR's kept family size K; `mean_admitted_n`: mean of DINA's admitted set (proposed candidates whose harm-oriented MD clears the admission floor 30). All three over declared replicates: DINA's selection object, which carries the proposed and admitted counts, exists only when a subgroup is selected.",
"- `share_size_larger_than_fs`, `share_size_equal_than_fs`, `share_size_smaller_than_fs` and the matching `count_size_*_than_fs`: DINA's |Ĥ| against FS's (`mdsgnb20`) on the same `sim_id`, over replicates both identifiers declared.",
"- FS rows (`identifier = fs`): `declaration_rate`, `mean_n_sel`, `sensitivity_mean`, `ppv_mean`, and the three coverage rows of the three-identifier table (`H`/`fld`/`cov1_lower`, `Hc`/`fld_s`/`cov1_upper`, `joint`/`bonf_s`/`joint_coverage`). **FS's `mean_n_harm` equals its `mean_n_sel`** in `md_field_metrics.csv` (a labelling defect: `n_harm` is |Ĥ|, not a true-positive count); it is not copied.",
"- GRF rows (`identifier = grf`): `declaration_rate`, `mean_size_hhat`, `sensitivity_mean`, `ppv_mean`, `mean_n_family`, `mean_admitted_n`, and the same three coverage rows.",
"- `secs_fit_mr_mean`, `secs_field_mean`, `secs_complement_mean` (block `all`): mean seconds per replicate for the fit with MR (all replicates), the field pass and the complement field (declared replicates), as the regime table prints them; `fit_mr_secs` contains the other two, which are never summed.", "",
"## Scale convention", "",
"Every estimate, bound, bias and threshold is on the **harm-oriented mean-difference scale** (positive = harm; adverse_outcome = FALSE); `betaHhat_*` in the bundles are raw cd4_change and are oriented with -1 here. The null cell's truth is a homogeneous +26.255 on this scale for every Ĥ and Ĥᶜ. **Every DINA and GRF coverage figure is coverage of the estimand conditional on the proposed family**: their candidate families are generated from fitted surfaces, so the fixed-family condition does not hold. FS's family is the prespecified cut grid. Comparisons across the three identifiers are descriptive: the identifier, the family construction and the detected set all differ. No significance language: bounds are read by location against the ladder.")
writeLines(cols, file.path(out_dir, "COLUMNS_md_dina.md"))
cat(sprintf("extract: %d rows (%d dina, %d fs, %d grf) -> %s ; %s\n", nrow(X), sum(gi), sum(X$identifier == "fs"), sum(X$identifier == "grf"), file.path(out_dir, "md_dina_metrics.csv"), file.path(out_dir, "COLUMNS_md_dina.md"))) 
