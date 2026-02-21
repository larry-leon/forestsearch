# Subgroup Identification for Regional Consistency in Multi-Regional Clinical Trials
**Larry Leon**
---

### Abstract

#### Subgroup Identification for Regional Consistency in Multi-Regional Clinical Trials:

The evaluation of regional consistency in multi-regional clinical trials (MRCTs) can be challenging. 
When a pivotal oncology trial meets its primary endpoint but a smaller region — such as an Asia-Pacific (AP) population comprising 15–20% of the trial — produces a hazard ratio 
estimate near 1.0, the consistency criterion may appear to fail even when the treatment is actually effective in that region. 
This apparent inconsistency can delay or jeopardize regulatory approval in affected regions despite strong global evidence of benefit.

We present a principled framework for utilizing exploratory subgroup identification to enhance regional consistency evaluation. The central question is whether a 
subgroup *G* with enhanced treatment benefit can be identified in the non-AP population such that the corresponding AP subgroup, AP(*G*), meets the regional 
consistency criterion (HR < 1.0). This train/test paradigm — identifying subgroups in one region and validating in another — provides a built-in safeguard against overfitting.

Our approach employs ForestSearch (Leon et al., 2024), a subgroup identification algorithm that combines Generalized Random Forests (GRF), LASSO regularization, 
and exhaustive combinatorial search to enumerate candidate subgroups defined by clinical factors and/or biomarkers. Three complementary objectives can be 
evaluated: (A) identifying the subgroup with the largest detrimental effect, (B) identifying the smallest subgroup indicative of limited benefit, 
and (C) directly targeting the largest subgroup with strong benefit. Subgroup candidates are screened against hazard ratio thresholds (e.g., HR ≥ 0.90 for harm detection) 
and evaluated for split-sample consistency across repeated random partitions.

Operating characteristics are evaluated through simulations based on observed data (case-study) representing real clinical characteristics: outcome and censoring models are fit to the case-study dataset, 
and the entire workflow — from subgroup identification through regional consistency testing — is replicated across simulations under both 
the alternative hypothesis (heterogeneous treatment effects) and the null hypothesis (uniform benefit). This enables assessment of false discovery rates
and the ability to characterize the underlying treatment benefit in the region (estimation, CI properties).

The forestsearch R package provides a complete implementation of this framework, including flexible data-generating mechanisms, parallel simulation infrastructure, 
and publication-quality summary tables and visualizations. We illustrate the approach using a synthetic oncology dataset and demonstrate that biomarker-defined subgroups 
identified in non-AP data can meaningfully improve regional consistency evaluation in the AP population.

---

**Keywords:** multi-regional clinical trials, subgroup identification, regional consistency, treatment effect heterogeneity, causal inference, hazard ratio, ForestSearch, survival analysis
