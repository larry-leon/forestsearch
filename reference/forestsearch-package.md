# forestsearch: Exploratory Subgroup Identification in Survival and GLM Outcomes

Implements statistical methods for exploratory subgroup identification
in clinical trials. Provides tools for identifying patient subgroups
with differential treatment effects using machine learning approaches
including Generalized Random Forests (GRF), LASSO regularization, and
exhaustive combinatorial search algorithms. Supports survival endpoints
(Cox proportional hazards), binary outcomes (log odds ratio, log
relative risk, risk difference), continuous outcomes (mean difference),
and count / rate outcomes (log incidence rate ratio via Poisson,
quasi-Poisson, or negative-binomial GLMs with optional person-time
offset). Features bootstrap bias correction using infinitesimal
jackknife methods to address selection bias in post-hoc analyses.
Designed for clinical researchers conducting exploratory subgroup
analyses in randomized controlled trials, particularly for
multi-regional clinical trials (MRCT) requiring regional consistency
evaluation. Methods are described in Leon et al. (2024)
[doi:10.1002/sim.10163](https://doi.org/10.1002/sim.10163) .
Post-selection inference for the identified subgroup follows León and
Anderson (2026)
[doi:10.48550/arXiv.2609.38361](https://doi.org/10.48550/arXiv.2609.38361)
.

## See also

Useful links:

- <https://github.com/larry-leon/forestsearch>

- <https://larry-leon.github.io/forestsearch/>

- Report bugs at <https://github.com/larry-leon/forestsearch/issues>

## Author

**Maintainer**: Larry Leon <larry.leon.05@post.harvard.edu>

Authors:

- Larry Leon <larry.leon.05@post.harvard.edu>
