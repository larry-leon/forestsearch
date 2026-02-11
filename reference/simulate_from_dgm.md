# Simulate Survival Data from AFT Data Generating Mechanism

Generates simulated survival data from a previously created AFT data
generating mechanism (DGM). Samples from the super population and
generates survival times with specified censoring.

## Usage

``` r
simulate_from_dgm(
  dgm,
  n = NULL,
  rand_ratio = 1,
  entry_var = NULL,
  max_entry = 24,
  analysis_time = 48,
  cens_adjust = 0,
  draw_treatment = TRUE,
  seed = NULL
)
```

## Arguments

- dgm:

  An object of class "aft_dgm_flex" created by
  [`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)

- n:

  Integer specifying the sample size. If NULL (default), uses the entire
  super population

- rand_ratio:

  Numeric specifying the randomization ratio (treatment:control).
  Default is 1 (1:1 allocation)

- entry_var:

  Character string specifying the name of an entry time variable in the
  data. If NULL, entry times are simulated uniformly. Default NULL

- max_entry:

  Numeric specifying maximum entry time for staggered entry. Entry times
  are simulated as Uniform(0, max_entry). Default 24

- analysis_time:

  Numeric specifying the calendar time at which analysis occurs.
  Follow-up time is calculated as analysis_time - entry_time. Default 48

- cens_adjust:

  Numeric adjustment to censoring distribution on log scale. Positive
  values increase censoring, negative values decrease it. Default is 0

- draw_treatment:

  Logical indicating whether to redraw treatment assignment. If TRUE
  (default), reassigns treatment according to rand_ratio. If FALSE,
  keeps original treatment assignments from super population

- seed:

  Integer random seed for reproducibility. Default is NULL (no seed set)

## Value

A data.frame containing simulated survival data with columns:

- id:

  Subject identifier

- treat:

  Treatment assignment (0 or 1)

- treat_sim:

  Simulated treatment assignment (may differ from treat if
  draw_treatment = TRUE)

- flag_harm:

  Subgroup indicator (1 if all subgroup conditions met, 0 otherwise)

- z\_\*:

  Standardized covariate values

- y_sim:

  Observed survival time (minimum of true time and censoring time)

- event_sim:

  Event indicator (1 = event observed, 0 = censored)

- t_true:

  True underlying survival time (before censoring)

- c_time:

  Censoring time

## Details

### Simulation Process

1.  **Sampling**: Draws n observations with replacement from the super
    population

2.  **Treatment Assignment**:

    - If `draw_treatment = TRUE`: Reassigns treatment based on
      `rand_ratio`

    - If `draw_treatment = FALSE`: Keeps original treatment assignments

3.  **Survival Times**: Generates from Weibull AFT model: \$\$\log(T) =
    \mu + \sigma \epsilon + X'\gamma\$\$ where ε ~ extreme value
    distribution

4.  **Censoring**: Applies specified censoring distribution (Weibull or
    uniform)

5.  **Administrative Censoring**: Applies max_follow cutoff if specified

### Censoring Adjustment

The `cens_adjust` parameter modifies the censoring distribution:

- `cens_adjust = log(2)` doubles expected censoring times

- `cens_adjust = log(0.5)` halves expected censoring times

## See also

[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
for creating the DGM

## Examples

``` r
if (FALSE) { # \dontrun{
# Create DGM first
dgm <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(er = 20, meno = 0),
  model = "alt"
)

# Simulate data with 1:1 randomization
sim_data <- simulate_from_dgm(
  dgm = dgm,
  n = 1000,
  rand_ratio = 1,
  max_follow = 84,
  cens_adjust = log(1.5),
  seed = 123
)

# Check results
table(sim_data$treat_sim)
mean(sim_data$event_sim)

# Simulate with 2:1 randomization
sim_data_2to1 <- simulate_from_dgm(
  dgm = dgm,
  n = 900,
  rand_ratio = 2,  # 2:1 treatment:control
  seed = 456
)
} # }
```
