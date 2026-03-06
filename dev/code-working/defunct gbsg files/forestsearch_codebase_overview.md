# ForestSearch Codebase Overview

**Package:** forestsearch v0.1.0  
**Total R source lines:** ~31,118 across 46 files  
**Reference:** León et al. (2024), *Statistics in Medicine*, doi:10.1002/sim.10163

---

## File Inventory by Module

### Core Algorithm
| File | Lines | Description |
|------|-------|-------------|
| `forestsearch_main.R` | 772 | Main `forestsearch()` entry point |
| `forestsearch_helpers.R` | 430 | Internal helpers for main search |
| `forestsearch_methods.R` | ~400 | S3 methods: print, summary, etc. |
| `subgroup_search.R` | ~500 | Exhaustive subgroup enumeration |
| `subgroup_consistency_main.R` | 804 | `subgroup.consistency()` — split-sample evaluation |
| `subgroup_consistency_helpers.R` | 1066 | Helpers for consistency evaluation |
| `find_k_inter_main.R` | ~400 | Interaction/k-factor search |

### Bootstrap Inference
| File | Lines | Description |
|------|-------|-------------|
| `bootstrap_dofuture_main.R` | ~600 | `forestsearch_bootstrap_dofuture()` — main bootstrap driver |
| `bootstrap_analysis_dofuture.R` | ~650 | Per-bootstrap-replicate analysis |
| `bootstrap_dofuture_setup_helpers.R` | ~300 | Setup/initialization helpers |
| `bootstrap_calculations_helpers.R` | ~450 | Bias correction calculations (IJ variance) |
| `bootstrap_summaries_helpers.R` | 1067 | Summary/diagnostics for bootstrap output |
| `summarize_bootstrap_results.R` | ~500 | `summarize_bootstrap_results()` — gt tables |
| `summarize_bootstrap_subgroups.R` | ~400 | Subgroup-level bootstrap summaries |

### Cross-Validation
| File | Lines | Description |
|------|-------|-------------|
| `forestsearch_cross-validation.R` | 1154 | `forestsearch_Kfold()`, `forestsearch_tenfold()` |
| `cv_summary_tables.R` | 820 | `cv_metrics_tables()` — sensCV, ppvCV |

### GRF / LASSO Factor Construction
| File | Lines | Description |
|------|-------|-------------|
| `grf_main.R` | ~500 | GRF-based CATE estimation wrapper |
| `grf_helpers.R` | ~450 | GRF utilities, quartile splits |

### Cox / Survival Estimation
| File | Lines | Description |
|------|-------|-------------|
| `Cox_estimation_helpers.R` | ~500 | Cox model fitting, HR extraction |
| `cox_ahr_cde_wrapper.R` | 849 | Average HR and CDE wrappers |
| `cox_spline_fit.R` | 780 | Spline-based treatment effect curves |
| `truefind_asymptotic.R` | 774 | Asymptotic power/type-1 error calculations |

### Data Preparation
| File | Lines | Description |
|------|-------|-------------|
| `get_fsdata.R` | ~400 | `get_fsdata()` — data prep for FS |
| `get_FSdata_helpers.R` | ~350 | Helpers for FSdata construction |

### Simulation Framework
| File | Lines | Description |
|------|-------|-------------|
| `generate_aft_dgm_main.R` | ~500 | `generate_aft_dgm_flex()` — AFT DGM |
| `generate_aft_dgm_helpers.R` | 1584 | AFT DGM helpers |
| `generate_aft_dgm_helpers_outcome.R` | ~600 | Outcome generation for AFT |
| `simulate_from_dgm.R` | 848 | `run_simulation_analysis()` |
| `sim_aft_gbsg.R` | 1280 | GBSG-based simulation |
| `oc_analyses_gbsg.R` | 2061 | Operating characteristics analyses |
| `simulation_tables.R` | 1120 | `summarize_simulation_results()` |
| `synthetic_data_perturbation.R` | ~400 | Bootstrap-based synthetic dataset generation |

### Multi-Regional Clinical Trials (MRCT)
| File | Lines | Description |
|------|-------|-------------|
| `mrct_simulation.R` | 1847 | `mrct_region_sims()`, `create_dgm_for_mrct()` |

### Visualization
| File | Lines | Description |
|------|-------|-------------|
| `plot_subgroup_results_forestplot.R` | 1050 | `plot_subgroup_results_forestplot()` |
| `plot_sg_weighted_km.R` | ~600 | `plot_sg_weighted_km()` |
| `plot_km_band_forestsearch.R` | 825 | `plot_km_band_forestsearch()` |
| `plot_sg_results.R` | 920 | Additional subgroup result plots |
| `plot_sg_distribution.R` | ~200 | Subgroup distribution plots |
| `gg_forest.R` | ~500 | Forest plot utilities |
| `render_forestplot.R` | ~500 | `render_forestplot()` |
| `cox_spline_fit.R` | 780 | `plot_spline_treatment_effect()` |

### Tables & Formatting
| File | Lines | Description |
|------|-------|-------------|
| `create_summary_table.R` | 963 | Publication-ready gt summary tables |
| `format_fs_details.R` | ~300 | Formatting FS result details |
| `format_subgroup_summary_tables.R` | ~450 | Subgroup summary table formatting |
| `summary_utility_functions.R` | ~400 | Shared table/summary utilities |

### Package Infrastructure
| File | Lines | Description |
|------|-------|-------------|
| `forestsearch-package.R` | ~50 | Package-level roxygen docs |
| `globals.R` | ~80 | Global variable declarations |

---

## Key Public Functions

```
forestsearch()                        # Main algorithm entry point
forestsearch_bootstrap_dofuture()     # Parallelized bootstrap bias correction
forestsearch_Kfold()                  # N-fold / LOO cross-validation
forestsearch_tenfold()                # Repeated K-fold CV
subgroup.consistency()                # Split-sample consistency evaluation
summarize_bootstrap_results()         # Bootstrap summary (gt tables)
cv_metrics_tables()                   # CV metrics (sensCV, ppvCV)
plot_subgroup_results_forestplot()    # Forest plot visualization
plot_sg_weighted_km()                 # Weighted KM curves
generate_aft_dgm_flex()              # AFT data-generating mechanism
run_simulation_analysis()             # Full simulation wrapper
```

## Dependencies
- **Core:** survival, data.table, glmnet, grf, doFuture, future, ggplot2, gt
- **Remote:** weightedsurv (larry-leon/weightedsurv)
- **Suggested:** forestploter, DiagrammeR, cubature
