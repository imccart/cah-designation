# Analysis

## Overview

The analysis has two entry points in the `analysis/` folder:

1. **`_build-estimation-data.r`**: Loads raw data from `data/output/`, merges, constructs all variables (financial, capacity, event indicators), applies winsorization and interpolation, computes IPW weights, and writes `est.dat` and `state.dat` CSVs to `data/output/`. Also sources `1-sum-stats.R` for summary statistics. Run this when underlying data inputs change.

2. **`_run-analysis.r`**: Reads saved CSVs, sets global parameters, builds three stacked datasets, and sequentially sources estimation scripts. Does not rebuild the data.

## Estimation Design

The core approach is **stacked difference-in-differences** to handle heterogeneous treatment timing. Treatment is CAH designation at the hospital level, with state-level program availability as the policy instrument.

Three stacking functions in `functions.R`:
- `stack_hosp()` — Hospital-level, controls from never-treated and not-yet-treated states
- `stack_hosp_elig()` — Hospital-level, controls are all non-CAH or not-yet-CAH hospitals (no state restriction)
- `stack_state()` — State-level, collapses to state-year counts

Primary estimators: Synthetic Difference-in-Differences (SDID) and Callaway-Sant'Anna (CS).

## Scripts

| Script | Level | Description |
|--------|-------|-------------|
| `0-setup.R` | — | Package loading (renv) |
| `1-sum-stats.R` | — | Summary statistics, anticipation figures |
| `2-hospital-dd.R` | Hospital | Main SDID/CS estimates, state-timing controls, cohorts 1999--2001 |
| `3-hospital-dd-alt.R` | Hospital | Eligibility-restricted controls, cohorts 1999--2005 |
| `4-changes-state-dd.R` | State | Closure/merger rates per 100 hospitals |
| `6-heterogeneity.R` | Hospital | Subgroup SDID across 4 dimensions |
| `7-gsynth.R` | Hospital | Factor-model estimation via fect (IFE) |
| `8-access-health.R` | Hospital/County | Displacement, net capacity Monte Carlo, mortality calibration |

### Appendix Scripts

| Script | Description |
|--------|-------------|
| `app-dd-diagnostics.R` | Pre-trends, sample comparison, weight concentration |
| `app-financial-preperiod.R` | Sensitivity to pre-period length for financial outcomes |
| `app-permutation.R` | Permutation inference (500 draws) |
| `app-statecut-sensitivity.R` | Sensitivity to state-level treatment timing cutoff |
| `app-bedcut-sensitivity.R` | Sensitivity to bed-size eligibility threshold |
| `app-anticipation.R` | Anticipation test (treatment shifted to t=-1) |

## Output

- `results/` — PNG figures and LaTeX tables for the manuscript
- `results/diagnostics/` — Pre-trend plots, CSVs, and LaTeX tables for the appendix
