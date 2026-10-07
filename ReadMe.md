# Survival versus Scale: The Effects of Critical Access Hospital Designation

This repo contains the code and supplementary materials for our analysis of the effects of Critical Access Hospital (CAH) designation on hospital closures, mergers, and operational outcomes.

## Data

See [data.md](data.md) for details on data sources and the cleaning/merging pipeline.

**External data repositories** (symlinked into `data/input/`):
- [imccart/aha-data](https://github.com/imccart/aha-data) — American Hospital Association Annual Survey
- [imccart/cah](https://github.com/imccart/cah) — Critical Access Hospital designations (Flex Monitoring Team)
- [imccart/HCRIS](https://github.com/imccart/HCRIS) — Healthcare Cost Report Information System
- [imccart/form-990s](https://github.com/imccart/form-990s) — IRS Form 990 nonprofit financials

A Census API key is required for geographic crosswalks (`data-code/api-keys.R`, not tracked). Request one [here](https://api.census.gov/data/key_signup.html).

## Analysis

See [analysis.md](analysis.md) for details on the estimation pipeline.

The analysis has two entry points:
- **`analysis/_build-estimation-data.r`** — Loads raw data, merges, constructs variables, writes estimation datasets. Run when data inputs change.
- **`analysis/_run-analysis.r`** — Reads estimation data, builds stacked datasets, and sequentially runs all estimation scripts.

Key scripts:
| Script | Description |
|--------|-------------|
| `1-sum-stats.R` | Summary statistics and anticipation figures |
| `2-hospital-dd.R` | Hospital-level SDID/CS with state-timing controls |
| `3-hospital-dd-alt.R` | Hospital-level SDID/CS with eligibility-restricted controls |
| `4-changes-state-dd.R` | State-level closure/merger rates |
| `6-heterogeneity.R` | Subgroup SDID across ownership, isolation, margin, system |
| `7-gsynth.R` | Factor-model estimation (fect/gsynth) |
| `8-access-health.R` | Displacement analysis, net capacity, mortality calibration |

Appendix scripts: `app-dd-diagnostics.R`, `app-financial-preperiod.R`, `app-permutation.R`, `app-statecut-sensitivity.R`, `app-bedcut-sensitivity.R`, `app-anticipation.R`

R packages are managed via `renv`. Run `source("renv/activate.R")` then `library()` calls in `analysis/0-setup.R`.

## Paper

The manuscript is in `paper/` and syncs with [Overleaf](https://www.overleaf.com/project/689f2b10a1bb823dcdec304d).

## Research Team

**Ian McCarthy**
Associate Professor, Department of Economics, Emory University
Research Associate, National Bureau of Economic Research
[ianmccarthyecon.com](https://www.ianmccarthyecon.com)

**Katie Leinenbach**
Senior Quantitative Analyst, Demand Side Analytics
[katieleinenbach.com](https://www.katieleinenbach.com/)
