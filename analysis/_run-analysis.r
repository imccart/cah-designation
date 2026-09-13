# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   9/12/2026
## Description:   Run Analysis Files
## Note:          Run _build-estimation-data.r first to generate the CSVs


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, haven, readxl, janitor, here, zoo, fedmatch, zipcodeR,
               fixest, did, did2s, BMisc, fect, glmnet, nnet, mlogit, survival,
               scales, plotly, panelView, dotwhisker, patchwork, sf,
               modelsummary, kableExtra, broom, synthdid)

source('analysis/functions.R')


# Read estimation data -----------------------------------------------------
est.dat   <- read_csv('data/output/estimation_data.csv', col_types = cols(ID = col_character()))
state.dat <- read_csv('data/output/state_estimation_data.csv', show_col_types = FALSE)
fin.vars  <- c('margin', 'current_ratio', 'net_fixed', 'capex', 'net_pat_rev', 'tot_operating_exp')


# Global parameters --------------------------------------------------------
set.seed(1234)       # the CS bootstrap is random; seed so runs reproduce
bed.cut   <- 50
post      <- 5
state.cut <- 0
financial.pre <- 3   # HCRIS financials begin in 1996-97; 3 pre-years gives every cohort its full window


# Stacked datasets ---------------------------------------------------------
stack.hosp  <- stack_hosp(pre.period=5, post.period=post, state.period=state.cut)
stack.state <- stack_state(pre.period=5, post.period=post, state.period=state.cut)
stack.elig  <- stack_hosp_elig(pre.period=5, post.period=post, cohort.years=1999:2005)


# Main analysis ------------------------------------------------------------
source('analysis/2-hospital-dd.R')        # -> hosp.results.table, hosp.cohort.results
source('analysis/3-hospital-dd-alt.R')     # -> elig.results, elig.cohort.results
source('analysis/4-changes-state-dd.R')    # -> state.results.table


# Combined results table (hospital + state) --------------------------------
# Reorder: financial (1-6), capacity (7-10), organizational (closures, mergers, system)
results.table <- bind_rows(hosp.results.table, state.results.table) %>%
  mutate(.order = case_when(
    outcome == "Operating margin"             ~ 1,
    outcome == "Current ratio"                ~ 2,
    outcome == "Net fixed assets"             ~ 3,
    outcome == "Capital expenditures per bed" ~ 4,
    outcome == "Net patient revenue per bed"  ~ 5,
    outcome == "Operating expenses per bed"   ~ 6,
    outcome == "Total beds"                   ~ 7,
    outcome == "OB beds"                      ~ 8,
    outcome == "FTE RNs"                      ~ 9,
    outcome == "Inpatient days per bed"       ~ 10,
    outcome == "Closures"                     ~ 11,
    outcome == "Mergers"                      ~ 12,
    outcome == "System membership"            ~ 13,
    TRUE                                      ~ 99
  )) %>%
  arrange(.order) %>%
  select(-.order)

fmt <- function(x) {
  ax <- abs(x)
  ifelse(ax < 5, sprintf("%.3f", x),
  ifelse(ax < 100, sprintf("%.2f", x),
         sprintf("%.1f", x)))
}
int <- function(lo, hi) {
  alo <- abs(lo); ahi <- abs(hi); ax <- max(alo, ahi)
  d <- ifelse(ax < 5, 2, ifelse(ax < 100, 1, 0))
  sprintf("[%s, %s]", formatC(lo, format="f", digits=d), formatC(hi, format="f", digits=d))
}

## Main table (SDID only) for paper — complete tabular block
## (LaTeX 2025 breaks \noalign/\omit inside \input within tabular)
tex.lines.sdid <- results.table %>%
  rowwise() %>%
  mutate(line = sprintf("%s & %s & %s & %d \\\\",
    outcome, fmt(sdid_att), int(sdid_ci_low, sdid_ci_high), sdid_ntr)) %>%
  ungroup() %>%
  pull(line)

# Insert group separators: financial (1-6), capacity (7-10), organizational (11-13)
tex.lines.sdid <- append(tex.lines.sdid, "\\addlinespace", after = 6)
tex.lines.sdid <- append(tex.lines.sdid, "\\addlinespace", after = 11)  # shifted by 1

writeLines(c(
  "\\begin{tabular}[t]{lccr}",
  "\\toprule",
  "Outcome & ATT & 95\\% CI & $N_{tr}$ \\\\",
  "\\midrule",
  tex.lines.sdid,
  "\\bottomrule",
  "\\end{tabular}"
), "results/att_overall.tex")

## CS table for appendix — complete tabular block
tex.lines.cs <- results.table %>%
  rowwise() %>%
  mutate(line = sprintf("%s & %s & %s \\\\",
    outcome, fmt(cs_att), int(cs_ci_low, cs_ci_high))) %>%
  ungroup() %>%
  pull(line)

# Insert group separators: financial (1-6), capacity (7-10), organizational (11-13)
tex.lines.cs <- append(tex.lines.cs, "\\addlinespace", after = 6)
tex.lines.cs <- append(tex.lines.cs, "\\addlinespace", after = 11)  # shifted by 1

writeLines(c(
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Outcome & CS ATT & CS 95\\% CI \\\\",
  "\\midrule",
  tex.lines.cs,
  "\\bottomrule",
  "\\end{tabular}"
), "results/att_cs_overall.tex")


# Forest plot: SDID results across control group constructions ---------------

forest_labels <- c(
  "Operating margin" = "Operating\nmargin",
  "Current ratio" = "Current ratio",
  "Net fixed assets" = "Net fixed assets",
  "Capital expenditures per bed" = "Capital expenditures\nper bed",
  "Net patient revenue per bed" = "Net patient revenue\nper bed",
  "Operating expenses per bed" = "Operating expenses\nper bed",
  "Total beds" = "Total beds",
  "OB beds" = "OB beds",
  "FTE RNs" = "FTE RNs",
  "Inpatient days per bed" = "Inpatient days\nper bed",
  "Closures" = "Closures",
  "Mergers" = "Mergers",
  "System membership" = "System\nmembership"
)

forest_dat <- bind_rows(
  results.table %>%
    transmute(outcome, controls = "State-timing",
              att = sdid_att, ci_low = sdid_ci_low, ci_high = sdid_ci_high),
  elig.results %>%
    transmute(outcome, controls = "Eligibility-restricted",
              att = sdid_att, ci_low = sdid_ci_low, ci_high = sdid_ci_high)
) %>%
  mutate(outcome_label = forest_labels[outcome],
         outcome_label = factor(outcome_label,
                                levels = forest_labels[unique(results.table$outcome)]),
         controls = factor(controls, levels = c("Eligibility-restricted", "State-timing")))

n_out <- length(unique(forest_dat$outcome_label))
p_forest <- ggplot(forest_dat, aes(x = att, y = controls, shape = controls)) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "gray40") +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high),
                  size = 0.4, linewidth = 0.4) +
  scale_shape_manual(values = c("State-timing" = 16, "Eligibility-restricted" = 1)) +
  facet_wrap(~ outcome_label, ncol = 1, scales = "free_x", strip.position = "left",
             dir = "v") +
  labs(x = "ATT (95% CI)", y = NULL, shape = "Controls") +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9,
                                     lineheight = 0.9),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )
ggsave("results/forest-sdid.png", p_forest,
       width = 7, height = 0.9 * n_out + 1.2, dpi = 300)


# Heterogeneity ------------------------------------------------------------
source('analysis/6-heterogeneity.R')         # -> het.results

# Factor-model estimation (fect) ------------------------------------------
source('analysis/7-gsynth.R')                # -> gsynth.results

# Access and health implications --------------------------------------------
source('analysis/8-access-health.R')          # -> map, histogram, PSA, county CSVs

# Diagnostics and sensitivity ----------------------------------------------
source('analysis/app-dd-diagnostics.R')
source('analysis/app-financial-preperiod.R')
source('analysis/app-permutation.R')         # -> perm.results, permutation-sdid.png
source('analysis/app-statecut-sensitivity.R') # -> statecut-sensitivity.png
source('analysis/app-bedcut-sensitivity.R')   # -> bedcut-sensitivity.png
source('analysis/app-anticipation.R')         # -> att_anticipation.tex
source('analysis/app-control-timing.R')       # -> control-timing-*.png

# --- R&R additional analyses (JPubE R1, JPUBE-D-26-00838) -------------------
source('analysis/R1_revision.R')

