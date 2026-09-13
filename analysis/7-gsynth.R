# Meta --------------------------------------------------------------------
## Factor-model estimation (fect) for financial and operational outcomes
## fect handles unbalanced panels natively via imputation, accommodating
## pre-trend heterogeneity and small balanced-panel samples
## Uses eligibility-restricted design (stack.elig), cohorts 1999-2005

# Expects from _run-analysis.r: stack.elig, est.dat, bed.cut, post

gsynth.cohorts <- 1999:2005

# Outcome map ------------------------------------------------------------------
gsynth_outcome_map <- list(
  ## Financial: fect drops units with fewer than 5 untreated periods, so the
  ## 5-year pre-period is kept here even though the SDID financial estimates
  ## use financial.pre (3). fect does not require a balanced panel, so the
  ## 2001-2005 cohorts retain their treated units on the HCRIS series.
  margin            = list(label = "Operating margin",             pre_period = 5),
  current_ratio     = list(label = "Current ratio",                pre_period = 5),
  net_fixed         = list(label = "Net fixed assets",             pre_period = 5),
  capex             = list(label = "Capital expenditures per bed", pre_period = 5),
  net_pat_rev       = list(label = "Net patient revenue per bed",  pre_period = 5),
  tot_operating_exp = list(label = "Operating expenses per bed",   pre_period = 5),
  ## Operational (5-year pre-period)
  BDTOT             = list(label = "Total beds",                   pre_period = 5),
  OBBD              = list(label = "OB beds",                      pre_period = 5),
  FTERN             = list(label = "FTE RNs",                      pre_period = 5),
  ip_per_bed        = list(label = "Inpatient days per bed",       pre_period = 5),
  system            = list(label = "System membership",            pre_period = 5)
)

# Results collector ------------------------------------------------------------
gsynth.results <- tibble(
  outcome = character(), att = numeric(), ci_low = numeric(), ci_high = numeric(),
  ntr = integer(), n_cohorts = integer()
)

gsynth.cohort.results <- tibble(
  outcome = character(), cohort = numeric(), att = numeric(),
  se = numeric(), Ntr = numeric()
)

# Main loop --------------------------------------------------------------------
for (oname in names(gsynth_outcome_map)) {
  o <- gsynth_outcome_map[[oname]]
  outcome_sym   <- sym(oname)
  outcome_label <- o$label
  pp            <- if (!is.null(o$pre_period)) o$pre_period else 5


  run_fect_cohort <- function(c) {
    fect.dat <- stack.elig %>%
      filter(stack_group == c) %>%
      group_by(ID) %>%
      mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(!is.na(!!outcome_sym), min_bedsize <= bed.cut,
             stacked_event_time >= -pp) %>%
      mutate(ID_num = as.numeric(factor(ID))) %>%
      select(ID_num, year, outcome = !!outcome_sym, post_treat, treated)

    ## Need at least 3 treated and 3 control units (counted on the unit
    ## indicator, not on post_treat, whose zeros include treated units' pre-years)
    n_tr <- fect.dat %>% filter(treated == 1) %>% distinct(ID_num) %>% nrow()
    n_co <- fect.dat %>% filter(treated == 0) %>% distinct(ID_num) %>% nrow()
    if (n_tr < 3 || n_co < 3) {
      message(sprintf("gsynth: %s, cohort %d dropped (n_tr=%d, n_co=%d, need >=3)", oname, c, n_tr, n_co))
      return(NULL)
    }

    out <- fect(
      data      = as.data.frame(fect.dat),
      Y         = "outcome",
      D         = "post_treat",
      index     = c("ID_num", "year"),
      method    = "ife",
      force     = "two-way",
      CV        = TRUE,
      r         = c(0, 5),
      se        = TRUE,
      vartype   = "bootstrap",
      nboots    = 150,
      parallel  = FALSE,
      seed      = 42
    )

    att_val <- out$att.avg
    se_val  <- sd(out$att.avg.boot)

    tibble(cohort = c, att = att_val, se = se_val, Ntr = n_tr)
  }

  out_all <- map(gsynth.cohorts, function(c) {
    tryCatch(run_fect_cohort(c),
             error = function(e) {
               message(sprintf("gsynth: %s, cohort %d dropped (%s)", oname, c, conditionMessage(e)))
               NULL })
  })
  out_all <- compact(out_all)
  atts_all <- bind_rows(out_all)

  if (nrow(atts_all) == 0) {
    message(sprintf("gsynth: %s skipped, fect failed for every cohort", oname))
    next
  }
  ## Transparency: report which cohorts entered the pooled estimate

  ## Collect cohort-specific results
  gsynth.cohort.results <- bind_rows(gsynth.cohort.results,
    atts_all %>% mutate(outcome = outcome_label))

  ## Aggregate across cohorts (weighted by Ntr)
  att_w  <- with(atts_all, sum(Ntr * att) / sum(Ntr))
  se_w   <- with(atts_all, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))
  ci_low <- att_w - 1.96 * se_w
  ci_high <- att_w + 1.96 * se_w

  gsynth.results <- bind_rows(gsynth.results, tibble(
    outcome   = outcome_label,
    att       = att_w,
    ci_low    = ci_low,
    ci_high   = ci_high,
    ntr       = as.integer(sum(atts_all$Ntr)),
    n_cohorts = as.integer(nrow(atts_all))
  ))


  rm(out_all, atts_all); gc()
}


# LaTeX table -----------------------------------------------------------------
fmt <- function(x) {
  ax <- abs(x)
  ifelse(ax < 5, sprintf("%.3f", x),
  ifelse(ax < 100, sprintf("%.2f", x),
         sprintf("%.1f", x)))
}
int <- function(lo, hi) {
  ax <- max(abs(lo), abs(hi))
  d <- ifelse(ax < 5, 2, ifelse(ax < 100, 1, 0))
  sprintf("[%s, %s]", formatC(lo, format="f", digits=d), formatC(hi, format="f", digits=d))
}

tex_lines <- gsynth.results %>%
  rowwise() %>%
  mutate(line = sprintf("%s & %s & %s & %d \\\\",
    outcome, fmt(att), int(ci_low, ci_high), ntr)) %>%
  ungroup() %>%
  pull(line)

# Insert group separators: financial (1-6), capacity (7-10), organizational (11)
tex_lines <- append(tex_lines, "\\addlinespace", after = 6)
tex_lines <- append(tex_lines, "\\addlinespace", after = 11)  # shifted by 1

writeLines(c(
  "\\begin{tabular}{lccr}",
  "\\toprule",
  "Outcome & IFE ATT & 95\\% CI & $N_{tr}$ \\\\",
  "\\midrule",
  tex_lines,
  "\\bottomrule",
  "\\end{tabular}"
), "results/att_gsynth.tex")

## Save cohort-level CSV for diagnostics
write_csv(gsynth.cohort.results, "results/diagnostics/gsynth_cohort_results.csv")

