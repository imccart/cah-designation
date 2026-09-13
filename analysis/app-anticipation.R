# app-anticipation.R
# Robustness to anticipatory downsizing: redefine treatment as starting at t-1
# Compares baseline SDID (treatment at t=0) with shifted SDID (treatment at t=-1)
# Uses stack.elig (eligibility-restricted design, cohorts 1999-2005)
# Expects from _run-analysis.r: stack.elig, bed.cut, financial.pre, elig.results


# Outcome map (matching 3-hospital-dd-alt.R) --------------------------------
antic_outcome_map <- list(
  margin            = list(label="Operating margin",             pre_period=financial.pre),
  current_ratio     = list(label="Current ratio",                pre_period=financial.pre),
  net_fixed         = list(label="Net fixed assets",             pre_period=financial.pre),
  capex             = list(label="Capital expenditures per bed", pre_period=financial.pre),
  net_pat_rev       = list(label="Net patient revenue per bed",  pre_period=financial.pre),
  tot_operating_exp = list(label="Operating expenses per bed",   pre_period=financial.pre),
  BDTOT             = list(label="Total beds"),
  OBBD              = list(label="OB beds"),
  FTERN             = list(label="FTE RNs"),
  ip_per_bed        = list(label="Inpatient days per bed"),
  system            = list(label="System membership")
)

antic_cohorts <- 1999:2005


# SDID helper with shifted treatment ----------------------------------------
run_sdid_antic <- function(c, outcome_sym, pp) {
  synth.c <- stack.elig %>%
    filter(stack_group == c) %>%
    group_by(ID) %>%
    mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
    ungroup() %>%
    # Shift treatment one year earlier: post_treat = 1 at stacked_event_time >= -1
    mutate(post_treat = ifelse(treated == 1 & stacked_event_time >= -1, 1, post_treat)) %>%
    filter(!is.na(!!outcome_sym),
           min_bedsize <= bed.cut,
           stacked_event_time >= -pp) %>%
    select(ID, year, outcome = !!outcome_sym, post_treat, treated)

  bal.c <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))

  n_yr <- length(unique(bal.c$year))
  n_tr <- sum(bal.c$treated == 1) / max(n_yr, 1)
  n_co <- sum(bal.c$treated == 0) / max(n_yr, 1)
  if (n_tr < 2 || n_co < 2) return(NULL)

  setup <- panel.matrices(as.data.frame(bal.c))
  Ntr   <- nrow(setup$Y) - setup$N0

  if (setup$T0 < 2 || Ntr < 1) return(NULL)

  est  <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
  se_c <- synthdid_se(est, method = "jackknife")

  tibble(cohort = c, att = as.numeric(est), se = as.numeric(se_c), Ntr = Ntr)
}


# Main loop -----------------------------------------------------------------
antic.results <- tibble(
  outcome = character(), sdid_att = numeric(),
  sdid_ci_low = numeric(), sdid_ci_high = numeric(),
  sdid_ntr = numeric()
)

for (oname in names(antic_outcome_map)) {
  o    <- antic_outcome_map[[oname]]
  osym <- sym(oname)
  pp   <- if (!is.null(o$pre_period)) o$pre_period else 5


  atts_all <- bind_rows(map(antic_cohorts, function(c) {
    tryCatch(
      run_sdid_antic(c, osym, pp),
      error = function(e) {
        message(sprintf("anticipation: cohort %d failed", c))
        NULL
      }
    )
  }))

  if (nrow(atts_all) == 0) {
    antic.results <- bind_rows(antic.results, tibble(
      outcome = o$label, sdid_att = NA_real_,
      sdid_ci_low = NA_real_, sdid_ci_high = NA_real_,
      sdid_ntr = NA_integer_
    ))
    next
  }

  att_w  <- with(atts_all, sum(Ntr * att) / sum(Ntr))
  se_w   <- with(atts_all, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))
  ntr_w  <- sum(atts_all$Ntr)


  antic.results <- bind_rows(antic.results, tibble(
    outcome     = o$label,
    sdid_att    = att_w,
    sdid_ci_low = att_w - 1.96 * se_w,
    sdid_ci_high= att_w + 1.96 * se_w,
    sdid_ntr    = ntr_w
  ))
}


# Comparison table -----------------------------------------------------------
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

# Join baseline (from elig.results) with anticipation-shifted
antic_table <- elig.results %>%
  transmute(outcome,
            base_att = sdid_att, base_ci_low = sdid_ci_low, base_ci_high = sdid_ci_high) %>%
  left_join(
    antic.results %>%
      transmute(outcome,
                antic_att = sdid_att, antic_ci_low = sdid_ci_low, antic_ci_high = sdid_ci_high),
    by = "outcome"
  )

tex_lines <- antic_table %>%
  rowwise() %>%
  mutate(line = sprintf("%s & %s & %s & %s & %s \\\\",
    outcome,
    fmt(base_att), int(base_ci_low, base_ci_high),
    fmt(antic_att), int(antic_ci_low, antic_ci_high))) %>%
  ungroup() %>%
  pull(line)

# Group separators: financial (1-6), capacity (7-10), organizational (11)
tex_lines <- append(tex_lines, "\\addlinespace", after = 6)
tex_lines <- append(tex_lines, "\\addlinespace", after = 11)  # shifted by 1

writeLines(c(
  "\\begin{tabular}[t]{lcccc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Baseline ($t = 0$)} & \\multicolumn{2}{c}{Anticipation ($t = -1$)} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "Outcome & ATT & 95\\% CI & ATT & 95\\% CI \\\\",
  "\\midrule",
  tex_lines,
  "\\bottomrule",
  "\\end{tabular}"
), "results/att_anticipation.tex")


# Save CSV ------------------------------------------------------------------
write_csv(antic_table, "results/diagnostics/anticipation-comparison.csv")

