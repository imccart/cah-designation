# app-financial-preperiod.R
# Pre-period sensitivity analysis for financial outcomes
# Runs SDID across pre-period lengths 2-5 for all financial outcomes
# Requires: stack.hosp, bed.cut, est.dat from _run-analysis.r
# Output: results/app-preperiod-financial.tex (LaTeX tabular innards)

financial_outcomes <- list(
  margin        = list(label = "Operating margin",             stub = "margin"),
  current_ratio = list(label = "Current ratio",                stub = "currentratio"),
  net_fixed     = list(label = "Net fixed assets",             stub = "netfixed"),
  capex         = list(label = "Capital expenditures per bed", stub = "capex"),
  net_pat_rev       = list(label = "Net patient revenue per bed",  stub = "netpatrev"),
  tot_operating_exp = list(label = "Operating expenses per bed",   stub = "totopexp")
)

pre_periods <- c(5, 4, 3, 2)
cohorts_fin <- 1999:2001

## SDID estimation for a single cohort × outcome × pre-period
run_sdid_cohort_pp <- function(c, outcome_sym, pp) {
  synth.c <- stack.hosp %>%
    filter(stack_group == c) %>%
    group_by(ID) %>%
    mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(!is.na(!!outcome_sym),
           min_bedsize <= bed.cut,
           stacked_event_time >= -pp) %>%
    select(ID, year, outcome = !!outcome_sym, post_treat)

  n_before <- n_distinct(synth.c$ID)
  bal.c <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
  n_after <- n_distinct(bal.c$ID)
  setup <- panel.matrices(as.data.frame(bal.c))
  Ntr   <- nrow(setup$Y) - setup$N0


  # Need at least 2 pre-periods and 1 treated unit
  if (setup$T0 < 2 || Ntr < 1) {
    return(NULL)
  }

  est  <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
  se_c <- synthdid_se(est, method = "jackknife")

  tibble(
    cohort  = c,
    att     = as.numeric(est),
    se      = as.numeric(se_c),
    Ntr     = Ntr,
    N_total = n_after
  )
}

## Collect results: outcome × pre-period
results <- tibble(
  outcome = character(), label = character(),
  pre_period = integer(),
  att = numeric(), se = numeric(),
  ci_low = numeric(), ci_high = numeric(),
  N_total = integer(), N_treated = integer()
)

for (oname in names(financial_outcomes)) {
  o <- financial_outcomes[[oname]]
  osym <- sym(oname)

  for (pp in pre_periods) {

    atts_all <- bind_rows(map(cohorts_fin, function(c) {
      tryCatch(
        run_sdid_cohort_pp(c, osym, pp),
        error = function(e) {
          message(sprintf("financial pre-period: cohort %d failed (%s)", c, conditionMessage(e)))
          NULL
        }
      )
    }))

    if (nrow(atts_all) == 0) {
      results <- bind_rows(results, tibble(
        outcome = oname, label = o$label, pre_period = as.integer(pp),
        att = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        N_total = NA_integer_, N_treated = NA_integer_
      ))
      next
    }

    att_w   <- with(atts_all, sum(Ntr * att) / sum(Ntr))
    se_w    <- with(atts_all, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))
    ci_low  <- att_w - 1.96 * se_w
    ci_high <- att_w + 1.96 * se_w
    n_tot   <- sum(atts_all$N_total)
    n_tr    <- sum(atts_all$Ntr)


    results <- bind_rows(results, tibble(
      outcome = oname, label = o$label, pre_period = as.integer(pp),
      att = att_w, se = se_w, ci_low = ci_low, ci_high = ci_high,
      N_total = as.integer(n_tot), N_treated = as.integer(n_tr)
    ))
  }
}

## Format and write complete LaTeX tabular block
## (LaTeX 2025 breaks \noalign/\omit inside \input within tabular)
# Columns: Outcome | Pre=5 ATT [CI] | Pre=4 | Pre=3 | Pre=2

format_cell <- function(att, ci_low, ci_high) {
  if (is.na(att)) return("--- & ---")
  ax <- abs(att); acx <- max(abs(ci_low), abs(ci_high))
  d_att <- ifelse(ax < 5, 3, ifelse(ax < 100, 2, 1))
  d_ci  <- ifelse(acx < 5, 2, ifelse(acx < 100, 1, 0))
  sprintf("%s & [%s, %s]",
    formatC(att, format="f", digits=d_att),
    formatC(ci_low, format="f", digits=d_ci),
    formatC(ci_high, format="f", digits=d_ci))
}

tex_lines <- c(
  "\\begin{tabular}{l cc cc cc cc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Pre $= 5$} & \\multicolumn{2}{c}{Pre $= 4$} & \\multicolumn{2}{c}{Pre $= 3$} & \\multicolumn{2}{c}{Pre $= 2$} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}",
  "Outcome & ATT & 95\\% CI & ATT & 95\\% CI & ATT & 95\\% CI & ATT & 95\\% CI \\\\",
  "\\midrule"
)

for (oname in names(financial_outcomes)) {
  o <- financial_outcomes[[oname]]
  row_data <- results %>% filter(outcome == oname)

  cells <- sapply(pre_periods, function(pp) {
    r <- row_data %>% filter(pre_period == pp)
    if (nrow(r) == 0 || is.na(r$att)) return("--- & ---")
    format_cell(r$att, r$ci_low, r$ci_high)
  })

  tex_lines <- c(tex_lines,
    paste0(o$label, " & ", paste(cells, collapse = " & "), " \\\\"))
}

# Add observation counts at the bottom (margin only, as representative)
tex_lines <- c(tex_lines, "\\midrule")
margin_data <- results %>% filter(outcome == "margin")

ctrl_cells <- sapply(pre_periods, function(pp) {
  r <- margin_data %>% filter(pre_period == pp)
  if (nrow(r) == 0 || is.na(r$N_total)) return("\\multicolumn{2}{c}{---}")
  sprintf("\\multicolumn{2}{c}{%d}", r$N_total - r$N_treated)
})
tex_lines <- c(tex_lines,
  paste0("Control obs. & ", paste(ctrl_cells, collapse = " & "), " \\\\"))

treat_cells <- sapply(pre_periods, function(pp) {
  r <- margin_data %>% filter(pre_period == pp)
  if (nrow(r) == 0 || is.na(r$N_treated)) return("\\multicolumn{2}{c}{---}")
  sprintf("\\multicolumn{2}{c}{%d}", r$N_treated)
})
tex_lines <- c(tex_lines,
  paste0("Treated obs. & ", paste(treat_cells, collapse = " & "), " \\\\"))

tex_lines <- c(tex_lines, "\\bottomrule", "\\end{tabular}")

writeLines(tex_lines, "results/app-preperiod-financial.tex")

## Also save CSV for reference
write_csv(results, "results/diagnostics/app-preperiod-financial.csv")
