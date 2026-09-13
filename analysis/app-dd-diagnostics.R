# Diagnostics for the SDID estimates
#   1. Sample size before and after the balanced-panel requirement
#   2. Pre-trend tests
#   3. Synthetic weight concentration
#
# Covers hospital-level outcomes (margin, current_ratio, etc.) and
# state-level outcomes (closures, mergers).
# Outputs CSVs, pre-trend PNGs, and LaTeX tabular innards.

cohorts <- 1999:2001
hosp_outcomes <- c("margin", "current_ratio", "net_fixed", "capex",
                   "net_pat_rev", "tot_operating_exp",
                   "BDTOT", "OBBD", "FTERN", "ip_per_bed", "system")
state_outcomes <- c("closures", "mergers")

outcome_labels <- c(
  margin = "Operating margin", current_ratio = "Current ratio",
  net_fixed = "Net fixed assets", capex = "Capital expenditures",
  net_pat_rev = "Net patient revenue per bed", tot_operating_exp = "Operating expenses per bed",
  BDTOT = "Total beds", OBBD = "OB beds", FTERN = "FTE RNs",
  ip_per_bed = "IP days per bed",
  system = "System membership", closures = "Closures", mergers = "Mergers"
)

# =============================================================================
# DIAGNOSTIC 1: Sample Size Comparison
# =============================================================================

sample_comparison <- tibble(
  outcome = character(),
  cohort = integer(),
  sdid_n_treated = integer(),
  sdid_n_control = integer(),
  sdid_n_total = integer(),
  cs_n_treated = integer(),
  cs_n_control = integer(),
  cs_n_total = integer(),
  pct_lost_to_balance = numeric()
)

# --- Hospital-level outcomes ---
for (outcome_var in hosp_outcomes) {
  outcome_sym <- sym(outcome_var)
  ## Financial outcomes use the same pre-period as the main SDID estimates
  pre_o <- if (outcome_var %in% fin.vars) financial.pre else 5

  for (c in cohorts) {
    # SDID sample (balanced panel)
    synth.c <- stack.hosp %>%
      filter(stack_group == c) %>%
      group_by(ID) %>%
      mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(!is.na(!!outcome_sym), min_bedsize <= bed.cut,
             stacked_event_time >= -pre_o) %>%
      select(ID, year, outcome = !!outcome_sym, treated = post_treat)

    # Before balancing. A unit is treated if it is ever post-treated; counting on
    # the row-level indicator would put treated units in both groups via their
    # own pre-period rows.
    pre_balance <- synth.c %>%
      group_by(ID) %>%
      summarize(ever_treated = max(treated), .groups = "drop")

    # After balancing
    bal.c <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
    post_balance <- bal.c %>%
      group_by(ID) %>%
      summarize(ever_treated = max(treated), .groups = "drop")

    sdid_treated <- sum(post_balance$ever_treated == 1)
    sdid_control <- sum(post_balance$ever_treated == 0)

    pre_treated <- sum(pre_balance$ever_treated == 1)
    pre_control <- sum(pre_balance$ever_treated == 0)

    # CS sample (unbalanced allowed, notyettreated control group)
    # Matches 2-hospital-dd.R: treat_group in c(0, 1999:2001), control_group="notyettreated"
    cs.c <- est.dat %>%
      group_by(ID) %>%
      mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(min_bedsize <= bed.cut) %>%
      mutate(
        y = !!outcome_sym,
        treat_group = case_when(
          !is.na(eff_year) ~ eff_year,
          is.na(eff_year) & state_treat_year > year + state.cut ~ 0,
          is.na(eff_year) & state_treat_year == 0 ~ 0,
          TRUE ~ NA
        )
      ) %>%
      filter(!is.na(y), !is.na(year), !is.na(BDTOT), !is.na(distance),
             treat_group %in% c(0, cohorts))

    cs_treated <- cs.c %>% filter(treat_group == c) %>% summarize(n = n_distinct(ID)) %>% pull(n)
    cs_control <- cs.c %>% filter(treat_group == 0 | treat_group > c) %>% summarize(n = n_distinct(ID)) %>% pull(n)

    # Calculate attrition
    pre_bal_total <- pre_treated + pre_control
    post_bal_total <- sdid_treated + sdid_control
    pct_lost <- ifelse(pre_bal_total > 0, 100 * (1 - post_bal_total / pre_bal_total), NA)

    sample_comparison <- bind_rows(sample_comparison, tibble(
      outcome = outcome_var,
      cohort = c,
      pre_n_treated = pre_treated,
      pre_n_control = pre_control,
      sdid_n_treated = sdid_treated,
      sdid_n_control = sdid_control,
      sdid_n_total = sdid_treated + sdid_control,
      cs_n_treated = cs_treated,
      cs_n_control = cs_control,
      cs_n_total = cs_treated + cs_control,
      pct_lost_to_balance = round(pct_lost, 1)
    ))
  }
}

# --- State-level outcomes ---
for (outcome_var in state_outcomes) {

  for (c in cohorts) {
    # SDID sample (balanced panel of states)
    synth.c <- stack.state %>%
      filter(stack_group == c) %>%
      transmute(ID = as.numeric(factor(MSTATE)),
                year,
                outcome = .data[[outcome_var]],
                treated = post_treat)

    pre_balance <- synth.c %>%
      group_by(ID) %>%
      summarize(ever_treated = max(treated), .groups = "drop")

    bal.c <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
    post_balance <- bal.c %>%
      group_by(ID) %>%
      summarize(ever_treated = max(treated), .groups = "drop")

    sdid_treated <- sum(post_balance$ever_treated == 1)
    sdid_control <- sum(post_balance$ever_treated == 0)

    pre_treated <- sum(pre_balance$ever_treated == 1)
    pre_control <- sum(pre_balance$ever_treated == 0)

    # CS sample: state.dat with notyettreated control group
    # CS uses control_group="notyettreated" on data filtered to state_treat_year %in% c(0, cohorts).
    # For cohort c, controls = never-treated (state_treat_year==0) + later cohorts.
    cs.c <- state.dat %>%
      filter(state_treat_year %in% c(0, cohorts))

    cs_treated <- cs.c %>% filter(state_treat_year == c) %>% summarize(n = n_distinct(state)) %>% pull(n)
    cs_control <- cs.c %>% filter(state_treat_year == 0 | state_treat_year > c) %>% summarize(n = n_distinct(state)) %>% pull(n)

    pre_bal_total <- pre_treated + pre_control
    post_bal_total <- sdid_treated + sdid_control
    pct_lost <- ifelse(pre_bal_total > 0, 100 * (1 - post_bal_total / pre_bal_total), NA)

    sample_comparison <- bind_rows(sample_comparison, tibble(
      outcome = outcome_var,
      cohort = c,
      pre_n_treated = pre_treated,
      pre_n_control = pre_control,
      sdid_n_treated = sdid_treated,
      sdid_n_control = sdid_control,
      sdid_n_total = sdid_treated + sdid_control,
      cs_n_treated = cs_treated,
      cs_n_control = cs_control,
      cs_n_total = cs_treated + cs_control,
      pct_lost_to_balance = round(pct_lost, 1)
    ))
  }
}

sample_summary <- sample_comparison %>%
  group_by(outcome) %>%
  summarize(
    sdid_treated_total = sum(sdid_n_treated),
    sdid_control_total = sum(sdid_n_control),
    cs_treated_total = sum(cs_n_treated),
    cs_control_total = sum(cs_n_control),
    avg_pct_lost = mean(pct_lost_to_balance, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================================
# DIAGNOSTIC 2: Pre-Trend Analysis
# =============================================================================

pretrend_results <- tibble(
  outcome = character(),
  cohort = integer(),
  trend_coef = numeric(),
  trend_se = numeric(),
  trend_pval = numeric(),
  has_pretrend = logical()
)

# --- Hospital-level outcomes ---
for (outcome_var in hosp_outcomes) {
  outcome_sym <- sym(outcome_var)
  ## Financial outcomes use the same pre-period as the main SDID estimates
  pre_o <- if (outcome_var %in% fin.vars) financial.pre else 5
  all_cohort_data <- tibble()

  for (c in cohorts) {
    cohort_dat <- stack.hosp %>%
      filter(stack_group == c) %>%
      group_by(ID) %>%
      mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(!is.na(!!outcome_sym), min_bedsize <= bed.cut,
             stacked_event_time >= -pre_o) %>%
      mutate(outcome = !!outcome_sym)

    pre_dat <- cohort_dat %>%
      filter(stacked_event_time < 0)

    if (nrow(pre_dat) > 100) {
      trend_reg <- tryCatch({
        feols(outcome ~ stacked_event_time * treated, data = pre_dat, cluster = ~ID)
      }, error = function(e) NULL)

      if (!is.null(trend_reg)) {
        coefs <- as.data.frame(coeftable(trend_reg))
        interaction_row <- grep("stacked_event_time:treated", rownames(coefs))

        if (length(interaction_row) > 0) {
          trend_coef <- coefs[interaction_row, "Estimate"]
          trend_se <- coefs[interaction_row, "Std. Error"]
          trend_pval <- coefs[interaction_row, "Pr(>|t|)"]

          pretrend_results <- bind_rows(pretrend_results, tibble(
            outcome = outcome_var,
            cohort = c,
            trend_coef = trend_coef,
            trend_se = trend_se,
            trend_pval = trend_pval,
            has_pretrend = trend_pval < 0.10
          ))
        }
      }
    }

    means_dat <- cohort_dat %>%
      group_by(stacked_event_time, treated) %>%
      summarize(mean_outcome = mean(outcome, na.rm = TRUE),
                se_outcome = sd(outcome, na.rm = TRUE) / sqrt(n()),
                n = n(),
                .groups = "drop") %>%
      mutate(cohort = c)

    all_cohort_data <- bind_rows(all_cohort_data, means_dat)
  }

  agg_means <- all_cohort_data %>%
    group_by(stacked_event_time, treated) %>%
    summarize(mean_outcome = mean(mean_outcome, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(group = ifelse(treated == 1, "Treated", "Control"))

  p <- ggplot(agg_means, aes(x = stacked_event_time, y = mean_outcome,
                              color = group, linetype = group)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Control" = "gray50", "Treated" = "black")) +
    scale_linetype_manual(values = c("Control" = "dashed", "Treated" = "solid")) +
    labs(
      x = "Event Time",
      y = "Mean Outcome",
      color = NULL,
      linetype = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(
    sprintf("results/diagnostics/diag-pretrend-%s.png", tolower(outcome_var)),
    p, width = 6, height = 4, dpi = 300
  )
}

# --- State-level outcomes ---
for (outcome_var in state_outcomes) {
  all_cohort_data <- tibble()

  for (c in cohorts) {
    cohort_dat <- stack.state %>%
      filter(stack_group == c) %>%
      mutate(outcome = .data[[outcome_var]],
             ID = as.numeric(factor(MSTATE)))

    pre_dat <- cohort_dat %>%
      filter(stacked_event_time < 0)

    if (nrow(pre_dat) > 10) {
      trend_reg <- tryCatch({
        feols(outcome ~ stacked_event_time * treated, data = pre_dat, cluster = ~MSTATE)
      }, error = function(e) NULL)

      if (!is.null(trend_reg)) {
        coefs <- as.data.frame(coeftable(trend_reg))
        interaction_row <- grep("stacked_event_time:treated", rownames(coefs))

        if (length(interaction_row) > 0) {
          trend_coef <- coefs[interaction_row, "Estimate"]
          trend_se <- coefs[interaction_row, "Std. Error"]
          trend_pval <- coefs[interaction_row, "Pr(>|t|)"]

          pretrend_results <- bind_rows(pretrend_results, tibble(
            outcome = outcome_var,
            cohort = c,
            trend_coef = trend_coef,
            trend_se = trend_se,
            trend_pval = trend_pval,
            has_pretrend = trend_pval < 0.10
          ))
        }
      }
    }

    means_dat <- cohort_dat %>%
      group_by(stacked_event_time, treated) %>%
      summarize(mean_outcome = mean(outcome, na.rm = TRUE),
                se_outcome = sd(outcome, na.rm = TRUE) / sqrt(n()),
                n = n(),
                .groups = "drop") %>%
      mutate(cohort = c)

    all_cohort_data <- bind_rows(all_cohort_data, means_dat)
  }

  agg_means <- all_cohort_data %>%
    group_by(stacked_event_time, treated) %>%
    summarize(mean_outcome = mean(mean_outcome, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(group = ifelse(treated == 1, "Treated", "Control"))

  p <- ggplot(agg_means, aes(x = stacked_event_time, y = mean_outcome,
                              color = group, linetype = group)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Control" = "gray50", "Treated" = "black")) +
    scale_linetype_manual(values = c("Control" = "dashed", "Treated" = "solid")) +
    labs(
      x = "Event Time",
      y = "Mean Count",
      color = NULL,
      linetype = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(
    sprintf("results/diagnostics/diag-pretrend-%s.png", tolower(outcome_var)),
    p, width = 6, height = 4, dpi = 300
  )
}

# =============================================================================
# DIAGNOSTIC 3: Synthetic Weight Concentration
# =============================================================================

weight_results <- tibble(
  outcome = character(),
  cohort = integer(),
  n_controls = integer(),
  n_nonzero_weights = integer(),
  max_weight = numeric(),
  top5_weight_share = numeric(),
  herfindahl = numeric()
)

# --- Hospital-level outcomes ---
for (outcome_var in hosp_outcomes) {
  outcome_sym <- sym(outcome_var)
  ## Financial outcomes use the same pre-period as the main SDID estimates
  pre_o <- if (outcome_var %in% fin.vars) financial.pre else 5

  for (c in cohorts) {
    synth.c <- stack.hosp %>%
      filter(stack_group == c) %>%
      group_by(ID) %>%
      mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(!is.na(!!outcome_sym), min_bedsize <= bed.cut,
             stacked_event_time >= -pre_o) %>%
      select(ID, year, outcome = !!outcome_sym, treated = post_treat)

    bal.c <- tryCatch({
      as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
    }, error = function(e) NULL)

    if (!is.null(bal.c) && nrow(bal.c) > 0) {
      setup <- tryCatch({
        panel.matrices(as.data.frame(bal.c))
      }, error = function(e) NULL)

      if (!is.null(setup) && setup$N0 > 0 && setup$T0 > 0) {
        est <- tryCatch({
          synthdid_estimate(setup$Y, setup$N0, setup$T0)
        }, error = function(e) NULL)

        if (!is.null(est)) {
          w <- attr(est, "weights")$omega

          n_controls <- length(w)
          n_nonzero <- sum(w > 0.001)
          max_w <- max(w)
          top5_share <- sum(sort(w, decreasing = TRUE)[1:min(5, length(w))])
          hhi <- sum(w^2)

          weight_results <- bind_rows(weight_results, tibble(
            outcome = outcome_var,
            cohort = c,
            n_controls = n_controls,
            n_nonzero_weights = n_nonzero,
            max_weight = round(max_w, 3),
            top5_weight_share = round(top5_share, 3),
            herfindahl = round(hhi, 4)
          ))
        }
      }
    }
  }
}

# --- State-level outcomes ---
for (outcome_var in state_outcomes) {

  for (c in cohorts) {
    synth.c <- stack.state %>%
      filter(stack_group == c) %>%
      transmute(ID = as.numeric(factor(MSTATE)),
                year,
                outcome = .data[[outcome_var]],
                treated = post_treat)

    bal.c <- tryCatch({
      as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
    }, error = function(e) NULL)

    if (!is.null(bal.c) && nrow(bal.c) > 0) {
      setup <- tryCatch({
        panel.matrices(as.data.frame(bal.c))
      }, error = function(e) NULL)

      if (!is.null(setup) && setup$N0 > 0 && setup$T0 > 0) {
        est <- tryCatch({
          synthdid_estimate(setup$Y, setup$N0, setup$T0)
        }, error = function(e) NULL)

        if (!is.null(est)) {
          w <- attr(est, "weights")$omega

          n_controls <- length(w)
          n_nonzero <- sum(w > 0.001)
          max_w <- max(w)
          top5_share <- sum(sort(w, decreasing = TRUE)[1:min(5, length(w))])
          hhi <- sum(w^2)

          weight_results <- bind_rows(weight_results, tibble(
            outcome = outcome_var,
            cohort = c,
            n_controls = n_controls,
            n_nonzero_weights = n_nonzero,
            max_weight = round(max_w, 3),
            top5_weight_share = round(top5_share, 3),
            herfindahl = round(hhi, 4)
          ))
        }
      }
    }
  }
}

# =============================================================================
# Save CSV outputs
# =============================================================================

write_csv(sample_comparison, "results/diagnostics/diag-sample-comparison.csv")
write_csv(sample_summary, "results/diagnostics/diag-sample-summary.csv")
write_csv(pretrend_results, "results/diagnostics/diag-pretrend-results.csv")
write_csv(weight_results, "results/diagnostics/diag-weight-results.csv")

# =============================================================================
# LaTeX table output (tabular innards only)
# =============================================================================

# Panel structure matching paper: A = Financial, B = Capacity, C = Organizational
panel_a <- c("margin", "current_ratio", "net_fixed", "capex", "net_pat_rev", "tot_operating_exp")
panel_b <- c("BDTOT", "OBBD", "FTERN", "ip_per_bed")
panel_c <- c("system", "closures", "mergers")

fmt_pct <- function(x) sprintf("%.1f\\%%", x)

# --- Table 1: Sample Comparison ---

tex_sample <- c(
  "\\begin{tabular}{l c cc cc c}",
  "\\toprule",
  " & & \\multicolumn{2}{c}{Before balancing} & \\multicolumn{2}{c}{SDID panel} & \\\\",
  "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6}",
  "Outcome & Cohort & Treated & Control & Treated & Control & \\% Lost \\\\",
  "\\midrule",
  "\\multicolumn{7}{l}{\\textit{Panel A: Financial Performance}} \\\\"
)

for (panel_outcomes in list(panel_a, panel_b, panel_c)) {
  if (identical(panel_outcomes, panel_b)) {
    tex_sample <- c(tex_sample,
      "\\midrule",
      "\\multicolumn{7}{l}{\\textit{Panel B: Capacity \\& Staffing}} \\\\")
  }
  if (identical(panel_outcomes, panel_c)) {
    tex_sample <- c(tex_sample,
      "\\midrule",
      "\\multicolumn{7}{l}{\\textit{Panel C: Organizational Changes}} \\\\")
  }

  for (ov in panel_outcomes) {
    rows <- sample_comparison %>% filter(outcome == ov)
    for (j in seq_len(nrow(rows))) {
      r <- rows[j, ]
      lab <- if (j == 1) outcome_labels[ov] else ""
      tex_sample <- c(tex_sample,
        sprintf("%s & %d & %d & %d & %d & %d & %s \\\\",
                lab, r$cohort,
                r$pre_n_treated, r$pre_n_control,
                r$sdid_n_treated, r$sdid_n_control,
                fmt_pct(r$pct_lost_to_balance)))
    }
  }
}
tex_sample <- c(tex_sample, "\\bottomrule", "\\end{tabular}")

writeLines(tex_sample, "results/diagnostics/diag-sample-comparison.tex")

# --- Table 2: Pre-Trend Tests ---

stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("$^{***}$")
  if (p < 0.05) return("$^{**}$")
  if (p < 0.10) return("$^{*}$")
  return("")
}

tex_pretrend <- c(
  "\\begin{tabular}{l c r r r}",
  "\\toprule",
  "Outcome & Cohort & Coefficient & Std.\\ Error & $p$-value \\\\",
  "\\midrule",
  "\\multicolumn{5}{l}{\\textit{Panel A: Financial Performance}} \\\\"
)

for (panel_outcomes in list(panel_a, panel_b, panel_c)) {
  if (identical(panel_outcomes, panel_b)) {
    tex_pretrend <- c(tex_pretrend,
      "\\midrule",
      "\\multicolumn{5}{l}{\\textit{Panel B: Capacity \\& Staffing}} \\\\")
  }
  if (identical(panel_outcomes, panel_c)) {
    tex_pretrend <- c(tex_pretrend,
      "\\midrule",
      "\\multicolumn{5}{l}{\\textit{Panel C: Organizational Changes}} \\\\")
  }

  for (ov in panel_outcomes) {
    rows <- pretrend_results %>% filter(outcome == ov)
    if (nrow(rows) == 0) {
      tex_pretrend <- c(tex_pretrend,
        sprintf("%s & --- & --- & --- & --- \\\\", outcome_labels[ov]))
      next
    }
    for (j in seq_len(nrow(rows))) {
      r <- rows[j, ]
      lab <- if (j == 1) outcome_labels[ov] else ""
      tex_pretrend <- c(tex_pretrend,
        sprintf("%s & %d & %.4f%s & %.4f & %.3f \\\\",
                lab, r$cohort, r$trend_coef, stars(r$trend_pval),
                r$trend_se, r$trend_pval))
    }
  }
}
tex_pretrend <- c(tex_pretrend, "\\bottomrule", "\\end{tabular}")

writeLines(tex_pretrend, "results/diagnostics/diag-pretrend.tex")

# --- Table 3: Weight Concentration ---

tex_weights <- c(
  "\\begin{tabular}{l c r r r r}",
  "\\toprule",
  "Outcome & Cohort & $N$ Controls & Max Weight & Top 5 Share & HHI \\\\",
  "\\midrule",
  "\\multicolumn{6}{l}{\\textit{Panel A: Financial Performance}} \\\\"
)

for (panel_outcomes in list(panel_a, panel_b, panel_c)) {
  if (identical(panel_outcomes, panel_b)) {
    tex_weights <- c(tex_weights,
      "\\midrule",
      "\\multicolumn{6}{l}{\\textit{Panel B: Capacity \\& Staffing}} \\\\")
  }
  if (identical(panel_outcomes, panel_c)) {
    tex_weights <- c(tex_weights,
      "\\midrule",
      "\\multicolumn{6}{l}{\\textit{Panel C: Organizational Changes}} \\\\")
  }

  for (ov in panel_outcomes) {
    rows <- weight_results %>% filter(outcome == ov)
    if (nrow(rows) == 0) {
      tex_weights <- c(tex_weights,
        sprintf("%s & --- & --- & --- & --- & --- \\\\", outcome_labels[ov]))
      next
    }
    for (j in seq_len(nrow(rows))) {
      r <- rows[j, ]
      lab <- if (j == 1) outcome_labels[ov] else ""
      tex_weights <- c(tex_weights,
        sprintf("%s & %d & %d & %.3f & %.3f & %.4f \\\\",
                lab, r$cohort, r$n_controls,
                r$max_weight, r$top5_weight_share, r$herfindahl))
    }
  }
}
tex_weights <- c(tex_weights, "\\bottomrule", "\\end{tabular}")

writeLines(tex_weights, "results/diagnostics/diag-weights.tex")
