# Hospital-level DD: state-timing identification
# Self-contained: owns outcome map, loop, results collection, and tex output
# Expects from _run-analysis.r: stack.hosp, est.dat, bed.cut, post, state.cut, financial.pre

# Outcome map ----------------------------------------------------------------
hosp_outcome_map <- list(
  margin        = list(label="Operating margin",             stub="margin",       cohorts=1999:2001, pre_period=financial.pre),
  current_ratio = list(label="Current ratio",                stub="currentratio", cohorts=1999:2001, pre_period=financial.pre),
  net_fixed     = list(label="Net fixed assets",             stub="netfixed",     cohorts=1999:2001, pre_period=financial.pre),
  capex         = list(label="Capital expenditures per bed", stub="capex",        cohorts=1999:2001, pre_period=financial.pre),
  net_pat_rev       = list(label="Net patient revenue per bed",  stub="netpatrev",  cohorts=1999:2001, pre_period=financial.pre),
  tot_operating_exp = list(label="Operating expenses per bed",   stub="totopexp",   cohorts=1999:2001, pre_period=financial.pre),
  BDTOT         = list(label="Total beds",                   stub="beds",         cohorts=1999:2001),
  OBBD          = list(label="OB beds",                      stub="beds_ob",      cohorts=1999:2001),
  FTERN         = list(label="FTE RNs",                      stub="ftern",        cohorts=1999:2001),
  ip_per_bed    = list(label="Inpatient days per bed",       stub="ipdays",       cohorts=1999:2001),
  system        = list(label="System membership",            stub="system",       cohorts=1999:2001)
)

# Results collectors ----------------------------------------------------------
hosp.results.table <- tibble(
  outcome = character(), sdid_att = numeric(), sdid_ci_low = numeric(), sdid_ci_high = numeric(),
  cs_att = numeric(), cs_ci_low = numeric(), cs_ci_high = numeric(),
  sdid_ntr = numeric()
)

hosp.cohort.results <- tibble(
  outcome = character(), cohort = numeric(), att = numeric(),
  se = numeric(), Ntr = numeric()
)

## The treated and synthetic paths plotted in the main-text figures. Section 4.4
## quotes the margin path directly, so it goes to a file as well as a plot.
hosp.event.paths <- tibble(
  outcome = character(), tau = numeric(), treated = numeric(),
  synthetic = numeric(), cohorts_n = numeric()
)

# Main loop -------------------------------------------------------------------
for (oname in names(hosp_outcome_map)) {
  o <- hosp_outcome_map[[oname]]
  outcome_var   <- oname
  outcome_sym   <- sym(outcome_var)
  outcome_label <- o$label
  file_stub     <- o$stub
  cohorts       <- o$cohorts
  pre_period    <- if (!is.null(o$pre_period)) o$pre_period else 5


  # SYNTH DD for single cohort ------------------------------------------------
  cohort.year <- 2000
  synth.year <- stack.hosp %>% group_by(ID) %>% mutate(min_bedsize=min(BDTOT, na.rm=TRUE)) %>% ungroup() %>%
              filter(stack_group==cohort.year, !is.na(!!outcome_sym), min_bedsize<=bed.cut,
                     stacked_event_time >= -pre_period) %>%
              select(ID, year, outcome=!!outcome_sym, post_treat, min_bedsize, treated)

  balance.year  <- as_tibble(makeBalancedPanel(synth.year, idname="ID", tname="year"))

  setup <- panel.matrices(as.data.frame(balance.year))
  synth.est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
  se  <- synthdid_se(synth.est, method="jackknife")

  w      <- attr(synth.est, "weights")$omega
  Y      <- setup$Y
  N0     <- setup$N0
  N      <- nrow(Y)
  years  <- as.numeric(colnames(Y))

  treated_path   <- colMeans(Y[(N0+1):N, , drop = FALSE])
  synthetic_path <- as.numeric(drop(t(w) %*% Y[1:N0, , drop = FALSE]))

  plot_df.year <- bind_rows(
    tibble(year = years, group = "Synthetic control", value = synthetic_path),
    tibble(year = years, group = "Treated",           value = treated_path)
  ) %>%
    mutate(group = factor(group, levels = c("Synthetic control","Treated")))

  att      <- as.numeric(synth.est)
  se_num   <- as.numeric(se)
  ci_low   <- att - 1.96 * se_num
  ci_high  <- att + 1.96 * se_num
  att_lab <- sprintf(
    "ATT and 95%% CI for %d: %.2f, [%.2f, %.2f]",
    cohort.year, att, ci_low, ci_high
  )

  plot.year <- ggplot(plot_df.year, aes(x = year, y = value, linetype = group)) +
    geom_line(linewidth = 0.8, color = "black") +
    geom_vline(xintercept = cohort.year - 0.5, linewidth = 1) +
    scale_linetype_manual(values = c("dashed","solid")) +
    scale_x_continuous(breaks = seq(min(plot_df.year$year), max(plot_df.year$year), by = 1)) +
    labs(x = "Year", y = outcome_label, linetype = NULL) +
    theme_bw() +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.88, 0.12),
      legend.background = element_blank(),
      legend.key.width = unit(2, "cm")
    ) +
    annotate("text", x = max(plot_df.year$year), y = Inf,
             label = att_lab, hjust = 1, vjust = 1.5, size = 3.5)

  ggsave(
    sprintf("results/%s-sdid-yearly.png", file_stub),
    plot.year, width = 6.5, height = 4.25, dpi = 300, scale = 1.5
  )


  # SYNTH DD across cohorts ------------------------------------------------
  run_sdid_cohort <- function(c){
    synth.c <- stack.hosp %>% filter (stack_group==c) %>%
      group_by(ID) %>% mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>% ungroup() %>%
      filter(!is.na(!!outcome_sym), min_bedsize <= bed.cut,
             stacked_event_time >= -pre_period) %>%
      select(ID, year, outcome=!!outcome_sym, post_treat, min_bedsize, treated)

    bal.c   <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
    setup   <- panel.matrices(as.data.frame(bal.c))
    est     <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
    se_c    <- synthdid_se(est, method = "jackknife")

    w       <- attr(est, "weights")$omega
    Y       <- setup$Y
    N0      <- setup$N0
    N       <- nrow(Y)
    years   <- as.numeric(colnames(Y))
    tau     <- years - c

    syn_path   <- as.numeric(drop(t(w) %*% Y[1:N0, , drop = FALSE]))
    treat_path <- colMeans(Y[(N0 + 1):N, , drop = FALSE])
    Ntr        <- N - N0

    list(
      cohort = c,
      paths  = tibble(cohort = c, tau = tau,
                      treated = treat_path, synthetic = syn_path,
                      Ntr = Ntr),
      att_se = tibble(cohort = c,
                      att = as.numeric(est),
                      se  = as.numeric(se_c),
                      Ntr = Ntr)
    )
  }

  out        <- map(cohorts, run_sdid_cohort)
  paths_all  <- bind_rows(map(out, "paths"))
  atts_all   <- bind_rows(map(out, "att_se"))

  agg_paths <- paths_all %>%
    group_by(tau) %>%
    summarise(
      treated   = weighted.mean(treated, w = Ntr, na.rm = TRUE),
      synthetic = weighted.mean(synthetic, w = Ntr, na.rm = TRUE),
      cohorts_n = n(),
      .groups = "drop"
    )

  att_w <- with(atts_all, sum(Ntr * att) / sum(Ntr))
  se_w  <- with(atts_all, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))

  ci_low  <- att_w - 1.96 * se_w
  ci_high <- att_w + 1.96 * se_w
  atts_table <- atts_all %>%
    mutate(ci_lo = att - 1.96 * se,
           ci_hi = att + 1.96 * se) %>%
    arrange(match(cohort, cohorts)) %>%
    transmute(line = sprintf("%d: %.2f [%.2f, %.2f]", cohort, att, ci_lo, ci_hi))

  att_lab <- paste(
    "ATT and 95%CI",
    sprintf("Overall: %.2f [%.2f, %.2f]", att_w, ci_low, ci_high),
    paste(atts_table$line, collapse = "\n"),
    sep = "\n"
  )

  att_header <- "ATT and 95%CI"
  att_body   <- sub("^ATT and 95%CI\\n", "", att_lab)

  yrange <- diff(range(c(agg_paths$treated, agg_paths$synthetic), na.rm = TRUE))
  y_top  <- max(c(agg_paths$treated, agg_paths$synthetic), na.rm = TRUE) + .25*yrange
  x_pos  <- max(agg_paths$tau) - 2

  p_evt <- ggplot(agg_paths, aes(x = tau)) +
    geom_line(aes(y = treated,  linetype = "Treated"),  linewidth = 0.8, color = "black") +
    geom_line(aes(y = synthetic,linetype = "Synthetic control"), linewidth = 0.8, color = "black") +
    geom_vline(xintercept = -0.5, linewidth = 1) +
    scale_linetype_manual(values = c("Treated" = "solid", "Synthetic control" = "dashed")) +
    scale_x_continuous(breaks = seq(min(agg_paths$tau), max(agg_paths$tau), by = 1)) +
    labs(x = "Event time", y = outcome_label, linetype = NULL) +
    theme_bw() +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.88, 0.12),
      legend.background = element_blank(),
      legend.key.width = unit(2, "cm")
    ) +
    annotate("label", x = x_pos, y = y_top,
             label = att_header, hjust = 0, vjust = 1,
             size = 4.0, fontface = "bold",
             fill = "white", alpha = 0.85, label.size = 0) +
    annotate("label", x = x_pos, y = y_top - 0.04*yrange,
             label = att_body, hjust = 0, vjust = 1,
             size = 3.8, lineheight = 1.05,
             fill = "white", alpha = 0.85, label.size = 0)

  ggsave(
    sprintf("results/%s-sdid.png", file_stub),
    p_evt, width = 6.5, height = 4.25, dpi = 300, scale = 1.5
  )


  # Callaway and Sant'Anna -----------------------------------------------------
  min.es <- -pre_period
  max.es <- 5

  cs.out <- cs_att(est.dat, oname, cohorts, bed.cut,
                   min.es = min.es, max.es = max.es)
  if (is.null(cs.out)) stop(sprintf("CS failed for every cohort on %s", oname))

  est.cs <- cs.out$es %>%
    mutate(
      conf.low  = if_else(event_time!= -1, estimate - 1.96 * se, 0),
      conf.high = if_else(event_time!= -1, estimate + 1.96 * se, 0),
      estimator = "CS"
    ) %>%
    select(-se)

  plot.cs <- ggplot(est.cs, aes(x = event_time, y = estimate)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(width = 0.5),
                  width = 0, linewidth = 0.5, alpha = 0.3, color = "black") +
    geom_point(position = position_dodge(width = 0.5),
               size = 2.5, color = "black", stroke = 0.1, fill="white") +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    scale_x_continuous(breaks = seq(-10, 10, by = 1)) +
    labs(x = "Event time", y = "Estimated Effects") +
    theme_bw() +
    theme(
      legend.position = "none"
    )

  ggsave(
    sprintf("results/%s-cs.png", file_stub),
    plot.cs, width = 6.5, height = 4.25, dpi = 300, scale = 1.5
  )

  # Collect results ------------------------------------------------------------
  cs_att_val <- cs.out$att
  cs_se_val  <- cs.out$se

  hosp.results.table <- bind_rows(hosp.results.table, tibble(
    outcome      = outcome_label,
    sdid_att     = as.numeric(att_w),
    sdid_ci_low  = as.numeric(ci_low),
    sdid_ci_high = as.numeric(ci_high),
    cs_att       = cs_att_val,
    cs_ci_low    = cs_att_val - 1.96 * cs_se_val,
    cs_ci_high   = cs_att_val + 1.96 * cs_se_val,
    sdid_ntr     = sum(atts_all$Ntr)
  ))

  hosp.cohort.results <- bind_rows(hosp.cohort.results,
    atts_all %>% mutate(outcome = outcome_label))

  hosp.event.paths <- bind_rows(hosp.event.paths,
    agg_paths %>% mutate(outcome = outcome_label) %>%
      select(outcome, tau, treated, synthetic, cohorts_n))
}

write_csv(hosp.event.paths, "results/diagnostics/sdid-event-paths.csv")

# Cohort-specific SDID LaTeX output -------------------------------------------
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

cohort_tex_lines <- c("Outcome & ATT & 95\\% CI & $N_{tr}$ \\\\")
for (c in sort(unique(hosp.cohort.results$cohort))) {
  cohort_tex_lines <- c(cohort_tex_lines,
    sprintf("\\midrule\n\\multicolumn{4}{l}{\\textit{Cohort %d}} \\\\", c))
  rows <- hosp.cohort.results %>% filter(cohort == c)
  for (i in 1:nrow(rows)) {
    r <- rows[i, ]
    cohort_tex_lines <- c(cohort_tex_lines,
      sprintf("%s & %s & %s & %d \\\\",
        r$outcome, fmt(r$att), int(r$att - 1.96*r$se, r$att + 1.96*r$se), r$Ntr))
  }
}
writeLines(cohort_tex_lines, "results/att_cohort.tex")

