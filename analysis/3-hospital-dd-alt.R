# Meta --------------------------------------------------------------------
## Alternative identification: eligibility-restricted, no state-level requirement
## Compares CAH converters to non-CAH hospitals meeting same bed size criteria
## Cohorts 1999-2005; controls = never-CAH + not-yet-CAH hospitals
## Hospital-level outcomes only (closures/mergers require state-level design)

# Expects from _run-analysis.r: stack.elig, est.dat, bed.cut, post, financial.pre
elig.cohorts <- 1999:2005

# Outcome map (hospital-level only) -----------------------------------------
elig_outcome_map <- list(
  margin        = list(label = "Operating margin",             stub = "elig-margin",       pre_period = financial.pre),
  current_ratio = list(label = "Current ratio",                stub = "elig-currentratio", pre_period = financial.pre),
  net_fixed     = list(label = "Net fixed assets",             stub = "elig-netfixed",     pre_period = financial.pre),
  capex         = list(label = "Capital expenditures per bed", stub = "elig-capex",        pre_period = financial.pre),
  net_pat_rev       = list(label = "Net patient revenue per bed",  stub = "elig-netpatrev",  pre_period = financial.pre),
  tot_operating_exp = list(label = "Operating expenses per bed",   stub = "elig-totopexp",   pre_period = financial.pre),
  BDTOT         = list(label = "Total beds",                   stub = "elig-beds"),
  OBBD          = list(label = "OB beds",                      stub = "elig-beds_ob"),
  FTERN         = list(label = "FTE RNs",                      stub = "elig-ftern"),
  ip_per_bed    = list(label = "Inpatient days per bed",       stub = "elig-ipdays"),
  system        = list(label = "System membership",            stub = "elig-system"),
  ## AHA visits and births, reported in their own table (Appendix Table 6)
  ed_per_bed    = list(label = "ED visits per bed",                stub = "edvisits-elig",       group = "visits"),
  op_per_bed    = list(label = "Outpatient visits per bed",        stub = "opvisits-elig",       group = "visits"),
  opnon_per_bed = list(label = "Non-ED outpatient visits per bed", stub = "opvisits-noned-elig", group = "visits"),
  BIRTHS        = list(label = "Births",                           stub = "births-elig",         group = "visits")
)

# Results collector ----------------------------------------------------------
elig.results <- tibble(
  outcome = character(), sdid_att = numeric(), sdid_ci_low = numeric(), sdid_ci_high = numeric(),
  cs_att = numeric(), cs_ci_low = numeric(), cs_ci_high = numeric(),
  sdid_ntr = numeric(), group = character(), pre_mean = numeric()
)

elig.cohort.results <- tibble(
  outcome = character(), cohort = numeric(), att = numeric(),
  se = numeric(), Ntr = numeric()
)

## The treated and synthetic paths plotted in the appendix figures. The treated
## pre-designation mean is the base for the percent effects quoted in the
## response letter, so the paths go to a file as well as a plot.
elig.event.paths <- tibble(
  outcome = character(), tau = numeric(), treated = numeric(), synthetic = numeric()
)

## Event-time effects behind the SDID/CS overlay figures, pooled across cohorts,
## and the cohort covariances the pre-trend sensitivity analysis needs.
elig.event.series <- tibble()
elig.event.cohort <- tibble()
elig.event.vcov   <- tibble()

## Size of the control pool this design draws on, which Section 4.6 quotes.
## Hospitals serve as controls in more than one cohort, so the distinct count
## across cohorts is the relevant figure; the per-cohort counts go alongside it.
elig.control.pool <- stack.elig %>%
  group_by(ID) %>% mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>% ungroup() %>%
  filter(min_bedsize <= bed.cut, treated == 0)

bind_rows(
  elig.control.pool %>%
    summarise(stack_group = NA_real_, n_control = n_distinct(ID)),
  elig.control.pool %>%
    group_by(stack_group) %>% summarise(n_control = n_distinct(ID), .groups = "drop")
) %>%
  write_csv("results/diagnostics/elig-control-pool.csv")

# Main loop ------------------------------------------------------------------
for (oname in names(elig_outcome_map)) {
  o <- elig_outcome_map[[oname]]
  outcome_var   <- oname
  outcome_sym   <- sym(outcome_var)
  outcome_label <- o$label
  file_stub     <- o$stub
  pp            <- if (!is.null(o$pre_period)) o$pre_period else 5


  # -----------------------------------------------------------------------
  # SDID across cohorts
  # -----------------------------------------------------------------------
  run_sdid_elig <- function(c) {
    synth.c <- stack.elig %>%
      filter(stack_group == c) %>%
      group_by(ID) %>% mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>% ungroup() %>%
      filter(!is.na(!!outcome_sym), min_bedsize <= bed.cut,
             stacked_event_time >= -pp) %>%
      select(ID, year, outcome = !!outcome_sym, post_treat, min_bedsize, treated)

    bal.c <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))

    ## guard against too few treated or control units
    n_yr <- length(unique(bal.c$year))
    n_tr <- sum(bal.c$treated == 1) / max(n_yr, 1)
    n_co <- sum(bal.c$treated == 0) / max(n_yr, 1)
    if (n_tr < 2 || n_co < 2) return(NULL)

    setup <- panel.matrices(as.data.frame(bal.c))
    est   <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
    se_c  <- synthdid_se(est, method = "jackknife")
    ev    <- sdid_event(est, setup, c)
    cs    <- cs_cohort(bal.c, c, min.es = -pp, max.es = 5)

    w      <- attr(est, "weights")$omega
    Y      <- setup$Y
    N0     <- setup$N0
    N      <- nrow(Y)
    years  <- as.numeric(colnames(Y))
    tau    <- years - c

    syn_path   <- as.numeric(drop(t(w) %*% Y[1:N0, , drop = FALSE]))
    treat_path <- colMeans(Y[(N0 + 1):N, , drop = FALSE])
    Ntr        <- N - N0

    list(
      cohort = c,
      paths  = tibble(cohort = c, tau = tau,
                       treated = treat_path, synthetic = syn_path, Ntr = Ntr),
      att_se = tibble(cohort = c, att = as.numeric(est),
                       se = as.numeric(se_c), Ntr = Ntr,
                       pre_mean = mean(Y[(N0 + 1):N, 1:setup$T0])),
      sdid_es   = ev$es,
      sdid_vcov = vcov_long(ev$vcov, c, ev$es$e),
      cs        = cs
    )
  }

  out       <- map(elig.cohorts, possibly(run_sdid_elig, otherwise = NULL))
  out       <- compact(out)
  paths_all <- bind_rows(map(out, "paths"))
  atts_all  <- bind_rows(map(out, "att_se"))

  if (nrow(atts_all) == 0) {
    message(sprintf("elig: %s skipped, SDID failed for every cohort", oname))
    next
  }

  # Collect cohort-specific SDID estimates
  elig.cohort.results <- bind_rows(elig.cohort.results,
    atts_all %>% mutate(outcome = outcome_label))

  ## aggregate paths
  agg_paths <- paths_all %>%
    group_by(tau) %>%
    summarise(treated   = weighted.mean(treated, w = Ntr, na.rm = TRUE),
              synthetic = weighted.mean(synthetic, w = Ntr, na.rm = TRUE),
              .groups = "drop")

  ## pooled ATT
  att_w   <- with(atts_all, sum(Ntr * att) / sum(Ntr))
  se_w    <- with(atts_all, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))
  ci_low  <- att_w - 1.96 * se_w
  ci_high <- att_w + 1.96 * se_w

  ## per-cohort labels
  atts_table <- atts_all %>%
    mutate(ci_lo = att - 1.96 * se, ci_hi = att + 1.96 * se) %>%
    arrange(cohort) %>%
    transmute(line = sprintf("%d: %.2f [%.2f, %.2f]", cohort, att, ci_lo, ci_hi))

  att_header <- "ATT and 95%CI"
  att_body   <- paste(
    sprintf("Overall: %.2f [%.2f, %.2f]", att_w, ci_low, ci_high),
    paste(atts_table$line, collapse = "\n"),
    sep = "\n"
  )

  ## SDID plot
  yrange <- diff(range(c(agg_paths$treated, agg_paths$synthetic), na.rm = TRUE))
  y_top  <- max(c(agg_paths$treated, agg_paths$synthetic), na.rm = TRUE) + 0.25 * yrange
  x_pos  <- max(agg_paths$tau) - 2

  p_sdid <- ggplot(agg_paths, aes(x = tau)) +
    geom_line(aes(y = treated,   linetype = "Treated"),           linewidth = 0.8, color = "black") +
    geom_line(aes(y = synthetic, linetype = "Synthetic control"), linewidth = 0.8, color = "black") +
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
    annotate("text", x = x_pos, y = y_top,
             label = att_header, hjust = 0, vjust = 1,
             size = 3.6, fontface = "bold") +
    annotate("text", x = x_pos, y = y_top - 0.04 * yrange,
             label = att_body, hjust = 0, vjust = 1,
             size = 3.0, lineheight = 1.05)

  ggsave(sprintf("results/%s-sdid.png", file_stub),
         p_sdid, width = 6.5, height = 4.25, dpi = 300, scale = 1.5)

  # -----------------------------------------------------------------------
  # Event-time effects: SDID per period and CS, cohort by cohort on the same
  # balanced panels, pooled with treated-unit-count weights
  # -----------------------------------------------------------------------
  sdid_es <- pool_event(bind_rows(map(out, "sdid_es"))) %>% mutate(estimator = "SDID")
  cs_list <- compact(map(out, "cs"))
  if (length(cs_list) == 0) stop(sprintf("CS failed for every cohort on %s", oname))
  cs_atts <- bind_rows(map(cs_list, function(x) tibble(att = x$att, se = x$se, Ntr = x$Ntr)))
  cs_es   <- pool_event(bind_rows(map(cs_list, "es"))) %>% mutate(estimator = "CS")
  cs_att_val <- with(cs_atts, sum(Ntr * att) / sum(Ntr))
  cs_se_val  <- with(cs_atts, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))

  series <- bind_rows(sdid_es, cs_es)
  ggsave(sprintf("results/%s-overlay.png", file_stub), overlay_plot(series, outcome_label),
         width = 6.5, height = 4.25, dpi = 300, scale = 1.5)

  ## CS event-study plot
  plot.cs <- ggplot(cs_es %>% mutate(se = replace_na(se, 0)), aes(x = e, y = estimate)) +
    geom_errorbar(aes(ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se),
                  width = 0, linewidth = 0.5, alpha = 0.3, color = "black") +
    geom_point(size = 2.5, color = "black", stroke = 0.1, fill = "white") +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    scale_x_continuous(breaks = seq(-10, 10, by = 1)) +
    labs(x = "Event time", y = "Estimated Effects") +
    theme_bw() +
    theme(legend.position = "none")

  ggsave(sprintf("results/%s-cs.png", file_stub),
         plot.cs, width = 6.5, height = 4.25, dpi = 300, scale = 1.5)

  # -----------------------------------------------------------------------
  # Collect results
  # -----------------------------------------------------------------------
  elig.results <- bind_rows(elig.results, tibble(
    outcome      = outcome_label,
    sdid_att     = as.numeric(att_w),
    sdid_ci_low  = as.numeric(ci_low),
    sdid_ci_high = as.numeric(ci_high),
    cs_att       = cs_att_val,
    cs_ci_low    = cs_att_val - 1.96 * cs_se_val,
    cs_ci_high   = cs_att_val + 1.96 * cs_se_val,
    sdid_ntr     = sum(atts_all$Ntr),
    group        = if (is.null(o$group)) "main" else o$group,
    pre_mean     = with(atts_all, sum(Ntr * pre_mean) / sum(Ntr))
  ))

  elig.event.paths <- bind_rows(elig.event.paths,
    agg_paths %>% mutate(outcome = outcome_label) %>%
      select(outcome, tau, treated, synthetic))

  elig.event.series <- bind_rows(elig.event.series,
    series %>% mutate(outcome = outcome_label))
  elig.event.cohort <- bind_rows(elig.event.cohort,
    bind_rows(map(out, "sdid_es")) %>% mutate(outcome = outcome_label, estimator = "SDID"),
    bind_rows(map(cs_list, "es")) %>% mutate(outcome = outcome_label, estimator = "CS"))
  elig.event.vcov <- bind_rows(elig.event.vcov,
    bind_rows(map(out, "sdid_vcov")) %>%
      mutate(outcome = outcome_label, estimator = "SDID"),
    bind_rows(map(cs_list, function(x) vcov_long(x$vcov, x$es$cohort[1], x$es$e))) %>%
      mutate(outcome = outcome_label, estimator = "CS"))
}

write_csv(elig.event.paths, "results/diagnostics/elig-event-paths.csv")
write_csv(elig.event.series, "results/diagnostics/event-series-elig.csv")
write_csv(elig.event.cohort, "results/diagnostics/event-cohort-elig.csv")
write_csv(elig.event.vcov, "results/diagnostics/event-vcov-elig.csv")

# LaTeX summary table --------------------------------------------------------
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

elig.results        <- elig.results %>% arrange(match(outcome, outcome.order))
elig.cohort.results <- elig.cohort.results %>% arrange(match(outcome, outcome.order))

tex.lines <- elig.results %>%
  rowwise() %>%
  mutate(line = sprintf("%s & %s & %s & %d & %s & %s \\\\",
    outcome, fmt(sdid_att), int(sdid_ci_low, sdid_ci_high), sdid_ntr,
    fmt(cs_att), int(cs_ci_low, cs_ci_high))) %>%
  ungroup() %>%
  pull(line)

## Complete tabular block (LaTeX 2025 breaks \noalign/\omit inside \input within tabular)
tex.lines <- with_group_breaks(tex.lines, elig.results$outcome)

writeLines(c(
  "\\begin{tabular}{lccrcc}",
  "\\toprule",
  "Outcome & ATT & 95\\% CI & $N_{tr}$ & CS ATT & CS 95\\% CI \\\\",
  "\\midrule",
  tex.lines,
  "\\bottomrule",
  "\\end{tabular}"
), "results/att_elig_overall.tex")

## Cohort-specific SDID LaTeX output
cohort_tex_lines <- c("Outcome & ATT & 95\\% CI & $N_{tr}$ \\\\")
for (c in sort(unique(elig.cohort.results$cohort))) {
  cohort_tex_lines <- c(cohort_tex_lines,
    sprintf("\\midrule\n\\multicolumn{4}{l}{\\textit{Cohort %d}} \\\\", c))
  rows <- elig.cohort.results %>% filter(cohort == c)
  for (i in 1:nrow(rows)) {
    r <- rows[i, ]
    cohort_tex_lines <- c(cohort_tex_lines,
      sprintf("%s & %s & %s & %d \\\\",
        r$outcome, fmt(r$att), int(r$att - 1.96*r$se, r$att + 1.96*r$se), r$Ntr))
  }
}
writeLines(cohort_tex_lines, "results/att_elig_cohort.tex")

