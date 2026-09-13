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
  system        = list(label = "System membership",            stub = "elig-system")
)

# Results collector ----------------------------------------------------------
elig.results <- tibble(
  outcome = character(), sdid_att = numeric(), sdid_ci_low = numeric(), sdid_ci_high = numeric(),
  cs_att = numeric(), cs_ci_low = numeric(), cs_ci_high = numeric(),
  sdid_ntr = numeric()
)

elig.cohort.results <- tibble(
  outcome = character(), cohort = numeric(), att = numeric(),
  se = numeric(), Ntr = numeric()
)

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
                       se = as.numeric(se_c), Ntr = Ntr)
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
  # Callaway and Sant'Anna
  # -----------------------------------------------------------------------
  cs.out <- cs_att(est.dat, oname, elig.cohorts, bed.cut,
                   min.es = -pp, max.es = 5)
  if (is.null(cs.out)) stop(sprintf("CS failed for every cohort on %s", oname))

  cs_att_val <- cs.out$att
  cs_se_val  <- cs.out$se

  ## CS event-study plot
  est.cs <- cs.out$es %>%
    transmute(
      event_time = event_time,
      estimate   = estimate,
      se         = se
  ) %>%
    mutate(
      conf.low  = if_else(event_time != -1, estimate - 1.96 * se, 0),
      conf.high = if_else(event_time != -1, estimate + 1.96 * se, 0),
      estimator = "CS"
    ) %>%
    select(-se)

  plot.cs <- ggplot(est.cs, aes(x = event_time, y = estimate)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
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
    sdid_ntr     = sum(atts_all$Ntr)
  ))
}

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

tex.lines <- elig.results %>%
  rowwise() %>%
  mutate(line = sprintf("%s & %s & %s & %d & %s & %s \\\\",
    outcome, fmt(sdid_att), int(sdid_ci_low, sdid_ci_high), sdid_ntr,
    fmt(cs_att), int(cs_ci_low, cs_ci_high))) %>%
  ungroup() %>%
  pull(line)

## Complete tabular block (LaTeX 2025 breaks \noalign/\omit inside \input within tabular)
# Insert group separators: financial (1-6), capacity (7-10), organizational (11)
tex.lines <- append(tex.lines, "\\addlinespace", after = 6)
tex.lines <- append(tex.lines, "\\addlinespace", after = 11)  # shifted by 1

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

