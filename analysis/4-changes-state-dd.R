# State-level DD: closures and mergers
# Self-contained: owns outcome map, loop, results collection
# Expects from _run-analysis.r: stack.state, state.dat, bed.cut, post, state.cut

# Outcome map ----------------------------------------------------------------
state_outcome_map <- list(
  closures = list(label="Closures", stub="closure-rate", cohorts=1999:2001),
  mergers  = list(label="Mergers",  stub="merger-rate",  cohorts=1999:2001)
)

# Results collector -----------------------------------------------------------
state.results.table <- tibble(
  outcome = character(), sdid_att = numeric(), sdid_ci_low = numeric(), sdid_ci_high = numeric(),
  cs_att = numeric(), cs_ci_low = numeric(), cs_ci_high = numeric(),
  sdid_ntr = numeric()
)

# Main loop -------------------------------------------------------------------
for (oname in names(state_outcome_map)) {
  o <- state_outcome_map[[oname]]
  outcome_var   <- oname
  outcome_sym   <- sym(outcome_var)
  outcome_label <- o$label
  file_stub     <- o$stub
  cohorts       <- o$cohorts


  # Synthetic DD ---------------------------------------------------------

  ## Single year
  cohort.year <- 2000
  denom.year <- stack.state %>%
    filter(stacked_event_time <= -1, treated==1, stack_group==cohort.year) %>%
    summarize(mean_hosp=mean(hospitals, na.rm=TRUE)) %>%
    pull(mean_hosp)

  synth.year <- stack.state %>%
    filter(stack_group == cohort.year) %>%
    transmute(ID = as.numeric(factor(MSTATE)),
              year,
              rate = .data[[outcome_var]],
              post_treat)

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
    tibble(year = years, group="Synthetic control",
           rate_out = (synthetic_path/denom.year)*100),
    tibble(year = years, group="Treated",
           rate_out = (treated_path/denom.year)*100)
  )

  att      <- (as.numeric(synth.est)/denom.year)*100
  se_num   <- (as.numeric(se)/denom.year)*100
  ci_low   <- att - 1.96 * se_num
  ci_high  <- att + 1.96 * se_num

  y_top <- max(plot_df.year$rate_out, na.rm = TRUE)
  x_pos <- max(plot_df.year$year)

  lab_txt <- sprintf("ATT and 95%% CI: %.2f [%.2f, %.2f]", att, ci_low, ci_high)

  plot.sdid.year <- ggplot(plot_df.year, aes(x = year, y = rate_out, linetype = group)) +
    geom_line(linewidth = 0.8, color = "black")   +
    geom_vline(xintercept = cohort.year - 1, linewidth = 1)  +
    scale_linetype_manual(values = c("Synthetic control" = "dashed", "Treated" = "solid")) +
    labs(x = "Year", y = paste0(outcome_label, " per 100 hospitals"), linetype = NULL) +
    theme_bw() +
    theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))  +
    annotate("text", x = x_pos, y = y_top, label = lab_txt, hjust = 1, vjust = 1, size = 3.4)

  ggsave(paste0("results/", file_stub, "-sdid-year.png"), plot.sdid.year,
         width = 6.5, height = 4.25, dpi = 300, scale=1.5)

  ## All cohorts
  denom.lag <- stack.state %>%
    filter(stacked_event_time <= -1, treated==1) %>%
    group_by(stack_group) %>%
    summarize(mean_hosp=mean(hospitals, na.rm=TRUE))

  run_sdid_state <- function(c){
    dat <- stack.state %>%
      filter(stack_group == c) %>%
      transmute(ID = as.numeric(factor(MSTATE)),
                year,
                rate = .data[[outcome_var]],
                post_treat)

    bal <- as_tibble(makeBalancedPanel(dat, idname = "ID", tname = "year"))
    setup <- panel.matrices(as.data.frame(bal))

    est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
    se  <- synthdid_se(est, method = "jackknife")

    w       <- attr(est, "weights")$omega
    Y       <- setup$Y
    N0      <- setup$N0
    N       <- nrow(Y)
    years   <- as.numeric(colnames(Y))
    tau     <- years - c

    syn_path   <- as.numeric(drop(t(w) %*% Y[1:N0, , drop = FALSE]))
    treat_path <- colMeans(Y[(N0 + 1):N, , drop = FALSE])
    Ntr        <- N - N0

    denom_c <- denom.lag %>% filter(stack_group == c) %>% pull(mean_hosp)
    list(
      cohort = c,
      paths  = tibble(cohort = c, tau = tau,
                      treated = (treat_path/denom_c)*100,
                      synthetic = (syn_path/denom_c)*100,
                      Ntr = Ntr),
      att_se = tibble(cohort = c,
                      att = (as.numeric(est)/denom_c)*100,
                      se  = (as.numeric(se)/denom_c)*100,
                      Ntr = Ntr)
    )
  }

  out        <- map(cohorts, run_sdid_state)
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

  att_w  <- with(atts_all, sum(Ntr * att) / sum(Ntr))
  se_w   <- with(atts_all, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))
  ci_low <- att_w - 1.96 * se_w
  ci_high<- att_w + 1.96 * se_w

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

  y_top  <- max(c(agg_paths$treated, agg_paths$synthetic), na.rm = TRUE)
  x_pos  <- max(agg_paths$tau) - 2
  yrange <- diff(range(c(agg_paths$treated, agg_paths$synthetic), na.rm = TRUE))

  plot.sdid <- ggplot(agg_paths, aes(x = tau)) +
    geom_line(aes(y = treated,   linetype = "Treated"),   linewidth = 0.8, color = "black") +
    geom_line(aes(y = synthetic, linetype = "Synthetic"), linewidth = 0.8, color = "black") +
    geom_vline(xintercept = -1, linewidth = 1) +
    scale_linetype_manual(values = c("Treated" = "solid", "Synthetic" = "dashed")) +
    scale_x_continuous(breaks = seq(min(agg_paths$tau), max(agg_paths$tau), by = 1)) +
    labs(x = "Event time", y = paste0(outcome_label, " per 100 hospitals"), linetype = NULL) +
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
    annotate("label", x = x_pos, y = y_top - 0.03*yrange,
             label = att_body, hjust = 0, vjust = 1,
             size = 3.8, lineheight = 1.05,
             fill = "white", alpha = 0.85, label.size = 0)

  ggsave(paste0("results/", file_stub, "-sdid.png"), plot.sdid,
         width = 6.5, height = 4.25, dpi = 300, scale = 1.5)


  # Callaway and Sant'Anna ------------------------------------------------
  min.es <- -5
  max.es <- 5

  ## Same denominator basis as the SDID estimates: the cohort-specific
  ## pre-period hospital counts, weighted by treated states per cohort.
  cs_weights <- stack.state %>%
    filter(treated == 1, stack_group %in% cohorts) %>%
    distinct(stack_group, MSTATE) %>%
    count(stack_group, name = "n_tr")

  denom_cs <- denom.lag %>%
    filter(stack_group %in% cohorts) %>%
    inner_join(cs_weights, by = "stack_group") %>%
    summarise(d = weighted.mean(mean_hosp, n_tr)) %>%
    pull(d)

  scale_cs <- 100 / denom_cs

  ## att_gt draws a multiplier bootstrap for its confidence bands. Seeding here
  ## keeps the intervals identical regardless of what ran earlier.
  set.seed(1234)

  csa.raw <- att_gt(yname=outcome_var,
                     gname="state_treat_year",
                     idname="state",
                     tname="year",
                     control_group="notyettreated",
                     panel=TRUE,
                     allow_unbalanced_panel=TRUE,
                     data = state.dat %>% filter(state_treat_year %in% c(0, cohorts)),
                     base_period="universal",
                     est_method="ipw",
                     clustervars="state",
                     bstrap=TRUE)
  csa.att <- aggte(csa.raw, type="simple", na.rm=TRUE)
  csa.es <- aggte(csa.raw, type="dynamic", na.rm=TRUE, min_e=-5, max_e=5)

  est.cs <- tibble(
    event_time = csa.es$egt,
    estimate   = scale_cs*csa.es$att,
    se         = scale_cs*csa.es$se
  ) %>%
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
    labs(x = "Event time", y = paste0(outcome_label, " per 100 hospitals"), linetype = NULL) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    )

  ggsave(paste0("results/", file_stub, "-cs.png"), plot.cs,
         width = 6.5, height = 4.25, dpi = 300, scale=1.5)

  # Collect results ----------------------------------------------------------
  cs_att_val <- csa.es$overall.att * scale_cs
  cs_se_val  <- csa.es$overall.se * scale_cs

  state.results.table <- bind_rows(state.results.table, tibble(
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

