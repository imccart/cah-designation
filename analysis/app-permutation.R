# Meta --------------------------------------------------------------------
## Permutation inference for SDID estimates
## Provides non-parametric inference for SDID with small N_tr
## Uses stack.hosp for the hospital-level outcomes and stack.state for
## closures and mergers, cohorts 1999-2001
## Produces: multi-panel histogram figure + p-value CSV

# Expects from _run-analysis.r: stack.hosp, stack.state, hosp.results.table, bed.cut, financial.pre

perm.cohorts <- c(1999, 2000, 2001)
n_perms <- 500

# Outcome map ------------------------------------------------------------------
perm_outcome_map <- list(
  ## Financial (use financial.pre for pre-period)
  margin            = list(label = "Operating margin",             pre_period = financial.pre),
  current_ratio     = list(label = "Current ratio",                pre_period = financial.pre),
  net_fixed         = list(label = "Net fixed assets",             pre_period = financial.pre),
  capex             = list(label = "Capital expenditures per bed", pre_period = financial.pre),
  net_pat_rev       = list(label = "Net patient revenue per bed",  pre_period = financial.pre),
  tot_operating_exp = list(label = "Operating expenses per bed",   pre_period = financial.pre),
  ## Operational (5-year pre-period)
  BDTOT             = list(label = "Total beds",                   pre_period = 5),
  OBBD              = list(label = "OB beds",                      pre_period = 5),
  FTERN             = list(label = "FTE RNs",                      pre_period = 5),
  ip_per_bed        = list(label = "Inpatient days per bed",       pre_period = 5),
  system            = list(label = "System membership",            pre_period = 5)
)

# Results collector ------------------------------------------------------------
perm.results <- tibble(
  outcome = character(), att_actual = numeric(), p_value = numeric()
)
perm.distributions <- tibble(
  outcome = character(), att_perm = numeric()
)

# Helper: SDID for one cohort --------------------------------------------------
run_sdid_cohort <- function(c, dat, outcome_sym, pp) {
  synth.c <- dat %>%
    filter(stack_group == c) %>%
    group_by(ID) %>%
    mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(!is.na(!!outcome_sym), min_bedsize <= bed.cut,
           stacked_event_time >= -pp) %>%
    select(ID, year, outcome = !!outcome_sym, post_treat)

  bal.c <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
  if (nrow(bal.c) == 0) return(NULL)

  setup <- panel.matrices(as.data.frame(bal.c))
  est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
  Ntr <- nrow(setup$Y) - setup$N0

  tibble(cohort = c, att = as.numeric(est), Ntr = Ntr)
}

# Helper: permuted SDID for one cohort ----------------------------------------
run_permuted_cohort <- function(c, dat, outcome_sym, pp) {
  synth.c <- dat %>%
    filter(stack_group == c) %>%
    group_by(ID) %>%
    mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(!is.na(!!outcome_sym), min_bedsize <= bed.cut,
           stacked_event_time >= -pp) %>%
    select(ID, year, outcome = !!outcome_sym, post_treat, treated)

  bal.c <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
  if (nrow(bal.c) == 0) return(NULL)

  units <- bal.c %>% distinct(ID, treated)
  n_tr <- sum(units$treated == 1)
  n_co <- sum(units$treated == 0)
  if (n_tr < 2 || n_co < 2) return(NULL)

  perm_treated <- sample(units$ID, n_tr)
  bal.perm <- bal.c %>%
    mutate(
      treated_perm = as.integer(ID %in% perm_treated),
      post = as.integer(year >= c),
      post_treat = treated_perm * post
    ) %>%
    select(ID, year, outcome, post_treat)

  setup <- tryCatch(panel.matrices(as.data.frame(bal.perm)), error = function(e) NULL)
  if (is.null(setup)) return(NULL)

  est <- tryCatch(synthdid_estimate(setup$Y, setup$N0, setup$T0), error = function(e) NULL)
  if (is.null(est)) return(NULL)

  tibble(cohort = c, att = as.numeric(est), Ntr = n_tr)
}

# Main loop --------------------------------------------------------------------
set.seed(42)
t0_total <- Sys.time()

for (oname in names(perm_outcome_map)) {
  o <- perm_outcome_map[[oname]]
  outcome_sym   <- sym(oname)
  outcome_label <- o$label
  pp            <- if (!is.null(o$pre_period)) o$pre_period else 5

  t0 <- Sys.time()

  # Actual SDID estimate
  actual <- bind_rows(map(perm.cohorts, ~run_sdid_cohort(.x, stack.hosp, outcome_sym, pp)))
  if (nrow(actual) == 0) {
    message(sprintf("permutation: %s skipped, SDID failed for every cohort", oname))
    next
  }
  att_actual <- with(actual, sum(Ntr * att) / sum(Ntr))

  # Permutation distribution
  perm_atts <- numeric(n_perms)
  for (i in seq_len(n_perms)) {
    perm_cohorts_out <- map(perm.cohorts, ~run_permuted_cohort(.x, stack.hosp, outcome_sym, pp))
    perm_cohorts_out <- compact(perm_cohorts_out)
    perm_all <- bind_rows(perm_cohorts_out)

    if (nrow(perm_all) == 0) {
      perm_atts[i] <- NA
      next
    }
    perm_atts[i] <- with(perm_all, sum(Ntr * att) / sum(Ntr))
  }

  perm_atts <- perm_atts[!is.na(perm_atts)]
  p_val <- mean(abs(perm_atts) >= abs(att_actual))

  t1 <- Sys.time()

  perm.results <- bind_rows(perm.results, tibble(
    outcome = outcome_label, att_actual = att_actual, p_value = p_val
  ))
  perm.distributions <- bind_rows(perm.distributions, tibble(
    outcome = outcome_label, att_perm = perm_atts
  ))
}

# State-level outcomes ---------------------------------------------------------
# Units are states rather than hospitals, and estimates are scaled to counts per
# 100 hospitals using each cohort's pre-period hospital count, as in
# 4-changes-state-dd.R. The scale is a units conversion, so the same cohort
# denominator is applied to the actual and permuted assignments.

state_denom <- stack.state %>%
  filter(stacked_event_time <= -1, treated == 1) %>%
  group_by(stack_group) %>%
  summarize(mean_hosp = mean(hospitals, na.rm = TRUE), .groups = "drop")

run_state_cohort <- function(c, outcome_var, permute) {
  dat <- stack.state %>%
    filter(stack_group == c) %>%
    transmute(ID = as.numeric(factor(MSTATE)), year,
              outcome = .data[[outcome_var]], treated, post_treat)

  bal <- as_tibble(makeBalancedPanel(dat, idname = "ID", tname = "year"))
  if (nrow(bal) == 0) return(NULL)

  units <- bal %>% distinct(ID, treated)
  n_tr <- sum(units$treated == 1)
  n_co <- sum(units$treated == 0)
  if (n_tr < 2 || n_co < 2) return(NULL)

  if (permute) {
    perm_treated <- sample(units$ID, n_tr)
    bal <- bal %>%
      mutate(post = as.integer(year >= c),
             post_treat = as.integer(ID %in% perm_treated) * post)
  }

  setup <- tryCatch(panel.matrices(as.data.frame(bal %>% select(ID, year, outcome, post_treat))),
                    error = function(e) NULL)
  if (is.null(setup)) return(NULL)
  est <- tryCatch(synthdid_estimate(setup$Y, setup$N0, setup$T0), error = function(e) NULL)
  if (is.null(est)) return(NULL)

  denom <- state_denom %>% filter(stack_group == c) %>% pull(mean_hosp)
  tibble(cohort = c, att = (as.numeric(est) / denom) * 100, Ntr = n_tr)
}

for (ov in c("closures", "mergers")) {
  outcome_label <- if (ov == "closures") "Closures" else "Mergers"
  t0 <- Sys.time()

  actual <- bind_rows(compact(map(perm.cohorts, ~run_state_cohort(.x, ov, FALSE))))
  if (nrow(actual) == 0) {
    message(sprintf("permutation: %s skipped, SDID failed for every cohort", ov))
    next
  }
  att_actual <- with(actual, sum(Ntr * att) / sum(Ntr))

  perm_atts <- numeric(n_perms)
  for (i in seq_len(n_perms)) {
    pa <- bind_rows(compact(map(perm.cohorts, ~run_state_cohort(.x, ov, TRUE))))
    perm_atts[i] <- if (nrow(pa) == 0) NA else with(pa, sum(Ntr * att) / sum(Ntr))
  }
  perm_atts <- perm_atts[!is.na(perm_atts)]
  p_val <- mean(abs(perm_atts) >= abs(att_actual))


  perm.results <- bind_rows(perm.results, tibble(
    outcome = outcome_label, att_actual = att_actual, p_value = p_val))
  perm.distributions <- bind_rows(perm.distributions, tibble(
    outcome = outcome_label, att_perm = perm_atts))
}


t1_total <- Sys.time()

# Multi-panel histogram --------------------------------------------------------
perm_plot_dat <- perm.distributions %>%
  left_join(perm.results, by = "outcome") %>%
  mutate(outcome = factor(outcome, levels = perm.results$outcome))

p_perm <- ggplot(perm_plot_dat, aes(x = att_perm)) +
  geom_histogram(bins = 40, fill = "gray70", color = "gray40", linewidth = 0.2) +
  geom_vline(aes(xintercept = att_actual), color = "red", linewidth = 0.6,
             linetype = "dashed") +
  geom_text(data = perm.results %>%
              mutate(outcome = factor(outcome, levels = perm.results$outcome)),
            aes(x = att_actual, y = Inf, label = sprintf("p = %.3f", p_value)),
            vjust = 1.5, hjust = -0.1, color = "red", size = 2.8) +
  facet_wrap(~ outcome, ncol = 3, scales = "free") +
  labs(x = "Placebo ATT", y = "Count") +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  )
ggsave("results/permutation-sdid.png", p_perm,
       width = 10, height = 12, dpi = 300)

# Save CSV ---------------------------------------------------------------------
write_csv(perm.results, "results/diagnostics/permutation_pvalues.csv")

