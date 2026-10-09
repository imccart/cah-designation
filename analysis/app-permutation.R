# Meta --------------------------------------------------------------------
## Permutation inference for SDID estimates
## Provides non-parametric inference for SDID with small N_tr
## Uses stack.hosp for the hospital-level outcomes and the clean-donor state
## SDID (adoption-timing permutation) for closures and mergers, cohorts 1999-2001
## Produces: multi-panel histogram figure + p-value CSV

# Expects from _run-analysis.r: stack.hosp, state.dat, state.timing, post, hosp.results.table, bed.cut, financial.pre

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
  system            = list(label = "System membership",            pre_period = 5),
  ## AHA visits and births
  ed_per_bed        = list(label = "ED visits per bed", pre_period = 5),
  op_per_bed        = list(label = "Outpatient visits per bed", pre_period = 5),
  opnon_per_bed     = list(label = "Non-ED outpatient visits per bed", pre_period = 5),
  BIRTHS            = list(label = "Births", pre_period = 5)
)

# Results collector ------------------------------------------------------------
perm.results <- tibble(
  outcome = character(), att_actual = numeric(), p_value = numeric()
)
perm.cohort.draws <- tibble()
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
  ## the visit outcomes were added after the others; the random-number state
  ## at this point is restored before the state-level loop below, so the
  ## closure and merger placebos are the same draws as without them
  if (oname == "ed_per_bed") rng_before_visits <- .Random.seed
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
    perm.cohort.draws <- bind_rows(perm.cohort.draws, perm_all %>% mutate(outcome = outcome_label, draw = i))
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

assign(".Random.seed", rng_before_visits, envir = globalenv())

# State-level outcomes ---------------------------------------------------------
# Units are states, and the estimator is the clean-donor SDID of
# 4-changes-state-dd.R (sdid_state_clean), whose donors in each post year are
# states not yet treated. Treated states cannot be reassigned to the control
# group under that design, so the placebo draws permute adoption timing: each
# draw reshuffles adoption years among the states adopting in 1999-2002 (cohort
# sizes kept), holds never-adopting and earlier-adopting states fixed, and reruns
# the estimator under the reshuffled timing.

state_att <- function(timing, outcome_var) {
  fits <- map(perm.cohorts, ~sdid_state_clean(state.dat, timing, .x, outcome_var,
                                              pre = 5, post = post, jack = FALSE))
  if (any(map_lgl(fits, is.null))) return(NULL)
  tibble(cohort = perm.cohorts, att = map_dbl(fits, "att"), Ntr = map_dbl(fits, "Ntr"))
}

timing_adopters <- names(state.timing)[state.timing %in% 1999:2002]

for (ov in c("closures", "mergers")) {
  outcome_label <- if (ov == "closures") "Closures" else "Mergers"

  actual <- state_att(state.timing, ov)
  if (is.null(actual)) {
    message(sprintf("permutation: %s skipped, SDID failed for a cohort", ov))
    next
  }
  att_actual <- with(actual, sum(Ntr * att) / sum(Ntr))

  perm_atts <- numeric(n_perms)
  for (i in seq_len(n_perms)) {
    g <- state.timing
    g[timing_adopters] <- sample(state.timing[timing_adopters])
    pa <- state_att(g, ov)
    perm_atts[i] <- if (is.null(pa)) NA else with(pa, sum(Ntr * att) / sum(Ntr))
    if (!is.null(pa)) perm.cohort.draws <- bind_rows(perm.cohort.draws, pa %>% mutate(outcome = outcome_label, draw = i))
    if (i %% 100 == 0) message(sprintf("permutation: %s draw %d of %d", ov, i, n_perms))
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
perm_levels <- outcome.order[outcome.order %in% perm.results$outcome]
perm_plot_dat <- perm.distributions %>%
  left_join(perm.results, by = "outcome", relationship = "many-to-one") %>%
  mutate(outcome = factor(outcome, levels = perm_levels))

p_perm <- ggplot(perm_plot_dat, aes(x = att_perm)) +
  geom_histogram(bins = 40, fill = "gray70", color = "gray40", linewidth = 0.2) +
  geom_vline(aes(xintercept = att_actual), color = "red", linewidth = 0.6,
             linetype = "dashed") +
  geom_text(data = perm.results %>%
              mutate(outcome = factor(outcome, levels = perm_levels)),
            aes(x = att_actual, y = Inf, label = sprintf("p = %.3f", p_value)),
            vjust = 1.5, hjust = "inward", color = "red", size = 2.8) +
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
write_csv(perm.distributions, "results/diagnostics/permutation_draws.csv")
write_csv(perm.cohort.draws, "results/diagnostics/permutation_cohort_draws.csv")

