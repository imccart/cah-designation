# Build stacked dataset at hospital level -----------------------------------------

stack_hosp <- function(pre.period, post.period, state.period) {

  stack.hosp <- tibble()
  for (i in unique(est.dat$eff_year)) {

    ## define treated group relative to eff_year i within event window
    treat.dat <- est.dat %>%
      filter(eff_year == i,
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "treated")

    ## define control group relative to eff_year i within event window
    ## never treated (+ never a CAH designation in the state)
    control.dat.never <- est.dat %>%
      filter(is.na(eff_year), state_treat_year == 0,
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "never")

    ## not yet treated (+ no CAH designation *yet* in the state)
    control.dat.notyet <- est.dat %>%
      filter((eff_year > (i + post.period) | is.na(eff_year)),
             state_treat_year > (i + state.period),
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "notyet")

    ## inner join back to est.dat
    stack.dat.group <- est.dat %>%
      inner_join(bind_rows(treat.dat, control.dat.never, control.dat.notyet),
                 by = c("ID", "year")) %>%
      mutate(state = as.numeric(factor(MSTATE)),
             stacked_event_time = year - i,
             stack_group = i)

    stack.hosp <- bind_rows(stack.hosp, stack.dat.group)
  }

  ## drop the eff_year==0 phantom treatment group
  stack.hosp <- stack.hosp %>% ungroup() %>%
    filter(stack_group > 0) %>%
    mutate(treated = ifelse(treat_type == "treated", 1, 0),
           control_any = ifelse(treat_type != "treated", 1, 0),
           control_notyet = ifelse(treat_type == "notyet", 1, 0),
           control_never = ifelse(treat_type == "never", 1, 0),
           post = ifelse(stacked_event_time >= 0, 1, 0),
           post_treat = post * treated,
           stacked_event_time_treat = stacked_event_time * treated,
           all_treated = sum(treated),
           all_control = sum(control_any),
           all_control_notyet = sum(control_notyet)) %>%
    group_by(stack_group) %>%
    mutate(group_treated = sum(treated),
           group_control_any = sum(control_any),
           group_control_notyet = sum(control_notyet)) %>%
    ungroup() %>%
    mutate(weight_any = case_when(
             treat_type == "treated" ~ 1,
             treat_type != "treated" ~ (group_treated / all_treated) / (group_control_any / all_control)),
           weight_notyet = case_when(
             treat_type == "treated" ~ 1,
             treat_type == "notyet" ~ (group_treated / all_treated) / (group_control_notyet / all_control_notyet))
    )

  stack.hosp
}


# Build stacked dataset at hospital level — no state restriction ----------------------
# Controls are all non-CAH or not-yet-CAH hospitals (no state-level timing requirement)

stack_hosp_elig <- function(pre.period, post.period, cohort.years = NULL) {

  eff_years <- if (!is.null(cohort.years)) cohort.years else unique(na.omit(est.dat$eff_year))

  stack.hosp <- tibble()
  for (i in eff_years) {

    ## treated: hospitals converting in year i
    treat.dat <- est.dat %>%
      filter(eff_year == i,
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "treated")

    ## controls (never): hospitals that never convert to CAH
    control.dat.never <- est.dat %>%
      filter(is.na(eff_year),
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "never")

    ## controls (not-yet): hospitals that convert after the event window
    control.dat.notyet <- est.dat %>%
      filter(eff_year > (i + post.period),
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "notyet")

    ## inner join back to est.dat
    stack.dat.group <- est.dat %>%
      inner_join(bind_rows(treat.dat, control.dat.never, control.dat.notyet),
                 by = c("ID", "year")) %>%
      mutate(state = as.numeric(factor(MSTATE)),
             stacked_event_time = year - i,
             stack_group = i)

    stack.hosp <- bind_rows(stack.hosp, stack.dat.group)
  }

  ## post-processing (same as stack_hosp)
  stack.hosp <- stack.hosp %>% ungroup() %>%
    filter(stack_group > 0) %>%
    mutate(treated = ifelse(treat_type == "treated", 1, 0),
           control_any = ifelse(treat_type != "treated", 1, 0),
           control_notyet = ifelse(treat_type == "notyet", 1, 0),
           control_never = ifelse(treat_type == "never", 1, 0),
           post = ifelse(stacked_event_time >= 0, 1, 0),
           post_treat = post * treated,
           stacked_event_time_treat = stacked_event_time * treated,
           all_treated = sum(treated),
           all_control = sum(control_any),
           all_control_notyet = sum(control_notyet)) %>%
    group_by(stack_group) %>%
    mutate(group_treated = sum(treated),
           group_control_any = sum(control_any),
           group_control_notyet = sum(control_notyet)) %>%
    ungroup() %>%
    mutate(weight_any = case_when(
             treat_type == "treated" ~ 1,
             treat_type != "treated" ~ (group_treated / all_treated) / (group_control_any / all_control)),
           weight_notyet = case_when(
             treat_type == "treated" ~ 1,
             treat_type == "notyet" ~ (group_treated / all_treated) / (group_control_notyet / all_control_notyet))
    )

  stack.hosp
}


# Build stacked dataset at state level ---------------------------------------------

stack_state <- function(pre.period, post.period, state.period) {

  ## loop over all possible values of state_treat_year (year of CAH availability)
  stack.state <- tibble()
  state.count <- est.dat %>%
    group_by(MSTATE, year) %>%
    summarize(hospitals = n(), .groups = "drop_last") %>%
    group_by(MSTATE) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(hospitals_lag = lag(hospitals)) %>%
    ungroup() %>%
    select(MSTATE, year, hospitals, hospitals_lag)

  for (i in unique(est.dat$state_treat_year)) {

    ## define treated group relative to state_treat_year i within event window
    treat.dat <- est.dat %>%
      filter(state_treat_year == i,
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "treated")

    ## define control group relative to state_treat_year i within event window
    ## never treated
    control.dat.never <- est.dat %>%
      filter(state_treat_year == 0,
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "never")

    ## not yet treated
    control.dat.notyet <- est.dat %>%
      filter(state_treat_year > (i + state.period),
             year >= (i - pre.period), year <= (i + post.period)) %>%
      select(ID, year) %>%
      mutate(treat_type = "notyet")

    ## inner join back to est.dat
    stack.dat.group <- est.dat %>%
      inner_join(bind_rows(treat.dat, control.dat.never, control.dat.notyet), by = c("ID", "year")) %>%
      left_join(state.count, by = c("MSTATE", "year")) %>%
      mutate(state = as.numeric(factor(MSTATE)),
             stacked_event_time = year - i,
             stack_group = i)

    ## collapse to state-year level
    stack.dat <- stack.dat.group %>%
      group_by(MSTATE, year) %>%
      summarize(changes = sum(closed_merged),
                sum_cah = sum(cah, na.rm = TRUE),
                group_type = first(treat_type),
                closures = sum(closed),
                mergers = sum(merged),
                hospitals = first(hospitals),
                hospitals_lag = first(hospitals_lag),
                .groups = "drop") %>%
      arrange(MSTATE, year) %>%
      group_by(MSTATE) %>%
      mutate(cumul_changes = cumsum(changes),
             rate_changes = changes / hospitals_lag,
             rate_closed  = closures / hospitals_lag,
             rate_merged  = mergers / hospitals_lag) %>%
      mutate(state = as.numeric(factor(MSTATE)),
             stacked_event_time = year - i,
             stack_group = i) %>%
      ungroup()

    stack.state <- bind_rows(stack.state, stack.dat)
  }

  ## drop the first_year_treat==0 phantom treatment group
  stack.state <- stack.state %>% ungroup() %>%
    filter(stack_group > 0) %>%
    mutate(treated = ifelse(group_type == "treated", 1, 0),
           control_any = ifelse(group_type != "treated", 1, 0),
           control_notyet = ifelse(group_type == "notyet", 1, 0),
           control_never = ifelse(group_type == "never", 1, 0),
           post = ifelse(stacked_event_time >= 0, 1, 0),
           post_treat = post * treated,
           stacked_event_time_treat = stacked_event_time * treated,
           all_treated = sum(treated),
           all_control = sum(control_any),
           all_control_notyet = sum(control_notyet)) %>%
    group_by(stack_group) %>%
    mutate(group_treated = sum(treated),
           group_control_any = sum(control_any),
           group_control_notyet = sum(control_notyet)) %>%
    ungroup() %>%
    mutate(weight_any = case_when(
             group_type == "treated" ~ 1,
             group_type != "treated" ~ (group_treated / all_treated) / (group_control_any / all_control)),
           weight_notyet = case_when(
             group_type == "treated" ~ 1,
             group_type == "notyet" ~ (group_treated / all_treated) / (group_control_notyet / all_control_notyet))
    )

  stack.state
}
# Event-time effects shared by the SDID scripts --------------------------------
# sdid_event(): per-period SDID effects for one cohort's balanced panel,
#   tau_t = (treated mean_t - omega'Y_t) - sum_s lambda_s (treated mean_s - omega'Y_s),
# whose post-period mean is the SDID estimate. Standard errors and the
# covariance across periods come from the jackknife synthdid itself uses: leave
# one unit out with the weights held fixed and omega renormalized over the
# remaining controls. Standard errors are NA in the two cases where that
# jackknife is undefined (one treated unit, one nonzero omega).
sdid_event <- function(est, setup, cohort) {
  w     <- attr(est, "weights")
  omega <- w$omega
  lambda <- w$lambda
  Y  <- setup$Y
  N0 <- setup$N0
  T0 <- setup$T0
  N  <- nrow(Y)
  e  <- as.numeric(colnames(Y)) - cohort

  per_period <- function(Y, N0, omega) {
    gap <- colMeans(Y[(N0 + 1):nrow(Y), , drop = FALSE]) -
           as.numeric(t(omega) %*% Y[1:N0, , drop = FALSE])
    gap - sum(lambda * gap[1:T0])
  }

  tau <- per_period(Y, N0, omega)
  if (N - N0 < 2 || sum(omega > 0) < 2) {
    return(list(es = tibble(cohort = cohort, e = e, estimate = tau, se = NA_real_, Ntr = N - N0),
                vcov = matrix(NA_real_, length(e), length(e))))
  }
  jk <- sapply(1:N, function(i) {
    if (i <= N0) {
      om <- omega[-i]
      om <- om / sum(om)
      per_period(Y[-i, , drop = FALSE], N0 - 1, om)
    } else {
      per_period(Y[-i, , drop = FALSE], N0, omega)
    }
  })
  dev <- jk - rowMeans(jk)
  V   <- (N - 1) / N * dev %*% t(dev)
  list(es = tibble(cohort = cohort, e = e, estimate = tau, se = sqrt(diag(V)), Ntr = N - N0),
       vcov = V)
}

# cs_cohort(): Callaway and Sant'Anna on one cohort's balanced panel, the same
# panel the SDID estimate uses. The stack has already fixed the treated and
# control hospitals, so the controls enter as the never-treated group of a
# single-group att_gt call: universal base period, no covariates, dynamic
# effects over the same window as SDID. Seeded per call so the bootstrap
# intervals do not depend on what ran earlier. `bal` carries ID, year, outcome,
# treated. Returns NULL when the estimator fails on the cohort.
cs_cohort <- function(bal, cohort, min.es, max.es, seed = 1234) {
  d <- bal %>%
    mutate(treat_group = ifelse(treated == 1, cohort, 0),
           ID2 = as.numeric(factor(ID))) %>%
    select(ID2, treat_group, year, outcome)
  set.seed(seed)
  raw <- tryCatch(
    att_gt(yname = "outcome", gname = "treat_group", idname = "ID2", tname = "year",
           control_group = "nevertreated", panel = TRUE, allow_unbalanced_panel = TRUE,
           data = as.data.frame(d), xformla = ~1,
           base_period = "universal", est_method = "reg"),
    error = function(err) {
      message(sprintf("cs_cohort: cohort %d failed (%s)", cohort, conditionMessage(err)))
      NULL })
  if (is.null(raw)) return(NULL)
  dyn <- tryCatch(aggte(raw, type = "dynamic", na.rm = TRUE, min_e = min.es, max_e = max.es),
                  error = function(err) NULL)
  if (is.null(dyn)) return(NULL)
  ## did reports no standard error when the bootstrap is degenerate (the treated
  ## outcome does not vary after treatment); the cohort then has no CS estimate.
  if (is.na(dyn$overall.se)) {
    message(sprintf("cs_cohort: cohort %d dropped from CS, degenerate bootstrap", cohort))
    return(NULL)
  }
  ntr <- n_distinct(d$ID2[d$treat_group > 0])
  inf <- dyn$inf.function$dynamic.inf.func.e
  V   <- if (is.null(inf) || ncol(inf) != length(dyn$egt)) {
    matrix(NA_real_, length(dyn$egt), length(dyn$egt))
  } else {
    t(inf) %*% inf / nrow(inf)^2
  }
  list(att = dyn$overall.att, se = dyn$overall.se,
       Ntr = ntr, N0 = n_distinct(d$ID2[d$treat_group == 0]),
       es  = tibble(cohort = cohort, e = dyn$egt, estimate = dyn$att, se = dyn$se, Ntr = ntr),
       vcov = V)
}

# pool_event(): pool cohort event-time series with treated-unit-count weights,
# the weighting the pooled ATTs use, treating cohorts as independent.
pool_event <- function(x) {
  x %>%
    group_by(e) %>%
    summarise(estimate = sum(Ntr * estimate) / sum(Ntr),
              se       = sqrt(sum(Ntr^2 * se^2)) / sum(Ntr),
              Ntr      = sum(Ntr),
              .groups  = "drop")
}

# vcov_long(): a covariance matrix as rows of (cohort, e1, e2, value), so the
# cohort covariances can go to CSV and be recombined by the sensitivity analysis.
vcov_long <- function(V, cohort, e) {
  expand_grid(i = seq_along(e), j = seq_along(e)) %>%
    transmute(cohort = cohort, e1 = e[i], e2 = e[j], value = V[cbind(i, j)])
}

# overlay_plot(): SDID (circles) and CS (triangles) event-time effects with 95%
# intervals on the same axes. `series` carries estimator, e, estimate, se.
overlay_plot <- function(series, ylab) {
  d <- series %>% mutate(estimator = factor(estimator, levels = c("SDID", "CS")))
  ggplot(d, aes(x = e, y = estimate, shape = estimator, linetype = estimator)) +
    geom_hline(yintercept = 0, linewidth = 0.8) +
    geom_vline(xintercept = -0.5, linetype = "dashed") +
    geom_errorbar(aes(ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se),
                  position = position_dodge(width = 0.5), width = 0, linewidth = 0.5, alpha = 0.5) +
    geom_point(position = position_dodge(width = 0.5), size = 2.6, fill = "white") +
    scale_shape_manual(values = c(SDID = 16, CS = 24)) +
    scale_linetype_manual(values = c(SDID = "solid", CS = "dashed")) +
    scale_x_continuous(breaks = seq(min(d$e), max(d$e), 1)) +
    labs(x = "Event time", y = ylab, shape = NULL, linetype = NULL) +
    theme_bw() +
    theme(legend.position = "inside", legend.position.inside = c(0.14, 0.14),
          legend.background = element_blank())
}
