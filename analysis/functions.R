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

# cs_state_cohort(): Callaway--Sant'Anna for one state-level cohort on its SDID
# sample, the balanced stack.state panel for cohort c (its treated states and
# SDID's control states). Control states keep their own adoption year, so each
# is a not-yet-treated comparison until it adopts. Only cohort c's group-time
# effects are kept: the event-time series, and the ATT as their mean over event
# times 0 to post, with an analytic standard error from the influence function.
cs_state_cohort <- function(stack, cohort, outcome_var, adopt, post) {
  bal <- stack %>%
    filter(stack_group == cohort) %>%
    transmute(MSTATE, year, outcome = .data[[outcome_var]], treated) %>%
    makeBalancedPanel(idname = "MSTATE", tname = "year") %>%
    as_tibble() %>%
    left_join(adopt, by = "MSTATE", relationship = "many-to-one") %>%
    mutate(g = if_else(treated == 1, cohort, coalesce(state_treat_year, 0)),
           ID2 = as.numeric(factor(MSTATE)))
  raw <- att_gt(yname = "outcome", gname = "g", idname = "ID2", tname = "year",
                control_group = "notyettreated", panel = TRUE,
                data = as.data.frame(bal), xformla = ~1,
                base_period = "universal", est_method = "reg", bstrap = FALSE)
  keep   <- which(raw$group == cohort)
  e      <- raw$t[keep] - cohort
  post_k <- which(e >= 0 & e <= post)
  inf    <- rowMeans(as.matrix(raw$inffunc[, keep[post_k], drop = FALSE]))
  inf_all <- as.matrix(raw$inffunc[, keep, drop = FALSE])
  list(att = mean(raw$att[keep][post_k]),
       se  = sqrt(sum(inf^2)) / length(inf),
       Ntr = n_distinct(bal$ID2[bal$treated == 1]),
       N0  = n_distinct(bal$ID2[bal$treated == 0]),
       es  = tibble(cohort = cohort, e = e, estimate = raw$att[keep], se = raw$se[keep]),
       vcov = t(inf_all) %*% inf_all / nrow(inf_all)^2)
}

# sdid_state_clean(): state-level SDID for one cohort with donors restricted to
# states still untreated in each post year. The cohort's panel is its treated
# states (adoption year equal to the cohort) and the states adopting later or
# never, over the pre and post window, balanced. For each post year t, SDID is
# fit on the pre years plus year t, using as donors the states whose program
# starts after t or never. The cohort ATT is the mean of these yearly effects,
# and `scale` converts counts to per 100 hospitals using the treated states'
# pre-period hospital count. The pre-period event-time effects come from the
# adoption-year fit (t = cohort), whose donors are all untreated in that year.
# Standard errors and the covariance across event times are a fixed-weight
# jackknife over all the fits at once (each fit's weights re-averaged when a
# state is left out, as in synthdid_se). `timing` is each state's
# adoption year (0 = never), named by state; placebo draws pass a reshuffled
# version. `lag` requires control states to adopt more than `lag` years after the
# cohort (the state-lag sensitivity). Returns NULL when a year has fewer than two
# donors.
sdid_state_clean <- function(sdat, timing, cohort, outcome_var, pre, post, jack = TRUE, lag = 0) {
  st <- names(timing)[timing == cohort | timing == 0 | timing > cohort + lag]
  d <- sdat %>% filter(MSTATE %in% st, year >= cohort - pre, year <= cohort + post)
  Yw <- d %>% select(MSTATE, year, v = all_of(outcome_var)) %>%
    pivot_wider(names_from = year, values_from = v) %>% arrange(MSTATE)
  Yw <- Yw[complete.cases(Yw), ]
  Y <- as.matrix(Yw[, -1]); years <- as.numeric(colnames(Yw)[-1]); g <- timing[Yw$MSTATE]
  tr <- which(g == cohort); pre_c <- which(years < cohort)
  if (length(tr) < 1) return(NULL)
  scale <- 100 / (d %>% filter(MSTATE %in% Yw$MSTATE[tr], year < cohort) %>% summarize(m = mean(hospitals)) %>% pull(m))
  ts <- cohort:(cohort + post)
  # leave-one-state-out replicate of a fit's effects, weights held fixed and
  # omega renormalized over the remaining donors (NA where undefined)
  jack_fit <- function(f, trs, dns, om) {
    sapply(seq_len(nrow(Y)), function(i) {
      if (i %in% trs) {
        if (length(trs) < 2) return(rep(NA_real_, length(f(trs, dns, om))))
        f(setdiff(trs, i), dns, om)
      } else if (i %in% dns) {
        k <- which(dns == i); o <- om[-k]
        if (sum(o) <= 0) return(rep(NA_real_, length(f(trs, dns, om))))
        f(trs, dns[-k], o / sum(o))
      } else f(trs, dns, om)
    })
  }
  fits <- lapply(ts, function(t) {
    tc <- which(years == t)
    dn <- which(g == 0 | (g > t & g > cohort + lag))
    if (length(dn) < 2) return(NULL)
    cols <- c(pre_c, tc)
    est <- tryCatch(synthdid_estimate(rbind(Y[dn, cols, drop = FALSE], Y[tr, cols, drop = FALSE]), length(dn), length(pre_c)),
                    error = function(e) NULL)
    if (is.null(est)) return(NULL)
    w <- attr(est, "weights")
    gaps <- function(trs, dns, om) {
      gap <- colMeans(Y[trs, cols, drop = FALSE]) - as.numeric(t(om) %*% Y[dns, cols, drop = FALSE])
      gap - sum(w$lambda * gap[seq_along(pre_c)])
    }
    tau_fixed <- function(trs, dns, om) { x <- gaps(trs, dns, om); x[length(x)] }
    pre_fixed <- function(trs, dns, om) gaps(trs, dns, om)[seq_along(pre_c)]
    list(att = as.numeric(est), donors = length(dn), weights = w$omega,
         synthetic = sum(w$omega * Y[dn, tc]),
         jk = if (jack) jack_fit(tau_fixed, tr, dn, w$omega) else NULL,
         # the adoption-year fit (t = cohort) also supplies the pre-period effects
         pre = if (t == cohort) pre_fixed(tr, dn, w$omega) else NULL,
         pre_jk = if (t == cohort && jack) matrix(jack_fit(pre_fixed, tr, dn, w$omega), nrow = length(pre_c)) else NULL,
         pre_synthetic = if (t == cohort) as.numeric(t(w$omega) %*% Y[dn, pre_c, drop = FALSE]) else NULL)
  })
  if (any(vapply(fits, is.null, TRUE))) return(NULL)
  tau <- vapply(fits, `[[`, 0, "att")
  e_pre <- years[pre_c] - cohort
  out <- list(att = mean(tau) * scale, Ntr = length(tr), scale = scale,
              donors = vapply(fits, `[[`, 0, "donors"),
              weights = lapply(fits, `[[`, "weights"),
              paths = tibble(cohort = cohort, tau = c(e_pre, ts - cohort),
                             treated = colMeans(Y[tr, c(pre_c, match(ts, years)), drop = FALSE]) * scale,
                             synthetic = c(fits[[1]]$pre_synthetic, vapply(fits, `[[`, 0, "synthetic")) * scale,
                             Ntr = length(tr)))
  if (jack) {
    J <- sapply(fits, `[[`, "jk")                      # states x post years
    R <- rbind(fits[[1]]$pre_jk, t(J)) * scale          # event times x states
    keep <- colSums(is.na(R)) == 0
    n <- sum(keep)
    dev <- R[, keep, drop = FALSE] - rowMeans(R[, keep, drop = FALSE])
    V <- (n - 1) / n * dev %*% t(dev)
    np <- length(pre_c)
    out$se <- sqrt(sum(V[-(1:np), -(1:np)]) / length(ts)^2)
    out$jk_dropped <- sum(!keep)
    out$es <- tibble(cohort = cohort, e = c(e_pre, ts - cohort),
                     estimate = c(fits[[1]]$pre, tau) * scale, se = sqrt(diag(V)), Ntr = length(tr))
    out$vcov <- V
  }
  out
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


## Table lines with \addlinespace between outcome groups (financial, capacity
## and staffing, visits and births, organizational). outcomes are the row
## labels in the order of lines, each a member of outcome.order.
with_group_breaks <- function(lines, outcomes) {
  grp <- findInterval(match(outcomes, outcome.order), c(1, 7, 11, 15))
  out <- character()
  for (i in seq_along(lines)) {
    if (i > 1 && grp[i] != grp[i - 1]) out <- c(out, "\\addlinespace")
    out <- c(out, lines[i])
  }
  out
}
