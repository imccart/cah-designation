# app-statecut-sensitivity.R
# Robustness to varying state-level CAH adoption lags (state.cut parameter)
# Runs SDID across state.cut = 1, 2, 3 for all hospital + state outcomes
# Baseline (state.cut = 0) results taken from hosp.results.table and state.results.table
# Expects from _run-analysis.r: est.dat, state.dat, bed.cut, post, financial.pre,
#   hosp.results.table, state.results.table


# Outcome maps (matching 2-hospital-dd.R and 4-changes-state-dd.R) ----------
sc_hosp_outcomes <- list(
  margin            = list(label="Operating margin",             pre_period=financial.pre),
  current_ratio     = list(label="Current ratio",                pre_period=financial.pre),
  net_fixed         = list(label="Net fixed assets",             pre_period=financial.pre),
  capex             = list(label="Capital expenditures per bed", pre_period=financial.pre),
  net_pat_rev       = list(label="Net patient revenue per bed",  pre_period=financial.pre),
  tot_operating_exp = list(label="Operating expenses per bed",   pre_period=financial.pre),
  BDTOT             = list(label="Total beds"),
  OBBD              = list(label="OB beds"),
  FTERN             = list(label="FTE RNs"),
  ip_per_bed        = list(label="Inpatient days per bed"),
  system            = list(label="System membership")
)

sc_state_outcomes <- list(
  closures = list(label="Closures"),
  mergers  = list(label="Mergers")
)

state_cuts    <- 1:2
cohorts_hosp  <- 1999:2001
cohorts_state <- 1999:2001


# Baseline results (state.cut = 0, already computed) ------------------------
sc_baseline <- bind_rows(
  hosp.results.table %>%
    transmute(outcome, att = sdid_att, ci_low = sdid_ci_low,
              ci_high = sdid_ci_high, state_cut = 0L),
  state.results.table %>%
    transmute(outcome, att = sdid_att, ci_low = sdid_ci_low,
              ci_high = sdid_ci_high, state_cut = 0L)
)


# SDID helper: hospital-level cohort ----------------------------------------
run_sdid_hosp_sc <- function(c, outcome_sym, pp, sh) {
  synth.c <- sh %>%
    filter(stack_group == c) %>%
    group_by(ID) %>%
    mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(!is.na(!!outcome_sym),
           min_bedsize <= bed.cut,
           stacked_event_time >= -pp) %>%
    select(ID, year, outcome = !!outcome_sym, post_treat)

  bal.c <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))
  setup <- panel.matrices(as.data.frame(bal.c))
  Ntr   <- nrow(setup$Y) - setup$N0
  N0    <- setup$N0

  if (setup$T0 < 2 || Ntr < 1 || N0 < 5) return(NULL)

  est  <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
  se_c <- synthdid_se(est, method = "jackknife")

  tibble(cohort = c, att = as.numeric(est), se = as.numeric(se_c), Ntr = Ntr)
}


# SDID helper: state-level cohort ------------------------------------------
run_sdid_state_sc <- function(c, outcome_var, ss) {
  denom_c <- ss %>%
    filter(stacked_event_time <= -1, treated == 1, stack_group == c) %>%
    summarize(mean_hosp = mean(hospitals, na.rm = TRUE)) %>%
    pull(mean_hosp)

  dat <- ss %>%
    filter(stack_group == c) %>%
    transmute(ID = as.numeric(factor(MSTATE)),
              year,
              rate = .data[[outcome_var]],
              post_treat)

  bal <- as_tibble(makeBalancedPanel(dat, idname = "ID", tname = "year"))
  setup <- panel.matrices(as.data.frame(bal))
  Ntr   <- nrow(setup$Y) - setup$N0

  if (setup$T0 < 2 || Ntr < 1) return(NULL)

  est <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
  se  <- synthdid_se(est, method = "jackknife")

  tibble(cohort = c,
         att = (as.numeric(est) / denom_c) * 100,
         se  = (as.numeric(se) / denom_c) * 100,
         Ntr = Ntr)
}


# Main loop: state.cut = 1, 2, 3 -------------------------------------------
sc_sensitivity <- tibble(
  outcome = character(), att = numeric(), ci_low = numeric(),
  ci_high = numeric(), state_cut = integer()
)

for (sc in state_cuts) {

  # Rebuild stacked datasets with new state.period
  sh <- stack_hosp(pre.period = 5, post.period = post, state.period = sc)
  ss <- stack_state(pre.period = 5, post.period = post, state.period = sc)

  # Hospital outcomes
  for (oname in names(sc_hosp_outcomes)) {
    o   <- sc_hosp_outcomes[[oname]]
    osym <- sym(oname)
    pp  <- if (!is.null(o$pre_period)) o$pre_period else 5


    atts_all <- bind_rows(map(cohorts_hosp, function(c) {
      tryCatch(
        run_sdid_hosp_sc(c, osym, pp, sh),
        error = function(e) {
          message(sprintf("statecut sensitivity: cohort %d failed", c))
          NULL
        }
      )
    }))

    if (nrow(atts_all) == 0) {
      sc_sensitivity <- bind_rows(sc_sensitivity, tibble(
        outcome = o$label, att = NA_real_, ci_low = NA_real_,
        ci_high = NA_real_, state_cut = as.integer(sc)
      ))
      next
    }

    att_w  <- with(atts_all, sum(Ntr * att) / sum(Ntr))
    se_w   <- with(atts_all, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))


    sc_sensitivity <- bind_rows(sc_sensitivity, tibble(
      outcome   = o$label,
      att       = att_w,
      ci_low    = att_w - 1.96 * se_w,
      ci_high   = att_w + 1.96 * se_w,
      state_cut = as.integer(sc)
    ))
  }

  # State outcomes
  for (oname in names(sc_state_outcomes)) {
    o <- sc_state_outcomes[[oname]]


    atts_all <- bind_rows(map(cohorts_state, function(c) {
      tryCatch(
        run_sdid_state_sc(c, oname, ss),
        error = function(e) {
          message(sprintf("statecut sensitivity: cohort %d failed", c))
          NULL
        }
      )
    }))

    if (nrow(atts_all) == 0) {
      sc_sensitivity <- bind_rows(sc_sensitivity, tibble(
        outcome = o$label, att = NA_real_, ci_low = NA_real_,
        ci_high = NA_real_, state_cut = as.integer(sc)
      ))
      next
    }

    att_w  <- with(atts_all, sum(Ntr * att) / sum(Ntr))
    se_w   <- with(atts_all, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))


    sc_sensitivity <- bind_rows(sc_sensitivity, tibble(
      outcome   = o$label,
      att       = att_w,
      ci_low    = att_w - 1.96 * se_w,
      ci_high   = att_w + 1.96 * se_w,
      state_cut = as.integer(sc)
    ))
  }
}


# Forest plot ---------------------------------------------------------------
sc_forest_labels <- c(
  "Operating margin"             = "Operating\nmargin",
  "Current ratio"                = "Current ratio",
  "Net fixed assets"             = "Net fixed assets",
  "Capital expenditures per bed" = "Capital expenditures\nper bed",
  "Net patient revenue per bed"  = "Net patient revenue\nper bed",
  "Operating expenses per bed"   = "Operating expenses\nper bed",
  "Total beds"                   = "Total beds",
  "OB beds"                      = "OB beds",
  "FTE RNs"                      = "FTE RNs",
  "Inpatient days per bed"       = "Inpatient days\nper bed",
  "System membership"            = "System\nmembership",
  "Closures"                     = "Closures",
  "Mergers"                      = "Mergers"
)

outcome_order <- c(
  sapply(sc_hosp_outcomes, function(o) o$label),
  sapply(sc_state_outcomes, function(o) o$label)
)

label_levels <- sc_forest_labels[outcome_order]

# Baseline band data (one row per outcome for geom_rect)
baseline_plot <- sc_baseline %>%
  filter(!is.na(att)) %>%
  mutate(outcome_label = sc_forest_labels[outcome],
         outcome_label = factor(outcome_label, levels = label_levels))

# Sensitivity points (sc = 1, 2, 3)
sens_plot <- sc_sensitivity %>%
  filter(!is.na(att)) %>%
  mutate(outcome_label = sc_forest_labels[outcome],
         outcome_label = factor(outcome_label, levels = label_levels),
         sc_label = factor(state_cut))

n_out <- length(levels(droplevels(sens_plot$outcome_label)))

p_sc <- ggplot() +
  # Shaded baseline CI band
  geom_rect(data = baseline_plot,
            aes(xmin = ci_low, xmax = ci_high, ymin = -Inf, ymax = Inf),
            fill = "gray85", alpha = 0.5) +
  # Baseline ATT vertical line
  geom_vline(data = baseline_plot,
             aes(xintercept = att),
             linetype = "dashed", linewidth = 0.5, color = "gray40") +
  # Zero line
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  # Sensitivity points + CIs
  geom_pointrange(data = sens_plot,
                  aes(x = att, y = sc_label, xmin = ci_low, xmax = ci_high,
                      shape = sc_label),
                  size = 0.5, linewidth = 0.5) +
  scale_shape_manual(values = c("1" = 16, "2" = 15),
                     labels = c("1" = "1 year", "2" = "2 years")) +
  facet_wrap(~ outcome_label, ncol = 1, scales = "free_x",
             strip.position = "left", dir = "v") +
  labs(x = "ATT (95% CI)", y = "State adoption lag", shape = "Lag") +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.placement  = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9,
                                      lineheight = 0.9),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

ggsave("results/statecut-sensitivity.png", p_sc,
       width = 7, height = 0.9 * n_out + 1.2, dpi = 300)


# Save CSV ------------------------------------------------------------------
write_csv(bind_rows(sc_baseline, sc_sensitivity),
          "results/diagnostics/statecut-sensitivity.csv")

