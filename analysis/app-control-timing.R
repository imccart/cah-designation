# app-control-timing.R
# Sensitivity to control group timing concerns (Carroll comment)
# Part A: Drop "immediate adopters" (eff_year == state_treat_year) from the sample, hospital-level outcomes
# Part B: Vary post-period length (5, 4, 3, 2) to tighten control pool
# Expects from _run-analysis.r: est.dat, state.dat, bed.cut, post, state.cut,
#   financial.pre, hosp.results.table, state.results.table


# Outcome maps (matching 2-hospital-dd.R and 4-changes-state-dd.R) ----------
ct_hosp_outcomes <- list(
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
  system            = list(label="System membership"),
  ## AHA visits and births
  ed_per_bed        = list(label="ED visits per bed"),
  op_per_bed        = list(label="Outpatient visits per bed"),
  opnon_per_bed     = list(label="Non-ED outpatient visits per bed"),
  BIRTHS            = list(label="Births")
)

ct_state_outcomes <- list(
  closures = list(label="Closures"),
  mergers  = list(label="Mergers")
)

cohorts_hosp  <- 1999:2001
cohorts_state <- 1999:2001


# Baseline results (from main analysis) ------------------------------------
ct_baseline <- bind_rows(
  hosp.results.table %>%
    transmute(outcome, att = sdid_att, ci_low = sdid_ci_low,
              ci_high = sdid_ci_high),
  state.results.table %>%
    transmute(outcome, att = sdid_att, ci_low = sdid_ci_low,
              ci_high = sdid_ci_high)
)


# SDID helper: hospital-level cohort ----------------------------------------
run_sdid_hosp_ct <- function(c, outcome_sym, pp, sh) {
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
# Clean-donor SDID as in 4-changes-state-dd.R; `timing` is the adoption-year
# vector (Part C drops the never-adopting states from it) and `pp_post` the
# post-period length. NULL when a post year has fewer than two donors.
run_sdid_state_ct <- function(c, outcome_var, timing, pp_post) {
  fit <- sdid_state_clean(state.dat, timing, c, outcome_var, pre = 5, post = pp_post, lag = state.cut)
  if (is.null(fit)) return(NULL)
  tibble(cohort = c, att = fit$att, se = fit$se, Ntr = fit$Ntr)
}


# Pooling helper -----------------------------------------------------------
pool_atts <- function(atts_all) {
  if (is.null(atts_all) || nrow(atts_all) == 0) return(NULL)
  atts_all <- atts_all %>% filter(!is.na(se))
  if (nrow(atts_all) == 0) return(NULL)
  att_w  <- with(atts_all, sum(Ntr * att) / sum(Ntr))
  se_w   <- with(atts_all, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))
  tibble(att = att_w, ci_low = att_w - 1.96 * se_w, ci_high = att_w + 1.96 * se_w)
}


# =========================================================================
# PART A: Drop immediate adopters
# =========================================================================


# Count how many treated hospitals are immediate adopters
n_immediate <- est.dat %>%
  filter(!is.na(eff_year), eff_year == state_treat_year) %>%
  summarize(n = n_distinct(ID)) %>%
  pull(n)
n_treated_total <- est.dat %>%
  filter(!is.na(eff_year)) %>%
  summarize(n = n_distinct(ID)) %>%
  pull(n)

# Drop immediate adopters entirely, so they are neither treated nor controls
immediate_ids <- est.dat %>%
  filter(!is.na(eff_year), eff_year == state_treat_year) %>%
  distinct(ID) %>%
  pull(ID)
est.dat.noimmed <- est.dat %>%
  filter(!ID %in% immediate_ids)

# Temporarily swap est.dat for stacking (stack_hosp reads from global)
est.dat.orig <- est.dat
est.dat <- est.dat.noimmed

sh_noimmed <- stack_hosp(pre.period = 5, post.period = post, state.period = state.cut)

est.dat <- est.dat.orig

partA_results <- tibble(
  outcome = character(), att = numeric(),
  ci_low = numeric(), ci_high = numeric(),
  spec = character()
)

# Hospital outcomes
for (oname in names(ct_hosp_outcomes)) {
  o   <- ct_hosp_outcomes[[oname]]
  osym <- sym(oname)
  pp  <- if (!is.null(o$pre_period)) o$pre_period else 5


  atts_all <- bind_rows(map(cohorts_hosp, function(c) {
    tryCatch(
      run_sdid_hosp_ct(c, osym, pp, sh_noimmed),
      error = function(e) { NULL }
    )
  }))

  pooled <- pool_atts(atts_all)
  if (is.null(pooled)) {
    partA_results <- bind_rows(partA_results, tibble(
      outcome = o$label, att = NA_real_, ci_low = NA_real_,
      ci_high = NA_real_, spec = "no_immediate"))
    next
  }

  partA_results <- bind_rows(partA_results, tibble(
    outcome = o$label, att = pooled$att, ci_low = pooled$ci_low,
    ci_high = pooled$ci_high, spec = "no_immediate"))
}

# State outcomes are not part of this check: state treatment is dated by the
# state's first designation, not by which hospitals convert.

rm(sh_noimmed, est.dat.noimmed); gc()


# =========================================================================
# PART B: Post-period sensitivity
# =========================================================================


post_periods <- c(5, 4, 3, 2)

partB_results <- tibble(
  outcome = character(), att = numeric(),
  ci_low = numeric(), ci_high = numeric(),
  post_period = integer()
)

for (pp_post in post_periods) {

  sh <- stack_hosp(pre.period = 5, post.period = pp_post, state.period = state.cut)

  # Hospital outcomes
  for (oname in names(ct_hosp_outcomes)) {
    o   <- ct_hosp_outcomes[[oname]]
    osym <- sym(oname)
    pp  <- if (!is.null(o$pre_period)) o$pre_period else 5


    atts_all <- bind_rows(map(cohorts_hosp, function(c) {
      tryCatch(
        run_sdid_hosp_ct(c, osym, pp, sh),
        error = function(e) { NULL }
      )
    }))

    pooled <- pool_atts(atts_all)
    if (is.null(pooled)) {
      partB_results <- bind_rows(partB_results, tibble(
        outcome = o$label, att = NA_real_, ci_low = NA_real_,
        ci_high = NA_real_, post_period = as.integer(pp_post)))
      next
    }

    partB_results <- bind_rows(partB_results, tibble(
      outcome = o$label, att = pooled$att, ci_low = pooled$ci_low,
      ci_high = pooled$ci_high, post_period = as.integer(pp_post)))
  }

  # State outcomes
  for (oname in names(ct_state_outcomes)) {
    o <- ct_state_outcomes[[oname]]


    atts_all <- bind_rows(map(cohorts_state, function(c) {
      tryCatch(
        run_sdid_state_ct(c, oname, state.timing, pp_post),
        error = function(e) { NULL }
      )
    }))

    pooled <- pool_atts(atts_all)
    if (is.null(pooled)) {
      partB_results <- bind_rows(partB_results, tibble(
        outcome = o$label, att = NA_real_, ci_low = NA_real_,
        ci_high = NA_real_, post_period = as.integer(pp_post)))
      next
    }

    partB_results <- bind_rows(partB_results, tibble(
      outcome = o$label, att = pooled$att, ci_low = pooled$ci_low,
      ci_high = pooled$ci_high, post_period = as.integer(pp_post)))
  }

  rm(sh); gc()
}


# =========================================================================
# Forest plot: Part A (drop immediate adopters)
# =========================================================================

ct_forest_labels <- c(
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
  "ED visits per bed"            = "ED visits\nper bed",
  "Outpatient visits per bed"    = "Outpatient visits\nper bed",
  "Non-ED outpatient visits per bed"= "Non-ED outpatient\nvisits per bed",
  "Births"                       = "Births",
  "Closures"                     = "Closures",
  "Mergers"                      = "Mergers"
)

outcome_order <- c(
  sapply(ct_hosp_outcomes, function(o) o$label),
  sapply(ct_state_outcomes, function(o) o$label)
)

label_levels <- ct_forest_labels[outcome_order]

# Baseline band
baseline_plot <- ct_baseline %>%
  filter(!is.na(att), outcome %in% partA_results$outcome) %>%
  mutate(outcome_label = ct_forest_labels[outcome],
         outcome_label = factor(outcome_label, levels = label_levels))

# Part A points
partA_plot <- partA_results %>%
  filter(!is.na(att)) %>%
  mutate(outcome_label = ct_forest_labels[outcome],
         outcome_label = factor(outcome_label, levels = label_levels))

n_out_A <- length(levels(droplevels(partA_plot$outcome_label)))

p_partA <- ggplot() +
  geom_rect(data = baseline_plot,
            aes(xmin = ci_low, xmax = ci_high, ymin = -Inf, ymax = Inf),
            fill = "gray85", alpha = 0.5) +
  geom_vline(data = baseline_plot,
             aes(xintercept = att),
             linetype = "dashed", linewidth = 0.5, color = "gray40") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  geom_pointrange(data = partA_plot,
                  aes(x = att, y = "",
                      xmin = ci_low, xmax = ci_high),
                  size = 0.5, linewidth = 0.5) +
  facet_wrap(~ outcome_label, ncol = 1, scales = "free_x",
             strip.position = "left", dir = "v") +
  labs(x = "ATT (95% CI)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.placement  = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9,
                                      lineheight = 0.9),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank()
  )

ggsave("results/control-timing-noimmediate.png", p_partA,
       width = 7, height = 0.9 * n_out_A + 1.2, dpi = 300)


# =========================================================================
# Forest plot: Part B (post-period sensitivity)
# =========================================================================

partB_plot <- partB_results %>%
  filter(!is.na(att)) %>%
  mutate(outcome_label = ct_forest_labels[outcome],
         outcome_label = factor(outcome_label, levels = label_levels),
         pp_label = factor(post_period))

n_out_B <- length(levels(droplevels(partB_plot$outcome_label)))

p_partB <- ggplot() +
  geom_rect(data = baseline_plot,
            aes(xmin = ci_low, xmax = ci_high, ymin = -Inf, ymax = Inf),
            fill = "gray85", alpha = 0.5) +
  geom_vline(data = baseline_plot,
             aes(xintercept = att),
             linetype = "dashed", linewidth = 0.5, color = "gray40") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  geom_pointrange(data = partB_plot,
                  aes(x = att, y = pp_label, xmin = ci_low, xmax = ci_high,
                      shape = pp_label),
                  size = 0.5, linewidth = 0.5) +
  scale_shape_manual(values = c("5" = 16, "4" = 17, "3" = 15, "2" = 18),
                     labels = c("5" = "5 years", "4" = "4 years",
                                "3" = "3 years", "2" = "2 years")) +
  facet_wrap(~ outcome_label, ncol = 1, scales = "free_x",
             strip.position = "left", dir = "v") +
  labs(x = "ATT (95% CI)", y = "Post-period length", shape = "Post") +
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

ggsave("results/control-timing-postperiod.png", p_partB,
       width = 7, height = 0.9 * n_out_B + 1.2, dpi = 300)


# =========================================================================
# Save CSVs
# =========================================================================

write_csv(partA_results, "results/diagnostics/control-timing-noimmediate.csv")
write_csv(partB_results, "results/diagnostics/control-timing-postperiod.csv")


# =========================================================================
# PART C: Drop never-treated states
# =========================================================================


# Count never-treated states and their hospitals
never_states <- est.dat %>%
  filter(state_treat_year == 0) %>%
  summarize(n_states = n_distinct(MSTATE),
            n_hosp = n_distinct(ID),
            states = paste(sort(unique(MSTATE)), collapse = ", "))

# Count how many of those pass the bed-size filter
never_small <- est.dat %>%
  filter(state_treat_year == 0) %>%
  group_by(ID) %>%
  summarize(min_beds = min(BDTOT, na.rm = TRUE)) %>%
  filter(min_beds <= bed.cut) %>%
  nrow()

## The appendix quotes the immediate-adopter share and the number of small
## hospitals in never-treated states, so both go to a file.
tibble(n_immediate = n_immediate,
       n_treated_total = n_treated_total,
       share_immediate = n_immediate / n_treated_total,
       n_never_states = never_states$n_states,
       n_never_hosp = never_states$n_hosp,
       n_never_small = never_small,
       never_state_names = never_states$states) %>%
  write_csv("results/diagnostics/control-timing-counts.csv")

# Create modified est.dat excluding never-treated states
est.dat.noNever <- est.dat %>% filter(state_treat_year > 0)
state.dat.noNever <- state.dat %>% filter(state_treat_year > 0)

# Temporarily swap est.dat for stacking
est.dat.orig <- est.dat
est.dat <- est.dat.noNever

sh_noNever <- stack_hosp(pre.period = 5, post.period = post, state.period = state.cut)

est.dat <- est.dat.orig

partC_results <- tibble(
  outcome = character(), att = numeric(),
  ci_low = numeric(), ci_high = numeric(),
  spec = character()
)

# Hospital outcomes
for (oname in names(ct_hosp_outcomes)) {
  o   <- ct_hosp_outcomes[[oname]]
  osym <- sym(oname)
  pp  <- if (!is.null(o$pre_period)) o$pre_period else 5


  atts_all <- bind_rows(map(cohorts_hosp, function(c) {
    tryCatch(
      run_sdid_hosp_ct(c, osym, pp, sh_noNever),
      error = function(e) { NULL }
    )
  }))

  pooled <- pool_atts(atts_all)
  if (is.null(pooled)) {
    partC_results <- bind_rows(partC_results, tibble(
      outcome = o$label, att = NA_real_, ci_low = NA_real_,
      ci_high = NA_real_, spec = "no_never_states"))
    next
  }

  partC_results <- bind_rows(partC_results, tibble(
    outcome = o$label, att = pooled$att, ci_low = pooled$ci_low,
    ci_high = pooled$ci_high, spec = "no_never_states"))
}

# State outcomes are not part of this check: without the never-adopting states,
# no state is still untreated after 2002 to serve as a donor.

rm(sh_noNever, est.dat.noNever, state.dat.noNever); gc()


# =========================================================================
# Forest plot: Part C (drop never-treated states)
# =========================================================================

partC_plot <- partC_results %>%
  filter(!is.na(att)) %>%
  mutate(outcome_label = ct_forest_labels[outcome],
         outcome_label = factor(outcome_label, levels = label_levels))

n_out_C <- length(levels(droplevels(partC_plot$outcome_label)))

p_partC <- ggplot() +
  geom_rect(data = baseline_plot,
            aes(xmin = ci_low, xmax = ci_high, ymin = -Inf, ymax = Inf),
            fill = "gray85", alpha = 0.5) +
  geom_vline(data = baseline_plot,
             aes(xintercept = att),
             linetype = "dashed", linewidth = 0.5, color = "gray40") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  geom_pointrange(data = partC_plot,
                  aes(x = att, y = "",
                      xmin = ci_low, xmax = ci_high),
                  size = 0.5, linewidth = 0.5) +
  facet_wrap(~ outcome_label, ncol = 1, scales = "free_x",
             strip.position = "left", dir = "v") +
  labs(x = "ATT (95% CI)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.placement  = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9,
                                      lineheight = 0.9),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank()
  )

ggsave("results/control-timing-nonever.png", p_partC,
       width = 7, height = 0.9 * n_out_C + 1.2, dpi = 300)

write_csv(partC_results, "results/diagnostics/control-timing-nonever.csv")


