# app-bedcut-sensitivity.R
# Robustness to varying bed size cutoffs for treatment eligibility
# Runs SDID at bed.cut = 25 and 75; baseline (bed.cut = 50) from hosp.results.table
# Uses existing stack.hosp (state-timing, state.cut = 0) — no rebuild needed
# State outcomes excluded (bed.cut does not affect state-level analysis)
# Expects from _run-analysis.r: stack.hosp, bed.cut, post, financial.pre, hosp.results.table


# Outcome map (hospital-level only, matching 2-hospital-dd.R) ---------------
bc_hosp_outcomes <- list(
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

bed_cuts <- c(25, 75)
cohorts  <- 1999:2001


# Baseline results (bed.cut = 50, already computed) -------------------------
bc_baseline <- hosp.results.table %>%
  transmute(outcome, att = sdid_att, ci_low = sdid_ci_low,
            ci_high = sdid_ci_high, bed_cut = as.integer(bed.cut))


# SDID helper ---------------------------------------------------------------
run_sdid_hosp_bc <- function(c, outcome_sym, pp, bc) {
  synth.c <- stack.hosp %>%
    filter(stack_group == c) %>%
    group_by(ID) %>%
    mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(!is.na(!!outcome_sym),
           min_bedsize <= bc,
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


# Main loop: bed.cut = 25, 75 ----------------------------------------------
bc_sensitivity <- tibble(
  outcome = character(), att = numeric(), ci_low = numeric(),
  ci_high = numeric(), bed_cut = integer()
)

for (bc in bed_cuts) {

  for (oname in names(bc_hosp_outcomes)) {
    o    <- bc_hosp_outcomes[[oname]]
    osym <- sym(oname)
    pp   <- if (!is.null(o$pre_period)) o$pre_period else 5


    atts_all <- bind_rows(map(cohorts, function(c) {
      tryCatch(
        run_sdid_hosp_bc(c, osym, pp, bc),
        error = function(e) {
          message(sprintf("bedcut sensitivity: cohort %d failed", c))
          NULL
        }
      )
    }))

    if (nrow(atts_all) == 0) {
      bc_sensitivity <- bind_rows(bc_sensitivity, tibble(
        outcome = o$label, att = NA_real_, ci_low = NA_real_,
        ci_high = NA_real_, bed_cut = as.integer(bc)
      ))
      next
    }

    att_w  <- with(atts_all, sum(Ntr * att) / sum(Ntr))
    se_w   <- with(atts_all, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))


    bc_sensitivity <- bind_rows(bc_sensitivity, tibble(
      outcome   = o$label,
      att       = att_w,
      ci_low    = att_w - 1.96 * se_w,
      ci_high   = att_w + 1.96 * se_w,
      bed_cut   = as.integer(bc)
    ))
  }
}


# Forest plot ---------------------------------------------------------------
bc_forest_labels <- c(
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
  "System membership"            = "System\nmembership"
)

outcome_order <- sapply(bc_hosp_outcomes, function(o) o$label)
label_levels  <- bc_forest_labels[outcome_order]

# Baseline band data
baseline_plot <- bc_baseline %>%
  filter(!is.na(att)) %>%
  mutate(outcome_label = bc_forest_labels[outcome],
         outcome_label = factor(outcome_label, levels = label_levels))

# Sensitivity points (bc = 25, 75)
sens_plot <- bc_sensitivity %>%
  filter(!is.na(att)) %>%
  mutate(outcome_label = bc_forest_labels[outcome],
         outcome_label = factor(outcome_label, levels = label_levels),
         bc_label = factor(bed_cut))

n_out <- length(levels(droplevels(sens_plot$outcome_label)))

p_bc <- ggplot() +
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
                  aes(x = att, y = bc_label, xmin = ci_low, xmax = ci_high,
                      shape = bc_label),
                  size = 0.5, linewidth = 0.5) +
  scale_shape_manual(values = c("25" = 16, "75" = 15),
                     labels = c("25" = "25 beds", "75" = "75 beds")) +
  facet_wrap(~ outcome_label, ncol = 1, scales = "free_x",
             strip.position = "left", dir = "v") +
  labs(x = "ATT (95% CI)", y = "Bed size cutoff", shape = "Cutoff") +
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

ggsave("results/bedcut-sensitivity.png", p_bc,
       width = 7, height = 0.9 * n_out + 1.2, dpi = 300)


# Save CSV ------------------------------------------------------------------
write_csv(bind_rows(bc_baseline, bc_sensitivity),
          "results/diagnostics/bedcut-sensitivity.csv")

