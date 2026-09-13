# Meta --------------------------------------------------------------------
## Heterogeneity analysis: subgroup-specific SDID ATTs
## Eligibility-restricted design (stack.elig), cohorts 1999-2005
## Dimensions: ownership, isolation, pre-treatment margin, system membership
## Each dimension is an independent subgroup split; all use SDID only

# Expects from _run-analysis.r: stack.elig, est.dat, bed.cut, post

het.cohorts <- 1999:2005
het.fin.pre <- 3   # fixed pre-period for financial outcomes (power)


# Outcome map --------------------------------------------------------------
het_outcome_map <- list(
  margin            = list(label = "Operating margin",             pre_period = het.fin.pre),
  current_ratio     = list(label = "Current ratio",                pre_period = het.fin.pre),
  net_fixed         = list(label = "Net fixed assets",             pre_period = het.fin.pre),
  capex             = list(label = "Capital expenditures per bed", pre_period = het.fin.pre),
  net_pat_rev       = list(label = "Net patient revenue per bed",  pre_period = het.fin.pre),
  tot_operating_exp = list(label = "Operating expenses per bed",   pre_period = het.fin.pre),
  BDTOT             = list(label = "Total beds"),
  OBBD              = list(label = "OB beds"),
  FTERN             = list(label = "FTE RNs"),
  ip_per_bed        = list(label = "Inpatient days per bed"),
  system            = list(label = "System membership")
)


# Hospital-level baseline characteristics ----------------------------------
## Computed once from est.dat; used to split stack.elig by ID

hosp_base <- est.dat %>%
  group_by(ID) %>%
  summarise(
    eff_yr = first(na.omit(eff_year)),

    ## Ownership: modal value (essentially time-invariant)
    own_gov_base = {
      vals <- own_gov[!is.na(own_gov)]
      if (length(vals) > 0) as.numeric(mean(vals) >= 0.5) else NA_real_
    },
    own_nfp_base = {
      vals <- own_nfp[!is.na(own_nfp)]
      if (length(vals) > 0) as.numeric(mean(vals) >= 0.5) else NA_real_
    },

    ## Distance to nearest hospital: mean across years (hospital doesn't move)
    distance_base = mean(distance, na.rm = TRUE),

    ## Pre-treatment margin: mean in 1996-1998 (pre-treatment for all cohorts)
    margin_base = mean(margin[year >= 1996 & year <= 1998], na.rm = TRUE),

    ## System membership: year before treatment, or 1998 for never-treated
    system_base = {
      ey <- first(na.omit(eff_year))
      ref <- if (!is.na(ey)) ey - 1 else 1998
      vals <- system[year == ref]
      if (length(vals) > 0 && !is.na(vals[1])) vals[1] else NA_real_
    },

    .groups = "drop"
  )

## Medians for continuous splits (among treated hospitals with beds <= bed.cut)
treated_base <- hosp_base %>%
  filter(!is.na(eff_yr)) %>%
  semi_join(
    est.dat %>% group_by(ID) %>%
      summarise(min_beds = min(BDTOT, na.rm = TRUE), .groups = "drop") %>%
      filter(min_beds <= bed.cut),
    by = "ID"
  )

med_distance <- median(treated_base$distance_base, na.rm = TRUE)
med_margin   <- median(treated_base$margin_base, na.rm = TRUE)

## med_distance and med_margin define the isolation and pre-treatment margin
## splits below. Neither is reported in the paper.


# Heterogeneity dimensions ------------------------------------------------
## Each dimension: list of two named subgroup filters (functions on stack.elig)

het_dims <- list(
  ownership = list(
    label = "Ownership",
    groups = list(
      Government = function(df) df %>% semi_join(hosp_base %>% filter(own_gov_base == 1), by = "ID"),
      Nonprofit  = function(df) df %>% semi_join(hosp_base %>% filter(own_nfp_base == 1), by = "ID")
    )
  ),
  isolation = list(
    label = "Geographic isolation",
    groups = list(
      Isolated  = function(df) df %>% semi_join(hosp_base %>% filter(distance_base >= med_distance), by = "ID"),
      Proximate = function(df) df %>% semi_join(hosp_base %>% filter(distance_base < med_distance), by = "ID")
    )
  ),
  margin = list(
    label = "Pre-treatment margin",
    groups = list(
      `Above median` = function(df) df %>% semi_join(hosp_base %>% filter(margin_base >= med_margin), by = "ID"),
      `Below median` = function(df) df %>% semi_join(hosp_base %>% filter(margin_base < med_margin), by = "ID")
    )
  ),
  system = list(
    label = "System membership",
    groups = list(
      System      = function(df) df %>% semi_join(hosp_base %>% filter(system_base == 1), by = "ID"),
      Independent = function(df) df %>% semi_join(hosp_base %>% filter(system_base == 0), by = "ID")
    )
  )
)


# SDID estimation ----------------------------------------------------------
## Runs SDID across cohorts for a subgroup-filtered stack, returns weighted ATT

run_sdid_het <- function(stack_data, outcome_sym, pp, cohorts) {

  run_one <- function(c) {
    synth.c <- stack_data %>%
      filter(stack_group == c) %>%
      group_by(ID) %>% mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>% ungroup() %>%
      filter(!is.na(!!outcome_sym), min_bedsize <= bed.cut,
             stacked_event_time >= -pp) %>%
      select(ID, year, outcome = !!outcome_sym, post_treat, min_bedsize, treated)

    bal.c <- as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year"))

    n_yr <- length(unique(bal.c$year))
    n_tr <- sum(bal.c$treated == 1) / max(n_yr, 1)
    n_co <- sum(bal.c$treated == 0) / max(n_yr, 1)
    if (n_tr < 2 || n_co < 2) return(NULL)

    setup <- panel.matrices(as.data.frame(bal.c))
    est   <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
    se_c  <- synthdid_se(est, method = "jackknife")
    Ntr   <- nrow(setup$Y) - setup$N0

    tibble(cohort = c, att = as.numeric(est), se = as.numeric(se_c), Ntr = Ntr)
  }

  atts_all <- bind_rows(map(cohorts, possibly(run_one, otherwise = NULL)))
  atts_all <- atts_all %>% filter(!is.na(att))

  if (nrow(atts_all) == 0) return(NULL)

  att_w  <- with(atts_all, sum(Ntr * att) / sum(Ntr))
  se_w   <- with(atts_all, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))

  tibble(
    att       = att_w,
    ci_low    = att_w - 1.96 * se_w,
    ci_high   = att_w + 1.96 * se_w,
    ntr       = sum(atts_all$Ntr),
    n_cohorts = nrow(atts_all)
  )
}


# Main loop ----------------------------------------------------------------
het.results <- tibble(
  dimension = character(), subgroup = character(), outcome = character(),
  att = numeric(), ci_low = numeric(), ci_high = numeric(),
  ntr = integer(), n_cohorts = integer()
)

for (dim_name in names(het_dims)) {
  dim <- het_dims[[dim_name]]

  for (grp_name in names(dim$groups)) {
    grp_filter <- dim$groups[[grp_name]]
    stack_sub  <- grp_filter(stack.elig)

    for (oname in names(het_outcome_map)) {
      # Skip self-referential: estimating system outcome split by system dimension
      if (dim_name == "system" && oname == "system") next

      o  <- het_outcome_map[[oname]]
      pp <- if (!is.null(o$pre_period)) o$pre_period else 5

      result <- run_sdid_het(stack_sub, sym(oname), pp, het.cohorts)

      if (!is.null(result)) {
        het.results <- bind_rows(het.results, tibble(
          dimension = dim$label,
          subgroup  = grp_name,
          outcome   = o$label,
          att       = result$att,
          ci_low    = result$ci_low,
          ci_high   = result$ci_high,
          ntr       = as.integer(result$ntr),
          n_cohorts = as.integer(result$n_cohorts)
        ))
      } else {
        message(sprintf("heterogeneity: SDID failed for %s / %s / %s", dim_name, grp_name, oname))
      }
    }
  }
}


# Save raw results ---------------------------------------------------------
write_csv(het.results, "results/het_results.csv")


# Forest plots (one per dimension) -----------------------------------------
het_forest_labels <- c(
  "Operating margin"             = "Op. margin",
  "Current ratio"                = "Current ratio",
  "Net fixed assets"             = "Net fixed\nassets",
  "Capital expenditures per bed" = "CapEx/bed",
  "Net patient revenue per bed"  = "Revenue/bed",
  "Operating expenses per bed"   = "Expenses/bed",
  "Total beds"                   = "Total beds",
  "OB beds"                      = "OB beds",
  "FTE RNs"                      = "FTE RNs",
  "Inpatient days per bed"       = "IP days/bed",
  "System membership"            = "System\nmembership"
)

outcome_order <- names(het_forest_labels)

for (dim_name in unique(het.results$dimension)) {
  dim_dat <- het.results %>%
    filter(dimension == dim_name) %>%
    mutate(
      outcome_label = het_forest_labels[outcome],
      outcome_label = factor(outcome_label, levels = het_forest_labels[outcome_order])
    )

  n_out <- length(unique(dim_dat$outcome_label))
  dim_stub <- tolower(gsub("[^a-z]", "", tolower(dim_name)))

  p <- ggplot(dim_dat, aes(x = att, y = subgroup, shape = subgroup)) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "gray40") +
    geom_pointrange(aes(xmin = ci_low, xmax = ci_high), size = 0.4, linewidth = 0.4) +
    facet_wrap(~ outcome_label, ncol = 1, scales = "free_x",
               strip.position = "left", dir = "v") +
    labs(x = "ATT (95% CI)", y = NULL, shape = dim_name,
         title = sprintf("Heterogeneity: %s", dim_name)) +
    theme_bw(base_size = 11) +
    theme(
      strip.background  = element_blank(),
      strip.placement   = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9, lineheight = 0.9),
      panel.grid.minor  = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      legend.position   = "bottom",
      plot.title        = element_text(size = 12, face = "bold")
    )

  ggsave(sprintf("results/het-forest-%s.png", dim_stub),
         p, width = 7, height = 0.9 * n_out + 1.5, dpi = 300)
}


# LaTeX tables (one per dimension) -----------------------------------------
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

for (dim_name in unique(het.results$dimension)) {
  dim_dat <- het.results %>% filter(dimension == dim_name)
  grps <- unique(dim_dat$subgroup)

  wide <- dim_dat %>%
    select(subgroup, outcome, att, ci_low, ci_high, ntr) %>%
    pivot_wider(
      names_from = subgroup,
      values_from = c(att, ci_low, ci_high, ntr),
      names_glue = "{subgroup}_{.value}"
    )

  g1 <- grps[1]; g2 <- grps[2]
  tex_lines <- wide %>%
    rowwise() %>%
    mutate(line = sprintf("%s & %s & %s & %d & %s & %s & %d \\\\",
      outcome,
      fmt(.data[[paste0(g1, "_att")]]),
      int(.data[[paste0(g1, "_ci_low")]], .data[[paste0(g1, "_ci_high")]]),
      .data[[paste0(g1, "_ntr")]],
      fmt(.data[[paste0(g2, "_att")]]),
      int(.data[[paste0(g2, "_ci_low")]], .data[[paste0(g2, "_ci_high")]]),
      .data[[paste0(g2, "_ntr")]]
    )) %>%
    ungroup() %>%
    pull(line)

  dim_stub <- tolower(gsub("[^a-z]", "", tolower(dim_name)))

  writeLines(c(
    sprintf("\\begin{tabular}{lccrccr}"),
    sprintf(" & \\multicolumn{3}{c}{%s} & \\multicolumn{3}{c}{%s} \\\\", g1, g2),
    sprintf("Outcome & ATT & 95\\%% CI & $N_{tr}$ & ATT & 95\\%% CI & $N_{tr}$ \\\\"),
    tex_lines,
    "\\end{tabular}"
  ), sprintf("results/att_het_%s.tex", dim_stub))
}

# Combined forest plot (all dimensions side-by-side) -----------------------
# Each column is a stack of per-outcome panels. An outcome that cannot be
# estimated for a dimension (system membership as an outcome when splitting on
# system membership) gets a blank spacer, so rows stay aligned across columns
# without drawing an empty panel.

make_het_cell <- function(dim_name, olab, first_col, first_row) {
  d <- het.results %>%
    mutate(outcome_label = unname(het_forest_labels[outcome])) %>%
    filter(dimension == dim_name, outcome_label == olab)
  if (nrow(d) == 0) return(plot_spacer())
  ggplot(d, aes(x = att, y = subgroup, shape = subgroup)) +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", color = "gray40") +
    geom_pointrange(aes(xmin = ci_low, xmax = ci_high), size = 0.35, linewidth = 0.35) +
    labs(x = NULL, y = if (first_col) olab else NULL, shape = NULL,
         title = if (first_row) dim_name else NULL) +
    theme_bw(base_size = 12) +
    theme(
      axis.title.y       = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 11, lineheight = 0.9),
      axis.text.y        = element_text(size = 10),
      axis.text.x        = element_text(size = 9),
      axis.ticks.y       = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position    = "none",
      plot.title         = element_text(size = 11, hjust = 0.5, face = "bold"),
      plot.margin        = margin(2, 4, 2, 2)
    )
}

dim_names <- unique(het.results$dimension)
olabs     <- unname(het_forest_labels[outcome_order])
columns <- map(seq_along(dim_names), function(j) {
  cells <- map(seq_along(olabs), ~ make_het_cell(dim_names[j], olabs[.x], first_col = (j == 1), first_row = (.x == 1)))
  wrap_plots(cells, ncol = 1)
})
combined <- wrap_plots(columns, nrow = 1)
n_out <- length(olabs)

ggsave("results/het-forest-combined.png", combined,
       width = 11, height = 0.85 * n_out + 1.4, dpi = 300)


