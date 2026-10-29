## =====================================================================
## Revision analyses for the JPubE R&R (manuscript JPUBE-D-26-00838).
##
## One script, one delineated chunk per reviewer comment that requires
## computation. Some chunks feed more than one comment (the Editor bundles
## R1.1, R2.4, and R3.3 into one required change), which is why these live
## together rather than in per-comment files. Outputs feed both the
## response-to-reviewers letter and the revised paper.
##
## Sourced by _run-analysis.r after the main (0-8) and appendix (app-*)
## scripts, so packages, parameters, helpers (functions.R), the
## stacked panels, the main results tables, and the accounting parameters
## (delta_B, delta_C from 8-access-health.R) are already in scope.
##
## Comments handled without a chunk here (writing-only or promotion of
## existing results): E1/R3.1, E5, R2.2 (text), R2.3 (text; numbers come
## from the financial-series chunk), R2.5 (answered in the first round; the
## Form 990 series are no longer built), R3.2, R2-minor(a,b,c,d).
## =====================================================================


fmt3 <- function(x) ifelse(abs(x) < 5, sprintf("%.3f", x), ifelse(abs(x) < 100, sprintf("%.2f", x), sprintf("%.1f", x)))
ci3  <- function(lo, hi) sprintf("[%s, %s]", fmt3(lo), fmt3(hi))

## Pooled SDID over cohorts for one outcome on one stacked panel. Mirrors the
## estimation in 2-hospital-dd.R / 3-hospital-dd-alt.R; `shift` counts the
## year before designation as treated (app-anticipation.R); `keep_ids`
## restricts the panel to a hospital list.
rr_sdid <- function(oname, stack, cohorts, pre, shift = FALSE, keep_ids = NULL, with_cs = FALSE) {
  osym <- sym(oname)
  run_c <- function(c) {
    d <- stack %>% filter(stack_group == c) %>%
      group_by(ID) %>% mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>% ungroup()
    if (shift) d <- d %>% mutate(post_treat = ifelse(treated == 1 & stacked_event_time >= -1, 1, post_treat))
    if (!is.null(keep_ids)) d <- d %>% filter(ID %in% keep_ids)
    d <- d %>% filter(!is.na(!!osym), is.finite(!!osym), min_bedsize <= bed.cut, stacked_event_time >= -pre) %>%
      select(ID, year, outcome = !!osym, post_treat, treated)
    bal <- as_tibble(makeBalancedPanel(d, idname = "ID", tname = "year"))
    n_yr <- length(unique(bal$year))
    if (n_yr < 3 || sum(bal$treated == 1) / n_yr < 2 || sum(bal$treated == 0) / n_yr < 2) return(NULL)
    s <- panel.matrices(as.data.frame(bal))
    if (s$T0 < 2) return(NULL)
    e <- synthdid_estimate(s$Y, s$N0, s$T0)
    se <- synthdid_se(e, method = "jackknife")
    ev <- sdid_event(e, s, c)
    cs <- if (with_cs) cs_cohort(bal, c, min.es = -pre, max.es = post) else NULL
    N   <- nrow(s$Y)
    N0  <- s$N0
    w   <- attr(e, "weights")$omega
    yrs <- as.numeric(colnames(s$Y))
    list(att = tibble(cohort = c, att = as.numeric(e), se = as.numeric(se), Ntr = N - N0, Nco = N0,
                      pre_mean = mean(s$Y[(N0 + 1):N, 1:s$T0])),
         path = tibble(cohort = c, tau = yrs - c, treated = colMeans(s$Y[(N0 + 1):N, , drop = FALSE]),
                       synthetic = as.numeric(drop(t(w) %*% s$Y[1:N0, , drop = FALSE])), Ntr = N - N0),
         sdid_es = ev$es,
         sdid_vcov = vcov_long(ev$vcov, c, ev$es$e),
         cs = cs)
  }
  out <- compact(map(cohorts, function(c) tryCatch(run_c(c), error = function(e) NULL)))
  if (length(out) == 0) return(NULL)
  atts  <- bind_rows(map(out, "att"))
  paths <- bind_rows(map(out, "path"))

  ## Pooled across cohorts, weighted by treated units. pre_mean is the treated
  ## hospitals' pre-period mean in the balanced panels actually estimated, which
  ## is the base for percent effects.
  pooled <- atts %>%
    summarise(att      = sum(Ntr * att) / sum(Ntr),
              se       = sqrt(sum(Ntr^2 * se^2)) / sum(Ntr),
              pre_mean = sum(Ntr * pre_mean) / sum(Ntr),
              ntr      = sum(Ntr),
              nco      = sum(Nco))

  agg <- paths %>%
    group_by(tau) %>%
    summarise(treated   = weighted.mean(treated, Ntr),
              synthetic = weighted.mean(synthetic, Ntr),
              .groups = "drop")

  ## Event-time effects, SDID per period and (with_cs) CS on the same balanced
  ## panels, pooled with the same treated-unit-count weights.
  sdid_es <- pool_event(bind_rows(map(out, "sdid_es"))) %>% mutate(estimator = "SDID")
  cs_list <- compact(map(out, "cs"))
  cs <- NULL
  if (length(cs_list) > 0) {
    cs_atts <- bind_rows(map(cs_list, function(x) tibble(att = x$att, se = x$se, Ntr = x$Ntr)))
    cs_att  <- with(cs_atts, sum(Ntr * att) / sum(Ntr))
    cs_se   <- with(cs_atts, sqrt(sum(Ntr^2 * se^2)) / sum(Ntr))
    cs_es   <- pool_event(bind_rows(map(cs_list, "es")))
    cs <- list(att = cs_att, lo = cs_att - 1.96 * cs_se, hi = cs_att + 1.96 * cs_se,
               es = cs_es %>% transmute(event_time = e, estimate, se))
  }
  series <- bind_rows(sdid_es,
                      if (is.null(cs)) NULL else
                        cs$es %>% transmute(e = event_time, estimate, se, estimator = "CS"))
  cohort_series <- bind_rows(bind_rows(map(out, "sdid_es")) %>% mutate(estimator = "SDID"),
                             bind_rows(map(cs_list, "es")) %>% mutate(estimator = "CS"))
  cohort_vcov <- bind_rows(bind_rows(map(out, "sdid_vcov")) %>% mutate(estimator = "SDID"),
                           bind_rows(map(cs_list, function(x) vcov_long(x$vcov, x$es$cohort[1], x$es$e))) %>%
                             mutate(estimator = "CS"))

  list(att = pooled$att,
       lo  = pooled$att - 1.96 * pooled$se,
       hi  = pooled$att + 1.96 * pooled$se,
       ntr = pooled$ntr,
       nco = pooled$nco,
       cohorts = atts,
       pre_mean = pooled$pre_mean,
       agg = agg,
       series = series,
       cohort_series = cohort_series,
       cohort_vcov = cohort_vcov,
       cs = cs)
}

write_tabular <- function(lines, header, path, align) {
  writeLines(c(sprintf("\\begin{tabular}{%s}", align), "\\toprule", header, "\\midrule", lines, "\\bottomrule", "\\end{tabular}"), path)
}


## ===== E4 / R1.2 — Outpatient, ED visits, births ==============================
## The visit outcomes are estimated with the other outcomes in 2-hospital-dd.R
## (state-timing) and 3-hospital-dd-alt.R (eligibility-restricted), and appear
## in Table 2 and the appendix eligibility table with the other outcomes. This
## chunk writes the visit estimates and series used elsewhere. pre_mean is the
## treated pre-period mean, the base for the percent effects quoted in the
## abstract and Section 4.4.
visit_labels <- c("ED visits per bed", "Outpatient visits per bed",
                  "Non-ED outpatient visits per bed", "Births")
visit_res <- bind_rows(
  hosp.results.table %>% filter(group == "visits") %>% mutate(design = "state"),
  elig.results %>% filter(group == "visits") %>% mutate(design = "elig")) %>%
  transmute(outcome, design, sdid_att, sdid_lo = sdid_ci_low, sdid_hi = sdid_ci_high, sdid_ntr,
            cs_att, cs_lo = cs_ci_low, cs_hi = cs_ci_high,
            pre_mean, pct_of_pre_mean = sdid_att / pre_mean) %>%
  arrange(match(outcome, visit_labels), match(design, c("state", "elig")))
visit_series <- bind_rows(
  hosp.event.series %>% filter(outcome %in% visit_labels) %>% mutate(design = "state"),
  elig.event.series %>% filter(outcome %in% visit_labels) %>% mutate(design = "elig"))
write_csv(visit_res, "results/diagnostics/visits-att.csv")
write_csv(visit_series, "results/diagnostics/event-series-visits.csv")


## ===== R2.2 — Pre-designation restructuring, state-timing design ==============
## The appendix anticipation analysis (treatment dated one year before
## designation) on the main state-timing design, for the main-text table.
## Writes results/att_antic_state.tex (+ CSV).
antic_map <- c(margin = "Operating margin", current_ratio = "Current ratio", net_fixed = "Net fixed assets",
               capex = "Capital expenditures per bed", net_pat_rev = "Net patient revenue per bed",
               tot_operating_exp = "Operating expenses per bed",
               BDTOT = "Total beds", OBBD = "OB beds", FTERN = "FTE RNs", ip_per_bed = "Inpatient days per bed",
               ed_per_bed = "ED visits per bed", op_per_bed = "Outpatient visits per bed",
               opnon_per_bed = "Non-ED outpatient visits per bed", BIRTHS = "Births",
               system = "System membership")
antic_rows <- list()
for (oname in names(antic_map)) {
  pp <- if (oname %in% fin.vars) financial.pre else 5
  b <- rr_sdid(oname, stack.hosp, 1999:2001, pre = pp, shift = FALSE)
  a <- rr_sdid(oname, stack.hosp, 1999:2001, pre = pp, shift = TRUE)
  if (!is.null(b) && !is.null(a)) antic_rows[[length(antic_rows) + 1]] <- tibble(outcome = antic_map[[oname]], base_att = b$att, base_lo = b$lo, base_hi = b$hi,
                                                                                antic_att = a$att, antic_lo = a$lo, antic_hi = a$hi, ntr = b$ntr)
}
antic_state <- bind_rows(antic_rows)
write_csv(antic_state, "results/diagnostics/antic-state.csv")
antic_lines <- antic_state %>%
  rowwise() %>%
  mutate(line = sprintf("%s & %s & %s & %s & %s \\\\", outcome,
                        fmt3(base_att), ci3(base_lo, base_hi),
                        fmt3(antic_att), ci3(antic_lo, antic_hi))) %>%
  pull(line)
antic_lines <- with_group_breaks(antic_lines, antic_state$outcome)

antic_header <- c(
  " & \\multicolumn{2}{c}{Treatment at designation} & \\multicolumn{2}{c}{Treatment one year earlier} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "Outcome & ATT & 95\\% CI & ATT & 95\\% CI \\\\")

write_tabular(antic_lines, antic_header, "results/att_antic_state.tex", "lcccc")
## Mean beds among eventual converters by year relative to designation, for the
## statement that most of the reduction occurs at or after designation. The 25th
## percentile and mean inpatient days per bed are the other two series plotted in
## Figure 5 and quoted in Sections 2 and 4.4.
conv_beds <- est.dat %>% filter(!is.na(eff_year)) %>%
  group_by(ID) %>% mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>% ungroup() %>%
  filter(min_bedsize <= bed.cut, !is.na(BDTOT)) %>%
  mutate(event_time = year - eff_year) %>% filter(event_time >= -5, event_time <= 5)
write_csv(bind_rows(conv_beds %>% filter(eff_year %in% 1999:2001) %>% mutate(cohorts = "1999-2001"),
                    conv_beds %>% filter(eff_year %in% 1999:2005) %>% mutate(cohorts = "1999-2005")) %>%
            group_by(cohorts, event_time) %>%
            summarise(mean_beds = mean(BDTOT), p25_beds = quantile(BDTOT, 0.25),
                      mean_ip_per_bed = mean(ip_per_bed, na.rm = TRUE),
                      n = n(), .groups = "drop"),
          "results/diagnostics/converter-beds-eventtime.csv")


## ===== E2 / R2.1 — What predicts a state's adoption cohort ====================
## Two sets of state characteristics by adoption cohort, all measured over
## 1996-98: HCRIS margins begin in 1996, and 1998 is the last year before any of
## these states adopts. Levels are 1996 values and changes run 1996 to 1998.
## Economic conditions (unemployment and per-capita income) are shown by year in
## results/adoption-timing-econ.png and tested in a regression of adoption cohort on
## their levels and changes, controlling for rurality and the number of small
## hospitals (1996 features of a state's hospital market that the economic
## changes cannot have caused). Hospital-market characteristics are summarized in
## results/att_adoption_timing.tex (+ CSV) and tested in a regression of cohort
## year on the 1996-98 closure rate and the margin level and change, controlling
## for rurality and the number of small hospitals. Both regressions use the
## 1999-2001 states.
econ_state <- read_csv("data/input/state_econ_1990_2002.csv", show_col_types = FALSE) %>%
  rename(MSTATE = state) %>%
  filter(year %in% c(1996, 1998)) %>%
  select(MSTATE, year, unemp_rate, pc_income) %>%
  pivot_wider(names_from = year, values_from = c(unemp_rate, pc_income)) %>%
  transmute(MSTATE, u_96 = unemp_rate_1996, d_u = unemp_rate_1998 - unemp_rate_1996,
            inc_96 = pc_income_1996 / 1000, d_inc = 100 * (pc_income_1998 / pc_income_1996 - 1))
st96 <- est.dat %>% filter(year == 1996) %>% group_by(MSTATE) %>%
  summarise(first_obs = first(state_treat_year), n_hosp = n_distinct(ID), n_small = sum(BDTOT <= 50, na.rm = TRUE),
            share_rural = mean(ever_rural, na.rm = TRUE), mean_beds = mean(BDTOT, na.rm = TRUE), mean_dist = mean(distance, na.rm = TRUE),
            share_gov = mean(own_gov, na.rm = TRUE), .groups = "drop")
cl <- est.dat %>% filter(year >= 1996, year <= 1998) %>% group_by(MSTATE) %>%
  summarise(closure_rate = 100 * sum(closed) / n(), .groups = "drop")
mg <- est.dat %>% filter(year %in% c(1996, 1998)) %>% group_by(MSTATE, year) %>%
  summarise(margin = mean(margin, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = year, values_from = margin, names_prefix = "m_") %>%
  transmute(MSTATE, margin_96 = m_1996, d_margin = m_1998 - m_1996)
adopt <- st96 %>%
  left_join(cl, by = "MSTATE", relationship = "one-to-one") %>%
  left_join(mg, by = "MSTATE", relationship = "one-to-one") %>%
  left_join(econ_state, by = "MSTATE", relationship = "one-to-one") %>%
  mutate(cohort = case_when(first_obs == 0 ~ "Never", first_obs < 1999 ~ "Before 1999", first_obs > 2001 ~ "After 2001", TRUE ~ as.character(first_obs)))
write_csv(adopt, "results/diagnostics/adoption-timing-states.csv")

## State economic conditions by calendar year and adoption cohort: unweighted
## state means of two series that designation cannot affect, for the treated
## cohorts (1999-2001) in black and all other states in gray.
econ_dat <- read_csv("data/input/state_econ_1990_2002.csv", show_col_types = FALSE) %>%
  rename(MSTATE = state) %>%
  inner_join(adopt %>% select(MSTATE, first_obs), by = "MSTATE", relationship = "many-to-one") %>%
  mutate(group = if_else(first_obs %in% 1999:2001, as.character(first_obs), "Other states"),
         pc_income = pc_income / 1000) %>%
  pivot_longer(c(unemp_rate, pc_income), names_to = "series", values_to = "value") %>%
  group_by(series, group, year) %>%
  summarise(value = mean(value), n_states = n(), .groups = "drop") %>%
  mutate(series = factor(series, levels = c("unemp_rate", "pc_income"),
                         labels = c("Unemployment rate (%)", "Per-capita income ($000s)")))
write_csv(econ_dat, "results/diagnostics/adoption-timing-econ.csv")
p_econ <- ggplot(econ_dat, aes(year, value, linetype = group, shape = group, color = group)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.2) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "solid")) +
  scale_shape_manual(values = c(16, 17, 15, 16)) +
  scale_color_manual(values = c("black", "black", "black", "gray65")) +
  facet_wrap(~ series, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(1990, 2002, 2)) +
  labs(x = NULL, y = NULL, linetype = "Adoption cohort", shape = "Adoption cohort", color = "Adoption cohort") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.key.width = unit(1.6, "lines"))
ggsave("results/adoption-timing-econ.png", p_econ, width = 7, height = 3.4, dpi = 300)
## State counts quoted in the abstract and Section 3.3: adopting states, states
## with a designation by 1995, states whose program began in 1999-2002, and the
## 2002 adopters.
write_csv(adopt %>% summarise(n_states = n(),
                              n_adopting = sum(first_obs > 0),
                              n_by_1995 = sum(first_obs > 0 & first_obs <= 1995),
                              n_1999_2002 = sum(first_obs >= 1999 & first_obs <= 2002),
                              n_2002 = sum(first_obs == 2002)),
          "results/diagnostics/adoption-timing-summary.csv")
adopt_vars <- c(share_rural = "Share of hospitals rural", n_small = "Hospitals with 50 or fewer beds", mean_beds = "Mean beds", share_gov = "Share government-owned",
                mean_dist = "Miles to nearest hospital", closure_rate = "Closures per 100 hospital-years, 1996-98",
                margin_96 = "Operating margin", u_96 = "Unemployment rate", inc_96 = "Per-capita income (\\$000s)",
                d_margin = "Operating margin", d_u = "Unemployment rate (points)", d_inc = "Per-capita income (percent)")
level_vars  <- c("share_rural", "n_small", "mean_beds", "share_gov", "mean_dist", "closure_rate", "margin_96", "u_96", "inc_96")
change_vars <- c("d_margin", "d_u", "d_inc")
cohort_order <- c("1999", "2000", "2001", "Never")
means <- adopt %>% group_by(cohort) %>% summarise(n_states = n(), across(all_of(names(adopt_vars)), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(cohort = factor(cohort, levels = cohort_order)) %>% arrange(cohort)
core <- adopt %>% filter(first_obs %in% 1999:2001)
adopt_econ_ols <- lm(first_obs ~ u_96 + d_u + inc_96 + d_inc + share_rural + n_small, data = core)
write_csv(tidy(adopt_econ_ols) %>% mutate(n_states = nobs(adopt_econ_ols), r_squared = summary(adopt_econ_ols)$r.squared),
          "results/diagnostics/adoption-timing-ols-econ.csv")
adopt_market_ols <- lm(first_obs ~ closure_rate + margin_96 + d_margin + share_rural + n_small, data = core)
write_csv(tidy(adopt_market_ols) %>% mutate(n_states = nobs(adopt_market_ols), r_squared = summary(adopt_market_ols)$r.squared),
          "results/diagnostics/adoption-timing-ols-market.csv")
## The same two regressions with an indicator for early adoption (1999, against
## 2000-2001) as the outcome, a linear probability model.
core <- core %>% mutate(early = as.integer(first_obs == 1999))
adopt_econ_lpm <- lm(early ~ u_96 + d_u + inc_96 + d_inc + share_rural + n_small, data = core)
write_csv(tidy(adopt_econ_lpm, conf.int = TRUE) %>% mutate(n_states = nobs(adopt_econ_lpm), r_squared = summary(adopt_econ_lpm)$r.squared),
          "results/diagnostics/adoption-timing-lpm-econ.csv")
adopt_market_lpm <- lm(early ~ closure_rate + margin_96 + d_margin + share_rural + n_small, data = core)
write_csv(tidy(adopt_market_lpm, conf.int = TRUE) %>% mutate(n_states = nobs(adopt_market_lpm), r_squared = summary(adopt_market_lpm)$r.squared),
          "results/diagnostics/adoption-timing-lpm-market.csv")
## One row per state characteristic, one column per adoption cohort, levels in
## 1996 above changes from 1996 to 1998, with state counts in the bottom row.
## Counts, bed and distance means, and income print to one decimal, the
## unemployment rate to two, shares and rates to three.
cell <- function(v, k) {
  x <- means %>% filter(cohort == k) %>% pull(v)
  if (length(x) == 0 || is.na(x)) return("")
  fmt <- if (v %in% c("n_small", "mean_beds", "mean_dist", "inc_96", "d_inc")) "%.1f" else if (v %in% c("u_96", "d_u")) "%.2f" else "%.3f"
  sprintf(fmt, x)
}
row_line <- function(v) sprintf("%s & %s \\\\", adopt_vars[[v]], paste(map_chr(cohort_order, ~ cell(v, .x)), collapse = " & "))
n_states <- cohort_order %>% map_chr(function(k) {
  x <- means %>% filter(cohort == k) %>% pull(n_states)
  if (length(x) == 0) "0" else as.character(x)
})
write_tabular(c("\\multicolumn{5}{l}{\\textit{Levels, 1996}} \\\\",
                map_chr(level_vars, row_line),
                "\\addlinespace",
                "\\multicolumn{5}{l}{\\textit{Changes, 1996 to 1998}} \\\\",
                map_chr(change_vars, row_line),
                "\\midrule",
                sprintf("States & %s \\\\", paste(n_states, collapse = " & "))),
              c(" & \\multicolumn{4}{c}{Adoption cohort} \\\\", "\\cmidrule(lr){2-5}", "State characteristic & 1999 & 2000 & 2001 & Never \\\\"),
              "results/att_adoption_timing.tex", "lcccc")


## ===== E3 / R1.1 / R2.4 / R3.3 — Inframarginal converters, break-even, bed-size split =====

## (a) Converters in the year before designation vs. hospitals in the year before closing
fin3 <- est.dat %>% arrange(ID, year) %>% group_by(ID) %>%
  mutate(margin_3y = rollapply(margin, 3, mean, na.rm = TRUE, align = "right", partial = TRUE, fill = NA),
         cr_3y = rollapply(current_ratio, 3, mean, na.rm = TRUE, align = "right", partial = TRUE, fill = NA),
         occ = ifelse(BDTOT > 0, IPDTOT / (BDTOT * 365), NA_real_),
         mcare_sh = ifelse(ADMTOT > 0, MCRDC / ADMTOT, NA_real_),
         close_next = lead(closed), year_next = lead(year)) %>% ungroup()
conv_pre <- fin3 %>% filter(!is.na(eff_year), year == eff_year - 1, eff_year >= 1999, eff_year <= 2010)
closers  <- fin3 %>% filter(close_next == 1, year_next == year + 1, year_next >= 1995, year_next <= 2010)
closers_small <- closers %>% filter(BDTOT <= 50)
cc_vars <- c(BDTOT = "Beds", occ = "Occupancy", margin_3y = "Operating margin (3-year mean)", cr_3y = "Current ratio (3-year mean)",
             distance = "Miles to nearest hospital", ever_rural = "Rural", own_gov = "Government-owned", mcare_sh = "Medicare share of admissions")
cc_groups <- list(`Converters, year before designation` = conv_pre, `Closers, year before closing` = closers, `Closers with 50 or fewer beds` = closers_small)
## One row per characteristic, one column per group. Beds and distance print to
## one decimal, the rest to three.
cc_means <- cc_groups %>%
  imap(~ .x %>%
         summarise(across(all_of(names(cc_vars)), ~ mean(.x, na.rm = TRUE))) %>%
         mutate(group = .y)) %>%
  bind_rows() %>%
  pivot_longer(all_of(names(cc_vars)), names_to = "variable", values_to = "mean")

cc_lines <- names(cc_vars) %>%
  map_chr(function(v) {
    fmt <- if (v %in% c("BDTOT", "distance")) "%.1f" else "%.3f"
    vals <- cc_means %>%
      filter(variable == v) %>%
      arrange(match(group, names(cc_groups))) %>%
      mutate(cell = sprintf(fmt, mean)) %>%
      pull(cell)
    sprintf("%s & %s \\\\", cc_vars[[v]], paste(vals, collapse = " & "))
  })

cc_counts <- cc_groups %>%
  map_int(~ n_distinct(pull(.x, ID))) %>%
  as.character()
cc_n <- sprintf("Hospitals & %s \\\\", paste(cc_counts, collapse = " & "))

write_tabular(c(cc_lines, "\\addlinespace", cc_n),
              "Characteristic & Converters & Closers & Small closers \\\\",
              "results/att_converters_closers.tex", "lccc")
med_closer_margin <- median(closers_small$margin_3y, na.rm = TRUE)
cc_summary <- tibble(median_small_closer_margin = med_closer_margin,
                     share_converters_below = mean(conv_pre$margin_3y < med_closer_margin, na.rm = TRUE),
                     n_converters_with_margin = sum(!is.na(conv_pre$margin_3y)))
write_csv(cc_summary, "results/diagnostics/converters-closers-summary.csv")

## (b) Break-even conversion share under alternative sizes of the marginal closer
##     (delta_B, delta_C, B_close from 8-access-health.R; rho = 0.30 as in Section 5.1)
if (exists("delta_B") && exists("delta_C")) {
  small98 <- est.dat %>% filter(year == 1998, BDTOT <= 50, is.na(eff_year) | eff_year > 1998) %>% pull(BDTOT)
  closer_beds_lag <- est.dat %>% arrange(ID, year) %>% group_by(ID) %>%
    mutate(beds_lag = lag(BDTOT), year_lag = lag(year)) %>% ungroup() %>%
    filter(closed == 1, year >= 1995, year <= 2010, !is.na(beds_lag), year_lag == year - 1)
  ## The last scenario draws the marginal closer from every closing hospital
  ## (mean 123 beds, Table G1). It goes to the CSV but not the paper's table,
  ## which keeps to hospital sizes within the program's reach.
  scen <- tibble(scenario = c("Mean converter before designation (baseline)", "Closers with 50 or fewer beds, mean", "Rural closers, mean",
                              "Smallest hospitals: 25th percentile of small hospitals", "Smallest hospitals: 10th percentile of small hospitals",
                              "All closers, mean"),
                 B = c(B_close, mean(closer_beds_lag$beds_lag[closer_beds_lag$beds_lag <= 50]), mean(closer_beds_lag$beds_lag[closer_beds_lag$ever_rural == 1], na.rm = TRUE),
                       quantile(small98, 0.25, na.rm = TRUE), quantile(small98, 0.10, na.rm = TRUE),
                       mean(closer_beds_lag$beds_lag))) %>%
    mutate(rho_star = (abs(delta_C) / 100 * B) / abs(delta_B), net_beds = (-delta_C / 100) * B + 0.30 * delta_B)
  write_csv(scen, "results/diagnostics/breakeven-sensitivity.csv")
  scen_lines <- scen %>%
    filter(scenario != "All closers, mean") %>%
    rowwise() %>%
    mutate(line = sprintf("%s & %.0f & %.3f & %.2f \\\\", scenario, B, rho_star, net_beds)) %>%
    pull(line)

  write_tabular(scen_lines,
                "Size of the marginal hospital & Beds & Break-even share $\\rho^*$ & Net beds per hospital at $\\rho = 0.30$ \\\\",
                "results/att_breakeven_sens.tex", "lccc")
} else {
  message("delta_B / delta_C not in scope, break-even sensitivity skipped")
}

## (c) R2.4: bed change among converters by bed size three years before designation (raw AHA beds)
raw_beds <- read_csv("data/output/aha_final.csv", show_col_types = FALSE, col_select = c(ID, year, BDTOT, eff_year)) %>%
  mutate(ID = as.character(ID)) %>% rename(beds_raw = BDTOT) %>%
  distinct(ID, year, .keep_all = TRUE) %>%
  inner_join(est.dat %>% distinct(ID, year), by = c("ID", "year")) %>%
  filter(!is.na(eff_year), eff_year >= 1999, eff_year <= 2005) %>%
  mutate(rel = year - eff_year)
split_dat <- raw_beds %>% filter(rel %in% c(-3, 0, 3)) %>% select(ID, rel, beds_raw) %>%
  pivot_wider(names_from = rel, values_from = beds_raw, names_prefix = "b") %>% rename(pre = `b-3`, at = b0, post = b3) %>%
  filter(!is.na(pre), !is.na(at), !is.na(post)) %>%
  mutate(bin = cut(pre, c(0, 25, 35, 50, 75, Inf), labels = c("25 or fewer", "26 to 35", "36 to 50", "51 to 75", "More than 75")))
split_tab <- split_dat %>% group_by(bin) %>%
  summarise(n = n(), beds_pre = mean(pre), beds_at = mean(at), beds_post = mean(post), change = mean(post - pre),
            share_at25 = mean(post == 25), share_le25 = mean(post <= 25), .groups = "drop")
write_csv(split_tab, "results/diagnostics/bedsize-split.csv")
split_lines <- split_tab %>%
  rowwise() %>%
  mutate(line = sprintf("%s & %d & %.1f & %.1f & %.1f & %.1f & %.2f \\\\",
                        bin, n, beds_pre, beds_at, beds_post, change, share_at25)) %>%
  pull(line)

split_header <- c(
  "Beds three years & & \\multicolumn{3}{c}{Mean beds} & Change & Share at exactly \\\\",
  "before designation & Hospitals & $t-3$ & $t$ & $t+3$ & $t-3$ to $t+3$ & 25 beds at $t+3$ \\\\")

write_tabular(split_lines, split_header, "results/att_bedsize_split.tex", "lcccccc")

## (d) R3.4: dollar accounting of the margin effect (net income gain per converter, aggregate, per averted closure)
if (exists("delta_C")) {
  rev_pre <- est.dat %>% filter(!is.na(eff_year), eff_year >= 1999, eff_year <= 2005, year >= eff_year - 3, year <= eff_year - 1) %>%
    group_by(ID) %>% mutate(mb = min(BDTOT, na.rm = TRUE)) %>% ungroup() %>% filter(mb <= bed.cut) %>%
    summarise(npr_per_bed = mean(net_pat_rev, na.rm = TRUE), beds_base = mean(beds_base, na.rm = TRUE), npr_total = mean(net_pat_rev * beds_base, na.rm = TRUE))
  m_state <- hosp.results.table %>% filter(outcome == "Operating margin") %>% pull(sdid_att)
  m_elig  <- elig.results %>% filter(outcome == "Operating margin") %>% pull(sdid_att)
  n_conv  <- if (exists("n_converters")) n_converters else 1240
  n_hosp  <- if (exists("n_hosp_cah_states")) n_hosp_cah_states else 4400
  closures_pyr <- abs(delta_C) / 100 * n_hosp
  dollars <- tibble(design = c("State-timing", "Eligibility-restricted"), margin_effect = c(m_state, m_elig)) %>%
    mutate(mean_revenue_k2010 = rev_pre$npr_total, gain_per_hospital_k2010 = margin_effect * mean_revenue_k2010,
           converters = n_conv, aggregate_billion_2010 = gain_per_hospital_k2010 * converters / 1e6,
           averted_closures_per_year = closures_pyr, gain_per_averted_closure_million_2010 = gain_per_hospital_k2010 * converters / closures_pyr / 1e3)
  write_csv(dollars, "results/diagnostics/dollar-accounting.csv")
}


## ===== R1.3 — Bed distribution around the 25-bed limit (raw AHA beds) =========
notch <- read_csv("data/output/aha_final.csv", show_col_types = FALSE, col_select = c(ID, year, BDTOT)) %>%
  mutate(ID = as.character(ID)) %>% rename(beds_raw = BDTOT) %>%
  inner_join(est.dat %>% distinct(ID, year, state_treat_year, ever_cah), by = c("ID", "year")) %>%
  filter(beds_raw > 0, beds_raw <= 100, year >= 1995, year <= 2008) %>%
  mutate(rel = ifelse(state_treat_year > 0, year - state_treat_year, NA_real_))
adopting <- notch %>% filter(state_treat_year %in% 1999:2001)
before <- adopting %>% filter(rel >= -4, rel <= -1) %>% mutate(panel = "Before state adoption (t-4 to t-1)")
after  <- adopting %>% filter(rel >= 3, rel <= 6) %>% mutate(panel = "After state adoption (t+3 to t+6)")
p_notch <- ggplot(bind_rows(before, after) %>% mutate(panel = factor(panel, levels = c("Before state adoption (t-4 to t-1)", "After state adoption (t+3 to t+6)"))), aes(beds_raw)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0.5, fill = "gray50", color = "white") +
  geom_vline(xintercept = 25.5, linetype = "dashed") + facet_wrap(~ panel, ncol = 2) +
  coord_cartesian(xlim = c(5, 60)) + labs(x = "Total beds", y = "Density") + theme_bw(base_size = 13)
ggsave("results/rr-bed-notch-beforeafter.png", p_notch, width = 10, height = 4.2, dpi = 300)
same_year <- bind_rows(notch %>% filter(year == 1999, state_treat_year > 0, state_treat_year <= 1999) %>% mutate(group = "Adopted"),
                       notch %>% filter(year == 1999, state_treat_year > 1999) %>% mutate(group = "Not yet adopted"))
p_notch2 <- ggplot(same_year, aes(beds_raw, fill = group)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0.5, position = "identity", alpha = 0.55, color = NA) +
  scale_fill_manual(values = c("Adopted" = "black", "Not yet adopted" = "gray60")) +
  geom_vline(xintercept = 25.5, linetype = "dashed") + coord_cartesian(xlim = c(5, 60)) +
  labs(x = "Total beds, 1999", y = "Density", fill = NULL) + theme_bw(base_size = 13) + theme(legend.position = "bottom")
ggsave("results/rr-bed-notch-adopting.png", p_notch2, width = 6.5, height = 4.2, dpi = 300)
## Shares at and around 25 beds: adopting states before/after adoption, the 1999
## adopted vs not-yet-adopted comparison, not-yet-adopted states in each of
## 1999-2001, eventual converters vs never-converters within adopting states,
## and adopting states in 2006.
notch_shares <- bind_rows(before %>% mutate(g = "adopting states, before"), after %>% mutate(g = "adopting states, after"), same_year %>% mutate(g = paste0("1999, ", group)),
                          map_dfr(1999:2001, ~ notch %>% filter(year == .x, state_treat_year > .x) %>% mutate(g = sprintf("%d, Not yet adopted", .x))),
                          before %>% mutate(g = ifelse(ever_cah == 1, "adopting states, before: eventual converters", "adopting states, before: never-converters")),
                          after %>% mutate(g = ifelse(ever_cah == 1, "adopting states, after: eventual converters", "adopting states, after: never-converters")),
                          adopting %>% filter(year == 2006) %>% mutate(g = "adopting states, 2006")) %>%
  group_by(g) %>% summarise(n = n(), share_at25 = mean(beds_raw == 25), share_26_30 = mean(beds_raw >= 26 & beds_raw <= 30), share_le25 = mean(beds_raw <= 25),
                            ratio_25_to_neighbors = mean(beds_raw == 25) / pmax((mean(beds_raw == 24) + mean(beds_raw == 26)) / 2, 1e-9), .groups = "drop")
write_csv(notch_shares, "results/diagnostics/bed-notch-shares.csv")


## ===== R2.4 — Payer mix ======================================================
est.dat <- est.dat %>% mutate(mcaid_sh = ifelse(ADMTOT > 0 & MCDDC / ADMTOT <= 1, MCDDC / ADMTOT, NA_real_),
                              mcare_sh = ifelse(ADMTOT > 0 & MCRDC / ADMTOT <= 1, MCRDC / ADMTOT, NA_real_))
stack.hosp <- stack_hosp(pre.period = 5, post.period = post, state.period = state.cut)
stack.elig <- stack_hosp_elig(pre.period = 5, post.period = post, cohort.years = 1999:2005)
payer_labels <- c(mcaid_sh = "Medicaid share of admissions",
                  mcare_sh = "Medicare share of admissions")

payer_rows <- list()
for (oname in names(payer_labels)) {
  for (dsg in c("state", "elig")) {
    stack   <- if (dsg == "state") stack.hosp else stack.elig
    cohorts <- if (dsg == "state") 1999:2001 else 1999:2005
    s <- rr_sdid(oname, stack, cohorts, pre = 5)
    if (is.null(s)) next
    payer_rows[[length(payer_rows) + 1]] <- tibble(
      outcome = payer_labels[[oname]], design = dsg,
      att = s$att, lo = s$lo, hi = s$hi, ntr = s$ntr)
  }
}
write_csv(bind_rows(payer_rows), "results/diagnostics/payer-mix.csv")


## ===== Financial series: pre-designation means and percent effects =========
## SDID on the HCRIS series for each financial outcome, with the treated
## pre-period mean as the base for the percent effects quoted in Section 4.2,
## plus the same series restricted to converters with a strict crosswalk match
## (score >= 0.90, or zip-exact with score >= 0.75). Never-converters are kept
## in both. Writes results/diagnostics/source-comparison.csv.
strict_ids <- est.dat %>% filter(is.na(eff_year) | (!is.na(match_strict) & match_strict)) %>% distinct(ID) %>% pull(ID)
src_rows <- list()
for (dsg in c("state", "elig")) {
  stack <- if (dsg == "state") stack.hosp else stack.elig; cohorts <- if (dsg == "state") 1999:2001 else 1999:2005
  for (o in fin.vars) {
    for (variant in c("h", "hs")) {
      s <- rr_sdid(o, stack, cohorts, pre = financial.pre, keep_ids = if (variant == "hs") strict_ids else NULL)
      if (!is.null(s)) src_rows[[length(src_rows) + 1]] <- tibble(outcome = o, design = dsg, series = c(h = "HCRIS", hs = "HCRIS, strict matches")[[variant]],
                                                                  att = s$att, lo = s$lo, hi = s$hi, ntr = s$ntr, pre_mean = s$pre_mean, pct_of_pre_mean = s$att / s$pre_mean)
    }
  }
}
src_res <- bind_rows(src_rows)
write_csv(src_res, "results/diagnostics/source-comparison.csv")


## ===== R1-minor(a) — Main SDID table with permutation p-values ===============
perm_file <- "results/diagnostics/permutation_pvalues.csv"
if (file.exists(perm_file)) {
  perm <- read_csv(perm_file, show_col_types = FALSE) %>% select(outcome, p_value)
  main_perm <- results.table %>% left_join(perm, by = "outcome", relationship = "one-to-one")
  lines <- main_perm %>% rowwise() %>%
    mutate(line = sprintf("%s & %s & %s & %s & %s & %s & %d \\\\", outcome,
                          fmt3(sdid_att), ci3(sdid_ci_low, sdid_ci_high),
                          ifelse(is.na(p_value), "", sprintf("%.3f", p_value)),
                          ifelse(is.na(cs_att), "", fmt3(cs_att)),
                          ifelse(is.na(cs_att), "", ci3(cs_ci_low, cs_ci_high)),
                          sdid_ntr)) %>% pull(line)
  lines <- with_group_breaks(lines, main_perm$outcome)
  write_tabular(lines,
                c(" & \\multicolumn{3}{c}{SDID} & \\multicolumn{2}{c}{Callaway--Sant'Anna} & \\\\",
                  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-6}",
                  "Outcome & ATT & 95\\% CI & Permutation $p$ & ATT & 95\\% CI & $N_{tr}$ \\\\"),
                "results/att_overall_perm.tex", "lcccccr")
}


## ===== R2-minor(e) — Multiple-testing adjustment for heterogeneity ===========
if (exists("het.results")) {
  het_diff <- het.results %>% mutate(se = (ci_high - ci_low) / 3.92) %>%
    group_by(dimension, outcome) %>% filter(n() == 2) %>%
    summarise(diff = diff(att), se_diff = sqrt(sum(se^2)), groups = paste(subgroup, collapse = " vs "), .groups = "drop") %>%
    mutate(z = diff / se_diff, p = 2 * pnorm(-abs(z)), tests = n(), p_bonferroni = pmin(p * tests, 1), survives = p_bonferroni < 0.05) %>%
    arrange(p)
  write_csv(het_diff, "results/diagnostics/het-multipletest.csv")
}


## ===== E5 — Closure ATT by cohort, SDID against Callaway--Sant'Anna ==========
## Cohort-level closure ATTs behind the closure footnote. SDID is the
## clean-donor estimator of 4-changes-state-dd.R (sdid_state_clean), and CS runs
## through cs_state_cohort() on the cohort's balanced state panel, averaged over
## the cohort's five post-treatment years. Values are per state-year, before the
## per-100-hospital scaling. n_control_states is the donor pool in the adoption
## year and n_donors_final the pool in the last post year.

adopt.closure <- state.dat %>% distinct(MSTATE, state_treat_year)

closure_cohort <- tibble()
for (c in 1999:2001) {
  sdid.c <- sdid_state_clean(state.dat, state.timing, c, "closures", pre = 5, post = post, lag = state.cut)
  if (is.null(sdid.c)) next

  cs.c <- tryCatch(cs_state_cohort(stack.state, c, "closures", adopt.closure, post)$att,
                   error = function(e) NA_real_)

  closure_cohort <- bind_rows(closure_cohort, tibble(
    cohort           = c,
    n_treated_states = sdid.c$Ntr,
    n_control_states = as.integer(sdid.c$donors[1]),
    n_donors_final   = as.integer(sdid.c$donors[length(sdid.c$donors)]),
    sdid_att         = sdid.c$att / sdid.c$scale,
    sdid_se          = sdid.c$se / sdid.c$scale,
    cs_att           = cs.c))
}

if (nrow(closure_cohort) > 0) {
  write_csv(closure_cohort, "results/diagnostics/closure-cohort-cs-sdid.csv")
} else {
  message("E5 closure-by-cohort produced no rows, skipping")
}

