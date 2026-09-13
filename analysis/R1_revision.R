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
## from the R2.5 chunk), R3.2, R2-minor(a,b,c,d). See response-plan.md.
## =====================================================================


winsor_by_year <- function(x, probs = c(0.05, 0.95)) {
  qs <- quantile(x, probs = probs, na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}

fmt3 <- function(x) ifelse(abs(x) < 5, sprintf("%.3f", x), ifelse(abs(x) < 100, sprintf("%.2f", x), sprintf("%.1f", x)))
ci3  <- function(lo, hi) sprintf("[%s, %s]", fmt3(lo), fmt3(hi))

## Pooled SDID over cohorts for one outcome on one stacked panel. Mirrors the
## estimation in 2-hospital-dd.R / 3-hospital-dd-alt.R; `shift` counts the
## year before designation as treated (app-anticipation.R); `keep_ids`
## restricts the panel to a hospital list.
rr_sdid <- function(oname, stack, cohorts, pre, shift = FALSE, keep_ids = NULL) {
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
    N   <- nrow(s$Y)
    N0  <- s$N0
    w   <- attr(e, "weights")$omega
    yrs <- as.numeric(colnames(s$Y))
    list(att = tibble(cohort = c, att = as.numeric(e), se = as.numeric(se), Ntr = N - N0, Nco = N0,
                      pre_mean = mean(s$Y[(N0 + 1):N, 1:s$T0])),
         path = tibble(cohort = c, tau = yrs - c, treated = colMeans(s$Y[(N0 + 1):N, , drop = FALSE]),
                       synthetic = as.numeric(drop(t(w) %*% s$Y[1:N0, , drop = FALSE])), Ntr = N - N0))
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

  list(att = pooled$att,
       lo  = pooled$att - 1.96 * pooled$se,
       hi  = pooled$att + 1.96 * pooled$se,
       ntr = pooled$ntr,
       nco = pooled$nco,
       cohorts = atts,
       pre_mean = pooled$pre_mean,
       agg = agg)
}

rr_cs <- function(oname, design) {
  out <- cs_att(est.dat, oname,
                if (design == "state") 1999:2001 else 1999:2005,
                bed.cut = bed.cut)
  if (is.null(out)) return(NULL)
  list(att = out$att, lo = out$ci_low, hi = out$ci_high)
}

write_tabular <- function(lines, header, path, align) {
  writeLines(c(sprintf("\\begin{tabular}{%s}", align), "\\toprule", header, "\\midrule", lines, "\\bottomrule", "\\end{tabular}"), path)
}


## ===== E4 / R1.2 — Outpatient, ED visits, births ==============================
## Same SDID (and CS) specification as the main results, both designs. Visits
## are counts per pre-designation bed, as for revenue and expenses; births stay
## a count. Writes results/att_visits.tex (state-timing) and
## results/att_visits_elig.tex (eligibility-restricted), each with SDID and CS
## estimates, plus per-outcome SDID path figures and a CSV.
visits <- read_csv("data/output/aha_visits.csv", show_col_types = FALSE,
                   col_types = cols(ID = col_character()))

est.dat <- est.dat %>%
  select(-any_of(c("VEM", "VTOT", "VOTH", "BIRTHS", "ed_per_bed", "op_per_bed", "opnon_per_bed"))) %>%
  left_join(visits, by = c("ID", "year")) %>%
  group_by(year, ever_cah) %>%
  mutate(across(c(VEM, VTOT, VOTH, BIRTHS), winsor_by_year)) %>%
  ungroup() %>%
  mutate(ed_per_bed = ifelse(beds_base > 0, VEM / beds_base, NA_real_),
         op_per_bed = ifelse(beds_base > 0, VTOT / beds_base, NA_real_),
         opnon_per_bed = ifelse(beds_base > 0, VOTH / beds_base, NA_real_))
stack.hosp <- stack_hosp(pre.period = 5, post.period = post, state.period = state.cut)
stack.elig <- stack_hosp_elig(pre.period = 5, post.period = post, cohort.years = 1999:2005)

## Counts per pre-designation bed, as for revenue and expenses. Births stay a
## count, since dividing births by beds would not mean anything useful.
visit_map <- list(ed_per_bed = list(label = "ED visits per bed", stub = "edvisits"),
                  op_per_bed = list(label = "Outpatient visits per bed", stub = "opvisits"),
                  opnon_per_bed = list(label = "Non-ED outpatient visits per bed", stub = "opvisits-noned"),
                  BIRTHS = list(label = "Births", stub = "births"))
visit_rows <- list()
for (oname in names(visit_map)) {
  o <- visit_map[[oname]]
  for (dsg in c("state", "elig")) {
    s <- rr_sdid(oname, if (dsg == "state") stack.hosp else stack.elig, if (dsg == "state") 1999:2001 else 1999:2005, pre = 5)
    cs <- rr_cs(oname, dsg)
    if (!is.null(s)) {
      visit_rows[[length(visit_rows) + 1]] <- tibble(outcome = o$label, design = dsg, sdid_att = s$att, sdid_lo = s$lo, sdid_hi = s$hi, sdid_ntr = s$ntr,
                                                     cs_att = if (is.null(cs)) NA_real_ else cs$att, cs_lo = if (is.null(cs)) NA_real_ else cs$lo, cs_hi = if (is.null(cs)) NA_real_ else cs$hi)
      p <- ggplot(s$agg, aes(tau)) +
        geom_line(aes(y = treated, linetype = "Treated"), linewidth = 0.9) +
        geom_line(aes(y = synthetic, linetype = "Synthetic control"), linewidth = 0.9) +
        geom_vline(xintercept = -0.5, linewidth = 1) +
        scale_linetype_manual(values = c("Treated" = "solid", "Synthetic control" = "dashed")) +
        labs(x = "Event time", y = o$label, linetype = NULL) + theme_bw(base_size = 13) + theme(legend.position = "bottom")
      ggsave(sprintf("results/%s-%s-sdid.png", o$stub, dsg), p, width = 6.5, height = 4.25, dpi = 300)
    }
  }
}
visit_res <- bind_rows(visit_rows)
write_csv(visit_res, "results/diagnostics/visits-att.csv")
dsg_label <- c(state = "State-timing", elig = "Eligibility-restricted")
## Main text carries the state-timing design; the eligibility-restricted design
## goes to the appendix, matching how the other outcomes are presented.
visit_line <- function(d) {
  d %>% rowwise() %>%
    mutate(line = sprintf("%s & %s & %s & %s & %s & %d \\\\", outcome,
                          fmt3(sdid_att), ci3(sdid_lo, sdid_hi),
                          ifelse(is.na(cs_att), "", fmt3(cs_att)),
                          ifelse(is.na(cs_att), "", ci3(cs_lo, cs_hi)),
                          sdid_ntr)) %>%
    pull(line)
}
visit_header <- c(" & \\multicolumn{2}{c}{SDID} & \\multicolumn{2}{c}{Callaway--Sant'Anna} & \\\\",
                  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
                  "Outcome & ATT & 95\\% CI & ATT & 95\\% CI & $N_{tr}$ \\\\")
write_tabular(visit_line(visit_res %>% filter(design == "state")),
              visit_header, "results/att_visits.tex", "lccccr")
write_tabular(visit_line(visit_res %>% filter(design == "elig")),
              visit_header, "results/att_visits_elig.tex", "lccccr")


## ===== R2.2 — Pre-designation restructuring, state-timing design ==============
## The appendix anticipation analysis (treatment dated one year before
## designation) on the main state-timing design, for the main-text table.
## Writes results/att_antic_state.tex (+ CSV).
antic_map <- c(BDTOT = "Total beds", OBBD = "OB beds", FTERN = "FTE RNs", ip_per_bed = "Inpatient days per bed", system = "System membership")
antic_rows <- list()
for (oname in names(antic_map)) {
  b <- rr_sdid(oname, stack.hosp, 1999:2001, pre = 5, shift = FALSE)
  a <- rr_sdid(oname, stack.hosp, 1999:2001, pre = 5, shift = TRUE)
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

antic_header <- c(
  " & \\multicolumn{2}{c}{Treatment at designation} & \\multicolumn{2}{c}{Treatment one year earlier} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "Outcome & ATT & 95\\% CI & ATT & 95\\% CI \\\\")

write_tabular(antic_lines, antic_header, "results/att_antic_state.tex", "lcccc")
## Mean beds among eventual converters by year relative to designation, for the
## statement that most of the reduction occurs at or after designation.
conv_beds <- est.dat %>% filter(!is.na(eff_year)) %>%
  group_by(ID) %>% mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>% ungroup() %>%
  filter(min_bedsize <= bed.cut, !is.na(BDTOT)) %>%
  mutate(event_time = year - eff_year) %>% filter(event_time >= -5, event_time <= 5)
write_csv(bind_rows(conv_beds %>% filter(eff_year %in% 1999:2001) %>% mutate(cohorts = "1999-2001"),
                    conv_beds %>% filter(eff_year %in% 1999:2005) %>% mutate(cohorts = "1999-2005")) %>%
            group_by(cohorts, event_time) %>% summarise(mean_beds = mean(BDTOT), n = n(), .groups = "drop"),
          "results/diagnostics/converter-beds-eventtime.csv")


## ===== E2 / R2.1 — What predicts a state's adoption cohort ====================
## State characteristics over 1990-96 (hospital market structure, closures,
## margins, unemployment) by adoption cohort, with correlations with cohort
## year among the 1999-2001 states. Writes results/att_adoption_timing.tex (+ CSV).
laus <- read_csv("data/input/laus_state_unemployment_1990_2000.csv", show_col_types = FALSE) %>% rename(MSTATE = state)
u_state <- bind_rows(
  laus %>% filter(year >= 1994, year <= 1996) %>% group_by(MSTATE) %>% summarise(u_94_96 = mean(unemp_rate), .groups = "drop"),
  laus %>% filter(year >= 1997, year <= 1998) %>% group_by(MSTATE) %>% summarise(u_97_98 = mean(unemp_rate), .groups = "drop")) %>%
  group_by(MSTATE) %>% summarise(across(c(u_94_96, u_97_98), ~ first(na.omit(.x))), .groups = "drop") %>%
  mutate(d_u_late = u_97_98 - u_94_96)
st96 <- est.dat %>% filter(year == 1996) %>% group_by(MSTATE) %>%
  summarise(first_obs = first(state_treat_year), n_hosp = n_distinct(ID), n_small = sum(BDTOT <= 50, na.rm = TRUE),
            share_rural = mean(ever_rural, na.rm = TRUE), mean_beds = mean(BDTOT, na.rm = TRUE), mean_dist = mean(distance, na.rm = TRUE),
            share_gov = mean(own_gov, na.rm = TRUE), .groups = "drop")
cl <- est.dat %>% filter(year >= 1990, year <= 1996) %>% group_by(MSTATE) %>%
  summarise(closure_rate = 100 * sum(closed) / n(), .groups = "drop")
cl_change <- est.dat %>% mutate(period = case_when(year >= 1990 & year <= 1993 ~ "early", year >= 1994 & year <= 1996 ~ "late")) %>%
  filter(!is.na(period)) %>% group_by(MSTATE, period) %>% summarise(rate = 100 * sum(closed) / n(), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = rate) %>% mutate(d_closure = late - early) %>% select(MSTATE, d_closure)
mg <- est.dat %>% filter(year >= 1994, year <= 1996) %>% group_by(MSTATE) %>% summarise(margin_94_96 = mean(margin, na.rm = TRUE), .groups = "drop")
adopt <- st96 %>% left_join(cl, by = "MSTATE") %>% left_join(cl_change, by = "MSTATE") %>% left_join(mg, by = "MSTATE") %>% left_join(u_state, by = "MSTATE") %>%
  mutate(cohort = case_when(first_obs == 0 ~ "Never", first_obs < 1999 ~ "Before 1999", first_obs > 2001 ~ "After 2001", TRUE ~ as.character(first_obs)))
write_csv(adopt, "results/diagnostics/adoption-timing-states.csv")
adopt_vars <- c(share_rural = "Share of hospitals rural", n_small = "Hospitals with 50 or fewer beds", mean_beds = "Mean beds", share_gov = "Share government-owned",
                mean_dist = "Miles to nearest hospital", closure_rate = "Closures per 100 hospital-years, 1990-96", d_closure = "Change in closure rate, 1990-93 to 1994-96",
                margin_94_96 = "Operating margin, 1994-96", u_94_96 = "Unemployment rate, 1994-96", d_u_late = "Change in unemployment, 1994-96 to 1997-98")
cohort_order <- c("Before 1999", "1999", "2000", "2001", "After 2001", "Never")
means <- adopt %>% group_by(cohort) %>% summarise(n_states = n(), across(all_of(names(adopt_vars)), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(cohort = factor(cohort, levels = cohort_order)) %>% arrange(cohort)
core <- adopt %>% filter(first_obs %in% 1999:2001)
corrs <- core %>%
  summarise(across(all_of(names(adopt_vars)), ~ cor(.x, first_obs, use = "complete.obs"))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "corr")
## Joint regression of cohort year on closures and unemployment (levels and
## changes) with rurality and the number of small hospitals, 1999-2001 states.
adopt_ols <- lm(first_obs ~ closure_rate + d_closure + u_94_96 + d_u_late + share_rural + n_small, data = core)
write_csv(tidy(adopt_ols) %>% mutate(n_states = nobs(adopt_ols), r_squared = summary(adopt_ols)$r.squared),
          "results/diagnostics/adoption-timing-ols.csv")
## One row per state characteristic, one column per adoption cohort. Counts and
## bed/distance means print to one decimal, shares and rates to three.
cell <- function(v, k) {
  x <- means %>% filter(cohort == k) %>% pull(v)
  if (length(x) == 0 || is.na(x)) return("")
  fmt <- if (v %in% c("n_small", "mean_beds", "mean_dist")) "%.1f" else "%.3f"
  sprintf(fmt, x)
}

lines <- names(adopt_vars) %>%
  map_chr(function(v) {
    vals <- cohort_order %>% map_chr(~ cell(v, .x))
    corr <- corrs %>% filter(variable == v) %>% pull(corr)
    sprintf("%s & %s & %s \\\\", adopt_vars[[v]], paste(vals, collapse = " & "), sprintf("%.2f", corr))
  })

n_states <- cohort_order %>% map_chr(function(k) {
  x <- means %>% filter(cohort == k) %>% pull(n_states)
  if (length(x) == 0) "0" else as.character(x)
})
n_line <- sprintf("States & %s & \\\\", paste(n_states, collapse = " & "))
write_tabular(c(n_line, "\\addlinespace", lines),
              c(" & \\multicolumn{6}{c}{Adoption cohort} & Correlation with \\\\", "State characteristic, 1990--1996 & Before 1999 & 1999 & 2000 & 2001 & After 2001 & Never & cohort year, 1999--2001 \\\\"),
              "results/att_adoption_timing.tex", "lcccccc c")


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
  scen <- tibble(scenario = c("Mean converter before designation (baseline)", "Closers with 50 or fewer beds, mean", "Rural closers, mean",
                              "Smallest hospitals: 25th percentile of small hospitals", "Smallest hospitals: 10th percentile of small hospitals"),
                 B = c(B_close, mean(closer_beds_lag$beds_lag[closer_beds_lag$beds_lag <= 50]), mean(closer_beds_lag$beds_lag[closer_beds_lag$ever_rural == 1], na.rm = TRUE),
                       quantile(small98, 0.25, na.rm = TRUE), quantile(small98, 0.10, na.rm = TRUE))) %>%
    mutate(rho_star = (abs(delta_C) / 100 * B) / abs(delta_B), net_beds = (-delta_C / 100) * B + 0.30 * delta_B)
  write_csv(scen, "results/diagnostics/breakeven-sensitivity.csv")
  scen_lines <- scen %>%
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


## ===== R2.5 — Data source: shares, definitional gaps, and same-hospital comparison =====
## Source of each financial observation in the original combined series by
## converter status and period; the Form 990 vs HCRIS difference for the same
## hospital-year; and SDID on the combined series, the HCRIS-only series, and
## the Form-990-only series for the same hospitals. Writes
## results/att_source_comparison.tex and diagnostics CSVs.

## Whether a hospital-year has a cost report, taken from the linked HCRIS file
## rather than from est.dat, whose financial columns are interpolated.
hcris_present <- read_csv('data/output/hcris_financial.csv', show_col_types = FALSE,
                          col_types = cols(ID = col_character())) %>%
  transmute(ID, year,
            has_hcris = !is.na(net_pat_rev) & !is.na(tot_operating_exp),
            margin_hcris = ifelse(net_pat_rev > 0, (net_pat_rev - tot_operating_exp) / net_pat_rev, NA_real_))

src <- est.dat %>%
  filter(year >= 1995, year <= 2010) %>%
  group_by(ID) %>%
  mutate(mb = min(BDTOT, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(mb <= bed.cut) %>%
  left_join(hcris_present, by = c("ID", "year")) %>%
  mutate(grp = case_when(ever_cah == 1 & year >= eff_year ~ "Converters, after designation",
                         ever_cah == 1 ~ "Converters, before designation",
                         TRUE ~ "Never-converters"),
         source = case_when(has_hcris ~ "HCRIS",
                            !is.na(margin_990) ~ "Form 990",
                            TRUE ~ "None"))
shares <- src %>% count(grp, source) %>% group_by(grp) %>% mutate(share = n / sum(n)) %>% ungroup()
write_csv(shares, "results/diagnostics/source-shares.csv")
## Same hospital-year observed in both sources: margin gap, and the within
## hospital-year ratio of the Form 990 to the HCRIS measure for revenue,
## expenses, and net fixed assets on the identically winsorized and normalized
## per-bed series (*_9 vs HCRIS-only). The median ratio is reported because the
## ratio has a heavy right tail (Form 990 filings that cover more than one
## facility), summarized by the share of hospital-years with a ratio above 2.
overlap <- src %>% filter(!is.na(margin_hcris), !is.na(margin_990)) %>%
  group_by(year, ever_cah) %>%
  mutate(mh = winsor_by_year(margin_hcris), m9 = winsor_by_year(margin_990)) %>%
  ungroup()
src_ratio <- function(h, n) { ok <- !is.na(h) & !is.na(n) & h > 0; r <- n[ok] / h[ok]; c(median = median(r), share_above_2 = mean(r > 2), n = sum(ok)) }
gap <- tibble(hospital_years = nrow(overlap), hcris_margin = mean(overlap$mh), form990_margin = mean(overlap$m9), gap_mean = mean(overlap$m9 - overlap$mh), gap_median = median(overlap$m9 - overlap$mh)) %>%
  bind_cols(as_tibble_row(setNames(src_ratio(overlap$net_pat_rev, overlap$net_pat_rev_9), paste0("rev_ratio_", c("median", "share_above_2", "n")))),
            as_tibble_row(setNames(src_ratio(overlap$tot_operating_exp, overlap$tot_operating_exp_9), paste0("exp_ratio_", c("median", "share_above_2", "n")))),
            as_tibble_row(setNames(src_ratio(overlap$net_fixed, overlap$net_fixed_9), paste0("net_fixed_ratio_", c("median", "share_above_2", "n")))))
write_csv(gap, "results/diagnostics/source-gap.csv")

## Three series: the HCRIS series used throughout the paper, a Form-990-only
## version, and the HCRIS series restricted to converters with a strict
## crosswalk match (score >= 0.90, or zip-exact with score >= 0.75).
## Never-converters are kept in every series. pre_mean is the treated pre-period
## mean, the base for percent effects.
strict_ids <- est.dat %>% filter(is.na(eff_year) | (!is.na(match_strict) & match_strict)) %>% distinct(ID) %>% pull(ID)
src_rows <- list()
for (dsg in c("state", "elig")) {
  stack <- if (dsg == "state") stack.hosp else stack.elig; cohorts <- if (dsg == "state") 1999:2001 else 1999:2005
  for (o in fin.vars) {
    for (variant in c("h", "9", "hs")) {
      oname <- if (variant == "9") paste0(o, "_9") else o
      s <- rr_sdid(oname, stack, cohorts, pre = financial.pre, keep_ids = if (variant == "hs") strict_ids else NULL)
      if (!is.null(s)) src_rows[[length(src_rows) + 1]] <- tibble(outcome = o, design = dsg, series = c(h = "HCRIS", `9` = "Form 990 only", hs = "HCRIS, strict matches")[[variant]],
                                                                  att = s$att, lo = s$lo, hi = s$hi, ntr = s$ntr, pre_mean = s$pre_mean, pct_of_pre_mean = s$att / s$pre_mean)
    }
  }
}
## The per-source estimates support the R2.5 response only. The paper uses the
## HCRIS series throughout, so no table is written; the CSV carries every series.
src_res <- bind_rows(src_rows)
write_csv(src_res, "results/diagnostics/source-comparison.csv")


## ===== R1-minor(a) — Main SDID table with permutation p-values ===============
perm_file <- "results/diagnostics/permutation_pvalues.csv"
if (file.exists(perm_file)) {
  perm <- read_csv(perm_file, show_col_types = FALSE) %>% select(outcome, p_value)
  main_perm <- results.table %>% left_join(perm, by = "outcome")
  lines <- main_perm %>% rowwise() %>%
    mutate(line = sprintf("%s & %s & %s & %s & %s & %s & %d \\\\", outcome,
                          fmt3(sdid_att), ci3(sdid_ci_low, sdid_ci_high),
                          ifelse(is.na(p_value), "", sprintf("%.3f", p_value)),
                          ifelse(is.na(cs_att), "", fmt3(cs_att)),
                          ifelse(is.na(cs_att), "", ci3(cs_ci_low, cs_ci_high)),
                          sdid_ntr)) %>% pull(line)
  lines <- append(lines, "\\addlinespace", after = 6); lines <- append(lines, "\\addlinespace", after = 11)
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
## The pooled CS closure estimate is smaller in magnitude than SDID. This shows
## where that comes from: the 2000 and 2001 cohorts, where the not-yet-treated
## control pool is smallest. CS is estimated once over all three cohorts, as in
## cs_att(), and the group-time estimates are then averaged over each cohort's
## own five post-treatment years so the horizon matches SDID. No CS standard
## error is reported, since averaging ATT(g,t) by hand does not carry one.

cs.closure <- state.dat %>%
  mutate(y = closures,
         ID2 = as.numeric(factor(MSTATE)),
         treat_group = ifelse(is.na(state_treat_year), 0, state_treat_year)) %>%
  filter(!is.na(y), treat_group == 0 | treat_group %in% 1999:2001) %>%
  select(ID2, treat_group, year, y)

cs.closure.raw <- tryCatch(
  att_gt(yname = "y", gname = "treat_group", idname = "ID2", tname = "year",
         control_group = "notyettreated", panel = TRUE, allow_unbalanced_panel = TRUE,
         data = as.data.frame(cs.closure), xformla = ~1,
         base_period = "universal", est_method = "reg"),
  error = function(e) NULL)

closure_cohort <- tibble()
for (c in 1999:2001) {
  synth.c <- stack.state %>%
    filter(stack_group == c) %>%
    transmute(ID = as.numeric(factor(MSTATE)), year, outcome = closures, treated = post_treat)

  bal.c <- tryCatch(as_tibble(makeBalancedPanel(synth.c, idname = "ID", tname = "year")),
                    error = function(e) NULL)
  if (is.null(bal.c) || nrow(bal.c) == 0) next
  setup.c <- tryCatch(panel.matrices(as.data.frame(bal.c)), error = function(e) NULL)
  if (is.null(setup.c)) next
  sdid.c <- tryCatch(synthdid_estimate(setup.c$Y, setup.c$N0, setup.c$T0), error = function(e) NULL)
  if (is.null(sdid.c)) next

  cs.c <- NA_real_
  if (!is.null(cs.closure.raw)) {
    keep <- cs.closure.raw$group == c &
      cs.closure.raw$t >= c & cs.closure.raw$t <= (c + post)
    if (any(keep)) cs.c <- mean(cs.closure.raw$att[keep], na.rm = TRUE)
  }

  closure_cohort <- bind_rows(closure_cohort, tibble(
    cohort           = c,
    n_treated_states = nrow(setup.c$Y) - setup.c$N0,
    n_control_states = setup.c$N0,
    sdid_att         = as.numeric(sdid.c),
    sdid_se          = tryCatch(as.numeric(synthdid_se(sdid.c, method = "jackknife")),
                                error = function(e) NA_real_),
    cs_att           = cs.c))
}

if (nrow(closure_cohort) > 0) {
  write_csv(closure_cohort, "results/diagnostics/closure-cohort-cs-sdid.csv")
} else {
  message("E5 closure-by-cohort produced no rows, skipping")
}

