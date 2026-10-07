# Trend-adjusted estimates for the SDID and CS event-time effects, state-timing
# design, following Rambachan and Roth (2023). For each outcome and estimator
# the cohort series are pooled with treated-unit-count weights (cohorts treated
# as independent, as the pooled ATTs are) and recentered on the year before
# designation so both estimators share a reference period. The differential
# trend observed before designation is extended into the post-period as a
# straight line, and the average effect over the five post-designation years is
# reported net of that extension, with a confidence interval that carries the
# uncertainty in the extrapolated trend (the linear-extrapolation case of the
# smoothness restriction in HonestDiD).
#
# Expects in scope: hosp.event.cohort, hosp.event.vcov (2-hospital-dd.R, visit
# outcomes included), and results.table and visit_res (R1_revision.R) for the
# main-text estimates and intervals.
# Writes results/pretrend-adjusted.png, results/diagnostics/pretrend-adjusted.csv,
# results/att_pretrend.tex.

sens_outcomes <- c(
  "Operating margin"             = "margin",
  "Current ratio"                = "currentratio",
  "Net fixed assets"             = "netfixed",
  "Capital expenditures per bed" = "capex",
  "Net patient revenue per bed"  = "netpatrev",
  "Operating expenses per bed"   = "totopexp",
  "Total beds"                   = "beds",
  "OB beds"                      = "beds_ob",
  "FTE RNs"                      = "ftern",
  "Inpatient days per bed"       = "ipdays",
  "System membership"            = "system",
  "ED visits per bed"            = "edvisits",
  "Outpatient visits per bed"    = "opvisits",
  "Non-ED outpatient visits per bed" = "opvisits-noned",
  "Births"                       = "births"
)

sens.cohort <- hosp.event.cohort
sens.vcov   <- hosp.event.vcov

## Pooled, recentered coefficients and covariance for one outcome and estimator.
## The e grid is the intersection across cohorts (capex 1999 starts at -2).
pooled_recentered <- function(out, estr) {
  cs <- sens.cohort %>% filter(outcome == out, estimator == estr)
  vs <- sens.vcov   %>% filter(outcome == out, estimator == estr)
  if (nrow(cs) == 0 || nrow(vs) == 0) return(NULL)
  cohorts <- sort(unique(cs$cohort))
  e_grid  <- reduce(map(cohorts, function(c) cs$e[cs$cohort == c]), intersect) %>% sort()
  if (!(-1 %in% e_grid) || sum(e_grid < -1) < 1) return(NULL)
  w <- cs %>% distinct(cohort, Ntr) %>% mutate(w = Ntr / sum(Ntr))
  beta <- rep(0, length(e_grid)); V <- matrix(0, length(e_grid), length(e_grid))
  for (c in cohorts) {
    wc <- w$w[w$cohort == c]
    bc <- cs %>% filter(cohort == c, e %in% e_grid)
    beta <- beta + wc * bc$estimate[match(e_grid, bc$e)]
    Vc <- vs %>% filter(cohort == c, e1 %in% e_grid, e2 %in% e_grid)
    M  <- matrix(0, length(e_grid), length(e_grid))
    M[cbind(match(Vc$e1, e_grid), match(Vc$e2, e_grid))] <- Vc$value
    V <- V + wc^2 * M
  }
  ## recenter on e = -1 and drop it
  ref <- which(e_grid == -1)
  A <- diag(length(e_grid)); A[, ref] <- A[, ref] - 1
  A <- A[-ref, , drop = FALSE]
  list(e = e_grid[-ref], beta = as.numeric(A %*% beta), sigma = A %*% V %*% t(A),
       numPre = sum(e_grid < -1), numPost = sum(e_grid >= 0))
}

## The paper's own estimates and intervals (Table 2 and the visits table).
table2 <- bind_rows(
  results.table %>% transmute(outcome, SDID_est = sdid_att, SDID_lb = sdid_ci_low, SDID_ub = sdid_ci_high,
                              CS_est = cs_att, CS_lb = cs_ci_low, CS_ub = cs_ci_high),
  visit_res %>% filter(design == "state") %>%
    transmute(outcome, SDID_est = sdid_att, SDID_lb = sdid_lo, SDID_ub = sdid_hi,
              CS_est = cs_att, CS_lb = cs_lo, CS_ub = cs_hi)
)

## The fixed-length confidence interval under the linear extension is symmetric
## about its estimator, so its midpoint is the trend-adjusted estimate.
sens_rows <- list()
for (out in names(sens_outcomes)) {
  message(sprintf("trend-adjusted estimates: %s", out))
  for (estr in c("SDID", "CS")) {
    p <- pooled_recentered(out, estr)
    if (is.null(p)) next
    l_vec <- rep(1 / p$numPost, p$numPost)
    adj <- tryCatch(createSensitivityResults(p$beta, p$sigma, p$numPre, p$numPost, Mvec = 0, l_vec = l_vec),
                    error = function(err) { message(sprintf("failed for %s / %s: %s", out, estr, conditionMessage(err))); NULL })
    if (is.null(adj)) next
    ref <- table2 %>% filter(outcome == out)
    sens_rows[[length(sens_rows) + 1]] <- tibble(
      outcome = out, estimator = estr,
      main_est = ref[[paste0(estr, "_est")]], main_lb = ref[[paste0(estr, "_lb")]], main_ub = ref[[paste0(estr, "_ub")]],
      adj_lb = as.numeric(adj$lb), adj_ub = as.numeric(adj$ub)) %>%
      mutate(adj_est = (adj_lb + adj_ub) / 2)
  }
}
sens.results <- bind_rows(sens_rows)
write_csv(sens.results, "results/diagnostics/pretrend-adjusted.csv")

## Figure: every outcome in its own panel, the main-text estimate beside the
## trend-adjusted estimate for each estimator.
fig_dat <- sens.results %>%
  pivot_longer(c(main_est, adj_est), names_to = "version", values_to = "est") %>%
  mutate(lb = ifelse(version == "main_est", main_lb, adj_lb),
         ub = ifelse(version == "main_est", main_ub, adj_ub),
         version = factor(ifelse(version == "main_est", "Main text", "Trend-adjusted"),
                          levels = c("Main text", "Trend-adjusted")),
         estimator = factor(estimator, levels = c("SDID", "CS")),
         outcome = factor(outcome, levels = names(sens_outcomes)))
p <- ggplot(fig_dat, aes(x = estimator, y = est, shape = estimator, fill = version, group = version)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_errorbar(aes(ymin = lb, ymax = ub), position = position_dodge(width = 0.6), width = 0, linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 2.4) +
  scale_shape_manual(values = c(SDID = 21, CS = 24), guide = "none") +
  scale_fill_manual(values = c(`Main text` = "white", `Trend-adjusted` = "black")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  facet_wrap(~ outcome, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Average effect, years 0 to 5", fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom", strip.text = element_text(size = 8))
ggsave("results/pretrend-adjusted.png", p, width = 6.5, height = 8, dpi = 300, scale = 1.3)

## Table: the paper's interval and the trend-adjusted estimate with its
## interval, per outcome and estimator.
fmt_ci <- function(lo, hi) {
  ax <- max(abs(c(lo, hi)))
  dg <- ifelse(ax < 5, 2, ifelse(ax < 100, 1, 0))
  sprintf("[%s, %s]", formatC(lo, format = "f", digits = dg), formatC(hi, format = "f", digits = dg))
}
fmt_est <- function(x, lo, hi) {
  ax <- max(abs(c(lo, hi)))
  dg <- ifelse(ax < 5, 2, ifelse(ax < 100, 1, 0))
  sprintf("%s %s", formatC(x, format = "f", digits = dg), fmt_ci(lo, hi))
}
tab_lines <- c()
for (out in names(sens_outcomes)) {
  cells <- c()
  for (estr in c("SDID", "CS")) {
    d <- sens.results %>% filter(outcome == out, estimator == estr)
    cells <- c(cells,
               if (nrow(d) == 1) fmt_ci(d$main_lb, d$main_ub) else "",
               if (nrow(d) == 1) fmt_est(d$adj_est, d$adj_lb, d$adj_ub) else "")
  }
  if (all(cells == "")) next
  tab_lines <- c(tab_lines, sprintf("%s & %s \\\\", out, paste(cells, collapse = " & ")))
}
tab_lines <- append(tab_lines, "\\addlinespace", after = 6)
tab_lines <- append(tab_lines, "\\addlinespace", after = 11)
writeLines(c(
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{SDID} & \\multicolumn{2}{c}{Callaway--Sant'Anna} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "Outcome & 95\\% CI & Trend-adjusted & 95\\% CI & Trend-adjusted \\\\",
  "\\midrule",
  tab_lines,
  "\\bottomrule",
  "\\end{tabular}"
), "results/att_pretrend.tex")
