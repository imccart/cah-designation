# Pre-trend sensitivity for the SDID and CS event-time effects, state-timing
# design, following Rambachan and Roth (2023). For each outcome and estimator
# the cohort series are pooled with treated-unit-count weights (cohorts treated
# as independent, as the pooled ATTs are), recentered on the year before
# designation so both estimators share a reference period, and the average
# effect over the five post-designation years is bounded under the smoothness
# restriction Delta^SD(M): the differential trend observed before designation
# continues afterward with a slope that changes by at most M per period, so
# M = 0 is the straight-line extension of the pre-trend. M is scaled in
# multiples of the largest change in slope observed in the pre-period.
#
# Expects in scope: hosp.event.cohort, hosp.event.vcov (2-hospital-dd.R),
# visit.event.cohort, visit.event.vcov (R1_revision.R, state-timing rows), and
# results.table and visit_res for the main-text intervals.
# Writes results/pretrend-sens-<stub>.png, results/diagnostics/pretrend-sensitivity.csv,
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

sens.cohort <- bind_rows(hosp.event.cohort,
                         visit.event.cohort %>% filter(design == "state") %>% select(-design))
sens.vcov   <- bind_rows(hosp.event.vcov,
                         visit.event.vcov %>% filter(design == "state") %>% select(-design))

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

sens_rows <- list()
for (out in names(sens_outcomes)) {
  message(sprintf("pretrend sensitivity: %s", out))
  for (estr in c("SDID", "CS")) {
    p <- pooled_recentered(out, estr)
    if (is.null(p)) next
    l_vec <- rep(1 / p$numPost, p$numPost)
    pre_b <- c(p$beta[seq_len(p$numPre)], 0)
    second_diff <- if (length(pre_b) >= 3) max(abs(diff(diff(pre_b)))) else abs(diff(pre_b))
    M_ref <- if (is.finite(second_diff) && second_diff > 0) second_diff else sqrt(p$sigma[p$numPre + 1, p$numPre + 1])
    Mvec  <- M_ref * c(0, 0.1, 0.25, 0.5, 1)

    orig <- tryCatch(constructOriginalCS(p$beta, p$sigma, p$numPre, p$numPost, l_vec = l_vec),
                     error = function(err) NULL)
    sd   <- tryCatch(createSensitivityResults(p$beta, p$sigma, p$numPre, p$numPost,
                                              Mvec = Mvec, l_vec = l_vec),
                     error = function(err) { message(sprintf("SD failed for %s / %s: %s", out, estr, conditionMessage(err))); NULL })
    ## HonestDiD returns its bounds as one-column matrices; flatten them so the
    ## rows bind and write as plain numbers.
    if (!is.null(orig)) sens_rows[[length(sens_rows) + 1]] <- tibble(outcome = out, estimator = estr, restriction = "original",
                                                                     M = NA_real_, M_mult = NA_real_,
                                                                     lb = as.numeric(orig$lb), ub = as.numeric(orig$ub))
    if (!is.null(sd))   sens_rows[[length(sens_rows) + 1]] <- tibble(outcome = out, estimator = estr, restriction = "smoothness",
                                                                     M = as.numeric(sd$M), M_mult = as.numeric(sd$M) / M_ref,
                                                                     lb = as.numeric(sd$lb), ub = as.numeric(sd$ub))
  }
}
sens.results <- bind_rows(sens_rows)
write_csv(sens.results, "results/diagnostics/pretrend-sensitivity.csv")

## The paper's own estimates and intervals (Table 2 and the visits table) are
## the reference the robust intervals are shown against.
table2_ci <- bind_rows(
  results.table %>% transmute(outcome, SDID_lb = sdid_ci_low, SDID_ub = sdid_ci_high,
                              CS_lb = cs_ci_low, CS_ub = cs_ci_high),
  visit_res %>% filter(design == "state") %>%
    transmute(outcome, SDID_lb = sdid_lo, SDID_ub = sdid_hi, CS_lb = cs_lo, CS_ub = cs_hi)
)

## Figures: one per outcome, confidence sets as M grows, both estimators, with
## the paper's intervals as reference bands.
for (out in names(sens_outcomes)) {
  stub <- sens_outcomes[[out]]
  rob <- sens.results %>% filter(outcome == out, restriction == "smoothness") %>%
    mutate(estimator = factor(estimator, levels = c("SDID", "CS")))
  if (nrow(rob) == 0) next
  ref <- table2_ci %>% filter(outcome == out)
  p <- ggplot(rob, aes(x = M_mult, y = (lb + ub) / 2, shape = estimator, linetype = estimator)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = ref$SDID_lb, ymax = ref$SDID_ub, fill = "grey80", alpha = 0.5) +
    geom_hline(yintercept = c(ref$CS_lb, ref$CS_ub), linetype = "dotted", colour = "grey40") +
    geom_hline(yintercept = 0, linewidth = 0.8) +
    geom_errorbar(aes(ymin = lb, ymax = ub), position = position_dodge(width = 0.06), width = 0, linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.06), size = 2.4, fill = "white") +
    scale_shape_manual(values = c(SDID = 16, CS = 24)) +
    scale_linetype_manual(values = c(SDID = "solid", CS = "dashed")) +
    scale_x_continuous(breaks = c(0, 0.1, 0.25, 0.5, 1)) +
    labs(x = "M (multiples of the largest pre-designation change in slope)",
         y = paste0(out, ": average effect, years 0 to 5"), shape = NULL, linetype = NULL) +
    theme_bw() +
    theme(legend.position = "inside", legend.position.inside = c(0.12, 0.12),
          legend.background = element_blank())
  ggsave(sprintf("results/pretrend-sens-%s.png", stub), p, width = 6.5, height = 4.25, dpi = 300, scale = 1.5)
}

## Table: the paper's interval and the interval under the straight-line
## extension of the pre-designation trend (smoothness, M = 0), per outcome and
## estimator.
fmt_ci <- function(lo, hi) {
  ax <- max(abs(c(lo, hi)))
  dg <- ifelse(ax < 5, 2, ifelse(ax < 100, 1, 0))
  sprintf("[%s, %s]", formatC(lo, format = "f", digits = dg), formatC(hi, format = "f", digits = dg))
}
tab_lines <- c()
for (out in names(sens_outcomes)) {
  cells <- c()
  for (estr in c("SDID", "CS")) {
    d   <- sens.results %>% filter(outcome == out, estimator == estr)
    ref <- table2_ci %>% filter(outcome == out)
    m0  <- d %>% filter(restriction == "smoothness", M_mult == 0)
    cells <- c(cells,
               if (nrow(ref) == 1) fmt_ci(ref[[paste0(estr, "_lb")]], ref[[paste0(estr, "_ub")]]) else "",
               if (nrow(m0) == 1) fmt_ci(m0$lb, m0$ub) else "")
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
  "Outcome & 95\\% CI & Extended pre-trend & 95\\% CI & Extended pre-trend \\\\",
  "\\midrule",
  tab_lines,
  "\\bottomrule",
  "\\end{tabular}"
), "results/att_pretrend.tex")
