# Summary Statistics -------------------------------------------------------

## Count of CAH and Non-CAH over time
fig.hosp.type <- est.dat %>% group_by(year, cah) %>%
  summarize(hosp_count=n())  %>%
  ggplot(aes(x=year, y=hosp_count, group=cah)) +
  geom_line() + geom_point() + theme_bw() +
  geom_text(data = est.dat %>% filter(year==2016) %>% group_by(cah, year) %>% summarize(hosp_count=n()),
            aes(label = c("Non-CAH", "CAH"),
                x = year,
                y = hosp_count-110)) +
  scale_y_continuous(labels = comma,
                     breaks=seq(0, 8000, 500)) +
  scale_x_continuous(breaks=seq(1980, 2019, 5)) +
  guides(linetype="none") +
  labs(
    x="Year",
    y="Count of Hospitals"
  )
ggsave("results/desc-hosp-types.png", fig.hosp.type, width = 6.5, height = 4.25, dpi = 300, scale=1.5)

## Section 3.3 quotes the panel's hospital count and the level at which the
## count of CAHs settles, so the series behind the figure goes to a file.
est.dat %>%
  count(year, cah, name = "hosp_count") %>%
  mutate(panel_hospitals = n_distinct(est.dat$ID),
         panel_year_min = min(est.dat$year),
         panel_year_max = max(est.dat$year)) %>%
  write_csv("results/diagnostics/hospital-counts-by-year.csv")


## Treatment timing --------------------------------------------------------

## panelView 1.3.1 needs a real column on the left of the formula (a literal
## `1 ~ treat_post` drops the treatment indicator), so `one` is a constant.
panel.cah <- panelview(one~treat_post, data=state.dat %>% mutate(MSTATE=as.character(MSTATE), one=1),
                      index=c("MSTATE","year"), legendOff=TRUE,
          theme.bw=TRUE, by.timing=TRUE, xlab="Year", ylab="State",
          main="", color=c("white","gray"), axis.lab.angle=45)
ggsave("results/desc-panelview-treat.png", panel.cah, width = 6.5, height = 4.25, dpi = 300, scale=2)


## Characteristics of CAH and non-CAH Hospitals ---------------------------------------------------

cah.desc <- est.dat %>% filter(year>=1995 & year<=2010) %>%
  filter(ever_cah==1) %>%
  summarize(
    bed_mean   = mean(BDTOT, na.rm = TRUE),
    bed_p10    = quantile(BDTOT, 0.10, na.rm = TRUE),
    bed_p90    = quantile(BDTOT, 0.90, na.rm = TRUE),
    ob_mean    = mean(OBBD, na.rm = TRUE),
    ob_p10     = quantile(OBBD, 0.10, na.rm = TRUE),
    ob_p90     = quantile(OBBD, 0.90, na.rm = TRUE),
    ftern_mean = mean(FTERN, na.rm = TRUE),
    ftern_p10  = quantile(FTERN, 0.10, na.rm = TRUE),
    ftern_p90  = quantile(FTERN, 0.90, na.rm = TRUE),
    dist_mean  = mean(distance, na.rm = TRUE),
    dist_p10   = quantile(distance, 0.10, na.rm = TRUE),
    dist_p90   = quantile(distance, 0.90, na.rm = TRUE),
    ip_days    = mean(ip_per_bed, na.rm = TRUE),
    ip_days_p10= quantile(ip_per_bed, 0.10, na.rm = TRUE),
    ip_days_p90= quantile(ip_per_bed, 0.90, na.rm = TRUE),
    rural      = mean(ever_rural, na.rm=TRUE),
    margin_mean=mean(margin, na.rm=TRUE),
    margin_p10 =quantile(margin, 0.10, na.rm=TRUE),
    margin_p90 =quantile(margin, 0.90, na.rm=TRUE),
    cr_mean    = mean(current_ratio, na.rm=TRUE),
    cr_p10     = quantile(current_ratio, 0.10, na.rm=TRUE),
    cr_p90     = quantile(current_ratio, 0.90, na.rm=TRUE),
    system_mean= mean(system, na.rm=TRUE),
    closures   = sum(closed, na.rm=TRUE),
    mergers    = sum(merged, na.rm=TRUE),
    hospitals  = n_distinct(ID),
    obs        = n()
  ) %>%
  mutate(rural = rural * 100,
         system_mean = system_mean * 100) %>%
  mutate(across(c(ip_days, ip_days_p10, ip_days_p90), ~round(., 2))) %>%
  mutate(across(-c(closures, mergers, hospitals, obs), ~round(., 2))) %>%
  mutate(across(c(bed_p10, bed_p90, ob_p10, ob_p90), ~round(.,0))) %>%
  rename_with(~paste0(., "_cah"))


cah.pre <- est.dat %>% filter(year>=1995 & year<=2010) %>%
  filter(ever_cah==1, eff_year > year) %>%
  summarize(
    bed_mean   = mean(BDTOT, na.rm = TRUE),
    bed_p10    = quantile(BDTOT, 0.10, na.rm = TRUE),
    bed_p90    = quantile(BDTOT, 0.90, na.rm = TRUE),
    ob_mean    = mean(OBBD, na.rm = TRUE),
    ob_p10     = quantile(OBBD, 0.10, na.rm = TRUE),
    ob_p90     = quantile(OBBD, 0.90, na.rm = TRUE),
    ftern_mean = mean(FTERN, na.rm = TRUE),
    ftern_p10  = quantile(FTERN, 0.10, na.rm = TRUE),
    ftern_p90  = quantile(FTERN, 0.90, na.rm = TRUE),
    dist_mean  = mean(distance, na.rm = TRUE),
    dist_p10   = quantile(distance, 0.10, na.rm = TRUE),
    dist_p90   = quantile(distance, 0.90, na.rm = TRUE),
    ip_days    = mean(ip_per_bed, na.rm = TRUE),
    ip_days_p10= quantile(ip_per_bed, 0.10, na.rm = TRUE),
    ip_days_p90= quantile(ip_per_bed, 0.90, na.rm = TRUE),
    rural      = mean(ever_rural, na.rm=TRUE),
    margin_mean=mean(margin, na.rm=TRUE),
    margin_p10 =quantile(margin, 0.10, na.rm=TRUE),
    margin_p90 =quantile(margin, 0.90, na.rm=TRUE),
    cr_mean    = mean(current_ratio, na.rm=TRUE),
    cr_p10     = quantile(current_ratio, 0.10, na.rm=TRUE),
    cr_p90     = quantile(current_ratio, 0.90, na.rm=TRUE),
    system_mean= mean(system, na.rm=TRUE),
    closures   = sum(closed, na.rm=TRUE),
    mergers    = sum(merged, na.rm=TRUE),
    hospitals  = n_distinct(ID),
    obs        = n()
  ) %>%
  mutate(rural = rural * 100,
         system_mean = system_mean * 100) %>%
  mutate(across(c(ip_days, ip_days_p10, ip_days_p90), ~round(., 2))) %>%
  mutate(across(-c(closures, mergers, hospitals, obs), ~round(., 2))) %>%
  mutate(across(c(bed_p10, bed_p90, ob_p10, ob_p90), ~round(.,0))) %>%
  rename_with(~paste0(., "_cahpre"))


cah.post <- est.dat %>% filter(year>=1995 & year<=2010) %>%
  filter(ever_cah==1, eff_year <= year) %>%
  summarize(
    bed_mean   = mean(BDTOT, na.rm = TRUE),
    bed_p10    = quantile(BDTOT, 0.10, na.rm = TRUE),
    bed_p90    = quantile(BDTOT, 0.90, na.rm = TRUE),
    ob_mean    = mean(OBBD, na.rm = TRUE),
    ob_p10     = quantile(OBBD, 0.10, na.rm = TRUE),
    ob_p90     = quantile(OBBD, 0.90, na.rm = TRUE),
    ftern_mean = mean(FTERN, na.rm = TRUE),
    ftern_p10  = quantile(FTERN, 0.10, na.rm = TRUE),
    ftern_p90  = quantile(FTERN, 0.90, na.rm = TRUE),
    dist_mean  = mean(distance, na.rm = TRUE),
    dist_p10   = quantile(distance, 0.10, na.rm = TRUE),
    dist_p90   = quantile(distance, 0.90, na.rm = TRUE),
    ip_days    = mean(ip_per_bed, na.rm = TRUE),
    ip_days_p10= quantile(ip_per_bed, 0.10, na.rm = TRUE),
    ip_days_p90= quantile(ip_per_bed, 0.90, na.rm = TRUE),
    rural      = mean(ever_rural, na.rm=TRUE),
    margin_mean=mean(margin, na.rm=TRUE),
    margin_p10 =quantile(margin, 0.10, na.rm=TRUE),
    margin_p90 =quantile(margin, 0.90, na.rm=TRUE),
    cr_mean    = mean(current_ratio, na.rm=TRUE),
    cr_p10     = quantile(current_ratio, 0.10, na.rm=TRUE),
    cr_p90     = quantile(current_ratio, 0.90, na.rm=TRUE),
    system_mean= mean(system, na.rm=TRUE),
    closures   = sum(closed, na.rm=TRUE),
    mergers    = sum(merged, na.rm=TRUE),
    hospitals  = n_distinct(ID),
    obs        = n()
  ) %>%
  mutate(rural = rural * 100,
         system_mean = system_mean * 100) %>%
  mutate(across(c(ip_days, ip_days_p10, ip_days_p90), ~round(., 2))) %>%
  mutate(across(-c(closures, mergers, hospitals, obs), ~round(., 2))) %>%
  mutate(across(c(bed_p10, bed_p90, ob_p10, ob_p90), ~round(.,0))) %>%
  rename_with(~paste0(., "_cahpost"))


non.cah.desc <- est.dat %>% filter(year>=1995 & year<=2010) %>%
  filter(ever_cah==0) %>%
  summarize(
    bed_mean   = mean(BDTOT, na.rm = TRUE),
    bed_p10    = quantile(BDTOT, 0.10, na.rm = TRUE),
    bed_p90    = quantile(BDTOT, 0.90, na.rm = TRUE),
    ob_mean    = mean(OBBD, na.rm = TRUE),
    ob_p10     = quantile(OBBD, 0.10, na.rm = TRUE),
    ob_p90     = quantile(OBBD, 0.90, na.rm = TRUE),
    ftern_mean = mean(FTERN, na.rm = TRUE),
    ftern_p10  = quantile(FTERN, 0.10, na.rm = TRUE),
    ftern_p90  = quantile(FTERN, 0.90, na.rm = TRUE),
    dist_mean  = mean(distance, na.rm = TRUE),
    dist_p10   = quantile(distance, 0.10, na.rm = TRUE),
    dist_p90   = quantile(distance, 0.90, na.rm = TRUE),
    ip_days    = mean(ip_per_bed, na.rm = TRUE),
    ip_days_p10= quantile(ip_per_bed, 0.10, na.rm = TRUE),
    ip_days_p90= quantile(ip_per_bed, 0.90, na.rm = TRUE),
    rural      = mean(ever_rural, na.rm=TRUE),
    margin_mean=mean(margin, na.rm=TRUE),
    margin_p10 =quantile(margin, 0.10, na.rm=TRUE),
    margin_p90 =quantile(margin, 0.90, na.rm=TRUE),
    cr_mean    = mean(current_ratio, na.rm=TRUE),
    cr_p10     = quantile(current_ratio, 0.10, na.rm=TRUE),
    cr_p90     = quantile(current_ratio, 0.90, na.rm=TRUE),
    system_mean= mean(system, na.rm=TRUE),
    closures   = sum(closed, na.rm=TRUE),
    mergers    = sum(merged, na.rm=TRUE),
    hospitals  = n_distinct(ID),
    obs        = n()
  ) %>%
  mutate(rural = rural * 100,
         system_mean = system_mean * 100) %>%
  mutate(across(c(ip_days, ip_days_p10, ip_days_p90), ~round(., 2))) %>%
  mutate(across(-c(closures, mergers, hospitals, obs), ~round(., 2))) %>%
  mutate(across(c(bed_p10, bed_p90, ob_p10, ob_p90), ~round(.,0))) %>%
  rename_with(~paste0(., "_noncah"))

int <- function(lo, hi) sprintf("{}[%.2f, %.2f]", lo, hi)
int0 <- function(lo, hi) sprintf("{}[%.0f, %.0f]", lo, hi)

make_col <- function(x, suf) {
  c(
    x[[paste0("rural_",       suf)]],
    x[[paste0("bed_mean_",    suf)]],
    int0(x[[paste0("bed_p10_", suf)]],    x[[paste0("bed_p90_", suf)]]),
    x[[paste0("ob_mean_",     suf)]],
    int0(x[[paste0("ob_p10_", suf)]],     x[[paste0("ob_p90_", suf)]]),
    x[[paste0("ftern_mean_",  suf)]],
    int(x[[paste0("ftern_p10_",suf)]],    x[[paste0("ftern_p90_",suf)]]),
    x[[paste0("dist_mean_",   suf)]],
    int(x[[paste0("dist_p10_",suf)]],    x[[paste0("dist_p90_",suf)]]),
    x[[paste0("ip_days_",     suf)]],
    int(x[[paste0("ip_days_p10_",suf)]], x[[paste0("ip_days_p90_",suf)]]),
    x[[paste0("margin_mean_", suf)]],
    int(x[[paste0("margin_p10_",suf)]],  x[[paste0("margin_p90_",suf)]]),
    x[[paste0("cr_mean_",     suf)]],
    int(x[[paste0("cr_p10_",  suf)]],    x[[paste0("cr_p90_",  suf)]]),
    x[[paste0("system_mean_", suf)]],
    x[[paste0("closures_",    suf)]],
    x[[paste0("mergers_",     suf)]],
    format(x[[paste0("hospitals_",   suf)]], big.mark=",", trim=TRUE),
    format(x[[paste0("obs_",         suf)]], big.mark=",", trim=TRUE)
  )
}

desc_tab <- tibble(
  Variable = c(
    "Rural (\\%)",
    "Bed size", "",
    "OB beds", "",
    "FTE RNs", "",
    "Distance to nearest hospital", "",
    "Inpatient days per bed", "",
    "Operating margin", "",
    "Current ratio", "",
    "System member (\\%)",
    "Closures",
    "Mergers",
    "Hospitals",
    "Observations"
  ),
  `Pre`       = make_col(cah.pre,      "cahpre"),
  `Post`      = make_col(cah.post,     "cahpost"),
  `Ever CAH`  = make_col(cah.desc,     "cah"),
  `Never CAH` = make_col(non.cah.desc, "noncah")
)

# Complete LaTeX tabular block (LaTeX 2025 breaks \noalign/\omit inside \input within tabular)
tex_lines <- c(
  "\\begin{tabular}[t]{lcccc}",
  "\\toprule",
  "\\multicolumn{1}{c}{ } & \\multicolumn{3}{c}{CAHs} & \\multicolumn{1}{c}{ } \\\\",
  "\\cmidrule(l{3pt}r{3pt}){2-4}",
  " & Pre & Post & Ever CAH & Never CAH\\\\",
  "\\midrule"
)

for (i in seq_len(nrow(desc_tab))) {
  # Add spacing before Hospitals row to separate counts from variables
  if (desc_tab$Variable[i] == "Hospitals") tex_lines <- c(tex_lines, "\\addlinespace")
  row_vals <- c(desc_tab$Variable[i], desc_tab$Pre[i], desc_tab$Post[i],
                desc_tab$`Ever CAH`[i], desc_tab$`Never CAH`[i])
  tex_lines <- c(tex_lines, paste(row_vals, collapse = " & ") %>% paste0(" \\\\"))
}

tex_lines <- c(tex_lines, "\\bottomrule", "\\end{tabular}")

writeLines(tex_lines, "results/desc_compare.tex")


## Anticipation figure: capacity trajectories for treated hospitals ----------

antic_treated <- est.dat %>%
  filter(!is.na(eff_year), eff_year %in% 1999:2005) %>%
  group_by(ID) %>%
  mutate(min_bedsize = min(BDTOT, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(min_bedsize <= 50) %>%
  mutate(event_time = year - eff_year) %>%
  filter(event_time >= -5, event_time <= 5)

antic_beds <- antic_treated %>%
  filter(!is.na(BDTOT)) %>%
  group_by(event_time) %>%
  summarise(
    mean_val = mean(BDTOT, na.rm = TRUE),
    p25 = quantile(BDTOT, 0.25, na.rm = TRUE),
    p75 = quantile(BDTOT, 0.75, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

antic_ipd <- antic_treated %>%
  filter(!is.na(ip_per_bed)) %>%
  group_by(event_time) %>%
  summarise(
    mean_val = mean(ip_per_bed, na.rm = TRUE),
    p25 = quantile(ip_per_bed, 0.25, na.rm = TRUE),
    p75 = quantile(ip_per_bed, 0.75, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p_antic_beds <- ggplot(antic_beds, aes(x = event_time)) +
  geom_ribbon(aes(ymin = p25, ymax = p75), fill = "gray85", alpha = 0.6) +
  geom_hline(yintercept = 25, linetype = "dotted", linewidth = 0.5,
             color = "gray40") +
  geom_vline(xintercept = -0.5, linewidth = 1) +
  geom_line(aes(y = mean_val), linewidth = 0.8) +
  geom_point(aes(y = mean_val), size = 1.8) +
  annotate("text", x = 5, y = 27, label = "25-bed statutory limit",
           size = 3, color = "gray40", hjust = 1) +
  scale_x_continuous(breaks = -5:5) +
  scale_y_continuous(limits = c(20, NA)) +
  labs(x = "Years relative to CAH designation",
       y = "Total beds") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave("results/anticipation-beds.png", p_antic_beds,
       width = 6.5, height = 4.25, dpi = 300, scale = 1.5)

p_antic_ipd <- ggplot(antic_ipd, aes(x = event_time)) +
  geom_ribbon(aes(ymin = p25, ymax = p75), fill = "gray85", alpha = 0.6) +
  geom_vline(xintercept = -0.5, linewidth = 1) +
  geom_line(aes(y = mean_val), linewidth = 0.8) +
  geom_point(aes(y = mean_val), size = 1.8) +
  scale_x_continuous(breaks = -5:5) +
  labs(x = "Years relative to CAH designation",
       y = "Inpatient days per bed") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave("results/anticipation-ipdays.png", p_antic_ipd,
       width = 6.5, height = 4.25, dpi = 300, scale = 1.5)
