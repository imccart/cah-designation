# Displacement analysis and capacity accounting
# Self-contained: expects from _run-analysis.r: est.dat, hosp.results.table, state.results.table
# Outputs: CAH-exposure map, net capacity Monte Carlo figure, diagnostic CSVs



# 0. Haversine function (inlined from data-code/functions.R) -----------------
haversine_vec <- function(lon1, lat1, lon2, lat2, units = "miles") {
  R <- if (units == "miles") 3958.8 else 6371.01
  lat1_rad <- lat1 * pi / 180
  lat2_rad <- lat2 * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  dlat <- (lat2 - lat1) * pi / 180
  a <- sin(dlat/2)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)^2
  c_val <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R * c_val
}


# 1. Hospital locations and 2005 snapshot ------------------------------------
modal_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  as.numeric(names(sort(table(x), decreasing = TRUE))[1])
}

hosp_locs <- est.dat %>%
  group_by(ID) %>%
  summarize(
    lat  = modal_value(LAT),
    lon  = modal_value(LONG),
    .groups = "drop"
  ) %>%
  filter(!is.na(lat), !is.na(lon))

hosp_2005 <- est.dat %>%
  filter(year == 2005, closed == 0, merged == 0) %>%
  select(ID, BDTOT, ADMTOT, IPDTOT, cah, eff_year, ever_rural, MSTATE) %>%
  left_join(hosp_locs, by = "ID") %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  mutate(is_cah_2005 = ifelse(!is.na(eff_year) & eff_year <= 2005, 1, 0))


# 2. Next-nearest hospital distance -----------------------------------------

hosp_coords <- hosp_2005 %>% select(ID, lat, lon)
n_hosp <- nrow(hosp_coords)

dist_next <- numeric(n_hosp)

for (i in seq_len(n_hosp)) {
  d <- haversine_vec(hosp_coords$lon[i], hosp_coords$lat[i],
                     hosp_coords$lon, hosp_coords$lat)
  d[i] <- Inf
  dist_next[i] <- sort(d)[1]
}

hosp_2005 <- hosp_2005 %>%
  mutate(dist_next = dist_next)

## Distance to the next-nearest hospital, overall and among converters. Section
## 5.2 quotes the converter median and IQR.
hosp_2005 %>%
  summarise(median_all = median(dist_next),
            median_cah = median(dist_next[is_cah_2005 == 1]),
            p25_cah    = quantile(dist_next[is_cah_2005 == 1], 0.25),
            p75_cah    = quantile(dist_next[is_cah_2005 == 1], 0.75),
            n_hospitals = n(),
            n_cah = sum(is_cah_2005 == 1)) %>%
  write_csv("results/diagnostics/next-nearest-distance.csv")


# 3. Pull SDID estimates and compute key quantities --------------------------

delta_B <- hosp.results.table %>%
  filter(outcome == "Total beds") %>%
  pull(sdid_att)

delta_B_se <- hosp.results.table %>%
  filter(outcome == "Total beds") %>%
  mutate(se = (sdid_ci_high - sdid_ci_low) / (2 * 1.96)) %>%
  pull(se)

delta_C <- state.results.table %>%
  filter(outcome == "Closures") %>%
  pull(sdid_att)

delta_C_se <- state.results.table %>%
  filter(outcome == "Closures") %>%
  mutate(se = (sdid_ci_high - sdid_ci_low) / (2 * 1.96)) %>%
  pull(se)

delta_IPD <- hosp.results.table %>%
  filter(outcome == "Inpatient days per bed") %>%
  pull(sdid_att)

delta_IPD_se <- hosp.results.table %>%
  filter(outcome == "Inpatient days per bed") %>%
  mutate(se = (sdid_ci_high - sdid_ci_low) / (2 * 1.96)) %>%
  pull(se)

# Observed rho
eligible_2005 <- est.dat %>%
  filter(year == 2005, closed == 0, merged == 0,
         BDTOT <= 50, state_treat_year > 0) %>%
  summarize(n_cah = sum(is.finite(eff_year) & eff_year <= 2005),
            n_eligible = n())
rho_obs <- eligible_2005$n_cah / eligible_2005$n_eligible

# Mean beds among pre-CAH hospitals (no bed cut — these are eventual CAHs)
pre_cah <- est.dat %>%
  filter(!is.na(eff_year), year < eff_year, year >= 1995, year <= 2010)

B_close    <- mean(pre_cah$BDTOT, na.rm = TRUE)
B_close_se <- sd(pre_cah$BDTOT, na.rm = TRUE) / sqrt(sum(!is.na(pre_cah$BDTOT)))

# Mean pre-CAH inpatient days (for IPD net capacity calculation)
IPD_close    <- mean(pre_cah$IPDTOT, na.rm = TRUE)
IPD_close_se <- sd(pre_cah$IPDTOT, na.rm = TRUE) / sqrt(sum(!is.na(pre_cah$IPDTOT)))

# Total hospitals in CAH states
n_hosp_cah_states <- est.dat %>%
  filter(year == 2005, state_treat_year > 0) %>%
  nrow()

# Closures prevented and average admissions
n_converters <- sum(hosp_2005$is_cah_2005)
closures_prevented <- abs(delta_C) / 100 * n_hosp_cah_states
avg_adm_cah <- hosp_2005 %>%
  filter(is_cah_2005 == 1) %>%
  summarize(mean_adm = mean(ADMTOT, na.rm = TRUE)) %>%
  pull(mean_adm)

# Closure mortality calibration (Gujral & Basu 2019: 0.78pp)
lives_saved_closure <- closures_prevented * avg_adm_cah * 0.0078

rho_star_beds <- (abs(delta_C) / 100 * B_close) / abs(delta_B)
rho_star_ipd  <- (abs(delta_C) / 100 * IPD_close) / (abs(delta_IPD) * B_close)

## The capacity accounting in Section 5.1 and the lives-saved figure in Section
## 5.2 are quoted from these quantities, so they go to a file rather than the
## console. rho_star_beds is the paper's break-even share.
tibble(delta_B = delta_B, delta_C = delta_C, delta_IPD = delta_IPD,
       rho_obs = rho_obs,
       B_close = B_close, B_close_se = B_close_se,
       IPD_close = IPD_close, IPD_close_se = IPD_close_se,
       rho_star_beds = rho_star_beds, rho_star_ipd = rho_star_ipd,
       closures_prevented = closures_prevented,
       avg_adm_cah = avg_adm_cah,
       lives_saved_closure = lives_saved_closure) %>%
  write_csv("results/diagnostics/capacity-accounting.csv")


# 4. Monte Carlo: net capacity accounting ------------------------------------
# The marginal hospital's bed size is drawn from the observed bed sizes of
# small hospitals in the year before they closed (R&R, R1 comment 1), rather
# than from a normal centered on the mean converter.
closer_beds <- est.dat %>%
  arrange(ID, year) %>%
  group_by(ID) %>%
  mutate(beds_lag = lag(BDTOT), year_lag = lag(year)) %>%
  ungroup() %>%
  ## year_lag == year - 1 keeps only bed counts observed the year before closure
  filter(closed == 1, year >= 1995, year <= 2010, !is.na(beds_lag),
         year_lag == year - 1, beds_lag <= 50) %>%
  pull(beds_lag)
## Section 5.1 quotes the count and mean of this distribution as the basis for
## the closer-based break-even alternative.
tibble(n = length(closer_beds), mean = mean(closer_beds), median = median(closer_beds),
       p25 = quantile(closer_beds, 0.25), p10 = quantile(closer_beds, 0.10)) %>%
  write_csv("results/diagnostics/closer-bedsize.csv")

set.seed(42)
n_draws <- 10000

mc_draws <- tibble(
  delta_B_draw   = rnorm(n_draws, delta_B, delta_B_se),
  delta_C_draw   = rnorm(n_draws, delta_C, delta_C_se),
  delta_IPD_draw = rnorm(n_draws, delta_IPD, delta_IPD_se),
  rho_draw       = runif(n_draws, 0.15, 0.45),
  B_close_draw   = sample(closer_beds, n_draws, replace = TRUE),
  IPD_close_draw = rnorm(n_draws, IPD_close, IPD_close_se)
) %>%
  mutate(
    net_beds = (-delta_C_draw / 100) * B_close_draw + rho_draw * delta_B_draw,
    rho_star = ifelse(delta_B_draw != 0,
                      ((-delta_C_draw / 100) * B_close_draw) / (-delta_B_draw),
                      NA),
    net_ipd  = (-delta_C_draw / 100) * IPD_close_draw +
               rho_draw * (delta_IPD_draw * B_close_draw),
    rho_star_ipd = ifelse(delta_IPD_draw * B_close_draw != 0,
                          ((-delta_C_draw / 100) * IPD_close_draw) /
                          (-(delta_IPD_draw * B_close_draw)),
                          NA)
  )

mc_summary <- tibble(
  metric = c("Net beds/hospital", "Break-even rho* (beds)",
             "Net inpatient days/hospital", "Break-even rho* (IPD)"),
  median = c(median(mc_draws$net_beds, na.rm = TRUE),
             median(mc_draws$rho_star, na.rm = TRUE),
             median(mc_draws$net_ipd, na.rm = TRUE),
             median(mc_draws$rho_star_ipd, na.rm = TRUE)),
  ci_05  = c(quantile(mc_draws$net_beds, 0.05, na.rm = TRUE),
             quantile(mc_draws$rho_star, 0.05, na.rm = TRUE),
             quantile(mc_draws$net_ipd, 0.05, na.rm = TRUE),
             quantile(mc_draws$rho_star_ipd, 0.05, na.rm = TRUE)),
  ci_95  = c(quantile(mc_draws$net_beds, 0.95, na.rm = TRUE),
             quantile(mc_draws$rho_star, 0.95, na.rm = TRUE),
             quantile(mc_draws$net_ipd, 0.95, na.rm = TRUE),
             quantile(mc_draws$rho_star_ipd, 0.95, na.rm = TRUE)),
  prob_positive = c(mean(mc_draws$net_beds > 0, na.rm = TRUE), NA,
                    mean(mc_draws$net_ipd > 0, na.rm = TRUE), NA)
)

## mc_summary is written to results/diagnostics/psa_summary.csv below.


# 5. Outputs ----------------------------------------------------------------

## 5a. County map: CAH dependence
county_pop <- read_csv("data/input/CenPop2010_Mean_CO.txt",
                       col_types = cols(STATEFP = "c", COUNTYFP = "c",
                                        POPULATION = "d", LATITUDE = "d",
                                        LONGITUDE = "d", .default = "c")) %>%
  mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
  filter(!STATEFP %in% c("72", "78", "66", "60", "69"))


county_nearest <- county_pop %>%
  select(fips, county_lat = LATITUDE, county_lon = LONGITUDE) %>%
  rowwise() %>%
  mutate(
    dists = list({
      d <- haversine_vec(county_lon, county_lat, hosp_2005$lon, hosp_2005$lat)
      idx_sorted <- order(d)
      tibble(
        nearest_id = hosp_2005$ID[idx_sorted[1]],
        nearest_dist = d[idx_sorted[1]],
        nearest_is_cah = hosp_2005$is_cah_2005[idx_sorted[1]]
      )
    })
  ) %>%
  unnest(dists) %>%
  ungroup() %>%
  select(-county_lat, -county_lon)

non_cah_hosp <- hosp_2005 %>% filter(is_cah_2005 == 0)

county_cah_dep <- county_nearest %>%
  filter(nearest_is_cah == 1) %>%
  left_join(county_pop %>% select(fips, county_lat = LATITUDE, county_lon = LONGITUDE),
            by = "fips") %>%
  rowwise() %>%
  mutate(
    dist_to_non_cah = {
      d <- haversine_vec(county_lon, county_lat, non_cah_hosp$lon, non_cah_hosp$lat)
      min(d)
    }
  ) %>%
  ungroup() %>%
  select(fips, nearest_dist, dist_to_non_cah) %>%
  mutate(extra_travel = dist_to_non_cah - nearest_dist)

county_map_dat <- county_nearest %>%
  left_join(county_cah_dep %>% select(fips, dist_to_non_cah, extra_travel), by = "fips") %>%
  mutate(cah_dependent = nearest_is_cah == 1)

counties_sf <- st_read(
  "data/input/county-shapefiles/cb_2017_us_county_5m.shp",
  quiet = TRUE
) %>%
  mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
  filter(!STATEFP %in% c("72", "78", "66", "60", "69", "02", "15"))

map_dat <- counties_sf %>%
  left_join(county_map_dat %>% select(fips, cah_dependent, extra_travel), by = "fips") %>%
  mutate(
    fill_val = case_when(
      !cah_dependent | is.na(cah_dependent) ~ NA_real_,
      TRUE ~ pmin(extra_travel, 60)
    )
  )

p_map <- ggplot(map_dat) +
  geom_sf(aes(fill = fill_val), color = "gray80", linewidth = 0.05) +
  scale_fill_gradient(
    low = "gray85",
    high = "gray15",
    na.value = "gray95",
    name = "Extra travel\n(miles)",
    breaks = c(0, 15, 30, 45, 60),
    labels = c("0", "15", "30", "45", "\u226560"),
    limits = c(0, 60)
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = c(0.92, 0.3),
    legend.key.height = unit(0.8, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave("results/map-cah-exposure.png", p_map, width = 10, height = 6, dpi = 300)

## 5b. Monte Carlo figure: net beds per hospital
p_mc_beds <- ggplot(mc_draws, aes(x = net_beds)) +
  geom_histogram(bins = 60, fill = "gray50", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = median(mc_draws$net_beds), linetype = "solid", color = "black") +
  annotate("text", x = median(mc_draws$net_beds),
           y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("Median = %.2f", median(mc_draws$net_beds)),
           size = 3.2) +
  annotate("text", x = max(mc_draws$net_beds) * 0.7,
           y = Inf, vjust = 3,
           label = sprintf("P(net > 0) = %.1f%%", 100 * mean(mc_draws$net_beds > 0)),
           size = 3.2) +
  labs(x = "Net beds per hospital", y = "Count") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("results/psa-net-beds.png", p_mc_beds, width = 5.5, height = 4.5, dpi = 300)

## 5b2. Monte Carlo figure: net inpatient days per hospital
p_mc_ipd <- ggplot(mc_draws, aes(x = net_ipd)) +
  geom_histogram(bins = 60, fill = "gray50", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = median(mc_draws$net_ipd), linetype = "solid", color = "black") +
  annotate("text", x = median(mc_draws$net_ipd),
           y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("Median = %.0f", median(mc_draws$net_ipd)),
           size = 3.2) +
  annotate("text", x = max(mc_draws$net_ipd) * 0.7,
           y = Inf, vjust = 3,
           label = sprintf("P(net > 0) = %.1f%%", 100 * mean(mc_draws$net_ipd > 0)),
           size = 3.2) +
  labs(x = "Net inpatient days per hospital", y = "Count") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("results/psa-net-ipd.png", p_mc_ipd, width = 4, height = 3.5, dpi = 300)

## 5c. Diagnostic CSVs
write_csv(mc_summary, "results/diagnostics/psa_summary.csv")

write_csv(
  hosp_2005 %>%
    filter(is_cah_2005 == 1) %>%
    select(ID, BDTOT, ADMTOT, IPDTOT, dist_next),
  "results/diagnostics/cah_converter_distances.csv"
)

