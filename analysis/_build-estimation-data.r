# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   2/24/2026
## Description:   Build estimation datasets (est.dat, state.dat)


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, haven, readxl, janitor, here, zoo, fedmatch, zipcodeR,
               fixest, did, did2s, BMisc, fect, glmnet, nnet, mlogit, survival,
               scales, plotly, panelView, dotwhisker, patchwork, sf,
               modelsummary, kableExtra, broom, synthdid)

source('analysis/functions.R')

## Winsorization function
winsor_by_year <- function(x, probs = c(0.05, 0.95)) {
  qs <- quantile(x, probs = probs, na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}

# Read-in data ------------------------------------------------------------
# The HCRIS fields in aha_final.csv are joined on the AHA provider number, which
# misses converters' pre-designation cost reports. hcris_financial.csv carries
# the same fields on the year-varying crosswalk, so it replaces them here and is
# the only source of the six financial outcomes.
hcris.financial <- read_csv('data/output/hcris_financial.csv', show_col_types = FALSE,
                            col_types = cols(ID = col_character()))
aha.data <- read_csv('data/output/aha_final.csv', col_types = cols(ID = col_character())) %>%
  select(-any_of(c('net_pat_rev', 'tot_operating_exp', 'fixed_assets', 'accum_dep',
                   'current_assets', 'current_liabilities'))) %>%
  left_join(hcris.financial, by = c('ID', 'year'))
aha.neighbors <- read_csv('data/output/aha_neighbors.csv')
cah.dates <- read_csv('data/input/cah-states.csv')

state.xwalk <- tibble(
  state = state.name,
  MSTATE = state.abb
)

copa.dat <- read_csv('data/input/copa-states.csv', col_names=c("state", "empty", "copa_status", "lat_lon", "details")) %>%
  separate(details, sep=",", into=c("state2","statute","status","other")) %>%
  mutate(
    copa_start = str_extract(status, "(?<=Enacted )\\d{4}") %>% as.integer(),
    copa_end = str_extract(status, "(?<=Repealed )\\d{4}"),
    copa_end = replace_na(copa_end, "2025") %>% as.integer()
  ) %>%
  mutate(
    copa_start = case_when(
      state == "Maine" ~ 2005,
      state == "New York" ~ 2014,
      TRUE ~ copa_start
    )
  ) %>%
  select(state, copa_status, copa_start, copa_end) %>%
  filter(!is.na(state)) %>%
  left_join(state.xwalk, by="state")

write_csv(copa.dat,'data/output/copa_data.csv')

cpi.data <- read_xlsx("data/input/CPI_1913_2019.xlsx", skip = 11)
cpi.data <- pivot_longer(cpi.data,
                         cols=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                         names_to="month",
                         values_to="index") %>%
  mutate(year=Year) %>%
  group_by(year) %>%
  summarize(index=mean(index, na.rm=TRUE))

cpi.2010 <- cpi.data %>%
  filter(year==2010) %>%
  select(index) %>% as.numeric()

cpi.final <- cpi.data %>%
  mutate(cpi_2010 = cpi.2010,
         cpi_deflator = index/cpi_2010) %>%
  select(year, cpi_deflator)


# Merge and finalize hospital-level data -----------------------------------

data.merge <- aha.data %>%
    left_join(aha.neighbors, by=c("ID", "year")) %>%
    left_join(cah.dates %>% select(cah_date_law=cah_date, cah_year_law=cah_year, MSTATE="abb"), by="MSTATE") %>%
    filter(! MSTATE %in% c("AK","HI","PR","VI","GU","MP","AS", "N","0", "DC", "DE", "MH", "ML","MD"),
           !is.na(MSTATE), MSTATE!="NA") %>%
    filter(COMMTY=="Y", hosp_type=="General") %>%
    group_by(ID, year) %>%
    mutate(hosp_count=n()) %>%
    filter(hosp_count==1) %>% ungroup () %>%
    select(-hosp_count)

## construct new variables
final.dat <- data.merge %>%
  mutate(closed=case_when(
                !is.na(change_type) & change_type=="Closure" ~ 1,
                !is.na(change_type) & change_type!="Closure" ~ 0,
                is.na(change_type) ~ 0),
         merged=case_when(
                !is.na(change_type) & change_type=="Merger" ~ 1,
                !is.na(change_type) & change_type!="Merger" ~ 0,
                is.na(change_type) ~ 0),
         closed_merged=ifelse(closed==1 | merged==1, 1, 0),
         ) %>%
    group_by(ID) %>% mutate(ever_cah=max(cah, na.rm=TRUE)) %>% ungroup() %>%
    group_by(MSTATE) %>% mutate(state_first_obs=min(eff_year, na.rm=TRUE),
                                state_first_law=min(cah_year_law, na.rm=TRUE)) %>% ungroup() %>%
    mutate(state_first_obs=ifelse(is.infinite(state_first_obs), 0, state_first_obs),
           state_first_law=ifelse(is.infinite(state_first_law), 0, state_first_law),
           state_treat_year=state_first_obs)

## mean hospital beds pre-treatment
hosp.beds <- final.dat %>%
  filter(year < eff_year | is.na(eff_year), year>1995) %>%
  group_by(ID) %>%
  summarize(beds_base=mean(BDTOT, na.rm=TRUE))

est.dat <- final.dat %>%
  left_join(hosp.beds, by="ID") %>%
  left_join(cpi.final , by="year") %>%
  mutate(margin_990=case_when(
            !is.na(margin_1) ~ margin_1,
            is.na(margin_1) & !is.na(margin_2) ~ margin_2,
            is.na(margin_1) & is.na(margin_2) & !is.na(margin_3) ~ margin_3),
         net_pat_rev_990=case_when(
            !is.na(total_revenue_1) ~ total_revenue_1,
            is.na(total_revenue_1) & !is.na(total_revenue_2) ~ total_revenue_2,
            is.na(total_revenue_1) & is.na(total_revenue_2) & !is.na(total_revenue_3) ~ total_revenue_3),
         tot_operating_exp_990=case_when(
            !is.na(total_expenses_1) ~ total_expenses_1,
            is.na(total_expenses_1) & !is.na(total_expenses_2) ~ total_expenses_2,
            is.na(total_expenses_1) & is.na(total_expenses_2) & !is.na(total_expenses_3) ~ total_expenses_3),
         net_fixed_990=case_when(
            !is.na(net_fixed_1) ~ net_fixed_1,
            is.na(net_fixed_1) & !is.na(net_fixed_2) ~ net_fixed_2,
            is.na(net_fixed_1) & is.na(net_fixed_2) & !is.na(net_fixed_3) ~ net_fixed_3),
         depreciation_990=case_when(
            !is.na(depreciation_1) ~ depreciation_1,
            is.na(depreciation_1) & !is.na(depreciation_2) ~ depreciation_2,
            is.na(depreciation_1) & is.na(depreciation_2) & !is.na(depreciation_3) ~ depreciation_3),
         current_ratio_990=case_when(
            !is.na(current_ratio_1) ~ current_ratio_1,
            is.na(current_ratio_1) & !is.na(current_ratio_2) ~ current_ratio_2,
            is.na(current_ratio_1) & is.na(current_ratio_2) & !is.na(current_ratio_3) ~ current_ratio_3)
          ) %>%
  mutate(
    ip_per_bed=IPDTOT/beds_base,
    gross_fixed = fixed_assets + accum_dep,
    margin = ifelse(net_pat_rev > 0, (net_pat_rev - tot_operating_exp) / net_pat_rev, NA_real_),
    current_ratio = ifelse(current_liabilities > 0, current_assets / current_liabilities, NA_real_),
    net_fixed = fixed_assets,
    state_event_time=case_when(
      state_treat_year>0 ~ year - state_treat_year,
      state_treat_year==0 ~ -1),
    hosp_event_time=case_when(
      !is.na(eff_year) ~ year - eff_year,
      TRUE ~ -1),
    aha_id=as.numeric(ID),
    treat_state=ifelse(state_treat_year>0, 1, 0),
    treat_state_post=ifelse(year>=state_treat_year & state_treat_year>0, 1, 0)
  ) %>%
  group_by(ID) %>%
  arrange(year, .by_group=TRUE) %>%
  mutate(capex_990=net_fixed_990 - lag(net_fixed_990) + depreciation_990,
         capex=gross_fixed - lag(gross_fixed)) %>%
  ungroup() %>%
  group_by(year, ever_cah) %>%
  mutate(
    margin        = winsor_by_year(margin),
    net_fixed     = winsor_by_year(net_fixed),
    current_ratio = winsor_by_year(current_ratio),
    capex         = winsor_by_year(capex),
    net_pat_rev       = winsor_by_year(net_pat_rev),
    tot_operating_exp = winsor_by_year(tot_operating_exp),
    ## Form-990-only versions, processed identically. Used only by the R2.5
    ## source comparison in R1_revision.R.
    margin_9        = winsor_by_year(ifelse(margin_990 > -1 & margin_990 < 1, margin_990, NA_real_)),
    net_fixed_9     = winsor_by_year(net_fixed_990),
    current_ratio_9 = winsor_by_year(current_ratio_990),
    capex_9         = winsor_by_year(capex_990),
    net_pat_rev_9       = winsor_by_year(net_pat_rev_990),
    tot_operating_exp_9 = winsor_by_year(tot_operating_exp_990),
    BDTOT         = winsor_by_year(BDTOT),
    OBBD          = winsor_by_year(OBBD),
    FTERN         = winsor_by_year(FTERN),
    IPDTOT        = winsor_by_year(IPDTOT),
    ip_per_bed    = if_else(ip_per_bed > (365*.85), 365*.85, ip_per_bed)
  ) %>%
  ungroup() %>%
  mutate(
    net_fixed    = net_fixed / 1e6 / beds_base / cpi_deflator,
    capex        = capex / 1e4 / beds_base / cpi_deflator,
    net_pat_rev       = net_pat_rev / 1e3 / beds_base / cpi_deflator,
    tot_operating_exp = tot_operating_exp / 1e3 / beds_base / cpi_deflator,
    net_fixed_9    = net_fixed_9 / 1e6 / beds_base / cpi_deflator,
    capex_9        = capex_9 / 1e4 / beds_base / cpi_deflator,
    net_pat_rev_9       = net_pat_rev_9 / 1e3 / beds_base / cpi_deflator,
    tot_operating_exp_9 = tot_operating_exp_9 / 1e3 / beds_base / cpi_deflator
  ) %>%
  arrange(ID, year) %>%
  group_by(ID) %>%
  fill(state, .direction = "downup") %>%
  mutate(margin=na.approx(margin, x=year, na.rm=FALSE),
         net_fixed=na.approx(net_fixed, x=year, na.rm=FALSE),
         current_ratio=na.approx(current_ratio, x=year, na.rm=FALSE),
         capex=na.approx(capex, x=year, na.rm=FALSE),
         BDTOT = na.approx(BDTOT, x=year, na.rm=FALSE),
         OBBD = na.approx(OBBD, x=year, na.rm=FALSE),
         FTERN = na.approx(FTERN, x=year, na.rm=FALSE),
         IPDTOT = na.approx(IPDTOT, x=year, na.rm=FALSE),
         ip_per_bed = na.approx(ip_per_bed, x=year, na.rm=FALSE),
         net_pat_rev = na.approx(net_pat_rev, x=year, na.rm=FALSE),
         tot_operating_exp = na.approx(tot_operating_exp, x=year, na.rm=FALSE),
         across(c(margin_9, net_fixed_9, current_ratio_9, capex_9,
                  net_pat_rev_9, tot_operating_exp_9),
                ~ na.approx(.x, x = year, na.rm = FALSE))) %>%
  ungroup()


## ever_rural indicator (hospital-level)
est.dat <- est.dat %>%
  group_by(ID) %>%
  mutate(ever_rural = as.integer(any(CBSATYPE == "Rural", na.rm = TRUE))) %>%
  ungroup()

## check missing values by year
missing.dat <- est.dat %>%
  group_by(year) %>%
  summarize(
    n=n(),
    missing_margin=sum(is.na(margin)),
    missing_net_fixed=sum(is.na(net_fixed)),
    missing_current_ratio=sum(is.na(current_ratio)),
    missing_capex=sum(is.na(capex))
  )

# Aggregate to state level --------------------------------------------------

state.dat <- est.dat %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=min(state_treat_year),
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE),
    mergers=sum(merged), changes=sum(closed_merged),
    state_treat_year=first(state_treat_year),
    mean_beds=mean(BDTOT, na.rm=TRUE), mean_distance=mean(distance, na.rm=TRUE)) %>%
  group_by(MSTATE) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(hospitals_lag = lag(hospitals)) %>%
  ungroup() %>%
  mutate(
    event_time=case_when(
      state_treat_year>0 ~ year - state_treat_year,
      state_treat_year==0 ~ -1),
    treat=ifelse(state_treat_year>0, 1, 0),
    treat_post=ifelse(year>=state_treat_year & state_treat_year>0, 1, 0),
    state=as.numeric(factor(MSTATE)),
    rate_closed=closures/hospitals_lag,
    rate_merged=mergers/hospitals_lag,
    rate_changes=changes/hospitals_lag)

# Write output --------------------------------------------------------------
write_csv(est.dat, 'data/output/estimation_data.csv')
write_csv(state.dat, 'data/output/state_estimation_data.csv')

# Summary statistics --------------------------------------------------------
source('analysis/1-sum-stats.R')
