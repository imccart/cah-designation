# Meta --------------------------------------------------------------------

## Author:        Ian McCarthy
## Date Created:  5/17/2023
## Date Edited:   1/29/2026
## Description:   Build Analytic Data
##
## Note: Run fuzzy matching scripts first:
##   - data-code/fuzzy990.R
##   - data-code/fuzzyhcris.R
##   - data-code/fuzzycah.R


# Preliminaries -----------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyverse, lubridate, stringr, modelsummary, broom, janitor, here,
               scales, zipcodeR, purrr, zoo)

source('data-code/functions.R')
source('data-code/api-keys.R')


# Read-in data ------------------------------------------------------------
aha.combine <- read_csv('data/input/aha_data.csv', show_col_types = FALSE)

col_names <- names(read_csv('data/input/zcta-to-county.csv', n_max = 0, show_col_types = FALSE))
state.zip.xwalk <- read_csv('data/input/zcta-to-county.csv', col_names = col_names, skip = 2, show_col_types = FALSE) %>%
  group_by(zcta5, stab) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(zcta5) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  select(zcta5, state_xwalk = stab)

hcris.data <- read_tsv('data/input/hcris_data.txt', show_col_types = FALSE) %>%
  rename(MCRNUM = provider_number)

form990.data <- read_tsv('data/input/form990_ahaid.txt', show_col_types = FALSE) %>%
  mutate(total_revenue=abs(total_revenue),
         total_expenses=abs(total_expenses),
         total_assets=abs(total_assets),
         current_assets=abs(current_assets),
         current_liabilities=abs(current_liabilities),
         total_liabilities=abs(total_liabilities),
         depreciation=abs(depreciation),
         margin=if_else(
              !is.na(total_revenue) & total_revenue>0 & !is.na(total_expenses) & total_expenses>0,
              (total_revenue-total_expenses)/total_revenue,
              NA),
         current_ratio=if_else(
              !is.na(current_assets) & current_assets>0 & !is.na(current_liabilities) & current_liabilities>0,
              current_assets/current_liabilities,
              NA)) %>%
  filter(margin>-1, margin<1)

# Read crosswalk files from fuzzy matching scripts ------------------------
unique.990 <- read_csv('data/output/unique_990.csv', show_col_types = FALSE) %>%
  mutate(ID = as.character(ID),
         across(starts_with("ein_"), as.character))

# Ensure form990.data$ein is character for consistent joins
form990.data <- form990.data %>%
  mutate(ein = as.character(ein))

aha.crosswalk <- read_csv('data/output/unique_hcris.csv', show_col_types = FALSE) %>%
  mutate(ID = as.character(ID)) %>%
  rename(MCRNUM_xw = MCRNUM)

fuzzy.unique.cah <- read_csv('data/output/unique_cah.csv', show_col_types = FALSE) %>%
  mutate(ID = as.character(ID))


# AHA data cleaning (manual) -----------------------------------------------

## Convert ID to character for consistent joins with crosswalks
aha.combine <- aha.combine %>%
  mutate(ID = as.character(ID))

## IDs with incorrect state codes in one year
aha.combine <- aha.combine %>%
  mutate(MSTATE = if_else(ID == "6540810", "MS", MSTATE),
         MSTATE = if_else(ID == "6710190", "AR", MSTATE),
         MSTATE = if_else(ID == "6931186", "CA", MSTATE))


# Final data --------------------------------------------------------------

aha.cah.dates <- aha.combine %>%
  filter(critical_access==1) %>%
  group_by(ID) %>%
  summarize(eff_year_aha=min(year, na.rm=TRUE))

aha.final <- aha.combine %>% 
  left_join(fuzzy.unique.cah, by='ID') %>%
  left_join(unique.990, by='ID') %>%
  left_join(form990.data %>% 
    select(ein_1=ein, year, margin_1=margin, current_ratio_1=current_ratio, 
          net_fixed_1=fixed_assets, depreciation_1=depreciation,
          total_revenue_1=total_revenue, total_expenses_1=total_expenses), by=c('ein_1','year')) %>%
  left_join(form990.data %>% 
    select(ein_2=ein, year, margin_2=margin, current_ratio_2=current_ratio, 
           net_fixed_2=fixed_assets, depreciation_2=depreciation,
           total_revenue_2=total_revenue, total_expenses_2=total_expenses), by=c('ein_2','year')) %>%
  left_join(form990.data %>% 
    select(ein_3=ein, year, margin_3=margin, current_ratio_3=current_ratio, 
           net_fixed_3=fixed_assets, depreciation_3=depreciation,
           total_revenue_3=total_revenue, total_expenses_3=total_expenses), by=c('ein_3','year')) %>%
  left_join(state.zip.xwalk, by=c("MLOCZIP"="zcta5")) %>%
  mutate(MSTATE=case_when(
      !is.na(MSTATE) ~ MSTATE,
      is.na(MSTATE) & !is.na(state_xwalk) ~ state_xwalk
    )) %>%
  left_join(aha.crosswalk, by=c("ID","year")) %>%
  mutate(MCRNUM=case_when(
      !is.na(MCRNUM) ~ MCRNUM,
      is.na(MCRNUM) & !is.na(MCRNUM_xw) ~ MCRNUM_xw
    )) %>%
  left_join(hcris.data, by=c('MCRNUM','year')) %>%
  left_join(aha.cah.dates, by='ID') %>%
  mutate(eff_year=case_when(
           !is.na(first_date) ~ year(first_date),
           is.na(first_date) & !is.na(eff_year_aha) ~ eff_year_aha,
           TRUE ~ NA_real_),
         cah = case_when(
           is.na(critical_access) & cah_sup==1 & year>=eff_year ~ 1,
           critical_access==0 & cah_sup==1 & year>=eff_year ~ 1,
           critical_access==1 ~ 1,
           TRUE ~ 0),
         hosp_type = case_when(
           SERV==10 ~ 'General',
           SERV==13 ~ 'Surgical',
           SERV==33 ~ 'Specialty - Respiratory',
           SERV==41 ~ 'Specialty - Cancer',
           SERV==42 ~ 'Specialty - Heart',
           SERV==44 ~ 'Specialty - OBGYN',
           SERV==45 ~ 'Specialty - ENT',
           SERV==47 ~ 'Specialty - Orthopedic',
           SERV==49 ~ 'Specialty - Other',
           SERV %in% c(50, 53, 55, 56, 57, 58, 59) ~ 'Specialty - Childrens',
           SERV==80 ~ 'Long-term Care',
           TRUE ~ 'Other'
         )) %>%
    write_csv('data/output/aha_final.csv')         


## hospital closures and mergers
merge.close <- aha.final %>% 
  group_by(year, change_type, cah) %>%
  summarize(hosp_count=n()) %>% ungroup() %>%
  group_by(year, cah) %>%
  mutate(hosp_type_count=sum(hosp_count)) %>% ungroup() %>%
  group_by(year, change_type) %>%
  mutate(change_type_count=sum(hosp_count)) %>% ungroup() %>%
  write_csv('data/output/merge_close.csv')


## identify each hospital's nearest neighbor and the distance to that neighbor using the haversine formula
aha.geo <- aha.final %>%
  mutate(zip = substr(MLOCZIP, 1, 5)) %>%
  select(ID, LAT, LONG, year, zip) %>%
  distinct(ID, LAT, LONG, year, zip) %>%
  mutate(across(c(LAT, LONG), as.numeric),
         LAT = na_if(LAT, 0),
         LONG = na_if(LONG, 0))

calc_nearest_neighbor <- function(yr, geo_data) {
  yr_data <- geo_data %>% filter(year == yr)

  # Self-join for pairs (both sides filtered to same year)
  aha.pairs <- yr_data %>%
    cross_join(yr_data %>% select(ID2 = ID, LAT2 = LAT, LONG2 = LONG, zip2 = zip)) %>%
    filter(ID != ID2)

  # Vectorized haversine (returns miles to match zip_distance)
  aha.latlongdist <- aha.pairs %>%
    filter(!is.na(LAT), !is.na(LAT2)) %>%
    mutate(dist_latlong = haversine_vec(LONG, LAT, LONG2, LAT2, units = "miles")) %>%
    group_by(ID) %>%
    summarize(min_dist_latlong = min(dist_latlong, na.rm = TRUE), .groups = "drop")

  # Zip-based distance
  zip.pairs <- aha.pairs %>%
    filter(!is.na(zip), !is.na(zip2)) %>%
    distinct(zip, zip2)

  if (nrow(zip.pairs) > 0) {
    zip.dist <- zip_distance(zip.pairs$zip, zip.pairs$zip2) %>%
      rename(zip = zipcode_a, zip2 = zipcode_b, dist_zip = distance) %>%
      filter(!is.na(dist_zip))

    aha.zipdist <- aha.pairs %>%
      inner_join(zip.dist, by = c("zip", "zip2")) %>%
      group_by(ID) %>%
      summarize(min_dist_zip = min(dist_zip, na.rm = TRUE), .groups = "drop")
  } else {
    aha.zipdist <- tibble(ID = character(), min_dist_zip = numeric())
  }

  # Combine results
  yr_data %>%
    select(ID) %>%
    distinct() %>%
    left_join(aha.zipdist, by = "ID") %>%
    left_join(aha.latlongdist, by = "ID") %>%
    mutate(
      year = yr,
      distance = coalesce(pmin(min_dist_latlong, min_dist_zip, na.rm = TRUE),
                          min_dist_latlong, min_dist_zip)
    )
}

# Process all years (collect results in list, bind once)
nearest.neighbor <- unique(aha.geo$year) %>%
  map(\(yr) calc_nearest_neighbor(yr, aha.geo), .progress = TRUE) %>%
  list_rbind() %>%
  write_csv('data/output/aha_neighbors.csv')


# AHA visit and birth volumes ---------------------------------------------------
# Outpatient, ED, and birth counts, which the aha-data keep-list drops, read
# from the same WRDS extract that repo reads.
source('data-code/build-aha-visits.R')


# Year-varying HCRIS crosswalk for CAH converters -------------------------------
# CAH designation assigns a new provider number; this links converters' earlier
# cost reports (needs aha_final.csv and unique_hcris.csv, written above).
source('data-code/hcris-yearvarying-xw.R')

# HCRIS financial fields on that crosswalk, keyed on ID and year ----------------
# The estimation build derives the six financial outcomes from this file.
source('data-code/build-hcris-financial.R')
