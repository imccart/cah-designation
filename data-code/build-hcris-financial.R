# build-hcris-financial.R
#
# HCRIS financial fields keyed on AHA ID and year, linked through the
# year-varying provider-number crosswalk. CAH designation assigns a new Medicare
# provider number, so a converter's pre-designation cost reports are filed under
# its earlier number; joining on the AHA number alone misses them.
#
# This is the only source of the six financial outcomes. The estimation build
# merges this file and derives margin, current ratio, net fixed assets, capital
# expenditures, revenue and expenses from these fields.
#
# Inputs:  data/input/hcris_data.txt, data/output/hcris_xw_yearvarying.csv
# Output:  data/output/hcris_financial.csv

xw <- read_csv('data/output/hcris_xw_yearvarying.csv',
               col_types = cols(ID = col_character(),
                                MCRNUM_yv = col_character(),
                                MCRNUM_pre_desig = col_character()))

## A hospital can file more than one cost report in a fiscal year, so collapse to
## one row per provider-year before joining. Summing with na.rm returns zero when
## every report is missing, so guard that case back to missing.
hcris <- read_tsv('data/input/hcris_data.txt', show_col_types = FALSE,
                  col_select = c(provider_number, year, net_pat_rev, tot_operating_exp,
                                 fixed_assets, accum_dep, current_assets, current_liabilities)) %>%
  mutate(MCRNUM = str_pad(str_remove_all(as.character(provider_number), "\\D"),
                          6, side = "left", pad = "0")) %>%
  select(-provider_number) %>%
  group_by(MCRNUM, year) %>%
  summarize(across(everything(),
                   ~ if_else(all(is.na(.x)), NA_real_, sum(.x, na.rm = TRUE))),
            .groups = "drop")

hcris_old <- hcris %>% rename_with(~ paste0(.x, "_old"), -c(MCRNUM, year))

## In the designation year a converter may file partial-year reports under both
## numbers. Revenue and expenses are flows and are summed across the two; the
## balance sheet items are taken from the report under the newer number.
linked <- xw %>%
  left_join(hcris, by = c("MCRNUM_yv" = "MCRNUM", "year")) %>%
  left_join(hcris_old, by = c("MCRNUM_pre_desig" = "MCRNUM", "year")) %>%
  mutate(
    net_pat_rev = if_else(is.na(net_pat_rev_old),
                          net_pat_rev,
                          coalesce(net_pat_rev, 0) + net_pat_rev_old),
    tot_operating_exp = if_else(is.na(tot_operating_exp_old),
                                tot_operating_exp,
                                coalesce(tot_operating_exp, 0) + tot_operating_exp_old),
    fixed_assets        = coalesce(fixed_assets, fixed_assets_old),
    accum_dep           = coalesce(accum_dep, accum_dep_old),
    current_assets      = coalesce(current_assets, current_assets_old),
    current_liabilities = coalesce(current_liabilities, current_liabilities_old)
  ) %>%
  select(-ends_with("_old"))

## match_strict describes the crosswalk match, which is a property of the
## hospital rather than the hospital-year.
strict_by_id <- linked %>%
  filter(!is.na(match_strict)) %>%
  distinct(ID, match_strict)

hcris.financial <- linked %>%
  select(-match_strict) %>%
  left_join(strict_by_id, by = "ID") %>%
  select(ID, year, net_pat_rev, tot_operating_exp, fixed_assets, accum_dep,
         current_assets, current_liabilities, match_strict) %>%
  ## aha_final.csv carries 41 duplicated ID-year rows from the AHA source, which
  ## reach here through the crosswalk. The duplicates are identical, so keep one
  ## and leave this file keyed on ID and year.
  distinct(ID, year, .keep_all = TRUE)


write_csv(hcris.financial, 'data/output/hcris_financial.csv')

## What the crosswalk buys: the share of converters' pre-designation
## hospital-years carrying a cost report, under the old time-invariant link and
## under this one. HCRIS begins in 1996, so earlier years are excluded. The
## appendix and the response letter both quote these shares.
xw.old <- read_csv('data/output/unique_hcris.csv', show_col_types = FALSE,
                   col_types = cols(ID = col_character())) %>%
  transmute(ID, year,
            MCRNUM_old = str_pad(str_remove_all(as.character(MCRNUM), "\\D"), 6, "left", "0"))

hcris.reports <- hcris %>%
  filter(!is.na(net_pat_rev), !is.na(tot_operating_exp)) %>%
  distinct(MCRNUM, year) %>%
  mutate(old_report = TRUE)

read_csv('data/output/aha_final.csv', show_col_types = FALSE,
         col_types = cols(ID = col_character()),
         col_select = c(ID, year, eff_year)) %>%
  filter(!is.na(eff_year), year < eff_year, year >= 1996) %>%
  select(ID, year) %>%
  left_join(hcris.financial %>%
              transmute(ID, year,
                        new_report = !is.na(net_pat_rev) & !is.na(tot_operating_exp)),
            by = c("ID", "year")) %>%
  left_join(xw.old, by = c("ID", "year")) %>%
  left_join(hcris.reports, by = c("MCRNUM_old" = "MCRNUM", "year")) %>%
  summarise(hospital_years = n(),
            hospitals = n_distinct(ID),
            share_old_link = mean(coalesce(old_report, FALSE)),
            share_new_link = mean(coalesce(new_report, FALSE))) %>%
  write_csv('data/output/hcris_crosswalk_coverage.csv')
