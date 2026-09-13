# hcris-yearvarying-xw.R
#
# Year-varying HCRIS-to-AHA crosswalk for CAH converters.
#
# CAH designation assigns a new Medicare provider number (xx13xx), and the AHA
# survey backfills that number to every year of a hospital's panel. The
# time-invariant crosswalk in unique_hcris.csv therefore links converters only
# to the cost reports filed under the CAH number, and their pre-designation
# reports (filed under the earlier number) are never found. This script matches
# each converter's pre-designation years to HCRIS providers with a non-CAH
# number, using the same name preprocessing and word-Jaccard scoring as
# fuzzyhcris.R, and writes a year-varying crosswalk that uses the earlier
# number before designation and the existing crosswalk number from designation
# onward. Does not modify est.dat or unique_hcris.csv.
#
# Inputs:  data/input/hcris_data.txt, data/output/aha_final.csv,
#          data/output/unique_hcris.csv
# Outputs: data/output/hcris_xw_yearvarying.csv          (ID, year, MCRNUM_yv, source)
#          data/output/fuzzy/qa_hcris_yearvarying.csv     (per-converter match diagnostics)

source('data-code/functions.R')

# HCRIS providers with a non-CAH number, one row per provider with name variants ----
hcris.data <- read_tsv('data/input/hcris_data.txt', show_col_types = FALSE,
                       col_select = c(provider_number, year, name, city, state, zip, net_pat_rev, tot_operating_exp)) %>%
  rename(MCRNUM = provider_number) %>%
  mutate(MCRNUM = str_pad(str_remove_all(as.character(MCRNUM), "\\D"), 6, side = "left", pad = "0"),
         pn4 = as.numeric(str_sub(MCRNUM, -4)),
         has_fin = !is.na(net_pat_rev) & !is.na(tot_operating_exp))

hcris.noncah <- hcris.data %>%
  filter(!(pn4 >= 1300 & pn4 <= 1399), has_fin) %>%
  mutate(zip = str_pad(substr(as.character(zip), 1, 5), width = 5, side = "left", pad = "0"),
         zip = if_else(str_detect(zip, "^[0-9]{5}$") & zip != "00000", zip, NA_character_),
         name_clean = preprocess_hospital_name(name),
         state = str_to_lower(state), city = str_to_lower(city)) %>%
  filter(!is.na(name_clean), name_clean != "")

hcris.years <- hcris.noncah %>% group_by(MCRNUM) %>% summarize(y_min = min(year), y_max = max(year), .groups = "drop")
hcris.variants <- hcris.noncah %>% distinct(MCRNUM, name_clean, state, city, zip)

# Converters: pre-designation name/location variants ---------------------------------
aha <- read_csv('data/output/aha_final.csv', show_col_types = FALSE,
                col_select = c(ID, year, MNAME, MSTATE, MLOCCITY, MLOCZIP, eff_year)) %>%
  mutate(ID = as.character(ID))

conv.pre <- aha %>%
  filter(!is.na(eff_year), year < eff_year) %>%
  mutate(zip = str_pad(substr(as.character(MLOCZIP), 1, 5), width = 5, side = "left", pad = "0"),
         zip = if_else(str_detect(zip, "^[0-9]{5}$") & zip != "00000", zip, NA_character_),
         state = str_to_lower(MSTATE), city = str_to_lower(MLOCCITY),
         name_clean = preprocess_hospital_name(MNAME)) %>%
  filter(!is.na(name_clean), name_clean != "", !is.na(state))

conv.variants <- conv.pre %>% distinct(ID, eff_year, name_clean, state, city, zip)

# Candidate pairs: state + zip blocking, then state + city blocking ------------------
cand.zip <- match_jaccard_words(
  data1 = conv.variants %>% filter(!is.na(zip)) %>% select(ID, name_clean, state, zip),
  data2 = hcris.variants %>% filter(!is.na(zip)) %>% select(MCRNUM, name_clean, state, zip),
  name_col1 = "name_clean", name_col2 = "name_clean", id_col1 = "ID", id_col2 = "MCRNUM",
  block_vars = c("state", "zip"), threshold = 0.50) %>%
  mutate(strategy = "jaccard_statezip", zip_exact = 1)

cand.city <- match_jaccard_words(
  data1 = conv.variants %>% filter(!is.na(city)) %>% select(ID, name_clean, state, city),
  data2 = hcris.variants %>% filter(!is.na(city)) %>% select(MCRNUM, name_clean, state, city),
  name_col1 = "name_clean", name_col2 = "name_clean", id_col1 = "ID", id_col2 = "MCRNUM",
  block_vars = c("state", "city"), threshold = 0.75) %>%
  mutate(strategy = "jaccard_statecity", zip_exact = 0)

candidates <- bind_rows(cand.zip %>% select(ID = id1, MCRNUM = id2, score, strategy, zip_exact),
                        cand.city %>% select(ID = id1, MCRNUM = id2, score, strategy, zip_exact)) %>%
  left_join(conv.variants %>% distinct(ID, eff_year), by = "ID", relationship = "many-to-many") %>%
  distinct(ID, MCRNUM, score, strategy, zip_exact, eff_year) %>%
  # the earlier provider must be observed with financials in at least one pre-designation year
  left_join(hcris.years, by = "MCRNUM") %>%
  filter(y_min < eff_year) %>%
  # accept: strong name match, or a moderate name match at the same zip
  filter(score >= 0.75 | (zip_exact == 1 & score >= 0.50))

# One earlier number per converter: best score, then zip-exact, then the lowest provider number
best <- candidates %>%
  group_by(ID, MCRNUM) %>% summarize(score = max(score), zip_exact = max(zip_exact), strategy = first(strategy), eff_year = first(eff_year), y_min = first(y_min), y_max = first(y_max), .groups = "drop") %>%
  group_by(ID) %>% arrange(desc(score), desc(zip_exact), MCRNUM, .by_group = TRUE) %>%
  mutate(n_candidates = n(), runner_up = lead(score)) %>% slice(1) %>% ungroup()


# Existing (time-invariant) crosswalk, used from designation onward and to screen collisions
xw.existing <- read_csv('data/output/unique_hcris.csv', show_col_types = FALSE) %>%
  mutate(ID = as.character(ID), MCRNUM = str_pad(str_remove_all(as.character(MCRNUM), "\\D"), 6, side = "left", pad = "0"))

# Screen out matches that cannot be a converter's own earlier number:
#   an earlier number claimed by two converters, one that is another hospital's number in the
#   existing crosswalk, or one whose reports continue more than two years after designation
other.numbers <- xw.existing %>% filter(!ID %in% best$ID) %>% distinct(MCRNUM) %>% pull(MCRNUM)
best <- best %>%
  group_by(MCRNUM) %>% mutate(claimed_by = n()) %>% ungroup() %>%
  mutate(exclude = case_when(claimed_by > 1 ~ "claimed by two converters",
                             MCRNUM %in% other.numbers ~ "another hospital's number",
                             y_max > eff_year + 2 ~ "reports continue after designation",
                             TRUE ~ NA_character_),
         match_strict = score >= 0.90 | (zip_exact == 1 & score >= 0.75))
best <- best %>% filter(is.na(exclude))

# Year-varying crosswalk: earlier number before designation, existing crosswalk number from designation on

xw.yv <- aha %>% select(ID, year, eff_year) %>%
  left_join(xw.existing %>% select(ID, year, MCRNUM_post = MCRNUM), by = c("ID", "year")) %>%
  left_join(best %>% select(ID, MCRNUM_pre = MCRNUM, match_strict), by = "ID") %>%
  mutate(MCRNUM_yv = case_when(!is.na(eff_year) & year < eff_year & !is.na(MCRNUM_pre) ~ MCRNUM_pre, TRUE ~ MCRNUM_post),
         source = case_when(!is.na(eff_year) & year < eff_year & !is.na(MCRNUM_pre) ~ "pre-designation match", TRUE ~ "existing crosswalk"),
         ## in the designation year a converter may file partial-year reports under both numbers
         MCRNUM_pre_desig = if_else(!is.na(eff_year) & year == eff_year & !is.na(MCRNUM_pre), MCRNUM_pre, NA_character_),
         match_strict = if_else(!is.na(MCRNUM_pre), match_strict, NA)) %>%
  select(ID, year, MCRNUM_yv, MCRNUM_pre_desig, match_strict, source)

# How many converter pre-designation years now find a cost report with financials?
check <- xw.yv %>% filter(source == "pre-designation match") %>%
  left_join(hcris.noncah %>% distinct(MCRNUM, year) %>% mutate(found = 1), by = c("MCRNUM_yv" = "MCRNUM", "year")) %>%
  summarize(pre_years = n(), found = sum(found, na.rm = TRUE))

write_csv(xw.yv, 'data/output/hcris_xw_yearvarying.csv')
write_csv(best %>% select(ID, eff_year, MCRNUM_pre = MCRNUM, score, zip_exact, strategy, match_strict, n_candidates, runner_up, hcris_y_min = y_min, hcris_y_max = y_max),
          'data/output/fuzzy/qa_hcris_yearvarying.csv')
