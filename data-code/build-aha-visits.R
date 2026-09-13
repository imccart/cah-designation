# Meta --------------------------------------------------------------------
## Outpatient visits, emergency department visits, and births from the AHA
## annual survey. These four fields are not in the aha-data repo's shared
## keep-list, so they never reach data/input/aha_data.csv. We read them from the
## same 1994-2021 WRDS extract that aha-data's code/1-historic.R reads, which is
## symlinked into data/input like every other external source.
##
##   VEM    = emergency department visits
##   VTOT   = total outpatient visits (VEM + VOTH)
##   VOTH   = outpatient visits other than emergency department
##   BIRTHS = births
##
## Input:  data/input/aha_survey_1994_recent.csv
## Output: data/output/aha_visits.csv (ID, year, VEM, VTOT, VOTH, BIRTHS)

## A zero surrounded by positive volume is a reporting gap rather than a real
## zero, so those become missing and are then filled by linear interpolation
## within hospital, as _build-estimation-data.r does for the financial and
## capacity measures. A zero that begins and never recovers is a service closing
## (an obstetrics unit, say) and is kept.
aha.visits <- read_csv("data/input/aha_survey_1994_recent.csv",
                       col_select = c(ID, YEAR, VEM, VTOT, VOTH, BIRTHS),
                       col_types = cols(ID = col_character()), progress = FALSE) %>%
  rename(year = YEAR) %>%
  arrange(ID, year) %>%
  group_by(ID) %>%
  mutate(across(c(VEM, VTOT, VOTH, BIRTHS),
                ~ {
                  positive <- !is.na(.x) & .x > 0
                  before <- lag(cumany(positive), default = FALSE)
                  after  <- lead(rev(cumany(rev(positive))), default = FALSE)
                  if_else(!is.na(.x) & .x == 0 & before & after, NA_real_, .x)
                })) %>%
  mutate(across(c(VEM, VTOT, VOTH, BIRTHS),
                ~ if (sum(!is.na(.x)) >= 2) na.approx(.x, x = year, na.rm = FALSE) else .x)) %>%
  ungroup()

write_csv(aha.visits, "data/output/aha_visits.csv")
