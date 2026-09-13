# CAH Supplement to AHA Fuzzy Matching
#
# Called by fuzzy-match.R (expects data already loaded)
#
# Ensemble matching strategy:
#   1. Stringdist (Jaro-Winkler) on preprocessed names, state+zip blocking
#   2. Stringdist (Jaro-Winkler) on preprocessed names, state+city blocking
#   3. Word-level Jaccard on preprocessed names, state+zip blocking
#   4. State-only fallback (stricter threshold) for unmatched IDs
#
# Selection: sort by score (desc), then city_exact (desc), zip_exact (desc), eff_date (asc), row ID (asc)
#
# Required objects from fuzzy-match.R:
#   - cah.supplement, aha.collapsed
#
# Output:
#   - data/output/unique_cah.csv (main crosswalk)
#   - data/output/fuzzy/qa_cah_diagnostics.csv (near-ties for manual review)
#   - data/output/fuzzy/qa_cah_summary.csv (match quality summaries)
#   - data/output/fuzzy/qa_cah_by_state.csv (match rates by state)

# =============================================================================
# Prepare CAH Supplement Data
# =============================================================================

# Create unique row ID (NOT eff_date, which can have duplicates)
cah.small <- cah.supplement %>%
  mutate(
    cah_row_id = row_number(),  # Unique ID for each CAH record
    zip = str_pad(substr(as.character(zip), 1, 5), width = 5, side = "left", pad = "0"),
    zip = if_else(str_detect(zip, "^[0-9]{5}$") & zip != "00000", zip, NA_character_)
  ) %>%
  select(cah_row_id, name, city, state, zip, eff_date) %>%
  mutate(
    name_raw = str_to_lower(name),
    name_clean = preprocess_hospital_name(name),
    state = str_to_lower(state),
    city = str_to_lower(city)
  ) %>%
  filter(!is.na(name_clean), name_clean != "")

# Note: aha.collapsed is created in fuzzy-match.R

# =============================================================================
# Strategy 1: Stringdist on preprocessed names, state+zip blocking
# =============================================================================

fuzzy.merge.1 <- merge_plus(
  data1 = aha.collapsed %>% select(aha_id, ID, name = name_clean, state, city, zip),
  data2 = cah.small %>% select(cah_row_id, name = name_clean, state, city, zip, eff_date),
  by = c("name", "city"),
  unique_key_1 = "aha_id",
  unique_key_2 = "cah_row_id",
  match_type = "multivar",
  multivar_settings = build_multivar_settings(
    compare_type = c("stringdist", "stringdist"),
    wgts = c(0.8, 0.2),
    blocks = c("state", "zip")
  )
)

matches.1 <- as_tibble(fuzzy.merge.1$matches) %>%
  filter(!is.na(multivar_score)) %>%
  # zip is a blocking variable, so join back to get separate zips
  left_join(aha.collapsed %>% select(ID, zip_aha = zip), by = "ID") %>%
  left_join(cah.small %>% select(cah_row_id, zip_cah = zip), by = "cah_row_id") %>%
  select(ID, cah_row_id, state,
         name_aha = name_1, name_cah = name_2,
         city_aha = city_1, city_cah = city_2,
         zip_aha, zip_cah,
         eff_date,
         score = multivar_score) %>%
  mutate(strategy = "stringdist_statezip")

# =============================================================================
# Strategy 2: Stringdist on preprocessed names, state+city blocking
# =============================================================================

fuzzy.merge.2 <- merge_plus(
  data1 = aha.collapsed %>% select(aha_id, ID, name = name_clean, state, city, zip),
  data2 = cah.small %>% select(cah_row_id, name = name_clean, state, city, zip, eff_date),
  by = c("name"),
  unique_key_1 = "aha_id",
  unique_key_2 = "cah_row_id",
  match_type = "multivar",
  multivar_settings = build_multivar_settings(
    compare_type = c("stringdist"),
    wgts = 1,
    blocks = c("state", "city")
  )
)

matches.2 <- as_tibble(fuzzy.merge.2$matches) %>%
  filter(!is.na(multivar_score)) %>%
  left_join(aha.collapsed %>% select(ID, city_aha = city, zip_aha = zip), by = "ID") %>%
  left_join(cah.small %>% select(cah_row_id, city_cah = city, zip_cah = zip), by = "cah_row_id") %>%
  select(ID, cah_row_id, state,
         name_aha = name_1, name_cah = name_2,
         city_aha, city_cah, zip_aha, zip_cah,
         eff_date,
         score = multivar_score) %>%
  mutate(strategy = "stringdist_statecity")

# =============================================================================
# Strategy 3: Word-level Jaccard on preprocessed names, state+zip blocking
# =============================================================================

matches.3 <- match_jaccard_words(
  data1 = aha.collapsed %>% select(ID, name_clean, state, city, zip),
  data2 = cah.small %>% select(cah_row_id, name_clean, state, city, zip),
  name_col1 = "name_clean",
  name_col2 = "name_clean",
  id_col1 = "ID",
  id_col2 = "cah_row_id",
  block_vars = c("state", "zip"),
  threshold = 0.75
) %>%
  rename(ID = id1, cah_row_id = id2, name_aha = name1, name_cah = name2,
         zip_aha = zip) %>%  # zip comes from block_vars (same for both)
  left_join(aha.collapsed %>% select(ID, city_aha = city), by = "ID") %>%
  left_join(cah.small %>% select(cah_row_id, city_cah = city, zip_cah = zip, eff_date), by = "cah_row_id") %>%
  mutate(strategy = "word_jaccard")

# =============================================================================
# Combine Strategies (before selection)
# =============================================================================

all_candidates <- bind_rows(matches.1, matches.2, matches.3) %>%
  filter(score >= 0.75) %>%
  # Add tie-breaker columns
  mutate(
    city_exact = as.integer(city_aha == city_cah),
    zip_exact = as.integer(zip_aha == zip_cah)
  ) %>%
  # Remove duplicates (same ID-cah_row_id pair from different strategies)
  group_by(ID, cah_row_id) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup()

# =============================================================================
# Create Diagnostics Table (top 3 candidates per ID with score gaps)
# =============================================================================

diagnostics <- all_candidates %>%
  group_by(ID) %>%
  arrange(desc(score), desc(city_exact), desc(zip_exact), eff_date, cah_row_id) %>%
  mutate(
    rank = row_number(),
    score_gap = score - lead(score, default = 0)
  ) %>%
  filter(rank <= 3) %>%
  ungroup() %>%
  select(ID, rank, cah_row_id, eff_date, score, score_gap, city_exact, zip_exact,
         name_aha, name_cah, city_aha, city_cah, zip_aha, zip_cah, strategy)

# Flag ambiguous matches (score gap < 0.05 between top 2)
ambiguous_ids <- diagnostics %>%
  filter(rank == 1, score_gap < 0.05) %>%
  pull(ID)

diagnostics <- diagnostics %>%
  mutate(ambiguous = ID %in% ambiguous_ids)

# =============================================================================
# Select Best Match with Deterministic Tie-Breaking
# =============================================================================

# Tie-breaker order: score > city_exact > zip_exact > earliest eff_date > row_id
best_matches <- all_candidates %>%
  group_by(ID) %>%
  arrange(desc(score), desc(city_exact), desc(zip_exact), eff_date, cah_row_id) %>%
  slice(1) %>%
  ungroup()

# =============================================================================
# Strategy 4: State-only fallback for unmatched IDs (stricter threshold)
# =============================================================================

matched_ids <- best_matches$ID
unmatched_aha <- aha.collapsed %>% filter(!ID %in% matched_ids)

if (nrow(unmatched_aha) > 0) {
  # State-only blocking with stricter 0.85 threshold
  matches.fallback <- match_jaccard_words(
    data1 = unmatched_aha %>% select(ID, name_clean, state, city, zip),
    data2 = cah.small %>% select(cah_row_id, name_clean, state, city, zip),
    name_col1 = "name_clean",
    name_col2 = "name_clean",
    id_col1 = "ID",
    id_col2 = "cah_row_id",
    block_vars = c("state"),  # State-only
    threshold = 0.85  # Stricter threshold
  ) %>%
    rename(ID = id1, cah_row_id = id2, name_aha = name1, name_cah = name2) %>%
    # state comes from block_vars; need to add city/zip from both sides
    left_join(aha.collapsed %>% select(ID, city_aha = city, zip_aha = zip), by = "ID") %>%
    left_join(cah.small %>% select(cah_row_id, city_cah = city, zip_cah = zip, eff_date), by = "cah_row_id") %>%
    mutate(
      strategy = "fallback_state_only",
      city_exact = as.integer(city_aha == city_cah),
      zip_exact = as.integer(zip_aha == zip_cah)
    ) %>%
    group_by(ID) %>%
    arrange(desc(score), desc(city_exact), desc(zip_exact), eff_date, cah_row_id) %>%
    slice(1) %>%
    ungroup()

  # Add fallback matches to diagnostics
  if (nrow(matches.fallback) > 0) {
    diagnostics.fallback <- matches.fallback %>%
      mutate(rank = 1, score_gap = NA_real_, ambiguous = FALSE) %>%
      select(ID, rank, cah_row_id, eff_date, score, score_gap, city_exact, zip_exact,
             name_aha, name_cah, city_aha, city_cah, zip_aha, zip_cah, strategy, ambiguous)
    diagnostics <- bind_rows(diagnostics, diagnostics.fallback)
  }

  # Combine with main matches
  best_matches <- bind_rows(best_matches, matches.fallback)
}

# =============================================================================
# Create Final Output (one row per AHA ID with earliest CAH date)
# =============================================================================

# Note: An AHA hospital might match multiple CAH records if it received CAH
# designation multiple times (rare but possible). We take the earliest date.
unique_cah <- best_matches %>%
  group_by(ID) %>%
  summarize(
    first_date = min(eff_date, na.rm = TRUE),
    best_score = max(score),
    strategy = first(strategy),
    .groups = "drop"
  ) %>%
  mutate(cah_sup = 1)

# =============================================================================
# QA Summary Statistics
# =============================================================================

# Handle case where best_matches might be empty
cah_scores <- if (nrow(best_matches) > 0) best_matches$score else NA_real_

qa_summary <- tibble(
  metric = c(
    "unique_aha_ids",
    "unique_cah_records",
    "matched_aha_ids",
    "fallback_matches",
    "ambiguous_matches",
    "mean_score",
    "median_score",
    "p10_score",
    "p90_score"
  ),
  value = c(
    n_distinct(aha.collapsed$ID),
    n_distinct(cah.small$cah_row_id),
    nrow(unique_cah),
    sum(best_matches$strategy == "fallback_state_only", na.rm = TRUE),
    length(ambiguous_ids),
    mean(cah_scores, na.rm = TRUE),
    median(cah_scores, na.rm = TRUE),
    if (length(cah_scores) > 1 && !all(is.na(cah_scores))) quantile(cah_scores, 0.10, na.rm = TRUE) else NA_real_,
    if (length(cah_scores) > 1 && !all(is.na(cah_scores))) quantile(cah_scores, 0.90, na.rm = TRUE) else NA_real_
  )
)

# Match rate by state
qa_by_state <- aha.collapsed %>%
  left_join(unique_cah %>% select(ID, cah_sup), by = "ID") %>%
  group_by(state) %>%
  summarize(
    n_hospitals = n(),
    n_cah = sum(cah_sup == 1, na.rm = TRUE),
    cah_rate = n_cah / n_hospitals,
    .groups = "drop"
  ) %>%
  arrange(desc(cah_rate))

# =============================================================================
# Save Outputs
# =============================================================================

write_csv(unique_cah, "data/output/unique_cah.csv")
write_csv(diagnostics, "data/output/fuzzy/qa_cah_diagnostics.csv")
write_csv(qa_summary, "data/output/fuzzy/qa_cah_summary.csv")
write_csv(qa_by_state, "data/output/fuzzy/qa_cah_by_state.csv")

