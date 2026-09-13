# HCRIS to AHA Fuzzy Matching
#
# Called by fuzzy-match.R (expects data already loaded)
#
# Ensemble matching strategy (time-invariant, no year blocking):
#   1. Stringdist (Jaro-Winkler) on preprocessed names, state blocking
#   2. Stringdist (Jaro-Winkler) on preprocessed names, state+city blocking
#   3. Word-level Jaccard on preprocessed names, state+zip blocking
#   4. Word-level Jaccard with all historical name variants
#   5. State-only fallback (stricter threshold) for unmatched IDs
#
# Matching is done at the hospital level (not hospital-year) since AHA ID and
# MCRNUM are time-invariant identifiers. The resulting crosswalk is then
# applied to all ID-year pairs.
#
# Selection: sort by score (desc), then city_exact (desc), zip_exact (desc), alphabetical MCRNUM
#
# Priority: Direct MCRNUM from AHA > Fuzzy matches
#
# Required objects from fuzzy-match.R:
#   - hcris.data, aha.collapsed, aha.combine
#
# Output:
#   - data/output/unique_hcris.csv (main crosswalk)
#   - data/output/fuzzy/qa_hcris_diagnostics.csv (near-ties for manual review)
#   - data/output/fuzzy/qa_hcris_summary.csv (match quality summaries)
#   - data/output/fuzzy/qa_hcris_by_state.csv (match rates by state)

# =============================================================================
# Prepare HCRIS Data (collapse to unique hospitals)
# =============================================================================

# Collapse to unique MCRNUM with most common name/city/zip
hcris.collapsed <- hcris.data %>%
  mutate(
    zip = str_pad(substr(as.character(zip), 1, 5), width = 5, side = "left", pad = "0"),
    zip = if_else(str_detect(zip, "^[0-9]{5}$") & zip != "00000", zip, NA_character_)
  ) %>%
  select(name, city, state, zip, MCRNUM) %>%
  mutate(
    name_clean = preprocess_hospital_name(name),
    state = str_to_lower(state),
    city = str_to_lower(city)
  ) %>%
  filter(!is.na(name_clean), name_clean != "") %>%
  group_by(MCRNUM) %>%
  # Use most frequent name/city/zip for each MCRNUM
  summarize(
    name_clean = modal_value(name_clean),
    state = first(state),
    city = modal_value(city),
    zip = modal_value(zip),
    .groups = "drop"
  ) %>%
  mutate(hcris_id = row_number())

# Also keep all name variants for matching (increases match chances)
hcris.variants <- hcris.data %>%
  mutate(
    zip = str_pad(substr(as.character(zip), 1, 5), width = 5, side = "left", pad = "0"),
    zip = if_else(str_detect(zip, "^[0-9]{5}$") & zip != "00000", zip, NA_character_)
  ) %>%
  select(name, city, state, zip, MCRNUM) %>%
  mutate(
    name_clean = preprocess_hospital_name(name),
    state = str_to_lower(state),
    city = str_to_lower(city)
  ) %>%
  filter(!is.na(name_clean), name_clean != "") %>%
  distinct(MCRNUM, name_clean, state, city, zip)

# Note: aha.collapsed is created in fuzzy-match.R

# =============================================================================
# Strategy 1: Stringdist on preprocessed names, state blocking
# =============================================================================

fuzzy.merge.1 <- merge_plus(
  data1 = aha.collapsed %>% select(aha_id, ID, name = name_clean, state, city, zip),
  data2 = hcris.collapsed %>% select(hcris_id, MCRNUM, name = name_clean, state, city, zip),
  by = c("name", "city", "zip"),
  unique_key_1 = "aha_id",
  unique_key_2 = "hcris_id",
  match_type = "multivar",
  multivar_settings = build_multivar_settings(
    compare_type = c("stringdist", "stringdist", "indicator"),
    wgts = c(0.4, 0.4, 0.2),
    blocks = c("state")
  )
)

matches.1 <- as_tibble(fuzzy.merge.1$matches) %>%
  filter(!is.na(multivar_score)) %>%
  select(ID, MCRNUM, state,
         name_aha = name_1, name_hcris = name_2,
         city_aha = city_1, city_hcris = city_2,
         zip_aha = zip_1, zip_hcris = zip_2,
         score = multivar_score) %>%  # Use weighted combination of name/city/zip
  mutate(strategy = "stringdist_state")

# =============================================================================
# Strategy 2: Stringdist on preprocessed names, state+city blocking
# =============================================================================

fuzzy.merge.2 <- merge_plus(
  data1 = aha.collapsed %>% select(aha_id, ID, name = name_clean, state, city, zip),
  data2 = hcris.collapsed %>% select(hcris_id, MCRNUM, name = name_clean, state, city, zip),
  by = c("name"),
  unique_key_1 = "aha_id",
  unique_key_2 = "hcris_id",
  match_type = "multivar",
  multivar_settings = build_multivar_settings(
    compare_type = c("stringdist"),
    wgts = 1,
    blocks = c("state", "city")
  )
)

matches.2 <- as_tibble(fuzzy.merge.2$matches) %>%
  filter(!is.na(multivar_score)) %>%
  # Need to join back city/zip info

  left_join(aha.collapsed %>% select(ID, city_aha = city, zip_aha = zip), by = "ID") %>%
  left_join(hcris.collapsed %>% select(MCRNUM, city_hcris = city, zip_hcris = zip), by = "MCRNUM") %>%
  select(ID, MCRNUM, state,
         name_aha = name_1, name_hcris = name_2,
         city_aha, city_hcris, zip_aha, zip_hcris,
         score = multivar_score) %>%
  mutate(strategy = "stringdist_statecity")

# =============================================================================
# Strategy 3: Word-level Jaccard on preprocessed names, state+zip blocking
# =============================================================================

matches.3 <- match_jaccard_words(
  data1 = aha.collapsed %>% select(ID, name_clean, state, city, zip),
  data2 = hcris.collapsed %>% select(MCRNUM, name_clean, state, city, zip),
  name_col1 = "name_clean",
  name_col2 = "name_clean",
  id_col1 = "ID",
  id_col2 = "MCRNUM",
  block_vars = c("state", "zip"),
  threshold = 0.75
) %>%
  rename(ID = id1, MCRNUM = id2, name_aha = name1, name_hcris = name2,
         zip_aha = zip) %>%  # zip comes from block_vars (same for both)
  left_join(aha.collapsed %>% select(ID, city_aha = city), by = "ID") %>%
  left_join(hcris.collapsed %>% select(MCRNUM, city_hcris = city, zip_hcris = zip), by = "MCRNUM") %>%
  mutate(strategy = "word_jaccard")

# =============================================================================
# Strategy 4: Match using all name variants (catches name changes over time)
# =============================================================================

matches.4 <- match_jaccard_words(
  data1 = aha.collapsed %>% select(ID, name_clean, state, city, zip),
  data2 = hcris.variants %>% select(MCRNUM, name_clean, state, city, zip),
  name_col1 = "name_clean",
  name_col2 = "name_clean",
  id_col1 = "ID",
  id_col2 = "MCRNUM",
  block_vars = c("state", "zip"),
  threshold = 0.75
) %>%
  rename(ID = id1, MCRNUM = id2, name_aha = name1, name_hcris = name2,
         zip_aha = zip) %>%  # zip comes from block_vars (same for both)
  left_join(aha.collapsed %>% select(ID, city_aha = city), by = "ID") %>%
  left_join(hcris.collapsed %>% select(MCRNUM, city_hcris = city, zip_hcris = zip), by = "MCRNUM") %>%
  mutate(strategy = "word_jaccard_variants")

# =============================================================================
# Combine Strategies (before selection)
# =============================================================================

all_candidates <- bind_rows(matches.1, matches.2, matches.3, matches.4) %>%
  filter(score >= 0.75) %>%
  # Add tie-breaker columns

  mutate(
    city_exact = as.integer(city_aha == city_hcris),
    zip_exact = as.integer(zip_aha == zip_hcris)
  ) %>%
  # Remove duplicates (same ID-MCRNUM pair from different strategies)
  group_by(ID, MCRNUM) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup()

# =============================================================================
# Create Diagnostics Table (top 3 candidates per ID with score gaps)
# =============================================================================

diagnostics <- all_candidates %>%
  group_by(ID) %>%
  arrange(desc(score), desc(city_exact), desc(zip_exact), MCRNUM) %>%
  mutate(
    rank = row_number(),
    score_gap = score - lead(score, default = 0)
  ) %>%
  filter(rank <= 3) %>%
  ungroup() %>%
  select(ID, rank, MCRNUM, score, score_gap, city_exact, zip_exact,
         name_aha, name_hcris, city_aha, city_hcris, zip_aha, zip_hcris, strategy)

# Flag ambiguous matches (score gap < 0.05 between top 2)
ambiguous_ids <- diagnostics %>%
  filter(rank == 1, score_gap < 0.05) %>%
  pull(ID)

diagnostics <- diagnostics %>%
  mutate(ambiguous = ID %in% ambiguous_ids)

# =============================================================================
# Select Best Match with Deterministic Tie-Breaking
# =============================================================================

# Tie-breaker order: score > city_exact > zip_exact > alphabetical MCRNUM
best_matches <- all_candidates %>%
  group_by(ID) %>%
  arrange(desc(score), desc(city_exact), desc(zip_exact), MCRNUM) %>%
  slice(1) %>%
  ungroup() %>%
  # Now ensure 1:1 (best match per MCRNUM)
  group_by(MCRNUM) %>%
  arrange(desc(score), desc(city_exact), desc(zip_exact), ID) %>%
  slice(1) %>%
  ungroup()

# =============================================================================
# Strategy 5: State-only fallback for unmatched IDs (stricter threshold)
# =============================================================================

matched_ids <- best_matches$ID
unmatched_aha <- aha.collapsed %>% filter(!ID %in% matched_ids)

if (nrow(unmatched_aha) > 0) {
  # State-only blocking with stricter 0.85 threshold
  matches.fallback <- match_jaccard_words(
    data1 = unmatched_aha %>% select(ID, name_clean, state, city, zip),
    data2 = hcris.collapsed %>% filter(!MCRNUM %in% best_matches$MCRNUM) %>%
      select(MCRNUM, name_clean, state, city, zip),
    name_col1 = "name_clean",
    name_col2 = "name_clean",
    id_col1 = "ID",
    id_col2 = "MCRNUM",
    block_vars = c("state"),  # State-only
    threshold = 0.85  # Stricter threshold
  ) %>%
    rename(ID = id1, MCRNUM = id2, name_aha = name1, name_hcris = name2) %>%
    # state comes from block_vars; need to add city/zip from both sides
    left_join(aha.collapsed %>% select(ID, city_aha = city, zip_aha = zip), by = "ID") %>%
    left_join(hcris.collapsed %>% select(MCRNUM, city_hcris = city, zip_hcris = zip), by = "MCRNUM") %>%
    mutate(
      strategy = "fallback_state_only",
      city_exact = as.integer(city_aha == city_hcris),
      zip_exact = as.integer(zip_aha == zip_hcris)
    ) %>%
    group_by(ID) %>%
    arrange(desc(score), desc(city_exact), desc(zip_exact), MCRNUM) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(MCRNUM) %>%
    arrange(desc(score), desc(city_exact), desc(zip_exact), ID) %>%
    slice(1) %>%
    ungroup()

  # Add fallback matches to diagnostics
  if (nrow(matches.fallback) > 0) {
    diagnostics.fallback <- matches.fallback %>%
      mutate(rank = 1, score_gap = NA_real_, ambiguous = FALSE) %>%
      select(ID, rank, MCRNUM, score, score_gap, city_exact, zip_exact,
             name_aha, name_hcris, city_aha, city_hcris, zip_aha, zip_hcris,
             strategy, ambiguous)
    diagnostics <- bind_rows(diagnostics, diagnostics.fallback)
  }

  # Combine with main matches
  best_matches <- bind_rows(best_matches, matches.fallback)
}

# =============================================================================
# Combine with Direct MCRNUM from AHA (Priority)
# =============================================================================

# Direct MCRNUM from AHA data - collapse to hospital level
aha.crosswalk.direct <- aha.combine %>%
  mutate(ID = as.character(ID)) %>%
  select(ID, MCRNUM) %>%
  filter(!is.na(MCRNUM)) %>%
  group_by(ID, MCRNUM) %>%
  summarize(n_years = n(), .groups = "drop") %>%
  group_by(ID) %>%
  slice_max(n_years, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(ID, MCRNUM_direct = MCRNUM)

# Fuzzy matches (already at hospital level)
aha.crosswalk.fuzzy <- best_matches %>%
  distinct(ID, MCRNUM) %>%
  rename(MCRNUM_fuzzy = MCRNUM)

# Final hospital-level crosswalk: direct > fuzzy
aha.crosswalk <- aha.collapsed %>%
  select(ID) %>%
  left_join(aha.crosswalk.direct, by = "ID") %>%
  left_join(aha.crosswalk.fuzzy, by = "ID") %>%
  mutate(
    MCRNUM = coalesce(MCRNUM_direct, MCRNUM_fuzzy),
    source = case_when(
      !is.na(MCRNUM_direct) ~ "direct",
      !is.na(MCRNUM_fuzzy) ~ "fuzzy",
      TRUE ~ "unmatched"
    )
  ) %>%
  filter(!is.na(MCRNUM)) %>%
  select(ID, MCRNUM, source)

# Expand to all ID-year pairs
unique_hcris <- aha.combine %>%
  mutate(ID = as.character(ID)) %>%
  distinct(ID, year) %>%
  inner_join(aha.crosswalk %>% select(ID, MCRNUM), by = "ID") %>%
  arrange(ID, year)

# =============================================================================
# QA Summary Statistics
# =============================================================================

# Handle case where best_matches might be empty
fuzzy_scores <- if (nrow(best_matches) > 0) best_matches$score else NA_real_

qa_summary <- tibble(
  metric = c(
    "unique_aha_ids",
    "unique_mcrnums",
    "direct_matches",
    "fuzzy_matches",
    "fallback_matches",
    "total_matched_hospitals",
    "unmatched_aha_ids",
    "ambiguous_matches",
    "mean_fuzzy_score",
    "median_fuzzy_score",
    "p10_fuzzy_score",
    "p90_fuzzy_score"
  ),
  value = c(
    n_distinct(aha.collapsed$ID),
    n_distinct(hcris.collapsed$MCRNUM),
    sum(aha.crosswalk$source == "direct"),
    sum(aha.crosswalk$source == "fuzzy" & !aha.crosswalk$ID %in%
          (best_matches %>% filter(strategy == "fallback_state_only") %>% pull(ID))),
    sum(best_matches$strategy == "fallback_state_only", na.rm = TRUE),
    nrow(aha.crosswalk),
    n_distinct(aha.collapsed$ID) - nrow(aha.crosswalk),
    length(ambiguous_ids),
    mean(fuzzy_scores, na.rm = TRUE),
    median(fuzzy_scores, na.rm = TRUE),
    if (length(fuzzy_scores) > 1 && !all(is.na(fuzzy_scores))) quantile(fuzzy_scores, 0.10, na.rm = TRUE) else NA_real_,
    if (length(fuzzy_scores) > 1 && !all(is.na(fuzzy_scores))) quantile(fuzzy_scores, 0.90, na.rm = TRUE) else NA_real_
  )
)

# Match rate by state
qa_by_state <- aha.collapsed %>%
  left_join(aha.crosswalk, by = "ID") %>%
  group_by(state) %>%
  summarize(
    n_hospitals = n(),
    n_matched = sum(!is.na(MCRNUM)),
    match_rate = n_matched / n_hospitals,
    n_direct = sum(source == "direct", na.rm = TRUE),
    n_fuzzy = sum(source == "fuzzy", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match_rate)

# =============================================================================
# Save Outputs
# =============================================================================

write_csv(unique_hcris, "data/output/unique_hcris.csv")
write_csv(diagnostics, "data/output/fuzzy/qa_hcris_diagnostics.csv")
write_csv(qa_summary, "data/output/fuzzy/qa_hcris_summary.csv")
write_csv(qa_by_state, "data/output/fuzzy/qa_hcris_by_state.csv")

