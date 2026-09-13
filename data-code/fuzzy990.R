# Form 990 to AHA Fuzzy Matching
#
# Called by fuzzy-match.R (expects data already loaded)
#
# Matching approach:
#   - Uses ALL historical name/state/zip variants for each hospital
#   - This allows matching via old hospital names that may appear in Form 990
#   - Best match per ID is selected after trying all variants
#
# Ensemble matching strategy:
#   1. Stringdist (Jaro-Winkler) on preprocessed names, state+zip blocking
#   2. Word-level Jaccard on preprocessed names, state+zip blocking
#   3. State-only fallback (stricter threshold) for unmatched IDs
#
# Filtering approach:
#   - Keyword filter requires hospital-related term: hospital, medical, health
#   - Threshold: 0.80 minimum score
#
# EIN handling:
#   - Prior manual matches take priority over fuzzy matches
#   - Fuzzy matches fill in gaps for unmatched AHA IDs
#   - Multiple EINs per ID are ranked and stored (ein_1, ein_2, ein_3)
#
# Tie-breaking: exact zip match > score > alphabetical EIN
#
# Required objects from fuzzy-match.R:
#   - form990.data, aha.variants, form990.small
#
# Output:
#   - data/output/unique_990.csv (main crosswalk)
#   - data/output/fuzzy/qa_990_diagnostics.csv (near-ties for manual review)
#   - data/output/fuzzy/qa_990_summary.csv (match quality summaries)
#   - data/output/fuzzy/qa_990_by_state.csv (match rates by state)

# Unique IDs for QA stats
aha.unique_ids <- aha.variants %>% distinct(ID)

# =============================================================================
# Get Prior Manual Matches
# =============================================================================

manual_matches <- form990.data %>%
  filter(!is.na(ID_hospital_1)) %>%
  select(ein, ID_hospital_1, ID_hospital_2, ID_hospital_3) %>%
  pivot_longer(
    cols = c(ID_hospital_1, ID_hospital_2, ID_hospital_3),
    names_to = "slot",
    values_to = "ID"
  ) %>%
  filter(!is.na(ID)) %>%
  mutate(ID = as.character(ID)) %>%
  distinct(ein, ID)

# IDs that already have manual matches (manual takes priority; no fuzzy additions)
manual_matched_ids <- unique(manual_matches$ID)

# =============================================================================
# Strategy 1: Stringdist on preprocessed names, state+zip blocking
# =============================================================================

# Note: aha.variants and form990.small are pre-filtered in fuzzy-match.R
# to exclude NA state/zip, preventing NA-to-NA matching

fuzzy.merge.1 <- merge_plus(
  data1 = aha.variants %>% select(variant_id, ID, name = name_clean, state, zip),
  data2 = form990.small %>% select(form990_row_id, name = name_clean, state, zip),
  by = c("name"),
  unique_key_1 = "variant_id",
  unique_key_2 = "form990_row_id",
  match_type = "multivar",
  multivar_settings = build_multivar_settings(
    compare_type = c("stringdist"),
    wgts = 1,
    blocks = c("state", "zip")
  )
)

matches.1 <- as_tibble(fuzzy.merge.1$matches) %>%
  filter(!is.na(multivar_score)) %>%
  left_join(aha.variants %>% select(variant_id, zip_aha = zip), by = "variant_id") %>%
  left_join(form990.small %>% select(form990_row_id, ein, zip_990 = zip), by = "form990_row_id") %>%
  select(ID, ein, state,
         name_aha = name_1, name_990 = name_2,
         zip_aha, zip_990,
         score = multivar_score) %>%
  mutate(strategy = "stringdist_statezip")

# =============================================================================
# Strategy 2: Word-level Jaccard on preprocessed names, state+zip blocking
# =============================================================================

matches.2 <- match_jaccard_words(
  data1 = aha.variants %>% select(variant_id, ID, name_clean, state, zip),
  data2 = form990.small %>% select(ein, name_clean, state, zip),
  name_col1 = "name_clean",
  name_col2 = "name_clean",
  id_col1 = "variant_id",
  id_col2 = "ein",
  block_vars = c("state", "zip"),
  threshold = 0.75
) %>%
  rename(variant_id = id1, ein = id2, name_aha = name1, name_990 = name2,
         zip_aha = zip) %>%
  left_join(aha.variants %>% select(variant_id, ID), by = "variant_id") %>%
  left_join(form990.small %>% distinct(ein, zip_990 = zip), by = "ein", relationship = "many-to-many") %>%
  mutate(strategy = "word_jaccard")

# =============================================================================
# Combine Strategies and Select Best Match per ID
# =============================================================================

all_candidates <- bind_rows(matches.1, matches.2) %>%
  filter(score >= 0.80) %>%
  # Add tie-breaker columns
  mutate(
    zip_exact = as.integer(zip_aha == zip_990)
  ) %>%
  # Select best match per ID-ein pair (across all variants and strategies)
  group_by(ID, ein) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup()

# =============================================================================
# Create Diagnostics Table (top 3 candidates per ID with score gaps)
# =============================================================================

diagnostics <- all_candidates %>%
  group_by(ID) %>%
  arrange(desc(zip_exact), desc(score), ein) %>%
  mutate(
    rank = row_number(),
    score_gap = score - lead(score, default = 0)
  ) %>%
  filter(rank <= 3) %>%
  ungroup() %>%
  select(ID, rank, ein, score, score_gap, zip_exact,
         name_aha, name_990, zip_aha, zip_990, strategy)

# Flag ambiguous matches (score gap < 0.05 between top 2)
ambiguous_ids <- diagnostics %>%
  filter(rank == 1, score_gap < 0.05) %>%
  pull(ID)

diagnostics <- diagnostics %>%
  mutate(ambiguous = ID %in% ambiguous_ids)

# =============================================================================
# Strategy 3: State-only fallback for unmatched IDs (stricter threshold)
# =============================================================================

# IDs with at least one fuzzy candidate
matched_ids <- unique(all_candidates$ID)
unmatched_variants <- aha.variants %>%
  filter(!ID %in% matched_ids, !ID %in% manual_matched_ids)

if (nrow(unmatched_variants) > 0) {
  # State-only blocking with stricter 0.85 threshold
  matches.fallback <- match_jaccard_words(
    data1 = unmatched_variants %>% select(variant_id, ID, name_clean, state, zip),
    data2 = form990.small %>% select(ein, name_clean, state, zip),
    name_col1 = "name_clean",
    name_col2 = "name_clean",
    id_col1 = "variant_id",
    id_col2 = "ein",
    block_vars = c("state"),  # State-only
    threshold = 0.85  # Stricter threshold
  ) %>%
    rename(variant_id = id1, ein = id2, name_aha = name1, name_990 = name2) %>%
    left_join(aha.variants %>% select(variant_id, ID, zip_aha = zip), by = "variant_id") %>%
    left_join(form990.small %>% distinct(ein, zip_990 = zip), by = "ein", relationship = "many-to-many") %>%
    mutate(
      strategy = "fallback_state_only",
      zip_exact = as.integer(zip_aha == zip_990)
    ) %>%
    # Select best match per ID (across all variants)
    group_by(ID, ein) %>%
    slice_max(score, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Add fallback to candidates
  if (nrow(matches.fallback) > 0) {
    all_candidates <- bind_rows(all_candidates, matches.fallback)

    # Add to diagnostics
    diagnostics.fallback <- matches.fallback %>%
      group_by(ID) %>%
      arrange(desc(zip_exact), desc(score), ein) %>%
      mutate(
        rank = row_number(),
        score_gap = score - lead(score, default = 0)
      ) %>%
      filter(rank <= 3) %>%
      ungroup() %>%
      mutate(ambiguous = FALSE) %>%
      select(ID, rank, ein, score, score_gap, zip_exact,
             name_aha, name_990, zip_aha, zip_990, strategy, ambiguous)
    diagnostics <- bind_rows(diagnostics, diagnostics.fallback)
  }
}

# =============================================================================
# Select Top 3 Matches per ID (for fuzzy matches)
# =============================================================================

# Tie-breaker order: zip_exact > score > alphabetical EIN
fuzzy_ranked <- all_candidates %>%
  # Exclude EINs already in manual matches for this ID
  anti_join(manual_matches, by = c("ID", "ein")) %>%
  group_by(ID) %>%
  arrange(desc(zip_exact), desc(score), ein) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  ungroup()

# Pivot to wide format
fuzzy_wide <- fuzzy_ranked %>%
  select(ID, ein, score, rank) %>%
  pivot_wider(
    id_cols = ID,
    names_from = rank,
    values_from = c(ein, score),
    names_glue = "{.value}_{rank}"
  ) %>%
  select(ID, starts_with("ein_"), starts_with("score_"))

# =============================================================================
# Combine Manual and Fuzzy Matches
# =============================================================================

# Manual matches in wide format (one row per ID)
manual_wide <- manual_matches %>%
  group_by(ID) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3) %>%
  pivot_wider(
    id_cols = ID,
    names_from = rank,
    values_from = ein,
    names_glue = "ein_{rank}"
  ) %>%
  mutate(source = "manual")

# Fuzzy matches (only for IDs not in manual)
fuzzy_new <- fuzzy_wide %>%
  filter(!ID %in% manual_matched_ids) %>%
  mutate(source = "fuzzy")

# Combine: manual matches + fuzzy matches for new IDs
unique_990 <- bind_rows(manual_wide, fuzzy_new) %>%
  arrange(ID)

# =============================================================================
# QA Summary Statistics
# =============================================================================

# Handle case where fuzzy_ranked might be empty
fuzzy_scores <- if (nrow(fuzzy_ranked) > 0) fuzzy_ranked$score else NA_real_

qa_summary <- tibble(
  metric = c(
    "unique_aha_ids",
    "unique_aha_variants",
    "unique_990_eins",
    "manual_matched_ids",
    "fuzzy_matched_ids",
    "fallback_matches",
    "total_matched_ids",
    "unmatched_aha_ids",
    "ambiguous_matches",
    "ids_with_multiple_eins",
    "mean_fuzzy_score",
    "median_fuzzy_score",
    "p10_fuzzy_score",
    "p90_fuzzy_score"
  ),
  value = c(
    nrow(aha.unique_ids),
    nrow(aha.variants),
    n_distinct(form990.small$ein),
    n_distinct(manual_matches$ID),
    nrow(fuzzy_new),
    sum(all_candidates$strategy == "fallback_state_only", na.rm = TRUE),
    nrow(unique_990),
    nrow(aha.unique_ids) - nrow(unique_990),
    length(ambiguous_ids),
    sum(!is.na(unique_990$ein_2)),
    mean(fuzzy_scores, na.rm = TRUE),
    median(fuzzy_scores, na.rm = TRUE),
    if (length(fuzzy_scores) > 1 && !all(is.na(fuzzy_scores))) quantile(fuzzy_scores, 0.10, na.rm = TRUE) else NA_real_,
    if (length(fuzzy_scores) > 1 && !all(is.na(fuzzy_scores))) quantile(fuzzy_scores, 0.90, na.rm = TRUE) else NA_real_
  )
)

# Match rate by state
qa_by_state <- aha.unique_ids %>%
  left_join(unique_990 %>% select(ID, source), by = "ID") %>%
  left_join(aha.variants %>% distinct(ID, state), by = "ID") %>%
  group_by(state) %>%
  summarize(
    n_hospitals = n(),
    n_matched = sum(!is.na(source)),
    match_rate = n_matched / n_hospitals,
    n_manual = sum(source == "manual", na.rm = TRUE),
    n_fuzzy = sum(source == "fuzzy", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match_rate)

# =============================================================================
# Save Outputs
# =============================================================================

write_csv(unique_990, "data/output/unique_990.csv")
write_csv(diagnostics, "data/output/fuzzy/qa_990_diagnostics.csv")
write_csv(qa_summary, "data/output/fuzzy/qa_990_summary.csv")
write_csv(qa_by_state, "data/output/fuzzy/qa_990_by_state.csv")

