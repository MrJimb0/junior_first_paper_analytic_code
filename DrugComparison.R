library(tidyverse)
library(readxl)
library(writexl)
library(boot)

source("config.R")

set.seed(RANDOM_SEED)

df <- read_xlsx(FINAL_DRUG_DF)

# Restrict to patients present in LLM data (handles partial runs, e.g. 20 of 100)
llm_patients <- df %>% filter(version %in% ALL_LLMS) %>% pull(gsID) %>% unique()
df <- df %>% filter(gsID %in% llm_patients)

# CRC censor date: remove lines started after each source's review cutoff
# CRC rows have censor_date = CRC_CENSOR_DATE (2022-06-01); v1/v2/LLMs have
# censor_date = patient last_note_date. This removes CRC lines after their
# review cutoff and any lines beyond the patient's last chart note.
df <- df %>%
  mutate(Start       = as.Date(Start),
         censor_date = as.Date(censor_date)) %>%
  filter(is.na(Start) | is.na(censor_date) | Start <= censor_date)

# CRC-censored v1 reference
# For fair CRC comparisons all v1 lines are also restricted to CRC_CENSOR_DATE;
# CRC abstracted only through that date so post-censor v1 lines were never seen.
df_v1_crc <- df %>% filter(version == "v1", is.na(Start) | Start <= CRC_CENSOR_DATE)

# Collapse only the generic "AI" / "aromatase inhibitor" shorthand; specific drug
# names are kept as written, so anastrozole != exemestane != letrozole.
normalize_drug_llm <- function(x) {
  ai_re <- regex("^(ai|aromatase inhibitor)$", ignore_case = TRUE)
  if_else(str_detect(x, ai_re), "aromatase_inhibitor", str_trim(tolower(x)))
}

# Helper: bijective line matching within a date window
# test_df must contain gsID + Start (plus any extra columns to carry through).
# ref_df  must contain gsID + Start_v1 (plus any extra columns to carry through).
# Returns the matched pairs with all columns from both sides preserved.
match_lines_window <- function(test_df, ref_df, window_days = 90) {
  test_df %>%
    inner_join(ref_df, by = "gsID", relationship = "many-to-many") %>%
    filter(abs(as.numeric(difftime(Start, Start_v1, units = "days"))) <= window_days) %>%
    mutate(day_diff = abs(as.numeric(difftime(Start, Start_v1, units = "days")))) %>%
    # Each test line → closest ref line
    group_by(gsID, Start) %>%
    slice_min(day_diff, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # Each ref line matched at most once (bijective matching)
    group_by(gsID, Start_v1) %>%
    slice_min(day_diff, n = 1, with_ties = FALSE) %>%
    ungroup()
}

#First, compare the bag of drugs for each patient
bag_of_drugs_by_patient_version <- df %>%
  separate_rows(regimen_bag_of_drugs, sep = ",") %>%
  # Every source keeps its specific drug names; only the generic "AI"/"aromatase
  # inhibitor" shorthand is collapsed. Class-level (AI + taxane) collapsing is a
  # separate sensitivity analysis.
  mutate(drug = normalize_drug_llm(str_trim(regimen_bag_of_drugs))) %>%
  filter(!is.na(drug) & drug != "") %>%
  distinct(gsID, version, drug) %>%
  group_by(gsID, version) %>%
  summarise(
    drugs = list(sort(unique(drug))),
    n_drugs = length(drugs[[1]]),
    .groups = "drop"
  )

#This gives us a wide df where each row is an ID and a vector bag of drugs in each col 
wide_bag_of_drugs_by_patient_version <- bag_of_drugs_by_patient_version %>%
  select(gsID, version, drugs) %>%
  pivot_wider(names_from = version, values_from = drugs)

#First get the human variability
wide_bag_of_drugs_by_patient_version <-
  wide_bag_of_drugs_by_patient_version %>%
  mutate(
    jaccard_v1_v2  = map2_dbl(v1, v2, jaccard_set),
    jaccard_v2_vs_v1 = jaccard_v1_v2
  )

x <- wide_bag_of_drugs_by_patient_version$jaccard_v1_v2

human_results <- tibble(
  mean = mean(x, na.rm = TRUE),
  median = median(x, na.rm = TRUE),
  IQR = IQR(x, na.rm = TRUE),
  lower_CI = boot_lo(x),
  upper_CI = boot_hi(x)
)

# Now we want to compare the LLMs to each other and to the human variability. Do v1 then v2
llms <- c(ALL_LLMS, "crc")

# LLM vs v1
for (m in llms) {
  wide_bag_of_drugs_by_patient_version[[paste0("jaccard_", m, "_vs_v1")]] <-
    map2_dbl(
      wide_bag_of_drugs_by_patient_version[[m]],
      wide_bag_of_drugs_by_patient_version[["v1"]],
      jaccard_set
    )
}

# LLM vs v2
for (m in llms) {
  wide_bag_of_drugs_by_patient_version[[paste0("jaccard_", m, "_vs_v2")]] <-
    map2_dbl(
      wide_bag_of_drugs_by_patient_version[[m]],
      wide_bag_of_drugs_by_patient_version[["v2"]],
      jaccard_set
    )
}

# Fair CRC bag Jaccard: rebuild v1 drug sets censored at CRC_CENSOR_DATE
v1_bags_crc <- df_v1_crc %>%
  separate_rows(regimen_bag_of_drugs, sep = ",") %>%
  mutate(drug = normalize_drug_llm(str_trim(regimen_bag_of_drugs))) %>%
  filter(!is.na(drug), drug != "") %>%
  distinct(gsID, drug) %>%
  group_by(gsID) %>%
  summarise(drugs_v1_crc = list(sort(unique(drug))), .groups = "drop")

wide_bag_of_drugs_by_patient_version <- wide_bag_of_drugs_by_patient_version %>%
  left_join(v1_bags_crc, by = "gsID") %>%
  mutate(jaccard_crc_vs_v1 = map2_dbl(crc, drugs_v1_crc, jaccard_set)) %>%
  select(-drugs_v1_crc)

#Summarize vs v1
llm_vs_v1_results <- map_dfr(llms, function(m) {

  x <- wide_bag_of_drugs_by_patient_version[[paste0("jaccard_", m, "_vs_v1")]]
  ci <- boot_ci_95(x)

  tibble(
    model = m,
    reference = "v1",
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    IQR = IQR(x, na.rm = TRUE),
    lower_CI = ci[1],
    upper_CI = ci[2]
  )
})

#Summarize vs v2
llm_vs_v2_results <- map_dfr(llms, function(m) {

  x <- wide_bag_of_drugs_by_patient_version[[paste0("jaccard_", m, "_vs_v2")]]
  ci <- boot_ci_95(x)

  tibble(
    model = m,
    reference = "v2",
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    IQR = IQR(x, na.rm = TRUE),
    lower_CI = ci[1],
    upper_CI = ci[2]
  )
})

human_summary <- human_results %>%
  mutate(
    model = "v1_vs_v2",
    reference = "human_variability"
  ) %>%
  select(model, reference, everything())

llm_threshold_summary <- map_dfr(llms, function(m) {
  
  x <- wide_bag_of_drugs_by_patient_version[[paste0("jaccard_", m, "_vs_v1")]]

  tibble(
    model = m,
    pct_exact = mean(x == 1, na.rm = TRUE),
    pct_ge_0_9 = mean(x >= 0.9, na.rm = TRUE)
  )
})

# Human benchmark row
human_row <- tibble(
  model = "Second Oncologist (v2)",
  mean = human_results$mean,
  lower_CI = human_results$lower_CI,
  upper_CI = human_results$upper_CI,
  median = human_results$median,
  IQR = human_results$IQR,
  mean_difference_from_human = 0,
  diff_lower_CI = NA_real_,
  diff_upper_CI = NA_real_
)

# LLM rows (vs v1) — consistent with lines-of-therapy table
llm_rows <- map_dfr(llms, function(m) {

  x <- wide_bag_of_drugs_by_patient_version[[paste0("jaccard_", m, "_vs_v1")]]
  ci <- boot_ci_95(x)

  diff_vec <- wide_bag_of_drugs_by_patient_version$jaccard_v1_v2 - x
  diff_ci <- boot_ci_95(diff_vec)

  tibble(
    model = m,
    mean = mean(x, na.rm = TRUE),
    lower_CI = ci[1],
    upper_CI = ci[2],
    median = median(x, na.rm = TRUE),
    IQR = IQR(x, na.rm = TRUE),
    mean_difference_from_human = mean(diff_vec, na.rm = TRUE),
    diff_lower_CI = diff_ci[1],
    diff_upper_CI = diff_ci[2]
  )
})

final_table <- bind_rows(human_row, llm_rows)

#Lines of Therapy
checklist <- as.data.frame(table(df$regimen_standardized))

writexl::write_xlsx(
  checklist,
  file.path(LAB_MTG_DIR, "lines_of_therapy_checklist.xlsx")
)

# Lines comparison based on regimen_standardized

# 1) Build regimen_standardized set per gsID × version (full regimen string as token)
# distinct() collapses repeated regimens and drops chronology, so this is a
# unique-regimen-set Jaccard. The order- and repeat-aware line metrics reported in
# Table 3 come from line_metrics.R.
lines_by_patient_version <- df %>%
  mutate(line_sig = str_squish(regimen_standardized)) %>%
  filter(!is.na(line_sig) & line_sig != "") %>%
  distinct(gsID, version, line_sig) %>%
  group_by(gsID, version) %>%
  summarise(
    lines = list(sort(unique(line_sig))),
    n_lines = length(lines[[1]]),
    .groups = "drop"
  )

# 2) Wide df (one row per gsID; list-columns per model/version)
wide_lines_by_patient_version <- lines_by_patient_version %>%
  select(gsID, version, lines) %>%
  pivot_wider(names_from = version, values_from = lines)

# 3) Human variability (v1 vs v2)
wide_lines_by_patient_version <- wide_lines_by_patient_version %>%
  mutate(jaccard_v1_v2 = map2_dbl(v1, v2, jaccard_set))

x <- wide_lines_by_patient_version$jaccard_v1_v2

human_results_lines <- tibble(
  mean = mean(x, na.rm = TRUE),
  median = median(x, na.rm = TRUE),
  IQR = IQR(x, na.rm = TRUE),
  lower_CI = boot_lo(x),
  upper_CI = boot_hi(x)
)

# LLM vs v1
for (m in llms) {
  wide_lines_by_patient_version[[paste0("jaccard_", m, "_vs_v1")]] <-
    map2_dbl(
      wide_lines_by_patient_version[[m]],
      wide_lines_by_patient_version[["v1"]],
      jaccard_set
    )
}

# LLM vs v2
for (m in llms) {
  wide_lines_by_patient_version[[paste0("jaccard_", m, "_vs_v2")]] <-
    map2_dbl(
      wide_lines_by_patient_version[[m]],
      wide_lines_by_patient_version[["v2"]],
      jaccard_set
    )
}

# Fair CRC lines Jaccard: rebuild v1 line sets censored at CRC_CENSOR_DATE
v1_lines_crc_wide <- df_v1_crc %>%
  mutate(line_sig = str_squish(regimen_standardized)) %>%
  filter(!is.na(line_sig), line_sig != "") %>%
  distinct(gsID, line_sig) %>%
  group_by(gsID) %>%
  summarise(lines_v1_crc = list(sort(unique(line_sig))), .groups = "drop")

wide_lines_by_patient_version <- wide_lines_by_patient_version %>%
  left_join(v1_lines_crc_wide, by = "gsID") %>%
  mutate(jaccard_crc_vs_v1 = map2_dbl(crc, lines_v1_crc, jaccard_set)) %>%
  select(-lines_v1_crc)

# Summary vs v2 + difference-from-human table
human_row_lines <- tibble(
  model = "Second Oncologist (v2)",
  mean = human_results_lines$mean,
  lower_CI = human_results_lines$lower_CI,
  upper_CI = human_results_lines$upper_CI,
  median = human_results_lines$median,
  IQR = human_results_lines$IQR,
  mean_difference_from_human = 0,
  diff_lower_CI = NA_real_,
  diff_upper_CI = NA_real_
)

# LLM rows (vs v1)
llm_rows_lines <- map_dfr(llms, function(m) {

  x <- wide_lines_by_patient_version[[paste0("jaccard_", m, "_vs_v1")]]
  ci <- boot_ci_95(x)

  diff_vec <- wide_lines_by_patient_version$jaccard_v1_v2 - x
  diff_ci <- boot_ci_95(diff_vec)

  tibble(
    model = m,
    mean = mean(x, na.rm = TRUE),
    lower_CI = ci[1],
    upper_CI = ci[2],
    median = median(x, na.rm = TRUE),
    IQR = IQR(x, na.rm = TRUE),
    mean_difference_from_human = mean(diff_vec, na.rm = TRUE),
    diff_lower_CI = diff_ci[1],
    diff_upper_CI = diff_ci[2]
  )
})

final_table_lines <- bind_rows(human_row_lines, llm_rows_lines)

# Inline stats for manuscript text
cat('\n=== Systemic therapy inline stats (for manuscript text) ===\n')

# n matched lines between v1 and v2 (second oncologist); manuscript text: "n = 651"
if (exists("intent_summary")) {
  n_matched_v2 <- intent_summary %>% filter(version == "v2") %>% pull(n_matched)
  cat(sprintf('  Matched lines (v1 vs v2 = second oncologist): n = %d\n',
              if (length(n_matched_v2) > 0) n_matched_v2[1] else NA_integer_))
}

# Inter-oncologist Jaccard
cat(sprintf('  Inter-oncologist drug Jaccard:  mean = %.2f (CI: %.2f-%.2f)\n',
            human_results$mean, human_results$lower_CI, human_results$upper_CI))
cat(sprintf('  Inter-oncologist lines Jaccard: mean = %.2f (CI: %.2f-%.2f)\n',
            human_results_lines$mean, human_results_lines$lower_CI, human_results_lines$upper_CI))

# LLM Jaccard ranked: drug
non_human_drugs <- final_table %>%
  filter(model != "Second Oncologist (v2)") %>%
  arrange(desc(mean))
cat('\n  Drug Jaccard by LLM (ranked):\n')
for (i in seq_len(nrow(non_human_drugs))) {
  r <- non_human_drugs[i, ]
  cat(sprintf('    %-35s mean=%.2f (%.2f-%.2f)  diff_from_oncologist=%.2f (%.1f%%)\n',
              r$model, r$mean, r$lower_CI, r$upper_CI,
              r$mean_difference_from_human, r$mean_difference_from_human * 100))
}

# LLM Jaccard ranked: lines
non_human_lines <- final_table_lines %>%
  filter(model != "Second Oncologist (v2)") %>%
  arrange(desc(mean))
cat('\n  Lines Jaccard by LLM (ranked):\n')
for (i in seq_len(nrow(non_human_lines))) {
  r <- non_human_lines[i, ]
  cat(sprintf('    %-35s mean=%.2f (%.2f-%.2f)  diff_from_oncologist=%.2f (%.1f%%)\n',
              r$model, r$mean, r$lower_CI, r$upper_CI,
              r$mean_difference_from_human, r$mean_difference_from_human * 100))
}

# Best vs 2nd-best LLM Jaccard difference (manuscript: "0.03 for drugs, 0.03 for lines")
if (nrow(non_human_drugs) >= 2 && nrow(non_human_lines) >= 2) {
  cat(sprintf('\n  Best vs 2nd-best LLM drug Jaccard diff:  %.2f\n',
              non_human_drugs$mean[1] - non_human_drugs$mean[2]))
  cat(sprintf('  Best vs 2nd-best LLM lines Jaccard diff: %.2f\n',
              non_human_lines$mean[1] - non_human_lines$mean[2]))
}

# Save wide per-patient tables so downstream analyses can be rebuilt without
# re-running the full DrugComparison pipeline
writexl::write_xlsx(
  wide_bag_of_drugs_by_patient_version %>%
    select(gsID, starts_with("jaccard_")),
  file.path(INT_DIR, "wide_bag_jaccard_per_patient.xlsx")
)

writexl::write_xlsx(
  wide_lines_by_patient_version %>%
    select(gsID, starts_with("jaccard_")),
  file.path(INT_DIR, "wide_lines_jaccard_per_patient.xlsx")
)

# Additional Analyses: Lines Count, Treatment Intent, ReasonStop
# v1 is ground truth for all comparisons
# Lines matched across versions by Start date within ± 90 days

non_v1 <- setdiff(unique(df$version), "v1")

# 1. Number of Lines Called

lines_count <- df %>%
  filter(!is.na(regimen_standardized) & regimen_standardized != "") %>%
  group_by(gsID, version) %>%
  summarise(
    # Count actual therapy entries (rows), not annotated line numbers.
    # Tx_LineN is an oncologist convention that can lump multiple sequential
    # agents into one line number — row count is consistent across annotators.
    n_lines = n(),
    .groups = "drop"
  )

lines_count_wide <- lines_count %>%
  pivot_wider(names_from = version, values_from = n_lines)

lines_count_summary <- map_dfr(non_v1, function(v) {
  # For CRC, compare against v1 lines censored at CRC_CENSOR_DATE
  if (v == "crc") {
    v1_crc_counts <- df_v1_crc %>%
      filter(!is.na(regimen_standardized) & regimen_standardized != "") %>%
      count(gsID, name = "n_v1_crc")
    merged <- lines_count_wide %>%
      select(gsID, n_crc = crc) %>%
      left_join(v1_crc_counts, by = "gsID") %>%
      mutate(n_v1_crc = replace_na(n_v1_crc, 0L))
    diff_vec <- merged$n_crc - merged$n_v1_crc
  } else {
    diff_vec <- lines_count_wide[[v]] - lines_count_wide[["v1"]]
  }
  tibble(
    version         = v,
    mean_diff       = mean(diff_vec, na.rm = TRUE),
    median_diff     = median(diff_vec, na.rm = TRUE),
    pct_exact_match = mean(diff_vec == 0, na.rm = TRUE),
    pct_overcall    = mean(diff_vec > 0, na.rm = TRUE),
    pct_undercall   = mean(diff_vec < 0, na.rm = TRUE)
  )
})

# 2. Treatment Intent (Curative / Palliative)
# Match each version's lines to v1 lines by gsID + Start within ± 90 days

v1_intent <- df %>%
  filter(version == "v1") %>%
  select(gsID, Start_v1 = Start, TreatmentIntent_v1 = TreatmentIntent)

v1_intent_crc <- df_v1_crc %>%
  select(gsID, Start_v1 = Start, TreatmentIntent_v1 = TreatmentIntent)

intent_summary <- map_dfr(non_v1, function(v) {

  ref  <- if (v == "crc") v1_intent_crc else v1_intent
  test <- df %>%
    filter(version == v) %>%
    select(gsID, Start, TreatmentIntent)

  joined <- match_lines_window(test, ref)

  cm <- table(
    predicted = joined$TreatmentIntent,
    actual    = joined$TreatmentIntent_v1
  )
  tp <- if ("Palliative" %in% rownames(cm) && "Palliative" %in% colnames(cm)) cm["Palliative", "Palliative"] else 0L
  fp <- if ("Palliative" %in% rownames(cm) && "Curative"   %in% colnames(cm)) cm["Palliative", "Curative"]   else 0L
  fn <- if ("Curative"   %in% rownames(cm) && "Palliative" %in% colnames(cm)) cm["Curative",   "Palliative"] else 0L
  tn <- if ("Curative"   %in% rownames(cm) && "Curative"   %in% colnames(cm)) cm["Curative",   "Curative"]   else 0L
  n_matched_valid <- tp + fp + fn + tn
  stopifnot(n_matched_valid == sum(!is.na(joined$TreatmentIntent) & !is.na(joined$TreatmentIntent_v1)))

  tibble(
    version     = v,
    n_matched   = nrow(joined),
    accuracy    = mean(joined$TreatmentIntent == joined$TreatmentIntent_v1, na.rm = TRUE),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    ppv         = tp / (tp + fp),
    npv         = tn / (tn + fn)
  )
})

cat(sprintf('\n  Matched lines (v1 vs v2): n = %d\n',
            intent_summary %>% filter(version == "v2") %>% pull(n_matched)))

# 3. ReasonStop Correctness (Progression / Other)
# Same ± 90-day matching on Start date

v1_reason <- df %>%
  filter(version == "v1") %>%
  select(gsID, Start_v1 = Start, ReasonStop_v1 = ReasonStop)

v1_reason_crc <- df_v1_crc %>%
  select(gsID, Start_v1 = Start, ReasonStop_v1 = ReasonStop)

reason_summary <- map_dfr(non_v1, function(v) {

  ref  <- if (v == "crc") v1_reason_crc else v1_reason
  test <- df %>%
    filter(version == v) %>%
    select(gsID, Start, ReasonStop)

  joined <- match_lines_window(test, ref)

  cm <- table(
    predicted = joined$ReasonStop,
    actual    = joined$ReasonStop_v1
  )
  tp <- if ("Progression" %in% rownames(cm) && "Progression" %in% colnames(cm)) cm["Progression", "Progression"] else 0L
  fp <- if ("Progression" %in% rownames(cm) && "Other"       %in% colnames(cm)) cm["Progression", "Other"]       else 0L
  fn <- if ("Other"       %in% rownames(cm) && "Progression" %in% colnames(cm)) cm["Other",       "Progression"] else 0L
  tn <- if ("Other"       %in% rownames(cm) && "Other"       %in% colnames(cm)) cm["Other",       "Other"]       else 0L
  n_matched_valid <- tp + fp + fn + tn
  stopifnot(n_matched_valid == sum(!is.na(joined$ReasonStop) & !is.na(joined$ReasonStop_v1)))

  tibble(
    version     = v,
    n_matched   = nrow(joined),
    accuracy    = mean(joined$ReasonStop == joined$ReasonStop_v1, na.rm = TRUE),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    ppv         = tp / (tp + fp),
    npv         = tn / (tn + fn)
  )
})

# Sensitivity / Error Quantification
# 1. Extra and missing drugs per patient
# 2. Per-patient complete agreement on TreatmentIntent (overall, curative, palliative)
# 3. Line detection F1

# 1. Extra and Missing Drugs Per Patient
# Extra   = drugs called by version not in v1 (additional predicted; set
#           differences, not adjudicated hallucinations)
# Missing = drugs in v1 not called by version (false negatives / omissions)

v1_drugs <- bag_of_drugs_by_patient_version %>%
  filter(version == "v1") %>%
  select(gsID, drugs_v1 = drugs)

# CRC-censored v1 drug reference (bag_of_drugs_by_patient_version only has all-v1;
# rebuild from raw df filtered to CRC_CENSOR_DATE)
v1_drugs_crc <- df_v1_crc %>%
  separate_rows(regimen_bag_of_drugs, sep = ",") %>%
  mutate(drug = normalize_drug_llm(str_trim(regimen_bag_of_drugs))) %>%
  filter(!is.na(drug), drug != "") %>%
  distinct(gsID, drug) %>%
  group_by(gsID) %>%
  summarise(drugs_v1 = list(sort(unique(drug))), .groups = "drop")

drug_error_summary <- map_dfr(non_v1, function(v) {

  ref <- if (v == "crc") v1_drugs_crc else v1_drugs
  joined <- bag_of_drugs_by_patient_version %>%
    filter(version == v) %>%
    select(gsID, drugs) %>%
    inner_join(ref, by = "gsID") %>%
    mutate(
      n_extra   = map2_int(drugs, drugs_v1, ~ length(setdiff(.x, .y))),
      n_missing = map2_int(drugs, drugs_v1, ~ length(setdiff(.y, .x))),
      n_tp      = map2_int(drugs, drugs_v1, ~ length(intersect(.x, .y)))
    )

  tibble(
    version        = v,
    n_pts          = nrow(joined),
    total_drugs    = sum(map_int(joined$drugs, length)),
    total_tp       = sum(joined$n_tp),
    total_missing  = sum(joined$n_missing),
    total_extra    = sum(joined$n_extra),
    mean_extra     = mean(joined$n_extra,   na.rm = TRUE),
    sd_extra       = sd(joined$n_extra,     na.rm = TRUE),
    median_extra   = median(joined$n_extra, na.rm = TRUE),
    mean_missing   = mean(joined$n_missing,   na.rm = TRUE),
    sd_missing     = sd(joined$n_missing,     na.rm = TRUE),
    median_missing = median(joined$n_missing, na.rm = TRUE),
    pct_no_extra   = mean(joined$n_extra   == 0, na.rm = TRUE),
    pct_no_missing = mean(joined$n_missing == 0, na.rm = TRUE),
    n_exact        = sum(joined$n_extra == 0 & joined$n_missing == 0)
  )
})

# Most commonly extra / missing drugs per version
all_drug_errors <- map_dfr(non_v1, function(v) {

  ref <- if (v == "crc") v1_drugs_crc else v1_drugs
  joined <- bag_of_drugs_by_patient_version %>%
    filter(version == v) %>%
    select(gsID, drugs) %>%
    inner_join(ref, by = "gsID") %>%
    mutate(
      extra   = map2(drugs, drugs_v1, ~ setdiff(.x, .y)),
      missing = map2(drugs, drugs_v1, ~ setdiff(.y, .x))
    )

  bind_rows(
    joined %>% unnest(extra)   %>% transmute(version = v, type = "extra",   drug = extra),
    joined %>% unnest(missing) %>% transmute(version = v, type = "missing", drug = missing)
  )
})

top_extra_drugs <- all_drug_errors %>%
  filter(type == "extra") %>%
  count(version, drug, sort = TRUE) %>%
  group_by(version) %>%
  slice_max(n, n = 10) %>%
  ungroup()

top_missing_drugs <- all_drug_errors %>%
  filter(type == "missing") %>%
  count(version, drug, sort = TRUE) %>%
  group_by(version) %>%
  slice_max(n, n = 10) %>%
  ungroup()

# Wide-format missed drug table: rows = drugs, cols = versions
# Includes any drug that appears in the top 10 for at least one version
missed_drugs_wide <- all_drug_errors %>%
  filter(type == "missing") %>%
  count(version, drug) %>%
  semi_join(top_missing_drugs, by = c("version", "drug")) %>%
  pivot_wider(names_from = version, values_from = n, values_fill = 0) %>%
  mutate(total_missed = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(total_missed))

writexl::write_xlsx(
  list(
    by_version = top_missing_drugs,
    wide       = missed_drugs_wide
  ),
  file.path(LAB_MTG_DIR, "missed_drugs_summary.xlsx")
)

# 2. Line Detection F1
# A line is a TP if a v1 line exists within ± 90 days (nearest match, one-to-one)
# Precision = TP / lines called by version
# Recall    = TP / lines in v1
# F1        = harmonic mean of precision and recall
#
# Start dates only: regimen composition is not required here, and
# match_lines_window() groups on (gsID, Start), so distinct lines sharing a start
# date can collapse. Table 3 uses the chronological matcher in line_metrics.R.

v1_lines <- df %>%
  filter(version == "v1") %>%
  select(gsID, Start_v1 = Start)

v1_lines_crc <- df_v1_crc %>%
  select(gsID, Start_v1 = Start)

line_f1_summary <- map_dfr(non_v1, function(v) {

  ref  <- if (v == "crc") v1_lines_crc else v1_lines
  test <- df %>%
    filter(version == v) %>%
    select(gsID, Start)

  matched <- match_lines_window(test, ref)

  tp        <- nrow(matched)
  n_called  <- nrow(test)
  # Recall denominator is restricted to patients present in both test and v1 versions.
  # Patients with complete extraction failure are excluded. See study limitations.
  n_v1      <- nrow(ref %>% filter(gsID %in% unique(test$gsID)))
  precision <- tp / n_called
  recall    <- tp / n_v1
  f1        <- 2 * precision * recall / (precision + recall)

  tibble(
    version   = v,
    n_called  = n_called,
    n_v1      = n_v1,
    tp        = tp,
    precision = precision,
    recall    = recall,
    f1        = f1
  )
})

all_versions <- unique(df$version)
all_versions <- all_versions[all_versions != "v1"]

# Per-Patient Accuracy Export
# One row per gsID × version with all accuracy metrics + patient context
# (oncologist, second opinion flag, chart complexity)
# Exported to Analysis_Output/final_dfs/

# Load patient context from intermediate files

# Per-version context (LLMs): oncologist, visit count, dates, genetics
llm_context <- read_xlsx(LLM_NON_DRUG_CLEANED) %>%
  select(gsID, version,
         primary_oncologist_shc,
         date_initial_diagnosis, date_first_metastatic_diagnosis, birth_date,
         stage_at_diagnosis, hr_her2_status,
         germline_testing_done, brca_status,
         somatic_testing_done, pik3ca_status, esr1_status)

# CRC: only dates available from this intermediate
crc_context <- read_xlsx(CRC_NON_DRUG_CLEANED) %>%
  mutate(version = "crc") %>%
  select(gsID, version, date_first_metastatic_diagnosis)

# Patient-level GT context (not version-specific — join by gsID only)
# Includes chart span proxy (last_note_date), mortality, treating site, oncologist
gt_patient_context <- read_xlsx(GT_NON_DRUG_CLEANED) %>%
  select(gsID,
         gt_oncologist          = primary_oncologist_shc,
         treated_at_shc,
         last_note_date,
         deceased_clinician,
         death_date_clinician,
         stage_at_diagnosis,
         hr_her2_status,
         date_initial_diagnosis,
         date_first_metastatic_diagnosis)

# Version-level context (LLM + CRC rows); GT patient context joined separately below
so_context <- bind_rows(llm_context, crc_context)

# Per-patient Jaccard bag vs v1 (already computed)
jaccard_long <- wide_bag_of_drugs_by_patient_version %>%
  select(gsID, all_of(paste0("jaccard_", all_versions, "_vs_v1"))) %>%
  pivot_longer(
    cols = -gsID,
    names_to  = "version",
    values_to = "jaccard_bag_vs_v1"
  ) %>%
  mutate(version = str_remove(version, "^jaccard_") %>% str_remove("_vs_v1$"))

# Per-patient intent agreement (all lines correct)
intent_per_patient <- map_dfr(non_v1, function(v) {
  ref <- if (v == "crc") v1_intent_crc else v1_intent
  df %>%
    filter(version == v) %>%
    select(gsID, Start, TreatmentIntent) %>%
    match_lines_window(ref) %>%
    mutate(intent_correct = TreatmentIntent == TreatmentIntent_v1) %>%
    group_by(gsID) %>%
    summarise(
      intent_lines_matched  = n(),
      intent_pct_correct    = mean(intent_correct, na.rm = TRUE),
      intent_all_correct    = all(intent_correct),
      .groups = "drop"
    ) %>%
    mutate(version = v)
})

# Per-patient ReasonStop agreement
reason_per_patient <- map_dfr(non_v1, function(v) {
  ref <- if (v == "crc") v1_reason_crc else v1_reason
  df %>%
    filter(version == v) %>%
    select(gsID, Start, ReasonStop) %>%
    match_lines_window(ref) %>%
    mutate(reason_correct = ReasonStop == ReasonStop_v1) %>%
    group_by(gsID) %>%
    summarise(
      reason_lines_matched = n(),
      reason_pct_correct   = mean(reason_correct, na.rm = TRUE),
      reason_all_correct   = all(reason_correct),
      .groups = "drop"
    ) %>%
    mutate(version = v)
})

# Per-patient line counts
# Tertile breakpoints for GT v1 line count complexity grouping
tertile_cuts <- quantile(lines_count_wide$v1, probs = c(1/3, 2/3), na.rm = TRUE)

line_counts_long <- lines_count_wide %>%
  select(gsID, n_lines_v1 = v1) %>%
  left_join(
    lines_count_wide %>%
      pivot_longer(cols = all_of(all_versions), names_to = "version", values_to = "n_lines_called"),
    by = "gsID"
  ) %>%
  mutate(
    line_count_diff = n_lines_called - n_lines_v1,
    complexity      = case_when(
      n_lines_v1 <= tertile_cuts[1] ~ "Low complexity (bottom third)",
      n_lines_v1 <= tertile_cuts[2] ~ "Moderate complexity (middle third)",
      n_lines_v1 >  tertile_cuts[2] ~ "High complexity (top third)",
      TRUE                           ~ NA_character_
    )
  )

# Assemble final per-patient table
patient_accuracy_export <- jaccard_long %>%
  left_join(line_counts_long,      by = c("gsID", "version")) %>%
  left_join(intent_per_patient,    by = c("gsID", "version")) %>%
  left_join(reason_per_patient,    by = c("gsID", "version")) %>%
  left_join(so_context,            by = c("gsID", "version")) %>%
  left_join(gt_patient_context,    by = "gsID") %>%      # patient-level, not version-specific
  mutate(
    # Chart span in years as proxy for chart size when visits_total is unavailable
    chart_span_years = as.numeric(difftime(last_note_date, date_initial_diagnosis.y,
                                           units = "days")) / 365.25
  ) %>%
  select(
    gsID,
    version,
    # Accuracy metrics
    jaccard_bag_vs_v1,
    intent_lines_matched,
    intent_pct_correct,
    intent_all_correct,
    reason_lines_matched,
    reason_pct_correct,
    reason_all_correct,
    # Line counts
    complexity,
    n_lines_v1,
    n_lines_called,
    line_count_diff,
    # Chart size proxies
    chart_span_years,
    last_note_date,
    # Patient context
    gt_oncologist,
    primary_oncologist_shc,
    treated_at_shc,
    stage_at_diagnosis = stage_at_diagnosis.y,
    hr_her2_status = hr_her2_status.y,
    date_initial_diagnosis = date_initial_diagnosis.y,
    date_first_metastatic_diagnosis = date_first_metastatic_diagnosis.y,
    birth_date,
    deceased_clinician,
    death_date_clinician,
    # Genetics
    germline_testing_done,
    brca_status,
    somatic_testing_done,
    pik3ca_status,
    esr1_status
  ) %>%
  arrange(version, gsID)

writexl::write_xlsx(
  patient_accuracy_export,
  file.path(FINAL_DIR, "patient_level_accuracy_by_version.xlsx")
)

print(paste0("Exported per-patient accuracy for ",
             n_distinct(patient_accuracy_export$gsID), " patients × ",
             n_distinct(patient_accuracy_export$version), " versions."))

# Appendix Tables S2–S4: Wide-format tables matching Word doc layout

# Column order and version→label map
appendix_cols <- c("Second oncologist", "GPT-5", "Gemini 2.5 Pro",
                   "GPT-4o", "DeepSeek-R1", "Research coordinator")
ver_col <- c(v2 = "Second oncologist", gpt5 = "GPT-5",
             gemini_2_5 = "Gemini 2.5 Pro", gpt4o = "GPT-4o",
             deepseek_r1 = "DeepSeek-R1", crc = "Research coordinator")
all_src_app <- names(ver_col)

# Bag-of-drugs: per-patient P/R/F1/Jaccard vs v1
bag_prf_app <- map_dfr(all_src_app, function(v) {
  d   <- bag_of_drugs_by_patient_version %>% filter(version == v)
  ref <- if (v == "crc") v1_drugs_crc else v1_drugs
  if (nrow(d) == 0) return(tibble(version = v, bag_jaccard = "\u2014",
    bag_precision = "\u2014", bag_recall = "\u2014", bag_f1 = "\u2014"))
  j <- d %>%
    inner_join(ref, by = "gsID") %>%
    mutate(
      tp   = map2_int(drugs, drugs_v1, ~ length(intersect(.x, .y))),
      fp   = map2_int(drugs, drugs_v1, ~ length(setdiff(.x, .y))),
      fn   = map2_int(drugs, drugs_v1, ~ length(setdiff(.y, .x))),
      prec = if_else(tp + fp > 0, tp / (tp + fp), NA_real_),
      rec  = if_else(tp + fn > 0, tp / (tp + fn), NA_real_),
      f1   = if_else(!is.na(prec) & !is.na(rec) & (prec + rec) > 0,
                     2 * prec * rec / (prec + rec), NA_real_),
      jacc = map2_dbl(drugs, drugs_v1, jaccard_set)
    )
  p_ci <- boot_ci(j$prec, mean); r_ci <- boot_ci(j$rec, mean)
  f_ci <- boot_ci(j$f1, mean)
  tibble(version = v,
    bag_jaccard   = { jc <- boot_ci(j$jacc, mean); fmt_ci(jc["est"], jc["lo"], jc["hi"]) },
    bag_precision = fmt_pct_ci(p_ci["est"], p_ci["lo"], p_ci["hi"]),
    bag_recall    = fmt_pct_ci(r_ci["est"], r_ci["lo"], r_ci["hi"]),
    bag_f1        = fmt_pct_ci(f_ci["est"], f_ci["lo"], f_ci["hi"]))
})

# Lines-of-therapy: Jaccard vs v1
v1_line_sets_app <- lines_by_patient_version %>%
  filter(version == "v1") %>% select(gsID, lines_v1 = lines)

# CRC-censored v1 line sets for fair CRC comparison
v1_line_sets_crc <- df_v1_crc %>%
  mutate(line_sig = str_squish(regimen_standardized)) %>%
  filter(!is.na(line_sig), line_sig != "") %>%
  distinct(gsID, line_sig) %>%
  group_by(gsID) %>%
  summarise(lines_v1 = list(sort(unique(line_sig))), .groups = "drop")

line_jacc_app <- map_dfr(all_src_app, function(v) {
  d   <- lines_by_patient_version %>% filter(version == v)
  ref <- if (v == "crc") v1_line_sets_crc else v1_line_sets_app
  if (nrow(d) == 0) return(tibble(version = v, line_jaccard = "\u2014"))
  j_vals <- d %>%
    inner_join(ref, by = "gsID") %>%
    mutate(jacc = map2_dbl(lines, lines_v1, jaccard_set)) %>%
    pull(jacc)
  jc <- boot_ci(j_vals, mean)
  tibble(version = v, line_jaccard = fmt_ci(jc["est"], jc["lo"], jc["hi"]))
})

# Line detection P/R/F1: per-patient bootstrap
v1_ls_app     <- df      %>% filter(version == "v1") %>% select(gsID, Start_v1 = Start)
v1_ls_app_crc <- df_v1_crc %>% select(gsID, Start_v1 = Start)

line_f1_pp <- map_dfr(all_src_app, function(v) {
  test   <- df %>% filter(version == v) %>% select(gsID, Start)
  v1_ref <- if (v == "crc") v1_ls_app_crc else v1_ls_app
  if (nrow(test) == 0) return(NULL)
  matched_pp <- match_lines_window(test, v1_ref)
  pats <- unique(test$gsID)
  map_dfr(pats, function(pid) {
    nc <- sum(test$gsID == pid)
    nv <- sum(v1_ref$gsID == pid)
    nt <- sum(matched_pp$gsID == pid)
    pr <- if (nc > 0) nt / nc else NA_real_
    re <- if (nv > 0) nt / nv else NA_real_
    # Zero-TP evaluable patients score F1 = 0 rather than NA, which would drop them
    # from the mean; only a patient with no predicted and no reference lines is
    # non-evaluable.
    f1 <- if (nc == 0 && nv == 0) NA_real_
          else if (is.na(pr) || is.na(re) || (pr + re) == 0) 0
          else 2 * pr * re / (pr + re)
    tibble(version = v, gsID = pid, prec = pr, rec = re, f1 = f1)
  })
})

line_f1_ci_app <- map_dfr(all_src_app, function(v) {
  d <- line_f1_pp %>% filter(version == v)
  if (nrow(d) == 0) return(tibble(version = v,
    line_precision = "\u2014", line_recall = "\u2014", line_f1 = "\u2014"))
  p_ci <- boot_ci(d$prec, mean); r_ci <- boot_ci(d$rec, mean)
  f_ci <- boot_ci(d$f1, mean)
  tibble(version = v,
    line_precision = fmt_pct_ci(p_ci["est"], p_ci["lo"], p_ci["hi"]),
    line_recall    = fmt_pct_ci(r_ci["est"], r_ci["lo"], r_ci["hi"]),
    line_f1        = fmt_pct_ci(f_ci["est"], f_ci["lo"], f_ci["hi"]))
})

# Line count exact match
# CRC-censored v1 line counts (for fair CRC comparison)
line_count_app <- map_dfr(all_src_app, function(v) {
  nv <- lines_count_wide[[v]]
  if (is.null(nv)) return(tibble(version = v, line_count_total = "\u2014"))
  tibble(version = v,
    line_count_total = as.character(sum(nv, na.rm = TRUE)))
})

# Start/stop date ±30 d for matched lines
v1_dates_app <- df %>% filter(version == "v1") %>%
  select(gsID, Start_v1 = Start, Stop_v1 = Stop)

v1_dates_app_crc <- df_v1_crc %>%
  select(gsID, Start_v1 = Start, Stop_v1 = Stop)

date_agree_app <- map_dfr(all_src_app, function(v) {
  ref  <- if (v == "crc") v1_dates_app_crc else v1_dates_app
  test <- df %>% filter(version == v) %>% select(gsID, Start, Stop)
  if (nrow(test) == 0) return(tibble(version = v,
    start_30d = "\u2014", stop_30d = "\u2014"))
  m2 <- match_lines_window(test, ref)
  if (nrow(m2) == 0) return(tibble(version = v,
    start_30d = "\u2014", stop_30d = "\u2014"))
  sp <- m2 %>% group_by(gsID) %>%
    summarise(p30 = mean(day_diff <= 30), .groups = "drop")
  s_ci <- boot_ci(sp$p30, mean)
  stp <- m2 %>% filter(!is.na(Stop) & !is.na(Stop_v1)) %>%
    mutate(de = abs(as.numeric(difftime(Stop, Stop_v1, units = "days")))) %>%
    group_by(gsID) %>% summarise(p30 = mean(de <= 30), .groups = "drop")
  stop_fmt <- if (nrow(stp) > 0) {
    e_ci <- boot_ci(stp$p30, mean)
    fmt_pct_ci(e_ci["est"], e_ci["lo"], e_ci["hi"])
  } else "\u2014"
  tibble(version = v,
    start_30d = fmt_pct_ci(s_ci["est"], s_ci["lo"], s_ci["hi"]),
    stop_30d  = stop_fmt)
})

# Intent / reason stop agreement
intent_ci_app <- map_dfr(all_src_app, function(v) {
  d <- intent_per_patient %>% filter(version == v)
  if (nrow(d) == 0) return(tibble(version = v, intent_agree = "\u2014"))
  ci <- boot_ci(d$intent_pct_correct, mean)
  tibble(version = v, intent_agree = fmt_pct_ci(ci["est"], ci["lo"], ci["hi"]))
})

reason_ci_app <- map_dfr(all_src_app, function(v) {
  d <- reason_per_patient %>% filter(version == v)
  if (nrow(d) == 0) return(tibble(version = v, reason_agree = "\u2014"))
  ci <- boot_ci(d$reason_pct_correct, mean)
  tibble(version = v, reason_agree = fmt_pct_ci(ci["est"], ci["lo"], ci["hi"]))
})

# Number of matched lines per source
matched_lines_app2 <- map_dfr(all_src_app, function(v) {
  test   <- df %>% filter(version == v) %>% select(gsID, Start)
  v1_ref <- if (v == "crc") v1_ls_app_crc else v1_ls_app
  if (nrow(test) == 0) return(tibble(version = v, n_matched_lines = "0"))
  matched <- match_lines_window(test, v1_ref)
  tibble(version = v, n_matched_lines = as.character(nrow(matched)))
})

# Assemble drug rows for Appendix Table 2
drug_metrics_combined <- Reduce(function(a, b) full_join(a, b, by = "version"),
  list(bag_prf_app, line_jacc_app, line_f1_ci_app,
       line_count_app, matched_lines_app2, date_agree_app,
       intent_ci_app, reason_ci_app)) %>%
  mutate(col = coalesce(ver_col[version], version))

pivot_drug_row <- function(metric_col, row_label, metric_label) {
  drug_metrics_combined %>%
    select(col, value = all_of(metric_col)) %>%
    pivot_wider(names_from = col, values_from = value) %>%
    mutate(`Clinical domain / variable` = row_label,
           `Metric reported`            = metric_label) %>%
    {tbl <- .; for (cn in appendix_cols) if (!cn %in% names(tbl)) tbl[[cn]] <- "\u2014"; tbl} %>%
    select(`Clinical domain / variable`, `Metric reported`, all_of(appendix_cols))
}

tableS2_drug_rows <- bind_rows(
  pivot_drug_row("bag_jaccard",    "All prescribed anti-cancer drugs: Jaccard similarity", "Mean [95% CI]"),
  pivot_drug_row("bag_precision",  "All prescribed anti-cancer drugs: precision",          "Mean (95% CI)"),
  pivot_drug_row("bag_recall",     "All prescribed anti-cancer drugs: recall",             "Mean (95% CI)"),
  pivot_drug_row("bag_f1",         "All prescribed anti-cancer drugs: F1 score",           "Mean (95% CI)"),
  # ±90-day line-start and unique-regimen-set values from this script's matcher.
  # manuscript_tables.R replaces these rows with the line_metrics.R output for
  # Table S2 and Table 3.
  pivot_drug_row("line_jaccard",   "Systemic therapy lines: unique regimen-set Jaccard (legacy matcher)", "Mean [95% CI]"),
  pivot_drug_row("line_precision", "Systemic therapy lines: line-start precision (legacy matcher)",       "Mean (95% CI)"),
  pivot_drug_row("line_recall",    "Systemic therapy lines: line-start recall (legacy matcher)",          "Mean (95% CI)"),
  pivot_drug_row("line_f1",        "Systemic therapy lines: line-start F1 (legacy matcher)",              "Mean (95% CI)"),
  pivot_drug_row("line_count_total","Number of therapy lines per patient",                 "Total lines, No."),
  pivot_drug_row("n_matched_lines", "Number of matched lines",                             "No."),
  pivot_drug_row("start_30d",      "Start date of matched treatment line (\u00b130 d)",    "% agreement (95% CI)"),
  pivot_drug_row("stop_30d",       "Stop date of matched treatment line (\u00b130 d)",     "% agreement (95% CI)"),
  pivot_drug_row("intent_agree",   "Treatment intent for matched lines",                   "% agreement (95% CI)"),
  pivot_drug_row("reason_agree",   "Reason for discontinuation for matched lines",         "% agreement (95% CI)")
)

writexl::write_xlsx(tableS2_drug_rows, file.path(INT_DIR, "TableS2_drug_rows.xlsx"))
cat("  TableS2_drug_rows.xlsx written (", nrow(tableS2_drug_rows), "rows)\n")

# Appendix Table 3: Line-level error metrics (wide)
# For each source: total lines, matched (TP), missed (FN), extra (FP), mean/SD
line_error_app <- map_dfr(all_src_app, function(v) {
  test <- df %>% filter(version == v) %>% select(gsID, Start)
  if (nrow(test) == 0) return(tibble(version = v,
    total = NA_integer_, tp = NA_integer_, missed = NA_integer_, extra = NA_integer_,
    missed_mean = NA_real_, missed_sd = NA_real_,
    extra_mean  = NA_real_, extra_sd  = NA_real_,
    n_exact = NA_integer_, n_no_missed = NA_integer_, n_no_extra = NA_integer_,
    n_pts = NA_integer_))
  # Use CRC-censored v1 reference for CRC for fair comparison
  v1_s <- if (v == "crc") {
    df_v1_crc %>% filter(gsID %in% unique(test$gsID)) %>% select(gsID, Start_v1 = Start)
  } else {
    df %>% filter(version == "v1", gsID %in% unique(test$gsID)) %>% select(gsID, Start_v1 = Start)
  }
  matched_e <- match_lines_window(test, v1_s)
  n_tp    <- nrow(matched_e)
  n_total <- nrow(test)
  n_v1    <- nrow(v1_s)
  # Per-patient
  pats <- unique(test$gsID)
  pp <- map_dfr(pats, function(pid) {
    nc <- sum(test$gsID == pid)
    nv <- sum(v1_s$gsID == pid)
    nt <- sum(matched_e$gsID == pid)
    tibble(gsID = pid,
           missed = nv - nt, extra = nc - nt,
           exact  = as.integer(nc == nv),
           no_m   = as.integer(nv - nt == 0),
           no_e   = as.integer(nc - nt == 0))
  })
  tibble(
    version     = v,
    total       = n_total,
    tp          = n_tp,
    missed      = n_v1 - n_tp,
    extra       = n_total - n_tp,
    missed_mean = mean(pp$missed),
    missed_sd   = sd(pp$missed),
    extra_mean  = mean(pp$extra),
    extra_sd    = sd(pp$extra),
    n_exact     = sum(pp$exact),
    n_no_missed = sum(pp$no_m),
    n_no_extra  = sum(pp$no_e),
    n_pts       = nrow(pp)
  )
})

# Add v1 (first oncologist) as reference column for S3
v1_lines_s3 <- df %>% filter(version == "v1") %>% select(gsID, Start)
v1_pats_s3  <- unique(v1_lines_s3$gsID)

line_error_app_s3 <- bind_rows(
  tibble(
    version     = "v1",
    total       = nrow(v1_lines_s3),
    tp          = nrow(v1_lines_s3),   # v1 vs v1: all lines matched
    missed      = 0L,
    extra       = 0L,
    missed_mean = 0,
    missed_sd   = 0,
    extra_mean  = 0,
    extra_sd    = 0,
    n_exact     = length(v1_pats_s3),
    n_no_missed = length(v1_pats_s3),
    n_no_extra  = length(v1_pats_s3),
    n_pts       = length(v1_pats_s3)
  ),
  line_error_app
)

# S3-specific column map — adds First oncologist before Second oncologist
appendix_cols_s3 <- c("First oncologist", appendix_cols)
ver_col_s3       <- c(v1 = "First oncologist", ver_col)
all_src_s3       <- c("v1", all_src_app)

# Build wide table S3
make_s3_row <- function(metric, label, fmt_fn) {
  vals <- setNames(
    map_chr(all_src_s3, ~ fmt_fn(line_error_app_s3[line_error_app_s3$version == .x, ])),
    ver_col_s3[all_src_s3]
  )
  bind_cols(tibble(`Line-level error metric` = label,
                   `Metric reported`         = metric),
            as_tibble_row(vals)) %>%
    {tbl <- .; for (cn in appendix_cols_s3) if (!cn %in% names(tbl)) tbl[[cn]] <- "\u2014"; tbl} %>%
    select(`Line-level error metric`, `Metric reported`, all_of(appendix_cols_s3))
}

fmt_n   <- function(r) if (nrow(r) == 0 || is.na(r$total[1]))  "\u2014" else as.character(r$total[1])
fmt_tp  <- function(r) if (nrow(r) == 0 || is.na(r$tp[1]))     "\u2014" else as.character(r$tp[1])
fmt_ms  <- function(r) if (nrow(r) == 0 || is.na(r$missed[1])) "\u2014" else as.character(r$missed[1])
fmt_ex  <- function(r) if (nrow(r) == 0 || is.na(r$extra[1]))  "\u2014" else as.character(r$extra[1])
fmt_msd <- function(r) if (nrow(r) == 0 || is.na(r$missed_mean[1])) "\u2014" else
  sprintf("%.1f (%.1f)", r$missed_mean[1], r$missed_sd[1])
fmt_esd <- function(r) if (nrow(r) == 0 || is.na(r$extra_mean[1]))  "\u2014" else
  sprintf("%.1f (%.1f)", r$extra_mean[1], r$extra_sd[1])
fmt_exact <- function(r) if (nrow(r) == 0 || is.na(r$n_exact[1])) "\u2014" else
  sprintf("%d (%.1f%%)", r$n_exact[1], r$n_exact[1] / r$n_pts[1] * 100)
fmt_nom   <- function(r) if (nrow(r) == 0 || is.na(r$n_no_missed[1])) "\u2014" else
  sprintf("%d (%.1f%%)", r$n_no_missed[1], r$n_no_missed[1] / r$n_pts[1] * 100)
fmt_noe   <- function(r) if (nrow(r) == 0 || is.na(r$n_no_extra[1])) "\u2014" else
  sprintf("%d (%.1f%%)", r$n_no_extra[1], r$n_no_extra[1] / r$n_pts[1] * 100)

tableS3_app <- bind_rows(
  make_s3_row("No.",       "Total predicted lines",              fmt_n),
  make_s3_row("No.",       "Correctly matched lines",            fmt_tp),
  make_s3_row("No.",       "Missed lines",                       fmt_ms),
  make_s3_row("No.",       "Additional predicted lines",         fmt_ex),
  make_s3_row("Mean (SD)", "Missed lines per patient",           fmt_msd),
  make_s3_row("Mean (SD)", "Additional predicted lines per patient", fmt_esd),
  make_s3_row("No. (%)",   "Patients with exact line count match", fmt_exact),
  make_s3_row("No. (%)",   "Patients with no missed lines",        fmt_nom),
  make_s3_row("No. (%)",   "Patients with no additional predicted lines", fmt_noe)
)

# Drug-level error section (second sheet of TableS3)
# Builds the same metrics as line-level errors but for unique drugs (bag-of-drugs)
v1_drugs_total <- sum(map_int(v1_drugs$drugs_v1, length))
v1_drugs_pts   <- nrow(v1_drugs)

drug_error_s3 <- bind_rows(
  tibble(
    version        = "v1",
    n_pts          = v1_drugs_pts,
    total_drugs    = v1_drugs_total,
    total_tp       = v1_drugs_total,
    total_missing  = 0L,
    total_extra    = 0L,
    mean_extra     = 0, sd_extra     = 0,
    mean_missing   = 0, sd_missing   = 0,
    pct_no_extra   = 1, pct_no_missing = 1,
    n_exact        = v1_drugs_pts
  ),
  drug_error_summary
)

make_s3_drug_row <- function(metric, label, fmt_fn) {
  vals <- setNames(
    map_chr(all_src_s3, ~ fmt_fn(drug_error_s3[drug_error_s3$version == .x, ])),
    ver_col_s3[all_src_s3]
  )
  bind_cols(tibble(`Drug-level error metric` = label,
                   `Metric reported`         = metric),
            as_tibble_row(vals)) %>%
    {tbl <- .; for (cn in appendix_cols_s3) if (!cn %in% names(tbl)) tbl[[cn]] <- "\u2014"; tbl} %>%
    select(`Drug-level error metric`, `Metric reported`, all_of(appendix_cols_s3))
}

fmt_d_n    <- function(r) if (nrow(r) == 0 || is.na(r$total_drugs[1]))   "\u2014" else as.character(r$total_drugs[1])
fmt_d_tp   <- function(r) if (nrow(r) == 0 || is.na(r$total_tp[1]))      "\u2014" else as.character(r$total_tp[1])
fmt_d_ms   <- function(r) if (nrow(r) == 0 || is.na(r$total_missing[1])) "\u2014" else as.character(r$total_missing[1])
fmt_d_ex   <- function(r) if (nrow(r) == 0 || is.na(r$total_extra[1]))   "\u2014" else as.character(r$total_extra[1])
fmt_d_msd  <- function(r) if (nrow(r) == 0 || is.na(r$mean_missing[1]))  "\u2014" else
  sprintf("%.1f (%.1f)", r$mean_missing[1], r$sd_missing[1])
fmt_d_esd  <- function(r) if (nrow(r) == 0 || is.na(r$mean_extra[1]))    "\u2014" else
  sprintf("%.1f (%.1f)", r$mean_extra[1], r$sd_extra[1])
fmt_d_nom  <- function(r) if (nrow(r) == 0 || is.na(r$pct_no_missing[1])) "\u2014" else
  sprintf("%d (%.1f%%)", round(r$pct_no_missing[1] * r$n_pts[1]),
          r$pct_no_missing[1] * 100)
fmt_d_noe  <- function(r) if (nrow(r) == 0 || is.na(r$pct_no_extra[1]))  "\u2014" else
  sprintf("%d (%.1f%%)", round(r$pct_no_extra[1] * r$n_pts[1]),
          r$pct_no_extra[1] * 100)
fmt_d_exact <- function(r) if (nrow(r) == 0 || is.na(r$n_exact[1])) "\u2014" else
  sprintf("%d (%.1f%%)", r$n_exact[1], r$n_exact[1] / r$n_pts[1] * 100)

tableS3_drug <- bind_rows(
  make_s3_drug_row("No.",       "Total predicted unique drugs",          fmt_d_n),
  make_s3_drug_row("No.",       "Correctly identified unique drugs",     fmt_d_tp),
  make_s3_drug_row("No.",       "Missed unique drugs",                   fmt_d_ms),
  make_s3_drug_row("No.",       "Additional predicted unique drugs",     fmt_d_ex),
  make_s3_drug_row("Mean (SD)", "Missed drugs per patient",              fmt_d_msd),
  make_s3_drug_row("Mean (SD)", "Extra drugs per patient",               fmt_d_esd),
  make_s3_drug_row("No. (%)",   "Patients with exact drug set match",    fmt_d_exact),
  make_s3_drug_row("No. (%)",   "Patients with no missed drugs",         fmt_d_nom),
  make_s3_drug_row("No. (%)",   "Patients with no extra drugs",          fmt_d_noe)
)

writexl::write_xlsx(
  list(line_errors = tableS3_app, drug_errors = tableS3_drug),
  file.path(TBL_DIR, "TableS3_line_errors.xlsx")
)
cat("  TableS3_line_errors.xlsx written (2 sheets: line_errors, drug_errors)\n")

# Appendix Table 4: Missed / added drugs (wide)
# All missed/added drugs ranked by total count across all sources
drug_counts_app <- all_drug_errors %>%
  count(version, type, drug) %>%
  mutate(col = coalesce(ver_col[version], version)) %>%
  select(drug, type, col, n)

all_missed_drugs_app <- all_drug_errors %>%
  filter(type == "missing") %>%
  count(drug, name = "total") %>%
  arrange(desc(total)) %>%
  slice_head(n = 20) %>%
  pull(drug)

all_added_drugs_app <- all_drug_errors %>%
  filter(type == "extra") %>%
  count(drug, name = "total") %>%
  arrange(desc(total)) %>%
  slice_head(n = 20) %>%
  pull(drug)

make_drug_row <- function(drug_name, err_type) {
  row_vals <- drug_counts_app %>%
    filter(drug == drug_name, type == err_type) %>%
    select(col, n) %>%
    deframe()
  bind_cols(
    tibble(Drug = str_to_title(drug_name)),
    as_tibble_row(setNames(
      map_chr(appendix_cols, ~ if (.x %in% names(row_vals))
                as.character(row_vals[[.x]]) else "0"),
      appendix_cols
    ))
  )
}

make_section_row <- function(label) {
  bind_cols(tibble(Drug = label),
            setNames(as_tibble(matrix("", 1, length(appendix_cols))),
                     appendix_cols))
}

tableS4_app <- bind_rows(
  make_section_row("Missed Drugs"),
  map_dfr(all_missed_drugs_app, ~ make_drug_row(.x, "missing")),
  make_section_row("Added Drugs"),
  map_dfr(all_added_drugs_app,  ~ make_drug_row(.x, "extra"))
)

writexl::write_xlsx(tableS4_app, file.path(TBL_DIR, "TableS4_missed_added_drugs.xlsx"))
cat("  TableS4_missed_added_drugs.xlsx written\n")

# TABLES OUTPUT — summary tables for manuscript / quick reference
# All objects were computed above; these writes consolidate them into tables/.

# 1. Drug accuracy summary
#    bag_jaccard    — Jaccard (bag of drugs) vs v1, all models + human benchmark
#    lines_jaccard  — same for lines of therapy
#    bag_thresholds — % patients with Jaccard == 1 and >= 0.9
writexl::write_xlsx(
  list(
    bag_jaccard    = final_table,
    lines_jaccard  = final_table_lines,
    bag_thresholds = llm_threshold_summary
  ),
  file.path(LAB_MTG_DIR, "drug_accuracy_summary.xlsx")
)

# 2. Lines comparison summary. Internal comparison only — the reported line metrics
# come from line_metrics.R.
writexl::write_xlsx(
  list(
    lines_count     = lines_count_summary,
    intent_accuracy = intent_summary,
    reason_stop     = reason_summary,
    line_f1         = line_f1_summary
  ),
  file.path(LAB_MTG_DIR, "SUPERSEDED_old_matcher_lines_comparison.xlsx")
)

# 3. Drug error analysis
#    error_summary   — mean/median extra + missing drugs
#    top_extra_drugs — top 10 additional (extra) predicted drugs per version
#    (top_missing_drugs already in missed_drugs_summary.xlsx)
writexl::write_xlsx(
  list(
    error_summary   = drug_error_summary,
    top_extra_drugs = top_extra_drugs
  ),
  file.path(LAB_MTG_DIR, "drug_error_analysis.xlsx")
)

print("Lab meeting files written to Analysis_Output/lab_meeting_figs_tables/:")
print("  drug_accuracy_summary.xlsx      (3 sheets)")
print("  lines_comparison_summary.xlsx   (4 sheets)")
print("  drug_error_analysis.xlsx        (2 sheets)")
print("  missed_drugs_summary.xlsx       (2 sheets)")
print("  TableS2_drug_rows.xlsx          (Appendix Table 2 drug section — read by manuscript_tables.R)")
print("Final tables written to Analysis_Output/tables/:")
print("  TableS3_line_errors.xlsx        (Appendix Table 3 — 2 sheets: line_errors, drug_errors)")
print("  TableS4_missed_added_drugs.xlsx (Appendix Table 4 — missed/added drugs wide)")

# FigS5 (chart span vs Jaccard scatter) is in manuscript_figures.R

dir.create(LOGS_DIR, showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(LOGS_DIR, "session_info_drug_comparison.txt"))
