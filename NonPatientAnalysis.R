library(tidyverse)
library(readxl)
library(writexl)

source("config.R")

# All path and identifier constants come from config.R
# BASE_DIR, INT_DIR, FINAL_DIR, TBL_DIR, RUN_ID, ALL_LLMS defined there
LLM_DIR <- EXP_DIR

set.seed(RANDOM_SEED)

### Step 1: Read vector rank data
vr <- read_excel(file.path(FINAL_DIR, 'llm_vector_rank.xlsx'))

### Step 2: Chunk rank summary by task × LLM
rank_cols <- names(vr)[grepl("vector_rank$", names(vr))]

chunk_rank_detail <- map_dfr(rank_cols, function(rc) {
  task_label <- rc %>%
    str_remove("_vector_rank$") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
  vr %>%
    select(gsID, version, rank_val = all_of(rc)) %>%
    mutate(
      variable  = task_label,
      rank_val  = as.numeric(rank_val)
    ) %>%
    filter(!is.na(rank_val))
})

chunk_rank_summary <- chunk_rank_detail %>%
  group_by(variable, version) %>%
  summarise(
    n             = n(),
    rank_median   = median(rank_val, na.rm = TRUE),
    rank_mean     = round(mean(rank_val, na.rm = TRUE), 2),
    rank_sd       = round(sd(rank_val, na.rm = TRUE), 2),
    pct_rank_1    = round(sum(rank_val == 1) / n() * 100, 1),
    pct_rank_1_3  = round(sum(rank_val <= 3) / n() * 100, 1),
    rank_min      = min(rank_val, na.rm = TRUE),
    rank_max      = max(rank_val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(variable, version)

### Step 3: Note frequency / documentation richness per patient × LLM
# Chart-size fields (visits, timespan, chunks) are already in llm_vector_rank.xlsx
# from NonDrugCleanup_VectorRank.R — no need to re-read the raw run files.
note_freq <- vr %>%
  select(gsID, version,
         visits_total       = file_visits_total,
         timespan_days      = file_timespan_days,
         total_chunks       = file_total_chunks,
         n_years_with_notes = file_years_with_notes) %>%
  mutate(across(c(visits_total, timespan_days, total_chunks, n_years_with_notes),
                ~ suppressWarnings(as.numeric(.x))))

### Step 4: Run metadata per LLM
load_run_metadata <- function(llm_id) {
  meta_path <- file.path(LLM_DIR, RUN_ID, llm_id,
    sprintf('run_metadata_%s.xlsx', RUN_ID))
  if (!file.exists(meta_path)) return(tibble())
  read_excel(meta_path)
}

run_metadata <- map_dfr(ALL_LLMS, load_run_metadata)

### Step 5: Token usage
# Read from the pre-processed file written by NonDrugCleanup_VectorRank.R.
token_df <- read_excel(file.path(FINAL_DIR, 'llm_token_usage.xlsx')) %>%
  filter(llm %in% ALL_LLMS)

### 5a: Per-LLM summary
per_pat <- token_df %>%
  group_by(llm, patient_id) %>%
  summarise(calls = n(), tokens = sum(total_tokens), latency = sum(latency_s),
            .groups = "drop")

token_usage <- token_df %>%
  group_by(llm) %>%
  summarise(
    total_calls       = n(),
    patients          = n_distinct(patient_id),
    prompt_tokens     = sum(prompt_tokens),
    completion_tokens = sum(completion_tokens),
    reasoning_tokens  = sum(reasoning_tokens),
    total_tokens      = sum(total_tokens),
    total_latency_s   = round(sum(latency_s), 1),
    mean_lat_call_s   = round(sum(latency_s) / n(), 2),
    .groups = "drop"
  ) %>%
  left_join(
    per_pat %>%
      group_by(llm) %>%
      summarise(
        mean_tok_pat = round(mean(tokens)),
        mean_lat_pat = round(mean(latency), 1),
        .groups = "drop"
      ),
    by = "llm"
  )

### 5b: Per-LLM × task breakdown
token_by_task <- token_df %>%
  group_by(llm, task) %>%
  summarise(
    n              = n(),
    mean_prompt    = round(mean(prompt_tokens)),
    mean_compl     = round(mean(completion_tokens)),
    mean_reasoning = round(mean(reasoning_tokens)),
    mean_total     = round(mean(total_tokens)),
    mean_latency_s = round(mean(latency_s), 2),
    .groups = "drop"
  )

### Step 6: Documentation richness vs accuracy
non_drug_long <- read_excel(file.path(FINAL_DIR, 'non_drug_long.xlsx'))

gt_df   <- non_drug_long %>% filter(version == "gt") %>% select(-version)
comp_df <- non_drug_long %>% filter(version != "gt" & version != "crc")

merged <- comp_df %>%
  inner_join(gt_df, by = "gsID", suffix = c("", "_gt"))

date_vars <- c('date_initial_diagnosis', 'date_first_locoregional_recurrence',
               'date_first_metastatic_diagnosis', 'birth_date')
cat_vars  <- c('stage_at_diagnosis', 'hr_her2_status', 'primary_oncologist_shc')
gen_vars  <- c('brca_status', 'pik3ca_status', 'esr1_status')

### Per-patient accuracy: count of correct extractions / total compared
# Date accuracy: within 90 days = correct
date_correct_mat <- map(date_vars, function(v) {
  d1 <- as.Date(merged[[v]])
  d2 <- as.Date(merged[[paste0(v, "_gt")]])
  both_present <- !is.na(d1) & !is.na(d2)
  correct <- both_present & abs(as.numeric(d1 - d2)) <= 90
  tibble(compared = as.integer(both_present), correct = as.integer(correct))
})

n_date_compared <- reduce(map(date_correct_mat, ~ .x$compared), `+`)
n_date_correct  <- reduce(map(date_correct_mat, ~ .x$correct), `+`)

# Categorical accuracy: exact match = correct
all_cat <- c(cat_vars, gen_vars)
cat_correct_mat <- map(all_cat, function(v) {
  a <- as.character(merged[[v]])
  b <- as.character(merged[[paste0(v, "_gt")]])
  both_present <- !is.na(a) & !is.na(b)
  correct <- both_present & a == b
  tibble(compared = as.integer(both_present), correct = as.integer(correct))
})

n_cat_compared <- reduce(map(cat_correct_mat, ~ .x$compared), `+`)
n_cat_correct  <- reduce(map(cat_correct_mat, ~ .x$correct), `+`)

patient_accuracy <- merged %>%
  select(gsID, version) %>%
  mutate(
    n_date_correct   = n_date_correct,
    n_date_compared  = n_date_compared,
    n_cat_correct    = n_cat_correct,
    n_cat_compared   = n_cat_compared,
    n_total_correct  = n_date_correct + n_cat_correct,
    n_total_compared = n_date_compared + n_cat_compared,
    pct_correct      = ifelse(n_total_compared > 0,
                              n_total_correct / n_total_compared, NA_real_)
  )

### Join note frequency to patient accuracy
doc_richness_vs_accuracy <- patient_accuracy %>%
  inner_join(note_freq, by = c("gsID", "version")) %>%
  select(gsID, version, visits_total, timespan_days, total_chunks,
         n_years_with_notes, n_total_correct, n_total_compared, pct_correct)

### Step 6b: Extraction rate — per LLM × variable completeness
# n_patients: total patients in the run
# n_extracted: patients where LLM returned a non-NA value
# n_gt_nonmissing: patients where GT has a non-NA value (reference coverage)
# extraction_pct: n_extracted / n_patients * 100
nondrug_vars <- c('date_initial_diagnosis', 'date_first_locoregional_recurrence',
                  'date_first_metastatic_diagnosis', 'birth_date',
                  'stage_at_diagnosis', 'hr_her2_status', 'primary_oncologist_shc',
                  'brca_status', 'pik3ca_status', 'esr1_status')

avail_vars <- intersect(nondrug_vars, names(non_drug_long))

# GT non-missing counts per variable (reference — how many patients have the value)
n_gt_per_var <- non_drug_long %>%
  filter(version == "gt") %>%
  summarise(across(all_of(avail_vars), ~ sum(!is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_gt_nonmissing")

# Per-LLM extraction counts (non-NA responses out of all patients)
extraction_rate <- non_drug_long %>%
  filter(version %in% ALL_LLMS) %>%
  group_by(version) %>%
  summarise(n_patients = n(),
            across(all_of(avail_vars), ~ sum(!is.na(.x))), .groups = "drop") %>%
  pivot_longer(-c(version, n_patients), names_to = "variable", values_to = "n_extracted") %>%
  left_join(n_gt_per_var, by = "variable") %>%
  mutate(extraction_pct = round(n_extracted / n_patients * 100, 1)) %>%
  select(version, variable, n_patients, n_gt_nonmissing, n_extracted, extraction_pct) %>%
  arrange(variable, version)

### Step 7: Export
write_xlsx(
  list(
    chunk_rank_summary       = chunk_rank_summary,
    chunk_rank_detail        = chunk_rank_detail,
    note_frequency           = note_freq,
    run_metadata             = run_metadata,
    token_usage              = token_usage,
    token_by_task            = token_by_task,
    doc_richness_vs_accuracy = doc_richness_vs_accuracy,
    extraction_rate          = extraction_rate
  ),
  file.path(LAB_MTG_DIR, 'non_patient_analysis_results.xlsx')
)

# Session info (for reproducibility)
dir.create(LOGS_DIR, showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()),
           file.path(LOGS_DIR, "session_info_non_patient_analysis.txt"))
