library(tidyverse)
library(readxl)
library(writexl)
library(jsonlite)

source("../config.R")

# All path and identifier constants come from config.R (sourced above)
run_id    <- RUN_ID
llm_dir   <- EXP_DIR
final_dir <- FINAL_DIR
int_dir   <- INT_DIR
llms      <- ALL_LLMS

dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

parse_jsonl_line <- function(ln) {
  tryCatch({
    rec   <- fromJSON(ln, flatten = FALSE)
    meta  <- rec[["_meta"]]  %||% list()
    usage <- rec[["usage"]]  %||% list()
    ctd   <- usage[["completion_tokens_details"]] %||% list()

    reasoning <- ctd[["reasoning_tokens"]] %||% 0L
    thinking  <- usage[["thoughts_tokens"]] %||% usage[["thinking_tokens"]] %||% 0L
    if (thinking > 0 && reasoning == 0) reasoning <- thinking

    tibble(
      llm               = meta[["llm_key"]]   %||% NA_character_,
      # Randomly generated study identifier assigned upstream; the
      # study-ID-to-chart crosswalk is held outside this repository and
      # outside the analytic output directory.
      patient_id        = meta[["study_id"]]   %||% NA_character_,
      task              = meta[["task"]]       %||% NA_character_,
      latency_s         = as.numeric(meta[["latency_s"]]  %||% 0),
      prompt_tokens     = as.integer(usage[["prompt_tokens"]]     %||% usage[["input_tokens"]]  %||% 0L),
      completion_tokens = as.integer(usage[["completion_tokens"]] %||% usage[["output_tokens"]] %||% 0L),
      total_tokens      = as.integer(usage[["total_tokens"]]      %||% 0L),
      reasoning_tokens  = as.integer(reasoning)
    )
  }, error = function(e) tibble())
}

####

#First, pull date of diagnosis evidence — initial dx, locoregional recurrence, and metastatic
load_dx_evidence <- function(llm_id) {
  read_xlsx(file.path(llm_dir, run_id, llm_id,
      sprintf('date_of_diagnosis_results_%s.xlsx', run_id))) %>%
    mutate(across(everything(), as.character)) %>%
    select(patient_id,
           original_certainty, original_evidence,
           original_evidence_chunk_id, original_evidence_vector_rank,
           locoregional_recurrence_certainty, locoregional_recurrence_evidence,
           locoregional_recurrence_evidence_chunk_id, locoregional_recurrence_evidence_vector_rank,
           metastatic_certainty, metastatic_evidence,
           metastatic_evidence_chunk_id, metastatic_evidence_vector_rank,
           rationale) %>%
    rename(gsID                    = patient_id,
           dx_original_certainty   = original_certainty,
           dx_original_evidence    = original_evidence,
           dx_original_chunk_id    = original_evidence_chunk_id,
           dx_original_vector_rank = original_evidence_vector_rank,
           dx_lrr_certainty        = locoregional_recurrence_certainty,
           dx_lrr_evidence         = locoregional_recurrence_evidence,
           dx_lrr_chunk_id         = locoregional_recurrence_evidence_chunk_id,
           dx_lrr_vector_rank      = locoregional_recurrence_evidence_vector_rank,
           dx_met_certainty        = metastatic_certainty,
           dx_met_evidence         = metastatic_evidence,
           dx_met_chunk_id         = metastatic_evidence_chunk_id,
           dx_met_vector_rank      = metastatic_evidence_vector_rank,
           dx_rationale            = rationale) %>%
    mutate(version = llm_id)
}

dx_ev <- map_dfr(llms, load_dx_evidence)

#Date of birth
load_dob_evidence <- function(llm_id) {
  read_xlsx(file.path(llm_dir, run_id, llm_id,
      sprintf('date_of_birth_results_%s.xlsx', run_id))) %>%
    mutate(across(everything(), as.character)) %>%
    select(patient_id,
           dob_certainty, dob_evidence,
           dob_evidence_chunk_id, dob_evidence_chunk_rank,
           rationale) %>%
    rename(gsID            = patient_id,
           dob_vector_rank = dob_evidence_chunk_rank,
           dob_chunk_id    = dob_evidence_chunk_id,
           dob_rationale   = rationale) %>%
    mutate(version = llm_id)
}

dob_ev <- map_dfr(llms, load_dob_evidence)

#Pathology and clinical stage
load_pa_evidence <- function(llm_id) {
  read_xlsx(file.path(llm_dir, run_id, llm_id,
      sprintf('pathology_stage_results_%s.xlsx', run_id))) %>%
    mutate(across(everything(), as.character)) %>%
    select(patient_id,
           pathology_chunk_id, pathology_chunk_rank,
           clinical_stage_chunk_id, clinical_stage_chunk_rank,
           clinical_stage_verbatim_stage_quote) %>%
    rename(gsID                  = patient_id,
           pathology_vector_rank = pathology_chunk_rank,
           stage_chunk_id        = clinical_stage_chunk_id,
           stage_vector_rank     = clinical_stage_chunk_rank,
           stage_evidence        = clinical_stage_verbatim_stage_quote) %>%
    mutate(version = llm_id)
}

pa_ev <- map_dfr(llms, load_pa_evidence)

#Second opinion file gives us chart size — visits, timespan, and chunks indexed per year
load_sop_evidence <- function(llm_id) {
  df <- read_xlsx(file.path(llm_dir, run_id, llm_id,
      sprintf('second_opinion_or_not_results_%s.xlsx', run_id))) %>%
    mutate(across(everything(), as.character))

  chunk_yr_cols <- names(df)[grepl('visit_summary_evidence_chunks_per_year_', names(df))]

  df %>%
    select(patient_id,
           treated_rationale, second_opinion_evidence,
           visit_summary_visits_total,
           visit_summary_visit_timespan_days,
           oncologist_from_bm25,
           all_of(chunk_yr_cols)) %>%
    rename(gsID                 = patient_id,
           oncologist_rationale = treated_rationale,
           oncologist_evidence  = second_opinion_evidence,
           file_visits_total    = visit_summary_visits_total,
           file_timespan_days   = visit_summary_visit_timespan_days,
           oncologist_bm25      = oncologist_from_bm25) %>%
    rowwise() %>%
    mutate(
      file_total_chunks     = sum(as.numeric(c_across(all_of(chunk_yr_cols))), na.rm = TRUE),
      file_years_with_notes = sum(as.numeric(c_across(all_of(chunk_yr_cols))) > 0, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(-all_of(chunk_yr_cols)) %>%
    mutate(
      version               = llm_id,
      file_visits_total     = as.numeric(file_visits_total),
      file_timespan_days    = as.numeric(file_timespan_days),
      file_chunks_per_visit = ifelse(
        !is.na(file_visits_total) & file_visits_total > 0,
        file_total_chunks / file_visits_total, NA_real_
      )
    )
}

sop_ev <- map_dfr(llms, load_sop_evidence)

#Germline genetics
load_germ_evidence <- function(llm_id) {
  read_xlsx(file.path(llm_dir, run_id, llm_id,
      sprintf('genetics_germline_results_%s.xlsx', run_id))) %>%
    mutate(across(everything(), as.character)) %>%
    select(patient_id, genetic_testing_done, year_of_testing, company,
           BRCA1_result, BRCA1_variant, BRCA2_result, BRCA2_variant) %>%
    rename(gsID                  = patient_id,
           germline_testing_done = genetic_testing_done,
           germline_year         = year_of_testing,
           germline_company      = company) %>%
    mutate(version = llm_id)
}

germ_ev <- map_dfr(llms, load_germ_evidence)

#Somatic genetics
load_som_evidence <- function(llm_id) {
  read_xlsx(file.path(llm_dir, run_id, llm_id,
      sprintf('genetics_somatics_results_%s.xlsx', run_id))) %>%
    mutate(across(everything(), as.character)) %>%
    select(patient_id, genetic_testing_done, year_of_testing, company,
           PIK3CA_result, PIK3CA_variant, ESR1_result, ESR1_variant) %>%
    rename(gsID                 = patient_id,
           somatic_testing_done = genetic_testing_done,
           somatic_year         = year_of_testing,
           somatic_company      = company) %>%
    mutate(version = llm_id)
}

som_ev <- map_dfr(llms, load_som_evidence)

#Merge all evidence tables and write to final_dfs
vector_rank_df <- dx_ev %>%
  left_join(dob_ev,  by = c('gsID', 'version')) %>%
  left_join(pa_ev,   by = c('gsID', 'version')) %>%
  left_join(sop_ev,  by = c('gsID', 'version')) %>%
  left_join(germ_ev, by = c('gsID', 'version')) %>%
  left_join(som_ev,  by = c('gsID', 'version'))

# Drop every column that can carry chart source text (evidence snippets,
# rationales, verbatim quotes, chunk identifiers, BM25 hits) before writing.
# Only ranks, certainty labels, chart-size counts and genetics codes survive.
vector_rank_output <- vector_rank_df %>%
  select(
    -matches("evidence|rationale|verbatim|chunk_id|bm25",
             ignore.case = TRUE)
  )

write_xlsx(
  vector_rank_output,
  file.path(final_dir, "llm_vector_rank.xlsx")
)

head(vector_rank_output)
table(vector_rank_output$version)

#Now parse the JSONL files to get token usage per patient × task
jsonl_paths <- map_chr(llms, ~ file.path(
  llm_dir, run_id, .x,
  sprintf('detailed_run_info_%s_%s.jsonl', run_id, .x)
)) %>% keep(file.exists)

token_raw <- map_dfr(jsonl_paths, function(path) {
  lines <- readLines(path, warn = FALSE) %>% keep(nzchar)
  map_dfr(lines, parse_jsonl_line)
}) %>%
  filter(llm %in% llms)

write_xlsx(token_raw, file.path(final_dir, 'llm_token_usage.xlsx'))

table(token_raw$llm)
table(token_raw$task)

####

vr  <- read_xlsx(file.path(final_dir, 'llm_vector_rank.xlsx'))
tok <- read_xlsx(file.path(final_dir, 'llm_token_usage.xlsx'))

#First, tidy the rank columns into long form so we can summarise across variables
rank_cols <- tribble(
  ~col,                       ~variable,
  'dx_original_vector_rank',  'date_initial_diagnosis',
  'dx_lrr_vector_rank',       'date_locoregional_recurrence',
  'dx_met_vector_rank',       'date_metastatic_diagnosis',
  'dob_vector_rank',          'date_of_birth',
  'pathology_vector_rank',    'pathology_stage',
  'stage_vector_rank',        'clinical_stage'
)

rank_long <- rank_cols %>%
  pmap_dfr(function(col, variable) {
    vr %>%
      select(gsID, version, rank_val = all_of(col)) %>%
      mutate(variable = variable,
             rank_val = suppressWarnings(as.integer(rank_val)))
  })

rank_by_variable <- rank_long %>%
  group_by(variable, version) %>%
  summarise(
    n              = n(),
    n_missing      = sum(is.na(rank_val)),
    pct_missing    = round(n_missing / n * 100, 1),
    median_rank    = median(rank_val, na.rm = TRUE),
    mean_rank      = round(mean(rank_val, na.rm = TRUE), 2),
    sd_rank        = round(sd(rank_val, na.rm = TRUE), 2),
    pct_rank_1     = round(sum(rank_val == 1,  na.rm = TRUE) / sum(!is.na(rank_val)) * 100, 1),
    pct_rank_1_3   = round(sum(rank_val <= 3,  na.rm = TRUE) / sum(!is.na(rank_val)) * 100, 1),
    pct_rank_1_5   = round(sum(rank_val <= 5,  na.rm = TRUE) / sum(!is.na(rank_val)) * 100, 1),
    pct_rank_gt_10 = round(sum(rank_val >  10, na.rm = TRUE) / sum(!is.na(rank_val)) * 100, 1),
    rank_min       = min(rank_val, na.rm = TRUE),
    rank_max       = max(rank_val, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(variable, version)

rank_by_variable

#Full rank distribution by bucket — useful for plotting
rank_distribution <- rank_long %>%
  filter(!is.na(rank_val)) %>%
  mutate(rank_bucket = case_when(
    rank_val == 1  ~ '1',
    rank_val == 2  ~ '2',
    rank_val == 3  ~ '3',
    rank_val <= 5  ~ '4-5',
    rank_val <= 10 ~ '6-10',
    rank_val <= 20 ~ '11-20',
    TRUE           ~ '21+'
  )) %>%
  count(variable, version, rank_bucket) %>%
  group_by(variable, version) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup() %>%
  mutate(rank_bucket = factor(rank_bucket,
           levels = c('1','2','3','4-5','6-10','11-20','21+'))) %>%
  arrange(variable, version, rank_bucket)

#Now look at certainty labels per variable × version
cert_cols <- tribble(
  ~col,                     ~variable,
  'dx_original_certainty',  'date_initial_diagnosis',
  'dx_lrr_certainty',       'date_locoregional_recurrence',
  'dx_met_certainty',       'date_metastatic_diagnosis',
  'dob_certainty',          'date_of_birth'
)

certainty_summary <- cert_cols %>%
  pmap_dfr(function(col, variable) {
    vr %>%
      select(gsID, version, certainty = all_of(col)) %>%
      mutate(variable = variable)
  }) %>%
  mutate(certainty = str_to_lower(str_trim(certainty))) %>%
  count(variable, version, certainty) %>%
  group_by(variable, version) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup() %>%
  arrange(variable, version, desc(n))

certainty_summary

chart_size_cols <- c('gsID', 'version',
                     'file_visits_total', 'file_timespan_days',
                     'file_total_chunks', 'file_years_with_notes',
                     'file_chunks_per_visit')

patient_chart_size <- vr %>%
  select(any_of(chart_size_cols)) %>%
  mutate(
    across(c(file_visits_total, file_timespan_days, file_total_chunks,
             file_years_with_notes, file_chunks_per_visit),
           ~ suppressWarnings(as.numeric(.x))),
    file_timespan_yrs = round(file_timespan_days / 365.25, 2),
    chart_size_tier   = case_when(
      is.na(file_visits_total) ~ NA_character_,
      file_visits_total < 10   ~ 'Small (<10 visits)',
      file_visits_total < 30   ~ 'Medium (10-29)',
      file_visits_total < 60   ~ 'Large (30-59)',
      TRUE                     ~ 'Very large (60+)'
    )
  )

chart_size_summary <- patient_chart_size %>%
  group_by(version, chart_size_tier) %>%
  summarise(
    n                   = n(),
    median_visits       = median(file_visits_total,  na.rm = TRUE),
    median_timespan_yrs = median(file_timespan_yrs,  na.rm = TRUE),
    median_total_chunks = median(file_total_chunks,   na.rm = TRUE),
    .groups = 'drop'
  )

chart_size_summary

#Do patients with larger charts have harder-to-find evidence (higher rank)?
chart_vs_rank <- rank_long %>%
  left_join(
    patient_chart_size %>%
      select(gsID, version, file_visits_total, file_total_chunks,
             file_timespan_yrs, chart_size_tier) %>%
      distinct(),
    by = c('gsID', 'version')
  ) %>%
  filter(!is.na(rank_val) & !is.na(chart_size_tier)) %>%
  group_by(variable, chart_size_tier) %>%
  summarise(
    n            = n(),
    median_rank  = median(rank_val, na.rm = TRUE),
    mean_rank    = round(mean(rank_val, na.rm = TRUE), 2),
    pct_rank_1   = round(sum(rank_val == 1) / n() * 100, 1),
    pct_rank_1_3 = round(sum(rank_val <= 3) / n() * 100, 1),
    .groups = 'drop'
  ) %>%
  arrange(variable, chart_size_tier)

chart_vs_rank

#Token usage summary per LLM
per_pat_tok <- tok %>%
  group_by(llm, patient_id) %>%
  summarise(total_tokens = sum(total_tokens), total_latency = sum(latency_s),
            .groups = 'drop')

token_summary <- tok %>%
  group_by(llm) %>%
  summarise(
    total_calls         = n(),
    n_patients          = n_distinct(patient_id),
    prompt_tokens       = sum(prompt_tokens),
    completion_tokens   = sum(completion_tokens),
    reasoning_tokens    = sum(reasoning_tokens),
    total_tokens        = sum(total_tokens),
    total_latency_s     = round(sum(latency_s), 1),
    mean_latency_call_s = round(mean(latency_s), 2),
    .groups = 'drop'
  ) %>%
  left_join(
    per_pat_tok %>%
      group_by(llm) %>%
      summarise(mean_tokens_per_patient  = round(mean(total_tokens), 0),
                mean_latency_per_patient = round(mean(total_latency), 1),
                .groups = 'drop'),
    by = 'llm'
  )

token_summary

#Token breakdown by task — helps identify which tasks are most expensive
token_by_task <- tok %>%
  group_by(llm, task) %>%
  summarise(
    n_calls         = n(),
    mean_prompt     = round(mean(prompt_tokens)),
    mean_completion = round(mean(completion_tokens)),
    mean_reasoning  = round(mean(reasoning_tokens)),
    mean_total      = round(mean(total_tokens)),
    mean_latency_s  = round(mean(latency_s), 2),
    total_tokens    = sum(total_tokens),
    .groups = 'drop'
  ) %>%
  arrange(llm, desc(mean_total))

token_by_task

write_xlsx(
  list(
    rank_by_variable   = rank_by_variable,
    rank_distribution  = rank_distribution,
    certainty_summary  = certainty_summary,
    patient_chart_size = patient_chart_size,
    chart_size_summary = chart_size_summary,
    chart_size_vs_rank = chart_vs_rank,
    token_summary      = token_summary,
    token_by_task      = token_by_task
  ),
  file.path(LAB_MTG_DIR, 'vector_rank_analytics.xlsx')
)
