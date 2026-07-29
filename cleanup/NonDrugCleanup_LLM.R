library(tidyverse)
library(readxl)
library(writexl)

source("../config.R")

# Step 1: Dates of diagnosis — keep all cols for inspection
load_dx <- function(llm_id) {
  llm_dir <- file.path(EXP_DIR, RUN_ID, llm_id)
  read_excel(file.path(llm_dir, sprintf('date_of_diagnosis_results_%s.xlsx', RUN_ID))) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(version = llm_id)
}

dx_raw <- map_dfr(ALL_LLMS, load_dx)

dx_df <- dx_raw %>%
  select(patient_id, version,
         date_original_diagnosis, date_locoregional_recurrence_diagnosis, date_metastatic_diagnosis) %>%
  rename(gsID = patient_id,
         date_initial_diagnosis = date_original_diagnosis,
         date_first_locoregional_recurrence = date_locoregional_recurrence_diagnosis,
         date_first_metastatic_diagnosis = date_metastatic_diagnosis)

# Step 2: Date of birth — keep all cols for inspection
load_dob <- function(llm_id) {
  llm_dir <- file.path(EXP_DIR, RUN_ID, llm_id)
  read_excel(file.path(llm_dir, sprintf('date_of_birth_results_%s.xlsx', RUN_ID))) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(version = llm_id)
}

dob_raw <- map_dfr(ALL_LLMS, load_dob)

dob_df <- dob_raw %>%
  select(patient_id, version, date_of_birth) %>%
  rename(gsID = patient_id, birth_date = date_of_birth)

# Step 3: Pathology / stage / receptor status — keep all cols for inspection
load_pa <- function(llm_id) {
  llm_dir <- file.path(EXP_DIR, RUN_ID, llm_id)
  read_excel(file.path(llm_dir, sprintf('pathology_stage_results_%s.xlsx', RUN_ID))) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(version = llm_id)
}

pa_raw <- map_dfr(ALL_LLMS, load_pa)

pa_df <- pa_raw %>%
  select(patient_id, version, clinical_stage_original_clinical_stage,
         pathology_er_status, pathology_pr_status, pathology_her2_result_status) %>%
  rename(gsID = patient_id,
         stage_at_diagnosis = clinical_stage_original_clinical_stage,
         er = pathology_er_status,
         pr = pathology_pr_status,
         her2_result = pathology_her2_result_status) %>%
  select(gsID, version, stage_at_diagnosis, er, pr, her2_result)

#Replace unknown w NA for er/pr/her2
pa_df <- pa_df %>%
  mutate(across(
    c(er, pr, her2_result),
    ~ na_if(.x, "unknown")
  ))
  
# Step 4: Second opinion / primary oncologist — keep all cols for inspection
load_sop <- function(llm_id) {
  llm_dir <- file.path(EXP_DIR, RUN_ID, llm_id)
  read_excel(file.path(llm_dir, sprintf('second_opinion_or_not_results_%s.xlsx', RUN_ID))) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(version = llm_id)
}

sop_raw <- map_dfr(ALL_LLMS, load_sop)

sop_df <- sop_raw %>%
  select(patient_id, version, treating_oncologist) %>%
  rename(gsID = patient_id,
         primary_oncologist_shc = treating_oncologist) %>%
  mutate(
    primary_oncologist_shc = primary_oncologist_shc %>%
      str_extract("^[^,]+") %>%      # text before comma
      word(-1) %>%                  # last word before comma (last name)
      tolower()
  )

# Step 5: Genetics germline — keep all cols for inspection
load_germ <- function(llm_id) {
  llm_dir <- file.path(EXP_DIR, RUN_ID, llm_id)
  read_excel(file.path(llm_dir, sprintf('genetics_germline_results_%s.xlsx', RUN_ID))) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(version = llm_id)
}

germ_raw <- map_dfr(ALL_LLMS, load_germ)

germ_df <- germ_raw %>%
  select(patient_id, version, genetic_testing_done, BRCA1_result, BRCA1_variant, BRCA2_result, BRCA2_variant) %>%
  rename(
    gsID = patient_id,
    germline_testing_done = genetic_testing_done
  ) %>%
  mutate(
    germline_testing_done = case_when(
      germline_testing_done == TRUE  ~ "Yes",
      germline_testing_done == FALSE ~ "No",
      TRUE ~ "No"   # NA = not tested = No
    ),
    BRCA1_result = case_when(
      str_to_lower(str_trim(as.character(BRCA1_result))) %in% c("unknown", "") ~ NA_character_,
      TRUE ~ str_to_lower(str_trim(as.character(BRCA1_result)))
    ),
    BRCA2_result = case_when(
      str_to_lower(str_trim(as.character(BRCA2_result))) %in% c("unknown", "") ~ NA_character_,
      TRUE ~ str_to_lower(str_trim(as.character(BRCA2_result)))
    ),
    brca_status = case_when(
      BRCA1_result == "mutated" | BRCA2_result == "mutated" ~ "mutated",
      BRCA1_result == "negative" & BRCA2_result == "negative" ~ "negative",
      TRUE ~ NA_character_
    ),
    brca1_variant = if_else(BRCA1_result == "mutated", as.character(BRCA1_variant), NA_character_),
    brca2_variant = if_else(BRCA2_result == "mutated", as.character(BRCA2_variant), NA_character_)
  ) %>%
  select(gsID, version, germline_testing_done, brca_status, brca1_variant, brca2_variant)

# Step 6: Genetics somatics — keep all cols for inspection
load_som <- function(llm_id) {
  llm_dir <- file.path(EXP_DIR, RUN_ID, llm_id)
  read_excel(file.path(llm_dir, sprintf('genetics_somatics_results_%s.xlsx', RUN_ID))) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(version = llm_id)
}

som_raw <- map_dfr(ALL_LLMS, load_som)

som_df <- som_raw %>%
  select(patient_id, version, genetic_testing_done, PIK3CA_result, PIK3CA_variant, ESR1_result, ESR1_variant) %>%
  rename(
    gsID = patient_id,
    somatic_testing_done = genetic_testing_done
  ) %>%
  mutate(
    somatic_testing_done = case_when(
      somatic_testing_done == TRUE  ~ "Yes",
      somatic_testing_done == FALSE ~ "No",
      TRUE ~ "No"
    ),
    PIK3CA_result = case_when(
      str_to_lower(str_trim(as.character(PIK3CA_result))) %in% c("unknown", "") ~ NA_character_,
      TRUE ~ str_to_lower(str_trim(as.character(PIK3CA_result)))
    ),
    ESR1_result = case_when(
      str_to_lower(str_trim(as.character(ESR1_result))) %in% c("unknown", "") ~ NA_character_,
      TRUE ~ str_to_lower(str_trim(as.character(ESR1_result)))
    ),
    pik3ca_status = case_when(
      PIK3CA_result == "mutated" ~ "mutated",
      PIK3CA_result == "negative" ~ "negative",
      TRUE ~ NA_character_
    ),
    pik3ca_variant = if_else(PIK3CA_result == "mutated", as.character(PIK3CA_variant), NA_character_),
    esr1_status = case_when(
      ESR1_result == "mutated" ~ "mutated",
      ESR1_result == "negative" ~ "negative",
      TRUE ~ NA_character_
    ),
    esr1_variant = if_else(ESR1_result == "mutated", as.character(ESR1_variant), NA_character_)
  ) %>%
  select(gsID, version, somatic_testing_done, pik3ca_status, pik3ca_variant, esr1_status, esr1_variant)

# Step 7: Date of death / vital status
load_dod <- function(llm_id) {
  path <- file.path(EXP_DIR, RUN_ID, llm_id,
                    sprintf('date_of_death_results_%s.xlsx', RUN_ID))
  if (!file.exists(path)) {
    message('  date_of_death not found for ', llm_id, ' — skipping')
    return(tibble())
  }
  read_excel(path) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(version = llm_id)
}

dod_raw <- map_dfr(ALL_LLMS, load_dod)

if (nrow(dod_raw) > 0) {
  dod_df <- dod_raw %>%
    select(patient_id, version, vital_status, date_of_death) %>%
    rename(gsID = patient_id) %>%
    mutate(vital_status = tolower(trimws(vital_status)))
} else {
  dod_df <- tibble(gsID = character(), version = character(),
                   vital_status = character(), date_of_death = character())
}

# Step 8: Merge all pieces
llm_combined <- dx_df %>%
  left_join(dob_df, by = c('gsID', 'version')) %>%
  left_join(pa_df, by = c('gsID', 'version')) %>%
  left_join(sop_df, by = c('gsID', 'version')) %>%
  left_join(germ_df, by = c('gsID', 'version')) %>%
  left_join(som_df, by = c('gsID', 'version')) %>%
  left_join(dod_df, by = c('gsID', 'version'))

# Step 9: Final cleanup — dates, hr/her2 combined col
llm_combined <- llm_combined %>%
  mutate(
    date_initial_diagnosis = parse_date(date_initial_diagnosis),
    date_first_locoregional_recurrence = parse_date(date_first_locoregional_recurrence),
    date_first_metastatic_diagnosis = parse_date(date_first_metastatic_diagnosis),
    birth_date = parse_date(birth_date),
    date_of_death = parse_date(date_of_death),
    hr = case_when(
      er == "positive" | pr == "positive" ~ "HR+",
      er %in% c("negative", "low") & pr %in% c("negative", "low") ~ "HR-",
      TRUE ~ "unknown"
    ),
    her2 = case_when(
      her2_result == "positive" ~ "HER2+",
      her2_result == "negative" ~ "HER2-",
      TRUE ~ NA_character_
    ),
    hr_her2_status = case_when(
      hr == "unknown" ~ "unknown",
      TRUE ~ paste0(hr, "/", her2)
    )
  )

#Drop no longer needed cols
llm_combined <- llm_combined %>%
  select(-c(er, pr, her2_result, hr, her2))

#Make sure dates are dates
llm_combined <- llm_combined %>%
  mutate(
    date_initial_diagnosis = as.Date(date_initial_diagnosis),
    date_first_locoregional_recurrence = as.Date(date_first_locoregional_recurrence),
    date_first_metastatic_diagnosis = as.Date(date_first_metastatic_diagnosis),
    birth_date = as.Date(birth_date),
    date_of_death = as.Date(date_of_death)
  )

#make sure stage is numeric
llm_combined <- llm_combined %>%
  mutate(stage_at_diagnosis = as.numeric(stage_at_diagnosis))

# Coerce all variant columns to character (some are logical when all NA)
llm_combined <- llm_combined %>%
  mutate(across(c(brca1_variant, brca2_variant, pik3ca_variant, esr1_variant), as.character))

write_xlsx(llm_combined, LLM_NON_DRUG_CLEANED)
