library(tidyverse)
library(readxl)
library(writexl)

source("../config.R")

input_df <- read_excel(GT_NON_DRUG_XLSX)

input_df <- input_df %>%
  select(gsID, birth_date, deceased_clinician, deceased_STARR, last_note_date, exact_death_date_available_clinician,
         death_date_clinician, death_date_STARR, treated_at_shc, primary_oncologist_shc,
         date_initial_diagnosis, date_first_locoregional_recurrence,
         date_first_metastatic_diagnosis, stage_at_diagnosis, hr_her2_status)

#Make all dates date format not POSIXct
input_df <- input_df %>%
  mutate(across(where(lubridate::is.POSIXct), lubridate::as_date)) %>%
  mutate(across(any_of(c("date_initial_diagnosis", "last_note_date", "birth_date",
                          "death_date_clinician", "death_date_STARR",
                          "date_first_locoregional_recurrence", "date_first_metastatic_diagnosis")),
                ~ as.Date(format(.x, "%Y-%m-%d")))) %>%
  mutate(stage_at_diagnosis = as.numeric(stage_at_diagnosis))

# Remove estimated clinician dates
input_df <- input_df %>%
  mutate(
    death_date_clinician = if_else(
      exact_death_date_available_clinician == "No",
      as.Date(NA),
      death_date_clinician
    )
  )

# Rename deceased answers for clarity
input_df <- input_df %>%
  mutate(
    deceased_clinician = case_when(
      deceased_clinician == TRUE  ~ "dead",
      deceased_clinician == FALSE ~ "not_dead",
      is.na(deceased_clinician)   ~ "unknown"
    ),
    deceased_STARR = case_when(
      deceased_STARR == TRUE  ~ "dead",
      deceased_STARR == FALSE ~ "not_dead",
      is.na(deceased_STARR)   ~ "unknown"
    )
  )

# Lowercase oncologist names for matching
input_df <- input_df %>%
  mutate(primary_oncologist_shc = tolower(trimws(primary_oncologist_shc)))

# Add version column
input_df <- input_df %>%
  mutate(version = "gt")

# Drop exact clinician yes/no column
input_df <- input_df %>%
  select(-exact_death_date_available_clinician)

# Merge genetics
genetics_df <- read_excel(GT_GENETICS_XLSX)

#Clean up genetics
genetics_df <- genetics_df %>%
  select(gsID, ESR1_status, ESR1_variant, PIK3CA_status, PIK3CA_variant,
         BRCA_testing, BRCA_status, BRCA1_variant, BRCA2_variant) %>%
  rename(
    esr1_status = ESR1_status,
    esr1_variant = ESR1_variant,
    pik3ca_status = PIK3CA_status,
    pik3ca_variant = PIK3CA_variant,
    germline_testing_done = BRCA_testing,
    brca_status = BRCA_status,
    brca1_variant = BRCA1_variant,
    brca2_variant = BRCA2_variant
  ) %>%
  mutate(
    germline_testing_done = case_when(
      str_to_lower(trimws(as.character(germline_testing_done))) == "yes" ~ "Yes",
      !is.na(germline_testing_done) &
        str_to_lower(trimws(as.character(germline_testing_done))) != "no" ~ "Yes",
      TRUE ~ "No"   # explicit "No" or NA → "No" (not tested)
    ),
    brca1_variant = if_else(brca_status == "negative", NA_character_, brca1_variant),
    brca2_variant = if_else(brca_status == "negative", NA_character_, brca2_variant)
  )

input_df <- input_df %>%
  left_join(genetics_df, by = "gsID")

write_xlsx(input_df, GT_NON_DRUG_CLEANED)
