library(tidyverse)
library(readxl)
library(writexl)

source("../config.R")

input_df <- read_excel(CRC_DATES_XLSX)

input_df <- input_df %>%
  select(gsID, Original_Dx_Date_CRC, Metastatic_Dx_Date_CRC) %>%
  mutate(version = "crc") %>%
  rename(
    date_initial_diagnosis          = Original_Dx_Date_CRC,
    date_first_metastatic_diagnosis = Metastatic_Dx_Date_CRC
  ) %>%
  mutate(
    date_initial_diagnosis          = as.Date(date_initial_diagnosis),
    date_first_metastatic_diagnosis = as.Date(date_first_metastatic_diagnosis)
  ) %>%
  # 1970-01-01 is a sentinel for missing dates (Excel numeric 0 coerced to date)
  mutate(
    date_initial_diagnosis          = if_else(date_initial_diagnosis          == as.Date("1970-01-01"), NA_Date_, date_initial_diagnosis),
    date_first_metastatic_diagnosis = if_else(date_first_metastatic_diagnosis == as.Date("1970-01-01"), NA_Date_, date_first_metastatic_diagnosis)
  )

write_xlsx(input_df, CRC_NON_DRUG_CLEANED)
