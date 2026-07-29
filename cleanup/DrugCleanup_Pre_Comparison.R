library(tidyverse)
library(readxl)
library(writexl)

source("../config.R")

# Load systemic therapy from each LLM and combine
all_llms <- bind_rows(
  map(ALL_LLMS, function(llm_id) {
    tx_path <- file.path(EXP_DIR, RUN_ID, llm_id,
                         sprintf("systemic_therapy_%s.xlsx", RUN_ID))
    if (!file.exists(tx_path)) return(tibble())
    read_xlsx(tx_path)
  })
)

#Cleaning up all_llms
all_llms <- all_llms %>%
  rename(gsID = patient_id)

all_llms <- all_llms %>%
  rename(version = llm)

all_llms <- all_llms %>%
  rename(Tx_LineN = line_number)

all_llms <- all_llms %>%
  rename(ReasonStop = stop_reason)

#If XX replace with 01 for start_date_exact and start_date_estimated and stop_date_exact and stop_date_estimated
all_llms <- all_llms %>%
  mutate(
    across(
      c(start_date_exact,
        start_date_estimated,
        stop_date_exact,
        stop_date_estimated),
      ~ str_replace_all(.x, "XX", "01")
    )
  )

#Now merge start exact and start estimated into a single Start col. if both exist, preserve start exact
all_llms <- all_llms %>%
  mutate(
    Start = coalesce(start_date_exact, start_date_estimated),
    Stop = coalesce(stop_date_exact, stop_date_estimated)
  )

#Change Start/Stop from chr to date
all_llms <- all_llms %>%
  mutate(
    Start = as.Date(Start, format = "%Y-%m-%d"),
    Stop = as.Date(Stop, format = "%Y-%m-%d")
  )

#Show all the responses for reason stop
table(all_llms$ReasonStop)

#Summate Progression Events
Progression <- c(
  "hospice enrollment",
  "disease progression",
  "progression",
  "patient death"
)

Other <- c(
  "clinician choice",
  "completed planned adjuvant course",
  "completed planned course",
  "completed planned neoadjuvant course",
  "ongoing",
  "patient preference",
  "toxicity",
  "trial protocol ended",
  "unknown"
)

all_llms <- all_llms %>%
  mutate(
    ReasonStop = case_when(
      str_detect(str_to_lower(ReasonStop),
                 str_c(Progression, collapse = "|")) ~ "Progression",
      TRUE ~ "Other"
    )
  )

#Now merge on TreatmentIntent
#Date files 
llm_dates <- read_xlsx(LLM_NON_DRUG_CLEANED) %>%
  mutate(across(where(lubridate::is.POSIXct), lubridate::as_date))

#make sure formatting is good
llm_dates <- llm_dates %>%
  select(gsID, version, date_initial_diagnosis, date_first_locoregional_recurrence,
         date_first_metastatic_diagnosis, birth_date) %>%
  mutate(date_first_metastatic_diagnosis = as.Date(date_first_metastatic_diagnosis, format = "%Y-%m-%d"))

#Now for each gsID in all_llms we create a new col called TreatmentIntent
#If Start is after date_metastatic_diagnosis then Palliative. If before or no date_metastatic_diagnosis exists then Curative
all_llms <- all_llms %>%
  left_join(llm_dates, by = c("gsID", "version")) %>%
  mutate(
    TreatmentIntent = case_when(
      is.na(Start)                                                                       ~ NA_character_,
      !is.na(date_first_metastatic_diagnosis) & Start >= date_first_metastatic_diagnosis ~ "Palliative",
      TRUE                                                                               ~ "Curative"
    )
  ) %>%
  select(-date_first_metastatic_diagnosis)

table(all_llms$TreatmentIntent, useNA = "ifany")

#Set up the ground truth
read_gt_folder <- function(dir_path) {
  
  files <- list.files(
    dir_path,
    pattern = "^gs.*\\.xlsx$",
    full.names = TRUE
  )
  
  files <- files[!grepl("(^|/|\\\\)~\\$", files)]  # remove Excel temp files
  
  bind_rows(lapply(files, function(f) {
    read_excel(f, col_types = "text") %>%          # read every column as text to avoid readxl date coercion
      mutate(
        source_file = basename(f),
        gsID        = str_extract(basename(f), "gs\\d+"),
        version     = str_extract(basename(f), "v\\d+"),   # derive from filename, not cell
        across(where(is.character), trimws),               # strip whitespace from all text cols
        across(any_of(c("Start", "Stop", "lastnotedate")), parse_gt_date)  # parse dates robustly on the way in
      )
  }))
}

# Create separate v1 and v2 dataframes
gt_v1 <- read_gt_folder(GT_DIR_V1)
gt_v2 <- read_gt_folder(GT_DIR_V2)
#Put em togther
gt_all <- bind_rows(gt_v1, gt_v2)

#keep relevant cols
gt_all <- gt_all %>%
  select(gsID, version, Start, Stop, ReasonStop, TreatmentIntent, Tx_LineN, regimen_standardized, regimen_bag_of_drugs)

#Fix caps 
table(all_llms$TreatmentIntent)
table(gt_all$TreatmentIntent)

gt_all <- gt_all %>%
  mutate(
    TreatmentIntent = case_when(
      str_to_lower(TreatmentIntent) == "palliative" ~ "Palliative",
      TRUE ~ "Curative"
    )
  )

table(gt_all$ReasonStop)

#Clean up ReasonStop in gt_all to match all_llms
gt_all <- gt_all %>%
  mutate(
    rs = str_to_lower(ReasonStop),
    ReasonStop = case_when(
      rs %in% c("death", "deceased", "hospice") ~ "Progression",
      str_detect(rs, "progression") ~ "Progression",
      TRUE ~ "Other"
    )
  ) %>%
  select(-rs)

# Build patient-level censor date from GT non-drug data (last_note_date per patient).
# LLMs use the same censor date as the GT patient they were evaluated on.
gt_censor <- read_excel(GT_NON_DRUG_XLSX) %>%
  mutate(censor_date = as.Date(last_note_date)) %>%
  select(gsID, censor_date) %>%
  filter(!is.na(censor_date))

#Cut down to key cols for matching
all_llms_final <- all_llms %>%
  left_join(gt_censor, by = "gsID") %>%
  select(gsID, version, Start, Stop, ReasonStop, TreatmentIntent, Tx_LineN, regimen_standardized, regimen_bag_of_drugs, censor_date)

gt_all_final <- gt_all %>%
  left_join(gt_censor, by = "gsID") %>%
  select(gsID, version, Start, Stop, ReasonStop, TreatmentIntent, Tx_LineN, regimen_standardized, regimen_bag_of_drugs, censor_date)

# CRC (Research Coordinator)

crc_all <- read_xlsx(CRC_LINES_XLSX) %>%
  mutate(across(where(lubridate::is.POSIXct), lubridate::as_date))

# Rename to shared schema — adjust if your CRC file uses different column names
crc_all <- crc_all %>%
  rename(Tx_LineN = TreatmentLine) %>%
  select(gsID, Tx_LineN, Start, Stop, ReasonStop, regimen_standardized, regimen_bag_of_drugs)

# XX → 01 for any partial dates (same pattern as LLM block)
crc_all <- crc_all %>%
  mutate(
    across(
      c(Start, Stop),
      ~ str_replace_all(.x, "XX", "01")
    )
  )

crc_all <- crc_all %>%
  mutate(
    Start = as.Date(Start, format = "%Y-%m-%d"),
    Stop  = as.Date(Stop,  format = "%Y-%m-%d")
  ) %>%
  # 1970-01-01 is a sentinel for missing dates (Excel numeric 0 coerced to date)
  mutate(
    Start = if_else(Start == as.Date("1970-01-01"), NA_Date_, Start),
    Stop  = if_else(Stop  == as.Date("1970-01-01"), NA_Date_, Stop)
  )

#add version col w crc
crc_all <- crc_all %>%
  mutate(version = "crc")

# TreatmentIntent from CRC diagnosis dates
crc_dates <- read_xlsx(CRC_DATES_XLSX) %>%
  mutate(across(where(lubridate::is.POSIXct), lubridate::as_date))

crc_dates <- crc_dates %>%
  mutate(
    date_first_metastatic_diagnosis = as.Date(Metastatic_Dx_Date_CRC, format = "%Y-%m-%d"),
    date_first_metastatic_diagnosis = if_else(
      date_first_metastatic_diagnosis == as.Date("1970-01-01"), NA_Date_, date_first_metastatic_diagnosis
    )
  )

crc_all <- crc_all %>%
  left_join(crc_dates %>% select(gsID, date_first_metastatic_diagnosis), by = "gsID") %>%
  mutate(
    TreatmentIntent = case_when(
      is.na(Start)                                                                       ~ NA_character_,
      !is.na(date_first_metastatic_diagnosis) & Start >= date_first_metastatic_diagnosis ~ "Palliative",
      TRUE                                                                               ~ "Curative"
    )
  ) %>%
  select(-date_first_metastatic_diagnosis)

table(crc_all$TreatmentIntent, useNA = "ifany")

# Fixed censor date for CRC — defined in config.R
crc_all <- crc_all %>%
  mutate(censor_date = CRC_CENSOR_DATE)

table(crc_all$ReasonStop)

#Clean up ReasonStop in gt_all to match all_llms
crc_all <- crc_all %>%
  mutate(
    rs = str_to_lower(ReasonStop),
    ReasonStop = case_when(
      rs %in% c("death", "deceased", "disease porgression", "disease progression", "disease progression or death", "mixed response?", "progression", "progression of disease", "toxicity + disease progression") ~ "Progression",
      str_detect(rs, "progression") ~ "Progression",
      TRUE ~ "Other"
    )
  ) %>%
  select(-rs)

# Ensure Start/Stop are Date and Tx_LineN is numeric in all three frames before combining
all_llms_final <- all_llms_final %>%
  mutate(Start = as.Date(Start), Stop = as.Date(Stop), Tx_LineN = as.numeric(Tx_LineN))

gt_all_final <- gt_all_final %>%
  mutate(Start = as.Date(Start), Stop = as.Date(Stop), Tx_LineN = as.numeric(Tx_LineN))

crc_all <- crc_all %>%
  mutate(Start = as.Date(Start), Stop = as.Date(Stop), Tx_LineN = as.numeric(Tx_LineN)) %>%
  select(gsID, version, Start, Stop, ReasonStop, TreatmentIntent, Tx_LineN, regimen_standardized, regimen_bag_of_drugs, censor_date)

final_drug_df <- bind_rows(all_llms_final, gt_all_final, crc_all) %>%
  filter(regimen_standardized != "supportive_only_line")

writexl::write_xlsx(final_drug_df, FINAL_DRUG_DF)

dir.create(LOGS_DIR, showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(LOGS_DIR, "session_info_drug_cleanup.txt"))
