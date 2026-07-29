library(tidyverse)
library(lubridate)
library(readxl)
library(survival)
library(writexl)
library(boot)

source("config.R")

set.seed(RANDOM_SEED)

dir.create(TBL_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(LAB_MTG_DIR, showWarnings = FALSE, recursive = TRUE)

# Helpers
# Extract HR and HER2 status from a combined hr_her2_status string column.
# Returns a two-column tibble (hr, her2).
extract_hr_her2 <- function(subtype_col) {
  tibble(
    hr = case_when(
      str_detect(coalesce(subtype_col, ""), "HR\\+")  ~ "HR+",
      str_detect(coalesce(subtype_col, ""), "HR-")    ~ "HR-",
      TRUE ~ NA_character_),
    her2 = case_when(
      str_detect(coalesce(subtype_col, ""), "HER2\\+") ~ "HER2+",
      str_detect(coalesce(subtype_col, ""), "HER2-")   ~ "HER2-",
      TRUE ~ NA_character_)
  )
}

# DATA LOADING  (reads from cleaned non_drug_long.xlsx)

# Read cleaned non-drug long form
non_drug <- read_excel(FINAL_NON_DRUG_LONG) %>%
  mutate(across(where(lubridate::is.POSIXct), lubridate::as_date))

# Cleaned intermediates
llm_cleaned <- read_xlsx(LLM_NON_DRUG_CLEANED)

# Ground truth (from cleaned intermediate)
gt_raw <- read_xlsx(GT_NON_DRUG_CLEANED) %>%
  mutate(across(where(lubridate::is.POSIXct), lubridate::as_date))

gt <- non_drug %>%
  filter(version == "gt") %>%
  transmute(
    gsID,
    gt_birth     = as.Date(birth_date),
    gt_init_dx   = as.Date(date_initial_diagnosis),
    gt_lrr_dx    = as.Date(date_first_locoregional_recurrence),
    gt_met_dx    = as.Date(date_first_metastatic_diagnosis),
    gt_stage     = as.character(stage_at_diagnosis),
    gt_subtype   = hr_her2_status,
    gt_hr   = extract_hr_her2(hr_her2_status)$hr,
    gt_her2 = extract_hr_her2(hr_her2_status)$her2,
    gt_primary_onc = tolower(trimws(as.character(primary_oncologist_shc))),
    # Genetics: "mutated" → "positive", everything else (negative/NA/unknown) → "negative"
    # Mirrors binarize_mutation() in manuscript_figures.R so accuracy is consistent
    gt_brca_mut = if_else(!is.na(brca_status)   & tolower(brca_status)   == "mutated",
                          "positive", "negative"),
    gt_pik3ca   = if_else(!is.na(pik3ca_status) & tolower(pik3ca_status) == "mutated",
                          "positive", "negative"),
    gt_esr1     = if_else(!is.na(esr1_status)   & tolower(esr1_status)   == "mutated",
                          "positive", "negative")
  ) %>%
  left_join(
    gt_raw %>% transmute(
      gsID         = as.character(gsID),
      gt_death     = parse_date(death_date_clinician),
      gt_deceased  = tolower(as.character(deceased_clinician)) %in% c("dead", "true"),
      gt_last_note = parse_date(last_note_date)
    ),
    by = "gsID"
  )

# LLM data from non_drug_long (supplemented with death from raw)

all_llm_data <- map(setNames(ALL_LLMS, ALL_LLMS), function(llm_id) {
  d <- non_drug %>%
    filter(version == llm_id) %>%
    transmute(
      gsID,
      birth    = as.Date(birth_date),
      init_dx  = as.Date(date_initial_diagnosis),
      lrr_dx   = as.Date(date_first_locoregional_recurrence),
      met_dx   = as.Date(date_first_metastatic_diagnosis),
      stage    = as.character(stage_at_diagnosis),
      subtype  = hr_her2_status,
      cap_hr   = extract_hr_her2(hr_her2_status)$hr,
      cap_her2 = extract_hr_her2(hr_her2_status)$her2,
      # Genetics: mirrors binarize_mutation() in manuscript_figures.R
      brca_mut = if_else(!is.na(brca_status)   & tolower(brca_status)   == "mutated",
                         "positive", "negative"),
      pik3ca   = if_else(!is.na(pik3ca_status) & tolower(pik3ca_status) == "mutated",
                         "positive", "negative"),
      esr1     = if_else(!is.na(esr1_status)   & tolower(esr1_status)   == "mutated",
                         "positive", "negative")
    )

  # Death date from intermediate cleaned file (not raw extractor)
  llm_death <- llm_cleaned %>%
    filter(version == llm_id) %>%
    transmute(gsID = as.character(gsID),
              death_date = as.Date(date_of_death))
  d <- d %>% left_join(llm_death, by = "gsID")
  if (!"death_date" %in% names(d)) d$death_date <- as.Date(NA)

  # Primary oncologist: exact last-name match from non_drug_long
  d$primary_onc <- tolower(trimws(as.character(
    non_drug %>% filter(version == llm_id) %>% pull(primary_oncologist_shc)
  )))

  d
})

# TABLE 2 (main text): Patient cohort characteristics
# Two-column format: Characteristic | Value  (GT ground truth only)
# Matches Word doc jr_pip_v1.0.docx Table 2

# Helpers
fmt_n_pct <- function(n, total) {
  if (is.na(n) || is.na(total) || total == 0) return("\u2014")
  sprintf("%d (%d%%)", as.integer(n), round(as.numeric(n) / as.numeric(total) * 100))
}

fmt_med_iqr_t2 <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("\u2014")
  sprintf(paste0("%.", digits, "f (%.", digits, "f\u2013%.", digits, "f)"),
          median(x), quantile(x, 0.25), quantile(x, 0.75))
}

# Ground truth data
N <- nrow(gt)  # 100 patients

# Stage 1-3 subset for recurrence denominators
s13_gsids <- gt %>% filter(gt_stage %in% c("1", "2", "3")) %>% pull(gsID)
s13       <- gt %>% filter(gsID %in% s13_gsids)
n_s13     <- length(s13_gsids)

# Follow-up time (diagnosis → death if deceased, else last note)
fu_yrs_vec <- gt %>%
  transmute(
    fu_end  = case_when(
      coalesce(gt_deceased, FALSE) & !is.na(gt_death) ~ gt_death,
      TRUE                                            ~ gt_last_note),
    fu_yrs  = as.numeric(difftime(fu_end, gt_init_dx, units = "days")) / 365.25
  ) %>%
  filter(!is.na(fu_yrs), fu_yrs >= 0) %>%
  pull(fu_yrs)

# Age at diagnosis from raw GT
age_dx_vec <- gt_raw %>%
  transmute(age_dx = as.numeric(difftime(parse_date(date_initial_diagnosis),
                                          parse_date(birth_date),
                                          units = "days")) / 365.25) %>%
  pull(age_dx) %>% na.omit()

# Subtype counts using pattern matching (robust to em-dash vs hyphen variants)
get_sub_n <- function(pat) sum(str_detect(coalesce(gt$gt_subtype, ""), pat), na.rm = TRUE)

# Drug/line stats from GT v1
drug_df_t2 <- read_xlsx(FINAL_DRUG_DF)
gt_v1_lines  <- drug_df_t2 %>%
  filter(version == "v1") %>%
  mutate(gsID = as.character(gsID))
total_lines_v1  <- nrow(gt_v1_lines)
lines_per_pat_v1 <- gt_v1_lines %>% count(gsID) %>% pull(n)
drug_strs_v1 <- gt_v1_lines %>%
  filter(!is.na(regimen_bag_of_drugs), nchar(trimws(regimen_bag_of_drugs)) > 0) %>%
  pull(regimen_bag_of_drugs)
all_drugs_v1 <- unique(trimws(tolower(
  unlist(strsplit(paste(drug_strs_v1, collapse = "+"), "[+,;/]"))
)))
all_drugs_v1 <- all_drugs_v1[nchar(all_drugs_v1) > 0]
n_drugs_v1 <- length(all_drugs_v1)

# Distinct primary treating oncologists
n_oncs_distinct <- non_drug %>%
  filter(version == "gt",
         !is.na(primary_oncologist_shc),
         trimws(primary_oncologist_shc) != "") %>%
  distinct(primary_oncologist_shc) %>%
  nrow()

# Assemble rows — exact row order and labels matching Word doc (35 rows)
tbl2_rows <- list()
add_char  <- function(char, val) {
  tbl2_rows <<- c(tbl2_rows, list(tibble(Characteristic = char, Value = val)))
}
section_hdr <- function(title) add_char(title, "")

# Demographic and disease characteristics (n=100)
section_hdr("Demographic and disease characteristics (n=100)")
add_char("Age at diagnosis, median (IQR), y", fmt_med_iqr_t2(age_dx_vec))

# Stage (sub-header row, then bare numbers)
add_char("Stage at diagnosis, No. (%)", "")
for (stg in c("1", "2", "3", "4")) {
  add_char(stg, fmt_n_pct(sum(coalesce(gt$gt_stage == stg, FALSE)), N))
}

# Subtype (sub-header row, then values)
add_char("Tumor subtype at diagnosis, No. (%)", "")
add_char("HR+/HER2\u2212",  fmt_n_pct(get_sub_n("HR\\+.*HER2-"),  N))
add_char("HR+/HER2+",       fmt_n_pct(get_sub_n("HR\\+.*HER2\\+"), N))
add_char("HR\u2212/HER2+",  fmt_n_pct(get_sub_n("HR-.*HER2\\+"),  N))
add_char("Triple-negative", fmt_n_pct(get_sub_n("HR-.*HER2-"),     N))

# Genetics (n=100 denominator)
add_char("Germline BRCA1/2 mutation present, No. (%)",
         fmt_n_pct(sum(coalesce(gt$gt_brca_mut == "positive", FALSE)), N))
add_char("Somatic PIK3CA mutation present, No. (%)",
         fmt_n_pct(sum(coalesce(gt$gt_pik3ca == "positive", FALSE)), N))
add_char("Somatic ESR1 mutation present, No. (%)",
         fmt_n_pct(sum(coalesce(gt$gt_esr1 == "positive", FALSE)), N))

# Recurrence Events (n = stage 1-3 count)
section_hdr(sprintf("Recurrence Events (n=%d)", n_s13))
n_any_recur <- sum(!is.na(s13$gt_lrr_dx) | !is.na(s13$gt_met_dx))
add_char("Any recurrence (locoregional or distant), No. (%)",
         fmt_n_pct(n_any_recur, n_s13))
add_char("Locoregional recurrence, No. (%)",
         fmt_n_pct(sum(!is.na(s13$gt_lrr_dx)), n_s13))
add_char("Distant metastatic recurrence, No. (%)",
         fmt_n_pct(sum(!is.na(s13$gt_met_dx)), n_s13))

# Survival status
section_hdr("Survival status")
# Death count masked ("not ready for prime time") — update when ready
add_char("Documented death by censor date, No. (%)", "***")
add_char("Documented follow-up from diagnosis, median (IQR), y",
         fmt_med_iqr_t2(fu_yrs_vec))

# Systemic therapy complexity
section_hdr("Systemic therapy complexity")
add_char("Total lines of systemic therapy abstracted, No.",
         if (!is.na(total_lines_v1)) as.character(total_lines_v1) else "***")
add_char("Lines of systemic therapy per patient, median (IQR)",
         local({
           v <- lines_per_pat_v1[!is.na(lines_per_pat_v1)]
           if (length(v) == 0) return("\u2014")
           # Use %g to suppress trailing zeros (e.g., "4.75" not "4.75000")
           fmtq <- function(x) sprintf("%.4g", x)
           sprintf("%s (%s\u2013%s)", fmtq(median(v)), fmtq(quantile(v, 0.25)), fmtq(quantile(v, 0.75)))
         }))
add_char("Total anti-cancer drugs abstracted, No.",
         if (!is.na(n_drugs_v1)) as.character(n_drugs_v1) else "***")

# Chart complexity — computed from raw chart files
{
  chart_base <- CHART_RAW_DIR
  pat_dirs   <- list.dirs(chart_base, recursive = FALSE, full.names = TRUE)
  pat_dirs   <- pat_dirs[grepl("/gs\\d+$", pat_dirs)]

  chart_mb <- vapply(pat_dirs, function(d) {
    fls <- list.files(d, recursive = TRUE, full.names = TRUE)
    sum(file.info(fls)$size, na.rm = TRUE) / 1024 / 1024
  }, numeric(1))

  chart_tokens_M <- chart_mb * 1024 * 1024 / 4 / 1e6  # ~4 bytes per token
  chart_pages    <- chart_tokens_M * 1e6 / 750          # ~750 tokens per page

  count_csv_rows <- function(csv_name) {
    vapply(pat_dirs, function(d) {
      fp <- file.path(d, csv_name)
      if (!file.exists(fp)) return(0L)
      nrow(tryCatch(read.csv(fp, nrows = Inf, header = TRUE), error = function(e) data.frame()))
    }, integer(1))
  }

  clinical_notes_n <- count_csv_rows("clinical_notes.csv")
  pathology_n      <- count_csv_rows("pathology_report.csv")
  med_admin_n      <- count_csv_rows("med_admin.csv")

  fmt_tokens <- function(x) {
    sprintf("%.1f (%.1f\u2013%.1f) million",
            median(x), quantile(x, 0.25), quantile(x, 0.75))
  }
  fmt_int_comma <- function(n) format(round(n), big.mark = ",", scientific = FALSE)
  fmt_med_iqr_comma <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return("\u2014")
    sprintf("%s (%s\u2013%s)",
            fmt_int_comma(median(x)),
            fmt_int_comma(quantile(x, 0.25)),
            fmt_int_comma(quantile(x, 0.75)))
  }
}

section_hdr("Chart complexity")
# Chart size/tokens/pages: show median only, IQR masked pending data review
add_char("Estimated chart size, median (IQR), MB",
         if (length(chart_mb) == 0) "\u2014"
         else sprintf("%.1f (***)", median(chart_mb, na.rm = TRUE)))
add_char("Estimated tokens per chart, median (IQR)",
         if (length(chart_tokens_M) == 0) "\u2014"
         else sprintf("%.1f million (***)", median(chart_tokens_M, na.rm = TRUE)))
add_char("Estimated pages per chart, median (IQR)",
         if (length(chart_pages) == 0) "\u2014"
         else sprintf("%s (***)", fmt_int_comma(median(chart_pages, na.rm = TRUE))))
add_char("Clinical notes per patient, median (IQR), No.",
         fmt_med_iqr_comma(clinical_notes_n))
add_char("Pathology documents per patient, median (IQR), No.",
         fmt_med_iqr_comma(pathology_n))
# Medication admin records masked pending data review
add_char("Medication administration records per patient, median (IQR), No.", "***")
add_char("Distinct primary treating oncologists across cohort, No.",
         as.character(n_oncs_distinct))

table1 <- bind_rows(tbl2_rows)

# TABLE S2: Accuracy across variables by ingestion strategy

compute_accuracy_row <- function(exp_df, exp_name, gt_df, tol) {
  m <- inner_join(exp_df, gt_df, by = "gsID")
  rows <- list()

  # Date variables
  date_vars <- list(
    "Date of birth"                   = c("birth",      "gt_birth"),
    "Date of initial diagnosis"       = c("init_dx",    "gt_init_dx"),
    "Date of locoregional recurrence" = c("lrr_dx",     "gt_lrr_dx"),
    "Date of metastatic diagnosis"    = c("met_dx",     "gt_met_dx")
  )
  for (vname in names(date_vars)) {
    cols <- date_vars[[vname]]
    if (!cols[1] %in% names(m)) next
    n_gt        <- sum(!is.na(m[[cols[2]]]))
    vec <- date_agree_vec(m[[cols[1]]], m[[cols[2]]], tol)
    if (length(vec) == 0) next
    ci <- boot_ci(vec, mean)
    rows <- c(rows, list(tibble(
      Variable   = vname, Category = "Dates", LLM = exp_name,
      Tolerance  = paste0("\u00b1", tol, " days"),
      N_GT       = n_gt,
      N          = length(vec),
      Extraction_Rate = if (n_gt > 0) sprintf("%.0f%%", length(vec) / n_gt * 100) else "\u2014",
      Agreement  = fmt_pct_ci(ci["est"], ci["lo"], ci["hi"])
    )))
  }

  # Categorical variables
  cat_vars <- list(
    "Initial cancer stage" = c("stage",   "gt_stage"),
    "HR/HER2 subtype"     = c("subtype", "gt_subtype"),
    "HR status"           = c("cap_hr",  "gt_hr"),
    "HER2 status"         = c("cap_her2", "gt_her2"),
    "Primary oncologist"  = c("primary_onc", "gt_primary_onc"),
    "Germline BRCA"       = c("brca_mut", "gt_brca_mut"),
    "Somatic PIK3CA"      = c("pik3ca",   "gt_pik3ca"),
    "Somatic ESR1"        = c("esr1",     "gt_esr1")
  )

  cat_categories <- c(
    "Initial cancer stage" = "Staging",
    "HR/HER2 subtype"     = "Receptor Status",
    "HR status"           = "Receptor Status",
    "HER2 status"         = "Receptor Status",
    "Primary oncologist"  = "Clinical",
    "Germline BRCA"       = "Genetics",
    "Somatic PIK3CA"      = "Genetics",
    "Somatic ESR1"        = "Genetics"
  )

  na_wrong_vars <- c("Initial cancer stage", "HR status", "HER2 status")
  for (vname in names(cat_vars)) {
    cols <- cat_vars[[vname]]
    if (!cols[1] %in% names(m) || !cols[2] %in% names(m)) next
    n_gt         <- sum(!is.na(m[[cols[2]]]))
    use_strict   <- vname %in% na_wrong_vars
    vec          <- cat_agree_vec(m[[cols[1]]], m[[cols[2]]], na_as_wrong = use_strict)
    if (length(vec) == 0) next
    # Extraction rate: patients where LLM gave a value (among GT-non-NA)
    n_extracted  <- sum(!is.na(m[[cols[1]]]) & !is.na(m[[cols[2]]]))
    ci <- boot_ci(vec, mean)
    rows <- c(rows, list(tibble(
      Variable   = vname, Category = cat_categories[vname], LLM = exp_name,
      Tolerance  = "Exact match",
      N_GT       = n_gt,
      N          = length(vec),
      Extraction_Rate = if (n_gt > 0) sprintf("%.0f%%", n_extracted / n_gt * 100) else "\u2014",
      Agreement  = fmt_pct_ci(ci["est"], ci["lo"], ci["hi"])
    )))
  }

  # Death date documented (agreement on presence/absence of exact death date)
  if ("death_date" %in% names(m) && "gt_death" %in% names(m)) {
    llm_has <- !is.na(m$death_date)
    gt_has  <- !is.na(m$gt_death)
    vec     <- as.integer(llm_has == gt_has)
    if (length(vec) > 0) {
      ci <- boot_ci(vec, mean)
      rows <- c(rows, list(tibble(
        Variable = "Exact date of death documented",
        Category = "Dates", LLM = exp_name,
        Tolerance       = "Exact match",
        N_GT            = sum(gt_has),
        N               = length(vec),
        Extraction_Rate = "100%",
        Agreement       = fmt_pct_ci(ci["est"], ci["lo"], ci["hi"])
      )))
    }
  }

  # Recurrence occurred (binary: did patient have any LRR or metastatic recurrence?)
  if (all(c("lrr_dx", "met_dx") %in% names(m)) &&
      all(c("gt_lrr_dx", "gt_met_dx") %in% names(m))) {
    # de novo Stage 4 is metastatic at presentation, not a recurrence
    keep      <- !coalesce(m$gt_stage == "4", FALSE)
    llm_recur <- (!is.na(m$lrr_dx) | !is.na(m$met_dx))[keep]
    gt_recur  <- (!is.na(m$gt_lrr_dx) | !is.na(m$gt_met_dx))[keep]
    vec       <- as.integer(llm_recur == gt_recur)
    if (length(vec) > 0) {
      ci <- boot_ci(vec, mean)
      rows <- c(rows, list(tibble(
        Variable = "Recurrence occurred",
        Category = "Clinical", LLM = exp_name,
        Tolerance       = "Exact match",
        N_GT            = sum(gt_recur),
        N               = length(vec),
        Extraction_Rate = "100%",
        Agreement       = fmt_pct_ci(ci["est"], ci["lo"], ci["hi"])
      )))
    }
  }

  bind_rows(rows)
}

# Compute for ALL LLMs at primary tolerance
tbl2_all <- map2_dfr(all_llm_data, ALL_LLMS, function(d, llm_id) {
  compute_accuracy_row(d, LLM_LABELS[llm_id], gt, TOL_PRIMARY)
})

# Add CRC — only initial diagnosis and metastatic diagnosis dates were extracted.
# The coordinator abstracted in 2022, so the expert reference is censored at
# CRC_CENSOR_DATE here, as in every other coordinator comparison.
crc_nondrug <- read_xlsx(CRC_NON_DRUG_CLEANED)
crc_for_tbl2 <- crc_nondrug %>%
  transmute(
    gsID    = as.character(gsID),
    init_dx = as.Date(date_initial_diagnosis),
    met_dx  = as.Date(date_first_metastatic_diagnosis)
  )
gt_crc_censored <- gt %>%
  mutate(gt_met_dx = if_else(!is.na(gt_met_dx) & gt_met_dx >= CRC_CENSOR_DATE,
                             as.Date(NA), gt_met_dx))
tbl2_all <- bind_rows(
  tbl2_all,
  compute_accuracy_row(crc_for_tbl2, "Research Coordinator", gt_crc_censored, TOL_PRIMARY)
)

# Sensitivity analyses for all LLMs
tbl2_sens <- map_dfr(TOL_SENS, function(tol) {
  map2_dfr(all_llm_data, ALL_LLMS, function(d, llm_id) {
    compute_accuracy_row(d, LLM_LABELS[llm_id], gt, tol) %>%
      filter(Category == "Dates")
  })
})

table2 <- bind_rows(tbl2_all, tbl2_sens) %>%
  # Include "Research Coordinator" as a valid factor level
  mutate(LLM = factor(LLM, levels = c(as.character(LLM_LABELS), "Research Coordinator"))) %>%
  arrange(Variable, Tolerance, LLM)

# Wide format for Appendix Table 2 (non-drug rows only)
# Drug rows are added from DrugComparison output file in the export section.
var_to_row <- c(
  "Exact date of death documented"   = "Exact date of death documented",
  "Date of birth"                    = "Date of birth (\u00b190 d)",
  "Date of initial diagnosis"        = "Date of diagnosis (\u00b190 d)",
  "Date of locoregional recurrence"  = "Date of first locoregional recurrence (\u00b190 d)",
  "Date of metastatic diagnosis"     = "Date of first metastatic diagnosis (\u00b190 d)",
  "Recurrence occurred"              = "Recurrence occurred (yes/no)",
  "HR/HER2 subtype"                  = "Original HR/HER2 subtype",
  "HR status"                        = "Hormone receptor status",
  "HER2 status"                      = "HER2 status",
  "Initial cancer stage"             = "Original stage",
  "Primary oncologist"               = "Primary treating oncologist",
  "Germline BRCA"                    = "Germline BRCA1/2 mutation identified",
  "Somatic PIK3CA"                   = "Somatic PIK3CA mutation present",
  "Somatic ESR1"                     = "Somatic ESR1 mutation present"
)

src_to_col2 <- c(
  "GPT-5"                = "GPT-5",
  "Gemini 2.5"           = "Gemini 2.5 Pro",
  "GPT-4o"               = "GPT-4o",
  "DeepSeek-R1"          = "DeepSeek-R1",
  "Research Coordinator" = "Research coordinator"
)

all_col_names <- c("Second oncologist", "GPT-5", "Gemini 2.5 Pro",
                   "GPT-4o", "DeepSeek-R1", "Research coordinator")

nondrug_wide <- tbl2_all %>%
  filter(
    Variable %in% names(var_to_row),
    # primary tolerance for date vars; exact match for categorical and vital/death
    (Category != "Dates") |
    (Variable == "Exact date of death documented" &
     Tolerance == "Exact match") |
    (Variable %in% c("Date of birth",
                     "Date of initial diagnosis",
                     "Date of locoregional recurrence",
                     "Date of metastatic diagnosis") &
     Tolerance == paste0("\u00b1", TOL_PRIMARY, " days"))
  ) %>%
  mutate(
    row_label = var_to_row[Variable],
    col_label = coalesce(src_to_col2[as.character(LLM)], as.character(LLM))
  ) %>%
  select(row_label, col_label, Agreement) %>%
  pivot_wider(names_from = col_label, values_from = Agreement) %>%
  # Ensure all source columns exist (fill missing with em-dash)
  {tbl <- .; for (cn in all_col_names) if (!cn %in% names(tbl)) tbl[[cn]] <- "\u2014"; tbl} %>%
  mutate(
    `Clinical domain / variable` = row_label,
    `Metric reported`            = "% agreement (95% CI)",
    `Second oncologist`          = "\u2014",
    across(all_of(all_col_names), ~ coalesce(.x, "\u2014"))
  ) %>%
  select(`Clinical domain / variable`, `Metric reported`,
         all_of(all_col_names))

# Row order matches Word doc (section headers will be added at export time)
row_order_nondrug <- c(
  "Date of birth (\u00b190 d)",
  "Exact date of death documented",
  "Date of diagnosis (\u00b190 d)",
  "Date of first locoregional recurrence (\u00b190 d)",
  "Date of first metastatic diagnosis (\u00b190 d)",
  "Recurrence occurred (yes/no)",
  "Original HR/HER2 subtype",
  "Hormone receptor status",
  "HER2 status",
  "Original stage",
  "Primary treating oncologist",
  "Germline BRCA1/2 mutation identified",
  "Somatic PIK3CA mutation present",
  "Somatic ESR1 mutation present"
)
nondrug_wide <- nondrug_wide %>%
  mutate(`Clinical domain / variable` = factor(`Clinical domain / variable`,
                                                levels = row_order_nondrug)) %>%
  arrange(`Clinical domain / variable`) %>%
  mutate(`Clinical domain / variable` = as.character(`Clinical domain / variable`))

# TABLE S5: Cox proportional hazards model estimates

stage4_pids <- gt %>% filter(gt_stage == "4") %>% pull(gsID)

# make_cox_df — OS event and time for LLM vs GT sources
# For GT rows: event = gt_deceased (from string "dead"/"not_dead" fields),
#              end_date = death_date if deceased+known, else last_note_date
# For LLM rows: event = !is.na(death_date),
#               end_date = date_of_death if deceased+known, else last_note_date (from GT)
# These Cox models censor conventionally: patients not known to have died are
# censored at last follow-up. The solid Kaplan-Meier curves in manuscript_figures.R
# instead apply the worst-case death-at-last-follow-up rule, so the two survival
# displays are different analyses and are reported as such.
make_cox_df <- function(d, label, is_gt = FALSE) {
  # Only join gt_last_note if not already present in d
  need_cols <- setdiff(c("gt_last_note"), names(d))
  m <- if (length(need_cols) > 0)
    inner_join(d, gt %>% select(gsID, all_of(need_cols)), by = "gsID")
  else d

  if (is_gt) {
    # GT — event from gt_deceased (already a logical from "dead"/"not_dead")
    need_surv <- setdiff(c("gt_deceased", "gt_death"), names(m))
    if (length(need_surv) > 0)
      m <- inner_join(m, gt %>% select(gsID, all_of(need_surv)), by = "gsID")
    m <- m %>%
      mutate(
        # coalesce to FALSE: unknown vital status → censored (alive) at last note
        gt_dec   = coalesce(gt_deceased, FALSE),
        # Like-for-like with the LLM Cox: a deceased patient without a documented
        # death date is censored at last note, exactly as the LLM treats a patient
        # with no extracted death date.
        event    = as.integer(gt_dec & !is.na(gt_death)),
        end_date = if_else(gt_dec & !is.na(gt_death), gt_death, gt_last_note)
      )
  } else {
    # LLM — use death_date presence from intermediate file
    m <- m %>%
      mutate(
        llm_deceased = !is.na(death_date),
        event        = as.integer(llm_deceased),
        end_date     = if_else(llm_deceased & !is.na(death_date),
                               death_date, gt_last_note)
      )
  }

  m %>%
    mutate(
      T_months = as.numeric(difftime(end_date, init_dx, units = "days")) / 30.4375,
      stage_iv = coalesce(stage == "4", FALSE),
      hr_neg   = coalesce(str_detect(subtype, "HR-"), FALSE)
    ) %>%
    filter(!is.na(T_months), T_months >= 0) %>%
    mutate(source_label = label)
}

# GT Cox data (OS) — pass gt with renamed columns so make_cox_df can use them
gt_cox <- gt %>%
  mutate(
    init_dx  = gt_init_dx,
    stage    = gt_stage,
    subtype  = gt_subtype
  ) %>%
  make_cox_df("Oncologist (GT)", is_gt = TRUE)

# All LLM Cox data (OS)
llm_cox_list <- map(setNames(ALL_LLMS, ALL_LLMS), function(llm_id) {
  d <- all_llm_data[[llm_id]]
  if (!all(c("init_dx", "stage", "subtype") %in% names(d))) return(NULL)
  make_cox_df(d, LLM_LABELS[llm_id], is_gt = FALSE)
})
llm_cox_list <- compact(llm_cox_list)

# GT Cox data (RFS — exclude stage IV, use GT lrr/met/death)
# gt already has gt_last_note, gt_deceased, gt_death from the left_join above
gt_rfs_cox <- gt %>%
  mutate(
    init_dx  = gt_init_dx,
    lrr_dx   = gt_lrr_dx,
    met_dx   = gt_met_dx,
    stage    = gt_stage,
    subtype  = gt_subtype
  ) %>%
  filter(!(gsID %in% stage4_pids)) %>%
  rowwise() %>%
  mutate(
    recur    = {
      cands <- na.omit(c(lrr_dx, met_dx))
      if (length(cands) > 0) min(cands) else as.Date(NA)
    },
    event    = as.integer(!is.na(recur) & !is.na(init_dx) & recur >= init_dx),
    end_date = if (event == 1L) recur
               else if (gt_deceased & !is.na(gt_death)) gt_death
               else gt_last_note,
    T_days   = as.numeric(difftime(end_date, init_dx, units = "days")),
    T_months  = T_days / 30.4375,
    stage_iv  = coalesce(gt_stage == "4", FALSE),
    stage_iii = coalesce(gt_stage == "3", FALSE),
    hr_neg    = coalesce(str_detect(gt_subtype, "HR-"), FALSE)
  ) %>%
  ungroup() %>%
  filter(!is.na(T_months), T_months >= 0) %>%
  mutate(source_label = "Oncologist (GT)")

# LLM Cox data (RFS — each LLM uses its own lrr/met dates; death from death_date)
llm_rfs_cox_list <- map(setNames(ALL_LLMS, ALL_LLMS), function(llm_id) {
  d <- all_llm_data[[llm_id]]
  if (!all(c("init_dx", "lrr_dx", "met_dx", "stage", "subtype") %in% names(d))) return(NULL)
  m <- d %>%
    inner_join(gt %>% select(gsID, gt_last_note), by = "gsID") %>%
    filter(!(coalesce(stage == "4", FALSE))) %>%   # each LLM's own Stage 4, not the expert list
    rowwise() %>%
    mutate(
      recur        = {
        cands <- na.omit(c(lrr_dx, met_dx))
        if (length(cands) > 0) min(cands) else as.Date(NA)
      },
      llm_deceased = !is.na(death_date),
      event        = as.integer(!is.na(recur) & !is.na(init_dx) & recur >= init_dx),
      end_date     = if (event == 1L) recur
                     else if (llm_deceased & !is.na(death_date)) death_date
                     else gt_last_note,
      T_days       = as.numeric(difftime(end_date, init_dx, units = "days")),
      T_months     = T_days / 30.4375,
      stage_iv     = coalesce(stage == "4", FALSE),
      stage_iii    = coalesce(stage == "3", FALSE),
      hr_neg       = coalesce(str_detect(subtype, "HR-"), FALSE)
    ) %>%
    ungroup() %>%
    filter(!is.na(T_months), T_months >= 0) %>%
    mutate(source_label = LLM_LABELS[llm_id])
  m
})
llm_rfs_cox_list <- compact(llm_rfs_cox_list)

# Covariates and model equations — OS and RFS
models <- list(
  list(eq = "OS ~ Stage IV",    var = "stage_iv",  outcome = "OS"),
  list(eq = "OS ~ HR-negative", var = "hr_neg",    outcome = "OS"),
  list(eq = "RFS ~ Stage III",  var = "stage_iii", outcome = "RFS"),
  list(eq = "RFS ~ HR-negative",var = "hr_neg",    outcome = "RFS")
)

s1_rows <- list()

# Run models for both OS and RFS sources
os_sources  <- c(list(`Oncologist (GT)` = gt_cox),
                 setNames(llm_cox_list,     LLM_LABELS[names(llm_cox_list)]))
rfs_sources <- c(list(`Oncologist (GT)` = gt_rfs_cox),
                 setNames(llm_rfs_cox_list, LLM_LABELS[names(llm_rfs_cox_list)]))

for (mod in models) {
  all_sources <- if (mod$outcome == "OS") os_sources else rfs_sources
  for (src_label in names(all_sources)) {
    df <- all_sources[[src_label]] %>%
      mutate(x = as.integer(.data[[mod$var]])) %>%
      filter(!is.na(x))
    if (nrow(df) < 10 || sum(df$event, na.rm = TRUE) < 5) next
    cx <- tryCatch(coxph(Surv(T_months, event) ~ x, data = df),
                   error = function(e) NULL)
    if (is.null(cx)) next
    s  <- summary(cx)
    ci <- tryCatch(exp(confint(cx)), error = function(e) NULL)
    if (is.null(ci)) next
    s1_rows <- c(s1_rows, list(tibble(
      Outcome      = mod$outcome,
      Equation     = mod$eq,
      Source       = src_label,
      N            = nrow(df),
      Events       = sum(df$event),
      # Raw log-HR and SE retained so Cochran Q uses full precision, not the
      # rounded HR/CI display strings.
      loghr        = unname(s$coefficients["x", "coef"]),
      se_loghr     = unname(s$coefficients["x", "se(coef)"]),
      `HR [95% CI]` = sprintf("%.2f\n[%.2f - %.2f]",
                              exp(coef(cx)["x"]), ci["x", 1], ci["x", 2]),
      `p-value`    = fmt_p(s$coefficients["x", "Pr(>|z|)"])
    )))
  }
}

if (length(s1_rows) > 0) {
  cox_tbl <- bind_rows(s1_rows)
} else {
  cox_tbl <- tibble(Equation = character())
  cat("  No Cox models could be fit\n")
}

# Wide format for Appendix Table 5
# Sources displayed in the table — Cochran Q is computed across these 3 only
s5_sources <- c("Oncologist (GT)", "GPT-5", "Gemini 2.5")

if (nrow(cox_tbl) > 0) {
  # Pairwise Cochran Q p-values, GT vs each LLM, from the full-precision log-HR and
  # its SE rather than the rounded display strings. Q assumes independent estimates
  # while these HRs come from the same patients, so it understates shared sampling
  # variation and these p-values are heterogeneity screens only.
  pairwise_q <- function(sub, src_a, src_b) {
    pair <- filter(sub, Source %in% c(src_a, src_b))
    ok   <- !is.na(pair$loghr) & !is.na(pair$se_loghr) & pair$se_loghr > 0
    if (sum(ok) < 2) return("\u2014")
    w  <- 1 / pair$se_loghr[ok]^2
    lh <- pair$loghr[ok]
    Q  <- sum((lh - weighted.mean(lh, w))^2 * w)
    fmt_p(pchisq(Q, sum(ok) - 1, lower.tail = FALSE))
  }

  cq_wide <- map_dfr(unique(cox_tbl$Equation), function(eq) {
    sub <- filter(cox_tbl, Equation == eq, Source %in% s5_sources)
    p_gpt5  <- pairwise_q(sub, "Oncologist (GT)", "GPT-5")
    p_gemini <- pairwise_q(sub, "Oncologist (GT)", "Gemini 2.5")
    tibble(Equation = eq,
           `Cochran Q P value` = paste0(p_gpt5, "\n", p_gemini))
  })

  eq_labels <- tibble(
    Equation  = c("OS ~ Stage IV", "OS ~ HR-negative",
                  "RFS ~ Stage III", "RFS ~ HR-negative"),
    Outcome   = c("Overall survival", "Overall survival",
                  "Invasive recurrence-free survival",
                  "Invasive recurrence-free survival"),
    Covariate = c("Stage IV vs stages I\u2013III",
                  "Hormone receptor-negative vs hormone receptor-positive",
                  "Stage III vs stages I\u2013II",
                  "Hormone receptor-negative vs hormone receptor-positive")
  )

  cox_tbl_wide <- cox_tbl %>%
    filter(Source %in% s5_sources) %>%
    mutate(Source = recode(Source,
      "Oncologist (GT)" = "Expert-derived HR (95% CI)",
      "GPT-5"           = "GPT-5-derived HR (95% CI)",
      "Gemini 2.5"      = "Gemini 2.5 Pro-derived HR (95% CI)")) %>%
    select(Equation, Source, `HR [95% CI]`) %>%
    pivot_wider(names_from = Source, values_from = `HR [95% CI]`) %>%
    left_join(eq_labels, by = "Equation") %>%
    left_join(cq_wide,   by = "Equation") %>%
    mutate(Equation = factor(Equation, levels = eq_labels$Equation)) %>%
    arrange(Equation) %>%
    select(Outcome, Covariate,
           `Expert-derived HR (95% CI)`,
           `GPT-5-derived HR (95% CI)`,
           `Gemini 2.5 Pro-derived HR (95% CI)`,
           `Cochran Q P value`)

} else {
  cox_tbl_wide <- tibble(
    Outcome = character(), Covariate = character(),
    `Expert-derived HR (95% CI)`        = character(),
    `GPT-5-derived HR (95% CI)`         = character(),
    `Gemini 2.5 Pro-derived HR (95% CI)`= character(),
    `Cochran Q P value`                 = character()
  )
}

# TABLE S8: Computational cost and runtime by LLM

# Read from the pre-processed token file written by NonDrugCleanup_VectorRank.R.
tok_path <- file.path(FINAL_DIR, "llm_token_usage.xlsx")
if (file.exists(tok_path)) {
  timing_df <- read_xlsx(tok_path) %>% filter(llm %in% ALL_LLMS)

  # Per-patient totals
  per_pat <- timing_df %>%
    group_by(llm, patient_id) %>%
    summarise(
      tokens   = sum(total_tokens),
      latency  = sum(latency_s),
      n_calls  = n(),
      .groups  = "drop"
    )

  # Extraction failures: patients with any task returning 0 completion tokens
  failures <- timing_df %>%
    filter(completion_tokens == 0) %>%
    group_by(llm) %>%
    summarise(`Extraction Failures` = n_distinct(patient_id), .groups = "drop")

  cost_tbl <- per_pat %>%
    group_by(llm) %>%
    summarise(
      `Mean Time/Patient (minutes)`                              = round(mean(latency) / 60, 1),
      `Mean Tokens/Patient`                                      = round(mean(tokens)),
      `Average API Calls Per Patient (Including Repair Attempts)` = round(mean(n_calls), 1),
      .groups = "drop"
    ) %>%
    left_join(failures, by = "llm") %>%
    mutate(
      `Extraction Failures`  = coalesce(`Extraction Failures`, 0L),
      LLM = coalesce(LLM_LABELS[llm], llm)
    ) %>%
    select(LLM,
           `Mean Time/Patient (minutes)`,
           `Mean Tokens/Patient`,
           `Average API Calls Per Patient (Including Repair Attempts)`,
           `Extraction Failures`) %>%
    arrange(LLM)
} else {
  cost_tbl <- tibble(LLM = character())
  cat("  llm_token_usage.xlsx not found; skipping Table S3\n")
}

# EXPORT: Wide-format appendix tables matching Word doc layout

cat("\n=== Assembling and exporting tables ===\n")

# Systemic-therapy line metrics, built before the Table S2 and Table 3 assembly so
# the same matcher feeds the main-text table and the Table S2 line rows.
source("line_metrics.R")
.line_metrics    <- build_line_metrics()
table3_lines     <- .line_metrics$table3
s2_line_rows     <- .line_metrics$s2_line_rows

# Supplementary Table S3 (line-level errors), from the same matcher as Table 4.
.s3col <- c(v2 = "Second oncologist", gpt5 = "GPT-5", gemini_2_5 = "Gemini 2.5 Pro",
            gpt4o = "GPT-4o", deepseek_r1 = "DeepSeek-R1", crc = "Research coordinator")
.s3 <- .line_metrics$per_patient %>% group_by(model) %>% summarise(
  n_pts = n(), total = sum(n_pred), matched = sum(n_match),
  missed = sum(n_ref) - sum(n_match), extra = sum(n_pred) - sum(n_match),
  miss_m = mean(n_ref - n_match), miss_s = sd(n_ref - n_match),
  extra_m = mean(n_pred - n_match), extra_s = sd(n_pred - n_match),
  n_exact = sum(n_pred == n_ref), n_nomiss = sum(n_match == n_ref),
  n_noextra = sum(n_match == n_pred), .groups = "drop")
.s3row <- function(metric, label, f, first = "—") {
  vals <- setNames(vapply(names(.s3col), function(m) {
    r <- .s3[.s3$model == m, ]; if (nrow(r) == 0) "—" else f(r) },
    character(1)), unname(.s3col))
  as_tibble(c(list(`Line-level error metric` = label, `Metric reported` = metric,
                   `First oncologist` = first), as.list(vals)))
}
s3_line_errors <- bind_rows(
  .s3row("No.", "Total predicted lines", function(r) as.character(r$total), first = "716"),
  .s3row("No.", "Correctly matched lines", function(r) as.character(r$matched)),
  .s3row("No.", "Missed lines", function(r) as.character(r$missed)),
  .s3row("No.", "Additional predicted lines", function(r) as.character(r$extra)),
  .s3row("Mean (SD)", "Missed lines per patient", function(r) sprintf("%.1f (%.1f)", r$miss_m, r$miss_s)),
  .s3row("Mean (SD)", "Additional predicted lines per patient", function(r) sprintf("%.1f (%.1f)", r$extra_m, r$extra_s)),
  .s3row("No. (%)", "Patients with exact line count match", function(r) sprintf("%d (%.1f%%)", r$n_exact, r$n_exact/r$n_pts*100)),
  .s3row("No. (%)", "Patients with no missed lines", function(r) sprintf("%d (%.1f%%)", r$n_nomiss, r$n_nomiss/r$n_pts*100)),
  .s3row("No. (%)", "Patients with no additional predicted lines", function(r) sprintf("%d (%.1f%%)", r$n_noextra, r$n_noextra/r$n_pts*100))
)
.s3path <- file.path(TBL_DIR, "TableS3_line_errors.xlsx")
.s3drug <- tryCatch(read_xlsx(.s3path, sheet = "drug_errors"), error = function(e) NULL)
write_xlsx(if (is.null(.s3drug)) list(line_errors = s3_line_errors)
           else list(line_errors = s3_line_errors, drug_errors = .s3drug), .s3path)
cat("  TableS3_line_errors.xlsx line_errors written from the line matcher\n")

# Assemble Appendix Table 2 (TableS2) wide format
# Non-drug rows (computed above) + drug rows (from DrugComparison output)
drug_rows_path <- file.path(INT_DIR, "TableS2_drug_rows.xlsx")

# Section header helper: blank row with section title in first column
make_header_row <- function(title) {
  tibble(`Clinical domain / variable` = title,
         `Metric reported`            = "",
         `Second oncologist`          = "",
         `GPT-5`                      = "",
         `Gemini 2.5 Pro`             = "",
         `GPT-4o`                     = "",
         `DeepSeek-R1`                = "",
         `Research coordinator`       = "")
}

# Section-labelled non-drug block
nondrug_sectioned <- bind_rows(
  make_header_row("Timing of key events"),
  filter(nondrug_wide, `Clinical domain / variable` %in%
           c("Date of birth (\u00b190 d)",
             "Exact date of death documented",
             "Date of diagnosis (\u00b190 d)",
             "Date of first locoregional recurrence (\u00b190 d)",
             "Date of first metastatic diagnosis (\u00b190 d)",
             "Recurrence occurred (yes/no)")),
  make_header_row("Clinical information"),
  filter(nondrug_wide, `Clinical domain / variable` %in%
           c("Original HR/HER2 subtype",
             "Hormone receptor status",
             "HER2 status",
             "Original stage",
             "Primary treating oncologist")),
  make_header_row("Genetic testing"),
  filter(nondrug_wide, `Clinical domain / variable` %in%
           c("Germline BRCA1/2 mutation identified",
             "Somatic PIK3CA mutation present",
             "Somatic ESR1 mutation present"))
)

# Drug rows: read from DrugComparison output if available
if (file.exists(drug_rows_path)) {
  drug_wide <- read_xlsx(drug_rows_path)
  # Keep only the drug rows from DrugComparison; the systemic-therapy line rows come
  # from line_metrics.R so that Table S2 agrees with Table 3.
  drug_only <- drug_wide %>%
    filter(grepl("^All prescribed anti-cancer drugs", `Clinical domain / variable`))
  line_block <- if (exists("s2_line_rows")) s2_line_rows else
    filter(drug_wide, !grepl("^All prescribed anti-cancer drugs", `Clinical domain / variable`))
  table2_wide <- bind_rows(
    nondrug_sectioned,
    make_header_row("Systemic therapy"),
    drug_only,
    line_block
  )
} else {
  cat("  NOTE: DrugComparison drug rows not found; Table S2 will be non-drug only.\n")
  cat("  Run DrugComparison.R first to generate", drug_rows_path, "\n")
  table2_wide <- bind_rows(
    nondrug_sectioned,
    make_header_row("Systemic therapy"),
    tibble(`Clinical domain / variable` = "(Run DrugComparison.R to populate)",
           `Metric reported` = "", `Second oncologist` = "",
           `GPT-5` = "", `Gemini 2.5 Pro` = "",
           `GPT-4o` = "", `DeepSeek-R1` = "",
           `Research coordinator` = "")
  )
}

# Chunk rank table (Supp Table S8)
chunk_rank_path <- file.path(LAB_MTG_DIR, "vector_rank_analytics.xlsx")
if (file.exists(chunk_rank_path)) {
  chunk_rank_raw <- read_xlsx(chunk_rank_path, sheet = "rank_by_variable")
  write_xlsx(chunk_rank_raw, file.path(TBL_DIR, "TableS8_chunk_rank.xlsx"))
  cat("  TableS8_chunk_rank.xlsx generated from vector_rank_analytics.xlsx\n")
} else {
  cat("  NOTE: vector_rank_analytics.xlsx not found; run NonDrugCleanup_VectorRank.R first.\n")
}

# Table 3 is built earlier so its line rows can also feed Table S2; write it here.
if (exists("table3_lines"))
  write_xlsx(table3_lines, file.path(TBL_DIR, "Table3_systemic_line_metrics.xlsx"))

# Save manuscript tables
write_xlsx(table1,       file.path(TBL_DIR, "Table2_patient_characteristics.xlsx"))

# TableS2: keep the display-wide sheet AND a detail sheet retaining evaluable-N
# and the implied correct count dropped by the wide pivot.
s2_detail <- table2 %>%
  mutate(`Agreement %` = suppressWarnings(as.numeric(sub("%.*", "", Agreement))),
         `Evaluable N`  = N,
         `Correct (approx.)` = ifelse(is.na(`Agreement %`) | is.na(N), NA_integer_,
                                      round(`Agreement %` / 100 * N))) %>%
  select(Variable, Category, LLM, Tolerance, N_GT, `Evaluable N`,
         `Correct (approx.)`, Extraction_Rate, Agreement)
write_xlsx(list(wide = table2_wide, detail_with_counts = s2_detail),
           file.path(TBL_DIR, "TableS2_full_model_performance.xlsx"))

# TableS5: keep the display-wide sheet AND a detail sheet retaining N, event
# counts, and raw log-HR/SE per source per equation.
s5_detail <- if (nrow(cox_tbl) > 0)
  cox_tbl %>% select(Outcome, Equation, Source, N, Events, loghr, se_loghr,
                     `HR [95% CI]`, `p-value`) else cox_tbl
write_xlsx(list(wide = cox_tbl_wide, detail_with_counts = s5_detail),
           file.path(TBL_DIR, "TableS5_cox_model.xlsx"))
write_xlsx(cost_tbl,      file.path(TBL_DIR, "TableS7_computational_cost.xlsx"))

cat("Manuscript tables written to", TBL_DIR, ":\n")
cat("  Table2_patient_characteristics.xlsx   (main text Table 2)\n")
cat("  TableS2_full_model_performance.xlsx   (Supp Table S2 — wide, matches appendix)\n")
cat("  TableS5_cox_model.xlsx                (Supp Table S5 — wide, matches appendix)\n")
cat("  TableS7_computational_cost.xlsx       (Supp Table S7)\n")
cat("  TableS8_chunk_rank.xlsx               (Supp Table S8)\n")
if (exists("table3_lines"))
  cat("  Table3_systemic_line_metrics.xlsx     (main text Table 3)\n")

dir.create(LOGS_DIR, showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()),
           file.path(LOGS_DIR, "session_info_manuscript_tables.txt"))
