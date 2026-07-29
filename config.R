#Root for the analysis for the first paper
BASE_DIR <- Sys.getenv("JUNIOR_FIRST_PAPER_BASE_DIR")

if (!nzchar(BASE_DIR)) {
  stop("Set the JUNIOR_FIRST_PAPER_BASE_DIR environment variable.")
}

# Run identifier & LLM identifiers
RUN_ID   <- "20260311_120527"
ALL_LLMS <- c("gpt4o", "gpt5", "deepseek_r1", "gemini_2_5")

LLM_LABELS <- c(
  gpt4o       = "GPT-4o",
  gpt5        = "GPT-5",
  deepseek_r1 = "DeepSeek-R1",
  gemini_2_5  = "Gemini 2.5"
)

# Input directories
EXP_DIR   <- file.path(BASE_DIR, "Run_Outputs/exp5")
GT_DIR_V1 <- file.path(BASE_DIR, "Ground_Truth/v1")
GT_DIR_V2 <- file.path(BASE_DIR, "Ground_Truth/v2")
CRC_DIR   <- file.path(BASE_DIR, "Run_Outputs/CRC_Review")
CHART_RAW_DIR <- file.path(BASE_DIR, "Raw_Data/100_full_charts/cleanedfullraw")

# Output directories
OUT_DIR   <- file.path(BASE_DIR, "Analysis_Output")
INT_DIR   <- file.path(BASE_DIR, "Analysis_Output/intermediate")
FINAL_DIR <- file.path(BASE_DIR, "Analysis_Output/final_dfs")
TBL_DIR   <- file.path(BASE_DIR, "Analysis_Output/tables")
FIG_DIR   <- file.path(BASE_DIR, "Analysis_Output/figures")
LOGS_DIR  <- file.path(BASE_DIR, "Analysis_Output/logs")

# Input file paths
GT_NON_DRUG_XLSX <- file.path(BASE_DIR, "Ground_Truth/non_therapy_line_ground_truth/ground_truth_data.xlsx")
GT_GENETICS_XLSX <- file.path(BASE_DIR, "Ground_Truth/non_therapy_line_ground_truth/ground_truth_genetics.xlsx")
CRC_LINES_XLSX   <- file.path(CRC_DIR,  "CRC_Lines_of_Treatment.xlsx")
CRC_DATES_XLSX   <- file.path(CRC_DIR,  "CRC_dates_diagnosis.xlsx")

# Intermediate file paths
LLM_NON_DRUG_CLEANED <- file.path(INT_DIR, "llm_non_drug_cleaned.xlsx")
GT_NON_DRUG_CLEANED  <- file.path(INT_DIR, "gt_non_drug_cleaned.xlsx")
CRC_NON_DRUG_CLEANED <- file.path(INT_DIR, "crc_non_drug_cleaned.xlsx")

# Final dataframe paths
FINAL_DRUG_DF        <- file.path(FINAL_DIR, "final_drug_df_for_comparisonfull.xlsx")
FINAL_NON_DRUG_LONG  <- file.path(FINAL_DIR, "non_drug_long.xlsx")
FINAL_NON_DRUG_WIDE  <- file.path(FINAL_DIR, "non_drug_wide.xlsx")

# Analysis constants
CRC_CENSOR_DATE <- as.Date("2022-06-01")  # fixed censor date for CRC reviewer
N_BOOT          <- 1000                   # bootstrap resamples
RANDOM_SEED     <- 42                     # set.seed() value — call explicitly in each script
TOL_PRIMARY     <- 90                     # primary date-agreement tolerance (days)

# Utility: null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# Date parsers
# Simple ISO parser — for LLM output files (format is always YYYY-MM-DD or YYYY-XX-XX).
parse_date <- function(x) {
  x <- substr(as.character(x), 1, 10)
  x <- gsub("[Xx]{2}", "01", x)
  as.Date(x, format = "%Y-%m-%d")
}

# Date parser 
parse_gt_date <- function(x) {
  x <- trimws(as.character(x))
  x <- stringr::str_replace_all(x, "XX", "01")
  x[x == "" | x == "NA"] <- NA_character_

  # ISO format: extract YYYY-MM-DD substring (handles trailing time component)
  d <- as.Date(stringr::str_extract(x, "\\d{4}-\\d{2}-\\d{2}"))

  # US slash format: MM/DD/YYYY, M/D/YYYY, MM/DD/YY
  needs <- !is.na(x) & is.na(d) & grepl("^\\d{1,2}/\\d{1,2}/\\d{2,4}$", x)
  if (any(needs)) {
    d[needs] <- as.Date(x[needs], tryFormats = c("%m/%d/%Y", "%m/%d/%y"))
  }

  # Excel serial integer (4–5 digit bare number string)
  serial <- !is.na(x) & is.na(d) & grepl("^\\d{4,5}$", x)
  if (any(serial)) {
    d[serial] <- as.Date(as.numeric(x[serial]), origin = "1899-12-30")
  }

  d
}

# Statistical helpers

# Wilson score interval 
wilson_ci <- function(k, n, conf = 0.95) {
  if (is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2); p <- k / n; den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half   <- (z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den
  c(max(0, center - half), min(1, center + half))
}

# Bootstrap 95% CI for a mean (BCa with percentile fallback). For binary (0/1)
boot_ci_95 <- function(x, R = N_BOOT) {
  x  <- x[!is.na(x)]
  binary <- length(x) > 0 && all(x %in% c(0, 1))
  if (length(x) < 3) {
    if (binary) return(wilson_ci(sum(x), length(x)))
    return(c(NA_real_, NA_real_))
  }
  b  <- boot::boot(data = x, statistic = function(d, i) mean(d[i]), R = R)
  ci <- tryCatch(
    boot::boot.ci(b, type = "bca",  conf = 0.95),
    error = function(e) tryCatch(boot::boot.ci(b, type = "perc", conf = 0.95),
                                 error = function(e2) NULL)
  )
  lo <- hi <- NA_real_
  if (!is.null(ci$bca))          { lo <- ci$bca[4];     hi <- ci$bca[5] }
  else if (!is.null(ci$percent)) { lo <- ci$percent[4]; hi <- ci$percent[5] }
  if ((!is.finite(lo) || !is.finite(hi)) && binary) return(wilson_ci(sum(x), length(x)))
  c(lo, hi)
}

boot_lo <- function(x) boot_ci_95(x)[1]
boot_hi <- function(x) boot_ci_95(x)[2]

# Cohen's kappa 
cohens_kappa <- function(a, b) {
  lvls <- sort(unique(c(a, b)))
  tb   <- table(factor(a, levels = lvls), factor(b, levels = lvls))
  if (sum(tb) == 0) return(NA_real_)
  po <- sum(diag(tb)) / sum(tb)
  pe <- sum(rowSums(tb) * colSums(tb)) / sum(tb)^2
  if (pe >= 1) return(1)
  (po - pe) / (1 - pe)
}

# Jaccard similarity for two sets (character or numeric vectors).
jaccard_set <- function(a, b) {
  if (is.null(a) || is.null(b) || length(a) == 0 || length(b) == 0) return(NA_real_)
  inter <- length(intersect(a, b))
  un    <- length(union(a, b))
  if (un == 0) return(NA_real_)
  inter / un
}

# Additional output paths
LAB_MTG_DIR <- file.path(BASE_DIR, "Analysis_Output/lab_meeting_figs_tables")

# Additional constants
TOL_SENS    <- 30               # sensitivity analysis date tolerance (days; paper reports 90-day primary + 30-day sensitivity)
M_DRAWS     <- 200              # convergence plot subsamples (exploratory viz only)

# Color palette
PAL_LLM <- c(
  "Gemini 2.5"  = "#377EB8",    # ColorBrewer Set1 blue
  "GPT-5"       = "#E41A1C",    # ColorBrewer Set1 red
  "GPT-4o"      = "#CC79A7",    # pink/magenta
  "DeepSeek-R1" = "#009E73"     # teal/green
)
COL_RC <- "#999999"   # CRC / research coordinator (grey)
COL_GT <- "#000000"   # oncologist ground truth (black)

# Domain ordering for forest plots and heatmaps
DOMAIN_LEVELS <- c("Dates", "Staging", "Receptor Status",
                   "Clinical", "Genetics", "Systemic Therapy")

# Generalized bootstrap CI
# Returns c(est=, lo=, hi=) for any summary statistic (use stat_fn = mean, median, etc.).
# BCa method with percentile fallback. Call set.seed(RANDOM_SEED) before use.
# For mean-only CI bounds, boot_ci_95() / boot_lo() / boot_hi() are also available.
boot_ci <- function(x, stat_fn, R = N_BOOT, conf = 0.95) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(c(est = NA_real_, lo = NA_real_, hi = NA_real_))
  est <- stat_fn(x)
  binary <- all(x %in% c(0, 1))
  if (length(x) < 3) {
    if (binary) { w <- wilson_ci(sum(x), length(x), conf); return(c(est = est, lo = w[1], hi = w[2])) }
    return(c(est = est, lo = NA_real_, hi = NA_real_))
  }
  b   <- boot::boot(data = x, statistic = function(d, i) stat_fn(d[i]), R = R)
  ci  <- tryCatch(
    boot::boot.ci(b, conf = conf, type = "bca"),
    error = function(e) tryCatch(
      boot::boot.ci(b, conf = conf, type = "perc"),
      error = function(e2) NULL
    )
  )
  lo <- hi <- NA_real_
  if (!is.null(ci$bca))          { lo <- ci$bca[4];     hi <- ci$bca[5] }
  else if (!is.null(ci$percent)) { lo <- ci$percent[4]; hi <- ci$percent[5] }
  # binary proportion with a degenerate bootstrap CI -> Wilson fallback (Task 3)
  if ((!is.finite(lo) || !is.finite(hi)) && binary) {
    w <- wilson_ci(sum(x), length(x), conf); lo <- w[1]; hi <- w[2]
  }
  c(est = est, lo = lo, hi = hi)
}

# Agreement vectors
date_agree_vec <- function(llm, gt, tol, na_as_wrong = FALSE, include_tn = FALSE) {
  if (include_tn) {
    # Full-cohort metric: denominator = all patients.
    # TN (both NA = no event on either side) counts as correct.
    # FP (LLM has date, GT is NA), FN (GT has date, LLM is NA), wrong date all count as wrong.
    as.integer(
      (is.na(gt) & is.na(llm)) |
      (!is.na(gt) & !is.na(llm) & abs(as.numeric(llm - gt)) <= tol)
    )
  } else if (na_as_wrong) {
    # Denominator = GT-non-NA; LLM-NA counts as wrong (0) — requires detection + date accuracy
    gt_valid <- !is.na(gt)
    as.integer(!is.na(llm[gt_valid]) & abs(as.numeric(llm[gt_valid] - gt[gt_valid])) <= tol)
  } else {
    valid <- !is.na(llm) & !is.na(gt)
    as.integer(abs(as.numeric(llm[valid] - gt[valid])) <= tol)
  }
}

cat_agree_vec <- function(llm, gt, na_as_wrong = FALSE) {
  if (na_as_wrong) {
    # Denominator = all GT-non-NA patients; LLM-NA counts as wrong (0)
    gt_valid <- !is.na(gt)
    as.integer(!is.na(llm[gt_valid]) & llm[gt_valid] == gt[gt_valid])
  } else {
    valid <- !is.na(llm) & !is.na(gt)
    as.integer(llm[valid] == gt[valid])
  }
}

# Formatting things

# Unified p-value formatter — use for log-rank, Cox, Wilcoxon, McNemar, chi-sq.
fmt_p <- function(p, digits = 2) {
  if (is.na(p)) return(NA_character_)
  threshold <- 10^(-digits)
  if (p < threshold) sprintf("p < %.*f", digits, threshold)
  else               sprintf("p = %.*f", digits, p)
}

# "X%\n[lo - hi]" — for agreement rates, accuracy, sensitivity/specificity in tables/figures.
# Point estimate on first line, CI on second line (renders as two lines in Word table cells). so it looks like econ tables
fmt_pct_ci <- function(est, lo, hi) {
  if (is.na(est)) return("\u2014")
  sprintf("%.0f%%\n[%.0f - %.0f]", est * 100, lo * 100, hi * 100)
}

# "X.XX\n[lo - hi]" — for HR, kappa, Jaccard, Spearman rho, etc.
# Point estimate on first line, CI on second line.
fmt_ci <- function(est, lo, hi, digits = 2) {
  if (is.na(est)) return("\u2014")
  sprintf("%.*f\n[%.*f - %.*f]", digits, est, digits, lo, digits, hi)
}

# Utility functions
# Binary mutation status: "mutated" vs "other" (treats NA as "other").
binarize_mutation <- function(x) {
  ifelse(!is.na(x) & tolower(as.character(x)) == "mutated", "mutated", "other")
}

# ggplot theme
theme_nature <- function(base_size = 8) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      plot.title        = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      plot.subtitle     = element_text(size = rel(0.95), color = "grey40", hjust = 0.5),
      axis.title        = element_text(size = rel(1.0)),
      axis.text         = element_text(size = rel(0.9), color = "black"),
      strip.text        = element_text(face = "bold", size = rel(1.0)),
      strip.background  = element_blank(),
      legend.title      = element_text(size = rel(0.9)),
      legend.text       = element_text(size = rel(0.85)),
      legend.key.size   = unit(0.4, "cm"),
      legend.background = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "grey92", linewidth = 0.2),
      plot.margin       = margin(5.5, 5.5, 5.5, 5.5)
    )
}

# Output dir creation
for (.d in c(OUT_DIR, INT_DIR, FINAL_DIR, TBL_DIR, FIG_DIR, LOGS_DIR, LAB_MTG_DIR)) {
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)
}
rm(.d)
