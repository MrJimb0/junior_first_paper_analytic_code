library(tidyverse)
library(readxl)
library(writexl)
library(boot)

source("config.R")

dir.create(FINAL_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(LAB_MTG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(RANDOM_SEED)

# Step 1: Read cleaned intermediates
gt  <- read_excel(GT_NON_DRUG_CLEANED)
crc <- read_excel(CRC_NON_DRUG_CLEANED)
llm <- read_excel(LLM_NON_DRUG_CLEANED)

# Step 2: Concatenate into long form (non_drug_long)
shared_cols <- c('gsID', 'version', 'birth_date', 'primary_oncologist_shc',
                 'date_initial_diagnosis', 'date_first_locoregional_recurrence',
                 'date_first_metastatic_diagnosis', 'stage_at_diagnosis',
                 'hr_her2_status', 'germline_testing_done', 'brca_status',
                 'brca1_variant', 'brca2_variant', 'pik3ca_status',
                 'pik3ca_variant', 'esr1_status', 'esr1_variant')

non_drug_long <- bind_rows(
  gt  %>% select(any_of(shared_cols)),
  crc %>% select(any_of(shared_cols)),
  llm %>% select(any_of(shared_cols))
) %>%
  mutate(
    birth_date                        = as.Date(birth_date),
    date_initial_diagnosis            = as.Date(date_initial_diagnosis),
    date_first_locoregional_recurrence = as.Date(date_first_locoregional_recurrence),
    date_first_metastatic_diagnosis   = as.Date(date_first_metastatic_diagnosis),
    stage_at_diagnosis                = as.numeric(stage_at_diagnosis)
  )

# Step 3: Create wide form (one row per gsID)
date_vars <- c('date_initial_diagnosis', 'date_first_locoregional_recurrence',
               'date_first_metastatic_diagnosis', 'birth_date')
cat_vars  <- c('stage_at_diagnosis', 'hr_her2_status', 'primary_oncologist_shc')
gen_vars  <- c('brca_status', 'pik3ca_status', 'esr1_status')

all_pivot_vars <- c(date_vars, cat_vars, gen_vars)

non_drug_wide <- non_drug_long %>%
  select(gsID, version, all_of(all_pivot_vars)) %>%
  pivot_wider(
    id_cols     = gsID,
    names_from  = version,
    values_from = all_of(all_pivot_vars),
    names_glue  = "{.value}_{version}"
  )

# Step 4: Compute analysis columns on wide form
llm_versions <- ALL_LLMS   # from config.R
tolerances   <- c(0, 30, 90, 180)

# 4a: Date diffs and agreement flags
for (dv in date_vars) {
  gt_col <- paste0(dv, "_gt")
  for (lv in c(llm_versions, "crc")) {
    llm_col  <- paste0(dv, "_", lv)
    if (!llm_col %in% names(non_drug_wide)) next
    diff_col <- paste0("diff_days_", dv, "_", lv)
    non_drug_wide[[diff_col]] <- as.numeric(
      as.Date(non_drug_wide[[llm_col]]) - as.Date(non_drug_wide[[gt_col]])
    )
    for (tol in tolerances) {
      flag_col <- paste0("agree_", dv, "_", lv, "_", tol)
      non_drug_wide[[flag_col]] <- as.integer(
        abs(non_drug_wide[[diff_col]]) <= tol
      )
    }
  }
}

# 4b: Categorical match flags (genetics use mutated vs other)
for (cv in cat_vars) {
  gt_col <- paste0(cv, "_gt")
  for (lv in llm_versions) {
    llm_col   <- paste0(cv, "_", lv)
    if (!llm_col %in% names(non_drug_wide)) next
    match_col <- paste0("match_", cv, "_", lv)
    non_drug_wide[[match_col]] <- as.integer(
      as.character(non_drug_wide[[llm_col]]) == as.character(non_drug_wide[[gt_col]])
    )
  }
}

for (cv in gen_vars) {
  gt_col <- paste0(cv, "_gt")
  for (lv in llm_versions) {
    llm_col   <- paste0(cv, "_", lv)
    if (!llm_col %in% names(non_drug_wide)) next
    match_col <- paste0("match_", cv, "_", lv)
    a_bin <- case_when(as.character(non_drug_wide[[llm_col]]) == "mutated" ~ "mutated", TRUE ~ "negative")
    b_bin <- case_when(as.character(non_drug_wide[[gt_col]])  == "mutated" ~ "mutated", TRUE ~ "negative")
    non_drug_wide[[match_col]] <- as.integer(a_bin == b_bin)
  }
}

# 4c: HR/HER2 decomposition from hr_her2_status string
for (v in c("gt", llm_versions)) {
  src_col <- paste0("hr_her2_status_", v)
  if (!src_col %in% names(non_drug_wide)) next
  non_drug_wide <- non_drug_wide %>%
    mutate(
      !!paste0("hr_", v) := case_when(
        str_detect(.data[[src_col]], "^HR\\+") ~ "HR+",
        str_detect(.data[[src_col]], "^HR\\-") ~ "HR-",
        TRUE ~ NA_character_
      ),
      !!paste0("her2_", v) := case_when(
        str_detect(.data[[src_col]], "HER2\\-$") ~ "HER2-",
        str_detect(.data[[src_col]], "HER2\\+$") ~ "HER2+",
        TRUE ~ NA_character_
      )
    )
}

# 4d: Per-patient error count + n_compared per LLM
for (lv in llm_versions) {
  err_cols <- c()
  for (dv in date_vars) {
    diff_col <- paste0("diff_days_", dv, "_", lv)
    ecol     <- paste0("err_", dv, "_", lv)
    if (diff_col %in% names(non_drug_wide)) {
      non_drug_wide[[ecol]] <- ifelse(
        is.na(non_drug_wide[[diff_col]]), NA_integer_,
        as.integer(abs(non_drug_wide[[diff_col]]) > TOL_PRIMARY)
      )
      err_cols <- c(err_cols, ecol)
    }
  }
  for (cv in c(cat_vars, gen_vars)) {
    mcol <- paste0("match_", cv, "_", lv)
    ecol <- paste0("err_", cv, "_", lv)
    if (mcol %in% names(non_drug_wide)) {
      non_drug_wide[[ecol]] <- ifelse(
        is.na(non_drug_wide[[mcol]]), NA_integer_,
        as.integer(non_drug_wide[[mcol]] == 0L)
      )
      err_cols <- c(err_cols, ecol)
    }
  }
  n_err_col  <- paste0("n_errors_", lv)
  n_comp_col <- paste0("n_compared_", lv)
  non_drug_wide[[n_err_col]]  <- rowSums(
    non_drug_wide[err_cols] == 1L, na.rm = TRUE
  )
  non_drug_wide[[n_comp_col]] <- rowSums(
    !is.na(non_drug_wide[err_cols])
  )
}

# Step 5: Export long + wide to final_dfs/
write_xlsx(non_drug_long, FINAL_NON_DRUG_LONG)
write_xlsx(non_drug_wide, FINAL_NON_DRUG_WIDE)

# Step 6: Date agreement analysis (Table 2 dates)

#Split gt and comparators, merge for pairwise comparison
gt_df   <- non_drug_long %>% filter(version == "gt") %>% select(-version)
comp_df <- non_drug_long %>% filter(version != "gt")
merged  <- comp_df %>%
  inner_join(gt_df, by = "gsID", suffix = c("", "_gt"))

date_results <- merged %>%
  group_by(version) %>%
  group_modify(~ {
    map_dfr(date_vars, function(v) {
      d1 <- as.Date(.x[[v]])
      d2 <- as.Date(.x[[paste0(v, "_gt")]])
      valid <- !is.na(d1) & !is.na(d2)
      diff  <- abs(as.numeric(d1[valid] - d2[valid]))

      tol_rows <- map_dfr(tolerances, function(tol) {
        vec <- as.integer(diff <= tol)
        ci  <- boot_ci_95(vec)
        tibble(
          variable  = v,
          tolerance = tol,
          n         = length(vec),
          agreement = if (length(vec) >= 3) mean(vec, na.rm = TRUE) else NA_real_,
          ci_lo     = ci[1],
          ci_hi     = ci[2]
        )
      })

      mae_row <- tibble(
        variable        = v,
        tolerance       = NA_real_,
        n               = length(diff),
        agreement       = NA_real_,
        ci_lo           = NA_real_,
        ci_hi           = NA_real_,
        median_abs_error = if (length(diff) > 0) median(diff) else NA_real_,
        iqr_lo           = if (length(diff) > 0) unname(quantile(diff, 0.25)) else NA_real_,
        iqr_hi           = if (length(diff) > 0) unname(quantile(diff, 0.75)) else NA_real_
      )

      bind_rows(
        tol_rows %>% mutate(median_abs_error = NA_real_, iqr_lo = NA_real_, iqr_hi = NA_real_),
        mae_row
      )
    })
  }) %>%
  ungroup()

# ── CRC override: met_dx comparison gives credit for EITHER GT LRR or GT met_dx
# Denominator = patients with any pre-censor GT recurrence (LRR or met_dx).
# This replaces the generic CRC rows for date_first_metastatic_diagnosis.
crc_met_rows <- merged %>%
  filter(version == "crc") %>%
  transmute(
    crc_met  = as.Date(date_first_metastatic_diagnosis),
    gt_met   = as.Date(date_first_metastatic_diagnosis_gt),
    gt_lrr   = as.Date(date_first_locoregional_recurrence_gt),
    in_denom = (!is.na(gt_lrr) & gt_lrr < CRC_CENSOR_DATE) |
               (!is.na(gt_met) & gt_met < CRC_CENSOR_DATE)
  ) %>%
  filter(in_denom) %>%
  mutate(
    d_met = ifelse(!is.na(crc_met) & !is.na(gt_met),
                   abs(as.numeric(crc_met - gt_met)), NA_real_),
    d_lrr = ifelse(!is.na(crc_met) & !is.na(gt_lrr),
                   abs(as.numeric(crc_met - gt_lrr)), NA_real_),
    best  = pmin(d_met, d_lrr, na.rm = TRUE),
    best  = ifelse(is.finite(best), best, NA_real_)
  )

crc_met_tol_rows <- map_dfr(tolerances, function(tol) {
  vec <- as.integer(!is.na(crc_met_rows$best) & crc_met_rows$best <= tol)
  ci  <- boot_ci_95(vec)
  tibble(
    version          = "crc",
    variable         = "date_first_metastatic_diagnosis",
    tolerance        = as.numeric(tol),
    n                = length(vec),
    agreement        = if (length(vec) >= 3) mean(vec, na.rm = TRUE) else NA_real_,
    ci_lo            = ci[1],
    ci_hi            = ci[2],
    median_abs_error = NA_real_,
    iqr_lo           = NA_real_,
    iqr_hi           = NA_real_
  )
})

best_valid <- crc_met_rows$best[!is.na(crc_met_rows$crc_met)]
crc_met_mae_row <- tibble(
  version          = "crc",
  variable         = "date_first_metastatic_diagnosis",
  tolerance        = NA_real_,
  n                = length(best_valid),
  agreement        = NA_real_,
  ci_lo            = NA_real_,
  ci_hi            = NA_real_,
  median_abs_error = if (length(best_valid) > 0) median(best_valid, na.rm = TRUE) else NA_real_,
  iqr_lo           = if (length(best_valid) > 0) unname(quantile(best_valid, 0.25, na.rm = TRUE)) else NA_real_,
  iqr_hi           = if (length(best_valid) > 0) unname(quantile(best_valid, 0.75, na.rm = TRUE)) else NA_real_
)

date_results <- date_results %>%
  filter(!(version == "crc" & variable == "date_first_metastatic_diagnosis")) %>%
  bind_rows(crc_met_tol_rows, crc_met_mae_row) %>%
  arrange(variable, version, tolerance)

write_xlsx(date_results, file.path(LAB_MTG_DIR, 'non_drug_date_agreement.xlsx'))

# Helper: bootstrap kappa CI
# Returns a named list: kappa, kappa_ci_lo, kappa_ci_hi.
compute_kappa_ci <- function(a, b, R = N_BOOT) {
  kappa_val  <- cohens_kappa(a, b)
  kappa_reps <- replicate(R, {
    idx <- sample(length(a), replace = TRUE)
    cohens_kappa(a[idx], b[idx])
  })
  list(
    kappa       = kappa_val,
    kappa_ci_lo = unname(quantile(kappa_reps, 0.025, na.rm = TRUE)),
    kappa_ci_hi = unname(quantile(kappa_reps, 0.975, na.rm = TRUE))
  )
}

# Step 7: Categorical agreement (Table 2 categorical)
cat_results <- merged %>%
  group_by(version) %>%
  group_modify(~ {
    map_dfr(cat_vars, function(v) {
      a <- as.character(.x[[v]])
      b <- as.character(.x[[paste0(v, "_gt")]])
      valid <- !is.na(a) & !is.na(b)

      vec <- as.integer(a[valid] == b[valid])
      ci  <- boot_ci_95(vec)

      kci <- compute_kappa_ci(a[valid], b[valid])

      tibble(
        variable        = v,
        n               = sum(valid),
        pct_exact_match = if (length(vec) >= 3) mean(vec, na.rm = TRUE) else NA_real_,
        match_ci_lo     = ci[1],
        match_ci_hi     = ci[2],
        kappa           = kci$kappa,
        kappa_ci_lo     = kci$kappa_ci_lo,
        kappa_ci_hi     = kci$kappa_ci_hi
      )
    })
  }) %>%
  ungroup()

write_xlsx(cat_results, file.path(LAB_MTG_DIR, 'non_drug_categorical_agreement.xlsx'))

# Cochran Q helper (manual implementation)
cochran_q <- function(mat) {
  k  <- ncol(mat); n <- nrow(mat)
  Lj <- rowSums(mat); Cj <- colSums(mat)
  Q  <- (k - 1) * (k * sum(Cj^2) - sum(Lj)^2) /
        (k * sum(Lj) - sum(Lj^2))
  p  <- pchisq(Q, df = k - 1, lower.tail = FALSE)
  list(Q = Q, df = k - 1, p = p)
}

# Helper: Cochran Q + pairwise McNemar for a set of match columns
# vars:         character vector of variable names (e.g. cat_vars or gen_vars)
# wide_df:      wide data frame with match columns
# llm_versions: character vector of LLM version ids
# Returns a tibble with Cochran Q stats and pairwise McNemar p-values per variable.
compute_paired_tests <- function(vars, wide_df, llm_versions) {
  map_dfr(vars, function(v) {
    match_cols <- paste0("match_", v, "_", llm_versions)
    available  <- match_cols[match_cols %in% names(wide_df)]
    if (length(available) < 2) return(NULL)

    mat <- as.matrix(wide_df[, available])
    complete_rows <- complete.cases(mat)
    mat <- mat[complete_rows, , drop = FALSE]
    colnames(mat) <- llm_versions[match_cols %in% names(wide_df)]

    lvs <- colnames(mat)

    cq <- tryCatch(cochran_q(mat),
                    error = function(e) list(Q = NA_real_, df = NA_integer_, p = NA_real_))

    row_out <- tibble(
      variable   = v,
      cochran_q  = cq$Q,
      cochran_df = cq$df,
      cochran_p  = cq$p
    )

    pairs <- combn(lvs, 2, simplify = FALSE)
    for (pr in pairs) {
      col_name <- paste0("mcnemar_", pr[1], "_vs_", pr[2])
      p_val <- tryCatch(
        mcnemar.test(table(mat[, pr[1]], mat[, pr[2]]))$p.value,
        error = function(e) NA_real_
      )
      row_out[[col_name]] <- p_val
    }

    row_out
  })
}

# Paired McNemar / Cochran Q for categorical variables
cat_mcnemar <- compute_paired_tests(cat_vars, non_drug_wide, llm_versions)

# Step 8: Genetics agreement (Table 2 genetics)
# Outcome of interest: mutated vs negative (missing = negative by design)

# binarize_genetic: genetics-specific binarizer returning "mutated"/"negative".
# Distinct from config's binarize_mutation() which returns "mutated"/"other".
# "negative" is used here because the confusion matrix and sensitivity/specificity
# calculations in this step label the reference category as "negative".
binarize_genetic <- function(x) {
  case_when(
    tolower(as.character(x)) == "mutated" ~ "mutated",
    TRUE                                  ~ "negative"
  )
}

gen_results <- merged %>%
  group_by(version) %>%
  group_modify(~ {
    map_dfr(gen_vars, function(v) {
      a <- binarize_genetic(.x[[v]])
      b <- binarize_genetic(.x[[paste0(v, "_gt")]])
      # No NAs remain after binarization; valid is always TRUE
      valid <- !is.na(a) & !is.na(b)

      vec <- as.integer(a[valid] == b[valid])
      ci  <- boot_ci_95(vec)

      kci <- compute_kappa_ci(a[valid], b[valid])

      #Sensitivity/specificity for detecting mutations
      tp <- sum(a[valid] == "mutated"  & b[valid] == "mutated")
      fn <- sum(a[valid] == "negative" & b[valid] == "mutated")
      fp <- sum(a[valid] == "mutated"  & b[valid] == "negative")
      tn <- sum(a[valid] == "negative" & b[valid] == "negative")

      tibble(
        variable        = v,
        n               = sum(valid),
        pct_exact_match = if (length(vec) >= 3) mean(vec, na.rm = TRUE) else NA_real_,
        match_ci_lo     = ci[1],
        match_ci_hi     = ci[2],
        kappa           = kci$kappa,
        kappa_ci_lo     = kci$kappa_ci_lo,
        kappa_ci_hi     = kci$kappa_ci_hi,
        sensitivity     = if (tp + fn > 0) tp / (tp + fn) else NA_real_,
        specificity     = if (tn + fp > 0) tn / (tn + fp) else NA_real_,
        tp = tp, fn = fn, fp = fp, tn = tn
      )
    })
  }) %>%
  ungroup()

write_xlsx(gen_results, file.path(LAB_MTG_DIR, 'non_drug_genetics_agreement.xlsx'))

# Paired McNemar / Cochran Q for genetics variables
gen_mcnemar <- compute_paired_tests(gen_vars, non_drug_wide, llm_versions)

# Step 9: Missing data rates per variable per version
all_vars <- c(date_vars, cat_vars, gen_vars)

missing_data <- merged %>%
  group_by(version) %>%
  group_modify(~ {
    map_dfr(all_vars, function(v) {
      llm_val <- .x[[v]]
      gt_val  <- .x[[paste0(v, "_gt")]]
      tibble(
        variable       = v,
        n_total        = nrow(.x),
        n_both_present = sum(!is.na(llm_val) & !is.na(gt_val)),
        n_llm_na       = sum(is.na(llm_val) & !is.na(gt_val)),
        n_gt_na        = sum(!is.na(llm_val) & is.na(gt_val)),
        n_both_na      = sum(is.na(llm_val) & is.na(gt_val))
      )
    })
  }) %>%
  ungroup()

write_xlsx(missing_data, file.path(LAB_MTG_DIR, 'non_drug_missing_data.xlsx'))

# Step 10: Per-patient error profile
per_patient <- merged %>%
  rowwise() %>%
  mutate(
    err_init_dx = {
      d1 <- as.Date(date_initial_diagnosis); d2 <- as.Date(date_initial_diagnosis_gt)
      if (!is.na(d1) & !is.na(d2)) as.integer(abs(as.numeric(d1 - d2)) > TOL_PRIMARY) else NA_integer_
    },
    err_lrr_dx = {
      d1 <- as.Date(date_first_locoregional_recurrence); d2 <- as.Date(date_first_locoregional_recurrence_gt)
      if (!is.na(d1) & !is.na(d2)) as.integer(abs(as.numeric(d1 - d2)) > TOL_PRIMARY) else NA_integer_
    },
    err_met_dx = {
      d1 <- as.Date(date_first_metastatic_diagnosis); d2 <- as.Date(date_first_metastatic_diagnosis_gt)
      if (!is.na(d1) & !is.na(d2)) as.integer(abs(as.numeric(d1 - d2)) > TOL_PRIMARY) else NA_integer_
    },
    err_birth = {
      d1 <- as.Date(birth_date); d2 <- as.Date(birth_date_gt)
      if (!is.na(d1) & !is.na(d2)) as.integer(abs(as.numeric(d1 - d2)) > TOL_PRIMARY) else NA_integer_
    },
    err_stage = {
      a <- as.character(stage_at_diagnosis); b <- as.character(stage_at_diagnosis_gt)
      if (!is.na(a) & !is.na(b)) as.integer(a != b) else NA_integer_
    },
    err_subtype = {
      a <- as.character(hr_her2_status); b <- as.character(hr_her2_status_gt)
      if (!is.na(a) & !is.na(b)) as.integer(a != b) else NA_integer_
    },
    err_oncologist = {
      a <- as.character(primary_oncologist_shc); b <- as.character(primary_oncologist_shc_gt)
      if (!is.na(a) & !is.na(b)) as.integer(a != b) else NA_integer_
    },
    err_brca = {
      a <- binarize_genetic(as.character(brca_status))
      b <- binarize_genetic(as.character(brca_status_gt))
      if (!is.na(a) & !is.na(b)) as.integer(a != b) else NA_integer_
    },
    err_pik3ca = {
      a <- binarize_genetic(as.character(pik3ca_status))
      b <- binarize_genetic(as.character(pik3ca_status_gt))
      if (!is.na(a) & !is.na(b)) as.integer(a != b) else NA_integer_
    },
    err_esr1 = {
      a <- binarize_genetic(as.character(esr1_status))
      b <- binarize_genetic(as.character(esr1_status_gt))
      if (!is.na(a) & !is.na(b)) as.integer(a != b) else NA_integer_
    }
  ) %>%
  ungroup() %>%
  mutate(
    n_errors   = rowSums(across(starts_with("err_"), ~ coalesce(.x, 0L))),
    n_compared = rowSums(across(starts_with("err_"), ~ as.integer(!is.na(.x))))
  ) %>%
  select(gsID, version, starts_with("err_"), n_errors, n_compared)

write_xlsx(per_patient, file.path(LAB_MTG_DIR, 'non_drug_per_patient_errors.xlsx'))

# Step 11: Confusion matrices — stage and subtype per LLM
confusion_stage <- merged %>%
  filter(!is.na(stage_at_diagnosis) & !is.na(stage_at_diagnosis_gt)) %>%
  count(version,
        llm_value = as.character(stage_at_diagnosis),
        gt_value  = as.character(stage_at_diagnosis_gt))

confusion_subtype <- merged %>%
  filter(!is.na(hr_her2_status) & !is.na(hr_her2_status_gt)) %>%
  count(version,
        llm_value = hr_her2_status,
        gt_value  = hr_her2_status_gt)

write_xlsx(
  list(confusion_stage = confusion_stage, confusion_subtype = confusion_subtype),
  file.path(LAB_MTG_DIR, 'non_drug_confusion_matrices.xlsx')
)

# Step 12: Date dispersion summary (signed diffs)
date_dispersion <- merged %>%
  group_by(version) %>%
  group_modify(~ {
    map_dfr(date_vars, function(v) {
      d1 <- as.Date(.x[[v]])
      d2 <- as.Date(.x[[paste0(v, "_gt")]])
      valid <- !is.na(d1) & !is.na(d2)
      signed_diff <- as.numeric(d1[valid] - d2[valid])

      tibble(
        variable    = v,
        n           = length(signed_diff),
        median_diff = if (length(signed_diff) > 0) median(signed_diff) else NA_real_,
        iqr_lo      = if (length(signed_diff) > 0) unname(quantile(signed_diff, 0.25)) else NA_real_,
        iqr_hi      = if (length(signed_diff) > 0) unname(quantile(signed_diff, 0.75)) else NA_real_,
        mean_diff   = if (length(signed_diff) > 0) mean(signed_diff) else NA_real_,
        sd_diff     = if (length(signed_diff) > 0) sd(signed_diff) else NA_real_,
        min_diff    = if (length(signed_diff) > 0) min(signed_diff) else NA_real_,
        max_diff    = if (length(signed_diff) > 0) max(signed_diff) else NA_real_
      )
    })
  }) %>%
  ungroup()

write_xlsx(date_dispersion, file.path(LAB_MTG_DIR, 'non_drug_date_dispersion.xlsx'))

# Step 13: Export analysis results to non_drug_analysis_results.xlsx
# All analysis sheets in one workbook for manuscript quick reference.

version_error_rate <- per_patient %>%
  group_by(version) %>%
  summarise(
    n_patients    = n(),
    mean_errors   = round(mean(n_errors,   na.rm = TRUE), 3),
    mean_compared = round(mean(n_compared, na.rm = TRUE), 3),
    error_rate    = round(sum(n_errors, na.rm = TRUE) / sum(n_compared, na.rm = TRUE), 4),
    pct_perfect   = round(mean(n_errors == 0, na.rm = TRUE), 3),
    .groups = "drop"
  )

write_xlsx(
  list(
    date_agreement     = date_results,
    cat_agreement      = cat_results,
    gen_agreement      = gen_results,
    missing_data       = missing_data,
    per_patient_errors = per_patient,
    confusion_stage    = confusion_stage,
    confusion_subtype  = confusion_subtype,
    date_dispersion    = date_dispersion,
    overall_by_version = version_error_rate,
    `Paired McNemar`   = cat_mcnemar,
    gen_paired_mcnemar = gen_mcnemar
  ),
  file.path(LAB_MTG_DIR, 'non_drug_analysis_results.xlsx')
)

dir.create(LOGS_DIR, showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(LOGS_DIR, "session_info_non_drug_analysis.txt"))
