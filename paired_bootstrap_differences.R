###############################################################################
# paired_bootstrap_differences.R
#
# Paired patient-level bootstrap of differences between abstraction sources
# measured on the same 100 patients.
#   - Mean patient-level Jaccard similarity, drugs and exact lines:
#     each LLM minus the second oncologist (both scored against the expert).
#   - Median overall and recurrence-free survival (death at last follow-up, as
#     in Figure 4): each LLM-derived dataset minus the expert-derived dataset.
#   - Cox log hazard ratios (conventional censoring, as in Supplementary Table 9):
#     each LLM-derived minus expert-derived, reported as a ratio of hazard ratios.
# One bootstrap draw resamples patients with replacement and recomputes every
# statistic on that draw, so each patient's expert, oncologist and LLM values
# stay together. Intervals are BCa with a percentile fallback, as in config.R.
#
# Usage: Rscript paired_bootstrap_differences.R <folder holding this code>
# Run after manuscript_tables.R, which writes the per-patient Jaccard files.
# Prints aggregate values only.
###############################################################################

analysis_code_folder <- normalizePath(commandArgs(TRUE)[1])
setwd(analysis_code_folder)
invisible(capture.output(suppressMessages(suppressWarnings(source("manuscript_tables.R")))))

# Reuse the Figure 4 survival builders exactly as written in manuscript_figures.R.
figure_code <- parse("manuscript_figures.R")
for (expression in figure_code) {
  if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      deparse(expression[[2]]) %in% c("build_os_lb", "build_rfs_lb")) eval(expression)
}

best_llms <- c(gpt5 = "GPT-5", gemini_2_5 = "Gemini 2.5")

# Per-patient Jaccard values written by the pipeline
drug_jaccard <- readxl::read_xlsx(file.path(INT_DIR, "wide_bag_jaccard_per_patient.xlsx"))
line_jaccard <- readxl::read_xlsx(file.path(TBL_DIR, "Table3_line_metrics.xlsx"), sheet = "per_patient")
line_jaccard_for <- function(model) line_jaccard %>% filter(model == !!model) %>% select(gsID, value = jaccard_matched)
drug_jaccard_for <- function(column) drug_jaccard %>% select(gsID, value = all_of(column))

# Figure 4 survival datasets (one row per patient per source)
expert_km <- gt %>% transmute(gsID, init_dx = gt_init_dx, lrr_dx = gt_lrr_dx, met_dx = gt_met_dx,
                              death_date = gt_death, last_note = gt_last_note)
survival_sets <- list(OS = list(Expert = build_os_lb(expert_km, "init_dx", "Expert")),
                      RFS = list(Expert = build_rfs_lb(expert_km, "init_dx", "lrr_dx", "met_dx", "Expert", stage4_pids)))
for (llm in names(best_llms)) {
  llm_km <- all_llm_data[[llm]] %>% left_join(gt %>% select(gsID, last_note = gt_last_note), by = "gsID")
  llm_stage_iv <- llm_km$gsID[coalesce(llm_km$stage == "4", FALSE)]
  survival_sets$OS[[llm]]  <- build_os_lb(llm_km, "init_dx", llm)
  survival_sets$RFS[[llm]] <- build_rfs_lb(llm_km, "init_dx", "lrr_dx", "met_dx", llm, llm_stage_iv)
}

# Cox datasets from manuscript_tables.R, keyed by source label
cox_sets <- list(OS = os_sources, RFS = rfs_sources)

resample <- function(per_patient, ids) {
  rows <- match(ids, per_patient$gsID)
  per_patient[rows[!is.na(rows)], , drop = FALSE]
}
median_survival <- function(d) unname(summary(survfit(Surv(T_months, event) ~ 1, data = d))$table["median"])
log_hazard_ratio <- function(d, covariate) {
  d <- d %>% mutate(x = as.integer(.data[[covariate]])) %>% filter(!is.na(x))
  fit <- tryCatch(suppressWarnings(coxph(Surv(T_months, event) ~ x, data = d)), error = function(e) NULL)
  if (is.null(fit)) NA_real_ else unname(coef(fit)["x"])
}

statistics_for <- function(ids) {
  values <- c()
  for (llm in names(best_llms)) {
    drug_column <- paste0("jaccard_", llm, "_vs_v1")
    values[paste("Drug Jaccard", llm)] <-
      mean(resample(drug_jaccard_for(drug_column), ids)$value, na.rm = TRUE) -
      mean(resample(drug_jaccard_for("jaccard_v2_vs_v1"), ids)$value, na.rm = TRUE)
    values[paste("Line Jaccard", llm)] <-
      mean(resample(line_jaccard_for(llm), ids)$value, na.rm = TRUE) -
      mean(resample(line_jaccard_for("v2"), ids)$value, na.rm = TRUE)
    for (outcome in c("OS", "RFS")) {
      values[paste("Median", outcome, llm)] <-
        median_survival(resample(survival_sets[[outcome]][[llm]], ids)) -
        median_survival(resample(survival_sets[[outcome]]$Expert, ids))
    }
    for (model in models) {
      sources <- cox_sets[[model$outcome]]
      values[paste("Log HR", model$eq, llm)] <-
        log_hazard_ratio(resample(sources[[best_llms[[llm]]]], ids), model$var) -
        log_hazard_ratio(resample(sources[["Oncologist (GT)"]], ids), model$var)
    }
  }
  values
}

patient_ids <- sort(unique(gt$gsID))
set.seed(RANDOM_SEED)
draws <- boot::boot(patient_ids, function(ids, i) statistics_for(ids[i]), R = N_BOOT)

interval_for <- function(i) {
  replicates <- draws$t[, i]
  if (any(!is.finite(replicates))) return(c(quantile(replicates, c(.025, .975), na.rm = TRUE), method = NA))
  ci <- tryCatch(boot::boot.ci(draws, type = "bca", index = i), error = function(e) NULL)
  if (!is.null(ci$bca) && all(is.finite(ci$bca[4:5]))) return(c(ci$bca[4:5], method = 1))
  c(quantile(replicates, c(.025, .975)), method = 0)
}
intervals <- t(sapply(seq_along(draws$t0), interval_for))
results <- tibble(statistic = names(draws$t0), difference = draws$t0,
                  lower = intervals[, 1], upper = intervals[, 2],
                  method = c("percentile", "BCa")[intervals[, 3] + 1],
                  failed_draws = colSums(!is.finite(draws$t)))
log_rows <- startsWith(results$statistic, "Log HR")
results$ratio_of_hazard_ratios <- ifelse(log_rows, sprintf("%.2f [%.2f, %.2f]", exp(results$difference),
                                          exp(results$lower), exp(results$upper)), "")
options(width = 200)
cat(sprintf("\nPaired bootstrap: %d patients, %d draws, seed %d\n", length(patient_ids), N_BOOT, RANDOM_SEED))
print(as.data.frame(results %>% mutate(across(c(difference, lower, upper), ~ round(.x, 4)))))

cat("\nMedian survival by source, months [95% CI from survfit]\n")
for (outcome in c("OS", "RFS")) for (source in names(survival_sets[[outcome]])) {
  fit_table <- summary(survfit(Surv(T_months, event) ~ 1, data = survival_sets[[outcome]][[source]]))$table
  cat(sprintf("  %-4s %-11s n=%3d median %.1f [%.1f, %.1f]\n", outcome, source, nrow(survival_sets[[outcome]][[source]]),
              fit_table["median"], fit_table["0.95LCL"], fit_table["0.95UCL"]))
}
