# Line-level metrics for systemic therapy. Matches each comparator's treatment
# lines to the expert reference, then summarises line detection, regimen
# agreement, and attribute agreement for Table 3 and the Table S2 line rows.
source("config.R")
suppressMessages({library(dplyr); library(tidyr); library(readxl); library(writexl)
                  library(purrr); library(boot); library(stringr)})
set.seed(RANDOM_SEED)

START_MATCH_WINDOW_DAYS  <- TOL_PRIMARY   # two line starts this close count as the same line
ATTRIBUTE_TOLERANCE_DAYS <- TOL_SENS      # start/stop agreement tolerance on matched lines

# Comparators shown in Table 3, in display order. GPT-4o and DeepSeek-R1 appear in Table S2 only.
TABLE3_COMPARATORS <- c(v2 = "Second oncologist", gpt5 = "GPT-5",
                        gemini_2_5 = "Gemini 2.5 Pro", crc = "Research coordinator")
ALL_COMPARATORS    <- c("v2", "gpt5", "gemini_2_5", "gpt4o", "deepseek_r1", "crc")

# Match predicted lines to reference lines on start date, keeping both sequences in
# chronological order: if a predicted line matches a reference line, no earlier
# predicted line may match a later reference line. Among the alignments that match
# the most lines, the one with the smallest total start-date difference is chosen.
# Both inputs are sorted by Start, so row order is chronological order.
match_lines_in_order <- function(predicted_lines, reference_lines) {
  n_predicted <- nrow(predicted_lines)
  n_reference <- nrow(reference_lines)
  no_pairs <- matrix(integer(0), ncol = 2, dimnames = list(NULL, c("pred", "ref")))
  if (n_predicted == 0 || n_reference == 0) return(no_pairs)

  MATCHED        <- 1L
  SKIP_PREDICTED <- 2L
  SKIP_REFERENCE <- 3L

  match_count <- matrix(0L, n_predicted + 1, n_reference + 1)  # best matches using first i predicted, j reference
  total_gap   <- matrix(0,  n_predicted + 1, n_reference + 1)  # smallest total date gap achieving that many matches
  decision    <- matrix(0L, n_predicted + 1, n_reference + 1)

  for (i in seq_len(n_predicted)) for (j in seq_len(n_reference)) {
    count_without_predicted <- match_count[i, j + 1]
    gap_without_predicted   <- total_gap[i, j + 1]
    count_without_reference <- match_count[i + 1, j]
    gap_without_reference   <- total_gap[i + 1, j]

    if (count_without_predicted > count_without_reference ||
        (count_without_predicted == count_without_reference &&
         gap_without_predicted <= gap_without_reference)) {
      best_count <- count_without_predicted
      best_gap   <- gap_without_predicted
      best_step  <- SKIP_PREDICTED
    } else {
      best_count <- count_without_reference
      best_gap   <- gap_without_reference
      best_step  <- SKIP_REFERENCE
    }

    if (!is.na(predicted_lines$Start[i]) && !is.na(reference_lines$Start[j])) {
      date_gap <- abs(as.numeric(predicted_lines$Start[i] - reference_lines$Start[j]))
      if (date_gap <= START_MATCH_WINDOW_DAYS) {
        count_if_matched <- match_count[i, j] + 1L
        gap_if_matched   <- total_gap[i, j] + date_gap
        if (count_if_matched > best_count ||
            (count_if_matched == best_count && gap_if_matched < best_gap)) {
          best_count <- count_if_matched
          best_gap   <- gap_if_matched
          best_step  <- MATCHED
        }
      }
    }

    match_count[i + 1, j + 1] <- best_count
    total_gap[i + 1, j + 1]   <- best_gap
    decision[i + 1, j + 1]    <- best_step
  }

  # Walk the recorded decisions back to recover the matched pairs.
  i <- n_predicted; j <- n_reference; pairs <- list()
  while (i > 0 && j > 0) {
    step <- decision[i + 1, j + 1]
    if (step == MATCHED) { pairs[[length(pairs) + 1]] <- c(i, j); i <- i - 1; j <- j - 1 }
    else if (step == SKIP_PREDICTED) i <- i - 1
    else j <- j - 1
  }
  if (!length(pairs)) return(no_pairs)
  matched_pairs <- do.call(rbind, rev(pairs))
  colnames(matched_pairs) <- c("pred", "ref")
  matched_pairs
}

# One patient: the chronological match above, then a regimen fallback for the few
# lines with a missing start date, which have no position in the sequence.
match_patient_lines <- function(predicted_lines, reference_lines) {
  n_predicted <- nrow(predicted_lines)
  n_reference <- nrow(reference_lines)
  no_pairs <- matrix(integer(0), ncol = 2, dimnames = list(NULL, c("pred", "ref")))
  if (n_predicted == 0 || n_reference == 0) return(no_pairs)

  matches <- match_lines_in_order(predicted_lines, reference_lines)
  unmatched_predicted <- setdiff(seq_len(n_predicted), matches[, "pred"])
  unmatched_reference <- setdiff(seq_len(n_reference), matches[, "ref"])

  for (predicted_index in unmatched_predicted) {
    regimen <- predicted_lines$bag[predicted_index]
    candidates <- unmatched_reference[reference_lines$bag[unmatched_reference] == regimen & regimen != ""]
    candidates <- candidates[is.na(predicted_lines$Start[predicted_index]) |
                             is.na(reference_lines$Start[candidates])]
    if (length(candidates)) {
      matches <- rbind(matches, c(pred = predicted_index, ref = candidates[1]))
      unmatched_reference <- setdiff(unmatched_reference, candidates[1])
    }
  }
  matches
}

# Jaccard that counts a repeated regimen once per occurrence rather than collapsing it.
jaccard_counting_repeats <- function(first_regimens, second_regimens) {
  first_regimens  <- first_regimens[!is.na(first_regimens) & first_regimens != ""]
  second_regimens <- second_regimens[!is.na(second_regimens) & second_regimens != ""]
  if (length(first_regimens) == 0 && length(second_regimens) == 0) return(NA_real_)
  regimens <- union(names(table(first_regimens)), names(table(second_regimens)))
  first_counts  <- as.integer(table(first_regimens)[regimens]);  first_counts[is.na(first_counts)]  <- 0L
  second_counts <- as.integer(table(second_regimens)[regimens]); second_counts[is.na(second_counts)] <- 0L
  union_size <- sum(pmax(first_counts, second_counts))
  if (union_size == 0) return(NA_real_)
  sum(pmin(first_counts, second_counts)) / union_size
}

# A patient with lines on either side but no true positives scores 0; a patient with
# no lines on either side is not evaluable.
f1_score <- function(true_positives, n_predicted, n_reference) {
  if (n_predicted == 0 && n_reference == 0) return(NA_real_)
  precision <- if (n_predicted > 0) true_positives / n_predicted else 0
  recall    <- if (n_reference > 0) true_positives / n_reference else 0
  if (precision + recall == 0) return(0)
  2 * precision * recall / (precision + recall)
}

# Wilson interval for a proportion; stays finite at 0% and 100% agreement.
wilson_confidence_interval <- function(successes, trials, conf = 0.95) {
  if (is.na(trials) || trials == 0) return(c(NA_real_, NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  proportion  <- successes / trials
  denominator <- 1 + z^2 / trials
  center     <- (proportion + z^2 / (2 * trials)) / denominator
  half_width <- (z * sqrt(proportion * (1 - proportion) / trials + z^2 / (4 * trials^2))) / denominator
  c(proportion, max(0, center - half_width), min(1, center + half_width))
}

# Bootstrap confidence interval for a patient-level mean, falling back to a Wilson
# interval on the pooled counts when there are too few patients or the bootstrap
# interval is degenerate (for example at perfect agreement).
bootstrap_mean_ci <- function(values, successes = NA, trials = NA) {
  values <- values[!is.na(values)]
  if (length(values) < 3) {
    if (!is.na(successes) && !is.na(trials) && trials > 0)
      return(wilson_confidence_interval(successes, trials))
    return(c(if (length(values)) mean(values) else NA_real_, NA_real_, NA_real_))
  }
  resamples <- boot(values, function(d, i) mean(d[i]), R = N_BOOT)
  interval  <- tryCatch(boot.ci(resamples, type = "bca"),
                        error = function(e) tryCatch(boot.ci(resamples, type = "perc"),
                                                     error = function(e2) NULL))
  lower <- upper <- NA_real_
  if (!is.null(interval$bca))          { lower <- interval$bca[4];     upper <- interval$bca[5] }
  else if (!is.null(interval$percent)) { lower <- interval$percent[4]; upper <- interval$percent[5] }
  if ((!is.finite(lower) || !is.finite(upper)) && !is.na(successes) && !is.na(trials) && trials > 0)
    return(wilson_confidence_interval(successes, trials))
  c(mean(values), lower, upper)
}

format_percent <- function(v) if (is.na(v[1])) "—" else
  sprintf("%.0f%% [%.0f–%.0f]", 100 * v[1], 100 * v[2], 100 * v[3])
format_jaccard <- function(v) if (is.na(v[1])) "—" else
  sprintf("%.2f [%.2f–%.2f]", v[1], v[2], v[3])

# Treatment lines with each regimen reduced to a sorted drug key. Censoring matches
# DrugComparison.R: coordinator rows stop at the review cutoff, every source stops
# at the patient's last chart note.
load_treatment_lines <- function() {
  read_xlsx(FINAL_DRUG_DF) %>%
    mutate(Start = as.Date(Start), Stop = as.Date(Stop), censor_date = as.Date(censor_date),
           bag = map_chr(regimen_bag_of_drugs, function(regimen) {
             if (is.na(regimen)) return("")
             paste(sort(unique(trimws(tolower(strsplit(regimen, ",")[[1]])))), collapse = "|")
           })) %>%
    filter(is.na(Start) | is.na(censor_date) | Start <= censor_date)
}

# One row per patient for a single comparator, scored against the expert reference.
per_patient_metrics <- function(treatment_lines, comparator) {
  reference_all <- treatment_lines %>% filter(version == "v1")
  if (comparator == "crc")
    reference_all <- reference_all %>% filter(is.na(Start) | Start <= CRC_CENSOR_DATE)
  patient_ids <- sort(unique(treatment_lines$gsID))

  map_dfr(patient_ids, function(patient_id) {
    reference_lines <- reference_all %>% filter(gsID == patient_id) %>% arrange(Start)
    predicted_lines <- treatment_lines %>%
      filter(version == comparator, gsID == patient_id) %>% arrange(Start)
    n_predicted <- nrow(predicted_lines)
    n_reference <- nrow(reference_lines)
    matches     <- match_patient_lines(predicted_lines, reference_lines)
    n_matched   <- nrow(matches)

    exact_content_matches  <- 0L
    exact_sequence_matches <- 0L
    start_agreements <- start_comparisons <- 0L
    stop_agreements  <- stop_comparisons  <- 0L
    intent_agreements <- intent_comparisons <- 0L
    reason_agreements <- reason_comparisons <- 0L

    if (n_matched) {
      reference_order <- rank(matches[, "ref"],  ties.method = "first")
      predicted_order <- rank(matches[, "pred"], ties.method = "first")
      for (k in seq_len(n_matched)) {
        predicted_index <- matches[k, "pred"]
        reference_index <- matches[k, "ref"]

        if (predicted_lines$bag[predicted_index] == reference_lines$bag[reference_index] &&
            reference_lines$bag[reference_index] != "") {
          exact_content_matches <- exact_content_matches + 1L
          if (reference_order[k] == predicted_order[k])
            exact_sequence_matches <- exact_sequence_matches + 1L
        }

        # Attributes are only compared where both sides documented them.
        if (!is.na(predicted_lines$Start[predicted_index]) &&
            !is.na(reference_lines$Start[reference_index])) {
          start_comparisons <- start_comparisons + 1L
          if (abs(as.numeric(predicted_lines$Start[predicted_index] -
                             reference_lines$Start[reference_index])) <= ATTRIBUTE_TOLERANCE_DAYS)
            start_agreements <- start_agreements + 1L
        }
        if (!is.na(predicted_lines$Stop[predicted_index]) &&
            !is.na(reference_lines$Stop[reference_index])) {
          stop_comparisons <- stop_comparisons + 1L
          if (abs(as.numeric(predicted_lines$Stop[predicted_index] -
                             reference_lines$Stop[reference_index])) <= ATTRIBUTE_TOLERANCE_DAYS)
            stop_agreements <- stop_agreements + 1L
        }
        if (!is.na(predicted_lines$TreatmentIntent[predicted_index]) &&
            !is.na(reference_lines$TreatmentIntent[reference_index])) {
          intent_comparisons <- intent_comparisons + 1L
          if (str_squish(tolower(predicted_lines$TreatmentIntent[predicted_index])) ==
              str_squish(tolower(reference_lines$TreatmentIntent[reference_index])))
            intent_agreements <- intent_agreements + 1L
        }
        if (!is.na(predicted_lines$ReasonStop[predicted_index]) &&
            !is.na(reference_lines$ReasonStop[reference_index])) {
          reason_comparisons <- reason_comparisons + 1L
          if (str_squish(tolower(predicted_lines$ReasonStop[predicted_index])) ==
              str_squish(tolower(reference_lines$ReasonStop[reference_index])))
            reason_agreements <- reason_agreements + 1L
        }
      }
    }

    agreement_rate <- function(agreements, comparisons)
      if (comparisons > 0) agreements / comparisons else NA_real_

    tibble(
      model = comparator, gsID = patient_id,
      n_pred = n_predicted, n_ref = n_reference, n_match = n_matched,
      exact_tp = exact_content_matches, exseq_tp = exact_sequence_matches,
      start_prec = if (n_predicted > 0) n_matched / n_predicted else 0,
      start_rec  = if (n_reference > 0) n_matched / n_reference else 0,
      start_f1   = f1_score(n_matched, n_predicted, n_reference),
      exact_prec = if (n_predicted > 0) exact_content_matches / n_predicted else 0,
      exact_rec  = if (n_reference > 0) exact_content_matches / n_reference else 0,
      exact_f1   = f1_score(exact_content_matches, n_predicted, n_reference),
      exseq_f1   = f1_score(exact_sequence_matches, n_predicted, n_reference),
      jaccard_multiset = jaccard_counting_repeats(predicted_lines$regimen_standardized,
                                                  reference_lines$regimen_standardized),
      jaccard_uniqueset = {
        predicted_regimens <- unique(predicted_lines$regimen_standardized[
          !is.na(predicted_lines$regimen_standardized)])
        reference_regimens <- unique(reference_lines$regimen_standardized[
          !is.na(reference_lines$regimen_standardized)])
        if (length(predicted_regimens) == 0 && length(reference_regimens) == 0) NA_real_
        else length(intersect(predicted_regimens, reference_regimens)) /
             max(1, length(union(predicted_regimens, reference_regimens)))
      },
      # Shared lines are those the chronological match pairs with identical regimen content.
      jaccard_matched = {
        union_size <- n_predicted + n_reference - exact_content_matches
        if (union_size > 0) exact_content_matches / union_size else NA_real_
      },
      start_agree   = agreement_rate(start_agreements, start_comparisons),
      stop_agree    = agreement_rate(stop_agreements, stop_comparisons),
      intent_agree  = agreement_rate(intent_agreements, intent_comparisons),
      discont_agree = agreement_rate(reason_agreements, reason_comparisons),
      st_ok = start_agreements,  st_n = start_comparisons,
      sp_ok = stop_agreements,   sp_n = stop_comparisons,
      it_ok = intent_agreements, it_n = intent_comparisons,
      ds_ok = reason_agreements, ds_n = reason_comparisons)
  })
}

build_line_metrics <- function() {
  treatment_lines <- load_treatment_lines()
  per_patient     <- map_dfr(ALL_COMPARATORS, ~ per_patient_metrics(treatment_lines, .x))

  line_counts <- per_patient %>% group_by(model) %>% summarise(
    ref_lines = sum(n_ref), predicted_lines = sum(n_pred), matched = sum(n_match),
    exact_content_tp = sum(exact_tp), exact_sequence_tp = sum(exseq_tp),
    missed = sum(n_ref) - sum(n_match), additional_predicted = sum(n_pred) - sum(n_match),
    .groups = "drop") %>% arrange(match(model, ALL_COMPARATORS))

  metric_row <- function(comparator) {
    patient_rows <- per_patient %>% filter(model == !!comparator)
    tibble(
      Comparator = unname(c(TABLE3_COMPARATORS, gpt4o = "GPT-4o",
                            deepseek_r1 = "DeepSeek-R1")[comparator]),
      model = comparator,
      `Exact-line (content) Jaccard`    = format_jaccard(bootstrap_mean_ci(patient_rows$jaccard_matched)),
      `Unique-regimen-set Jaccard`      = format_jaccard(bootstrap_mean_ci(patient_rows$jaccard_uniqueset)),
      `Line-start precision`            = format_percent(bootstrap_mean_ci(patient_rows$start_prec)),
      `Line-start recall`               = format_percent(bootstrap_mean_ci(patient_rows$start_rec)),
      `Line-start F1`                   = format_percent(bootstrap_mean_ci(patient_rows$start_f1)),
      `Exact-line (content) precision`  = format_percent(bootstrap_mean_ci(patient_rows$exact_prec)),
      `Exact-line (content) recall`     = format_percent(bootstrap_mean_ci(patient_rows$exact_rec)),
      `Exact-line (content) F1`         = format_percent(bootstrap_mean_ci(patient_rows$exact_f1)),
      `Exact-line (content+order) F1`   = format_percent(bootstrap_mean_ci(patient_rows$exseq_f1)),
      `Start-date agreement (±30d)`     = format_percent(bootstrap_mean_ci(patient_rows$start_agree,
                                            sum(patient_rows$st_ok), sum(patient_rows$st_n))),
      `Stop-date agreement (±30d)`      = format_percent(bootstrap_mean_ci(patient_rows$stop_agree,
                                            sum(patient_rows$sp_ok), sum(patient_rows$sp_n))),
      `Treatment-intent agreement`      = format_percent(bootstrap_mean_ci(patient_rows$intent_agree,
                                            sum(patient_rows$it_ok), sum(patient_rows$it_n))),
      `Discontinuation-reason agreement`= format_percent(bootstrap_mean_ci(patient_rows$discont_agree,
                                            sum(patient_rows$ds_ok), sum(patient_rows$ds_n))))
  }
  metrics <- map_dfr(names(c(TABLE3_COMPARATORS, gpt4o = 1, deepseek_r1 = 1)), metric_row)

  table3 <- metrics %>% filter(model %in% names(TABLE3_COMPARATORS)) %>%
    mutate(model = factor(model, levels = names(TABLE3_COMPARATORS))) %>%
    arrange(model) %>% select(-model)

  # Table S2 line rows for all six comparators, wide, for manuscript_tables.R. The
  # same matcher feeds both, so the supplement and the main text agree.
  s2_columns <- c("Second oncologist", "GPT-5", "Gemini 2.5 Pro", "GPT-4o",
                  "DeepSeek-R1", "Research coordinator")
  s2_comparator <- c("Second oncologist" = "v2", "GPT-5" = "gpt5", "Gemini 2.5 Pro" = "gemini_2_5",
                     "GPT-4o" = "gpt4o", "DeepSeek-R1" = "deepseek_r1",
                     "Research coordinator" = "crc")
  metric_values <- function(column) setNames(vapply(s2_columns, function(label)
    as.character(metrics[[column]][metrics$Comparator == label]), character(1)), s2_columns)
  matched_counts <- function() setNames(vapply(s2_columns, function(label) {
    counts <- line_counts[line_counts$model == s2_comparator[label], ]
    sprintf("%d / %d / %d", counts$matched, counts$missed, counts$additional_predicted)
  }, character(1)), s2_columns)
  s2_row <- function(label, metric_reported, values)
    as_tibble(c(list(`Clinical domain / variable` = label, `Metric reported` = metric_reported),
                as.list(values)))

  s2_line_rows <- bind_rows(
    s2_row("Systemic therapy lines: exact-line Jaccard (order-aware matched lines)", "Mean [95% CI]", metric_values("Exact-line (content) Jaccard")),
    s2_row("Systemic therapy lines: regimen-set Jaccard (order-agnostic, reference)", "Mean [95% CI]", metric_values("Unique-regimen-set Jaccard")),
    s2_row("Systemic therapy lines: line-start detection precision", "Mean (95% CI)", metric_values("Line-start precision")),
    s2_row("Systemic therapy lines: line-start detection recall", "Mean (95% CI)", metric_values("Line-start recall")),
    s2_row("Systemic therapy lines: line-start detection F1", "Mean (95% CI)", metric_values("Line-start F1")),
    s2_row("Systemic therapy lines: exact reconstruction (regimen content) F1", "Mean (95% CI)", metric_values("Exact-line (content) F1")),
    s2_row("Systemic therapy lines: exact reconstruction (content + order) F1", "Mean (95% CI)", metric_values("Exact-line (content+order) F1")),
    s2_row("Line-start matched / missed / additional predicted, No.", "No.", matched_counts()),
    s2_row("Start date of matched treatment line (±30 d)", "% agreement (95% CI)", metric_values("Start-date agreement (±30d)")),
    s2_row("Stop date of matched treatment line (±30 d)", "% agreement (95% CI)", metric_values("Stop-date agreement (±30d)")),
    s2_row("Treatment intent for matched lines", "% agreement (95% CI)", metric_values("Treatment-intent agreement")),
    s2_row("Reason for discontinuation for matched lines", "% agreement (95% CI)", metric_values("Discontinuation-reason agreement"))
  )

  out_path <- file.path(TBL_DIR, "Table3_line_metrics.xlsx")
  write_xlsx(list(Table3 = table3, all_comparators = select(metrics, -model),
                  s2_line_rows = s2_line_rows, per_patient = per_patient,
                  reconciliation = line_counts), out_path)

  cat("\nLine-start detection counts\n")
  print(as.data.frame(line_counts))
  cat("\nTable 3\n")
  print(as.data.frame(t(table3)), right = FALSE)
  cat("\nwrote", out_path, "\n")

  invisible(list(per_patient = per_patient, line_counts = line_counts,
                 table3 = table3, s2_line_rows = s2_line_rows))
}
