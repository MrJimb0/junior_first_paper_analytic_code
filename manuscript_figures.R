library(tidyverse)
library(readxl)
library(writexl)
library(survival)
library(broom)
library(patchwork)
library(boot)

source("config.R")

# Figures display the leading model as "Gemini 2.5 Pro" (config keeps the short
# internal label, so the tables' own relabeling is unaffected).
LLM_LABELS[['gemini_2_5']] <- 'Gemini 2.5 Pro'
names(PAL_LLM)[names(PAL_LLM) == 'Gemini 2.5'] <- 'Gemini 2.5 Pro'

set.seed(RANDOM_SEED)

### Helper functions

# Event determined solely by presence of death_date.
build_os <- function(df, dx_col, source_label) {
  df %>%
    mutate(
      dx_date  = .data[[dx_col]],
      event    = as.integer(!is.na(death_date)),
      end_date = if_else(!is.na(death_date), death_date, last_note),
      T_days   = as.numeric(difftime(end_date, dx_date, units = 'days')),
      T_months = T_days / 30.4375
    ) %>%
    filter(!is.na(dx_date), !is.na(end_date), T_days >= 0) %>%
    mutate(source = source_label) %>%
    select(gsID, source, event, T_months)
}

# Recurrence event = earliest LRR or metastatic date from that source.
build_rfs <- function(df, dx_col, lrr_col, met_col,
                      source_label, stage4_pids) {
  df %>%
    filter(!(gsID %in% stage4_pids)) %>%
    rowwise() %>%
    mutate(
      dx_date  = .data[[dx_col]],
      recur    = {
        cands <- na.omit(c(.data[[lrr_col]], .data[[met_col]]))
        if (length(cands) > 0) min(cands) else as.Date(NA)
      },
      event    = as.integer(!is.na(recur) & !is.na(dx_date) & recur >= dx_date),
      end_date = if (event == 1L) recur else coalesce(death_date, last_note),
      T_days   = as.numeric(difftime(end_date, dx_date, units = 'days')),
      T_months = T_days / 30.4375
    ) %>%
    ungroup() %>%
    filter(!is.na(dx_date), !is.na(end_date), T_days >= 0) %>%
    mutate(source = source_label) %>%
    select(gsID, source, event, T_months)
}

plot_km <- function(data, title_str, pal,
                    y_label = 'Survival Probability', max_yr = 5,
                    use_months = FALSE, show_risk_labels = TRUE,
                    overlay_data = NULL) {

  max_mo <- max_yr * 12

  # Compute median survival from FULL (uncensored) data
  km_fit_full <- survfit(Surv(T_months, event) ~ source, data = data)

  # Administratively censor for plotting only
  data_plot <- data %>% mutate(
    event    = if_else(T_months > max_mo, 0L, event),
    T_months = pmin(T_months, max_mo)
  )

  km_fit <- survfit(Surv(T_months, event) ~ source, data = data_plot)
  km_df  <- tidy(km_fit) %>%
    mutate(source = factor(str_remove(strata, 'source='), levels = names(pal)))

  censor_df <- km_df %>% filter(n.censor > 0)
  src_names <- names(pal)
  n_src     <- length(src_names)

  if (use_months) {
    # x-axis in months
    x_breaks <- seq(0, max_mo, by = 6)
    km_df     <- km_df     %>% mutate(time_x = time)
    censor_df <- censor_df %>% mutate(time_x = time)
    data_plot <- data_plot %>% mutate(T_x = T_months)
    x_max  <- max_mo
    x_lab  <- 'Months from Diagnosis'
    # Origin
    origins <- tibble(
      source = factor(src_names, levels = names(pal)),
      time = 0, time_x = 0, estimate = 1, conf.low = 1, conf.high = 1,
      n.event = 0L, n.censor = 0L, n.risk = NA_integer_,
      std.error = 0, strata = paste0('source=', src_names))
    km_df <- bind_rows(origins, km_df)
    # Extend
    km_df <- km_df %>%
      group_by(source) %>%
      group_modify(~ {
        last_row <- .x %>% slice_tail(n = 1)
        if (last_row$time_x < x_max) {
          ext <- last_row; ext$time_x <- x_max; ext$time <- x_max
          ext$n.event <- 0L; ext$n.censor <- 0L
          bind_rows(.x, ext)
        } else .x
      }) %>% ungroup()
  } else {
    # x-axis in years
    x_breaks <- seq(0, max_yr, by = 1)
    km_df     <- km_df     %>% mutate(time_x = time / 12)
    censor_df <- censor_df %>% mutate(time_x = time / 12)
    data_plot <- data_plot %>% mutate(T_x = T_months / 12)
    x_max  <- max_yr
    x_lab  <- 'Years from Diagnosis'
    origins <- tibble(
      source = factor(src_names, levels = names(pal)),
      time = 0, time_x = 0, estimate = 1, conf.low = 1, conf.high = 1,
      n.event = 0L, n.censor = 0L, n.risk = NA_integer_,
      std.error = 0, strata = paste0('source=', src_names))
    km_df <- bind_rows(origins, km_df)
    km_df <- km_df %>%
      group_by(source) %>%
      group_modify(~ {
        last_row <- .x %>% slice_tail(n = 1)
        if (last_row$time_x < x_max) {
          ext <- last_row; ext$time_x <- x_max; ext$time <- x_max * 12
          ext$n.event <- 0L; ext$n.censor <- 0L
          bind_rows(.x, ext)
        } else .x
      }) %>% ungroup()
  }

  # Optional overlay: compute KM curves from reference data for dashed lines
  overlay_km_df <- NULL
  if (!is.null(overlay_data)) {
    ov_plot <- overlay_data %>% mutate(
      event    = if_else(T_months > max_mo, 0L, event),
      T_months = pmin(T_months, max_mo))
    ov_fit <- survfit(Surv(T_months, event) ~ source, data = ov_plot)
    overlay_km_df <- tidy(ov_fit) %>%
      mutate(source = factor(str_remove(strata, 'source='), levels = names(pal)))
    if (use_months) {
      overlay_km_df <- overlay_km_df %>% mutate(time_x = time)
      ov_origins <- tibble(
        source = factor(src_names, levels = names(pal)),
        time = 0, time_x = 0, estimate = 1, conf.low = 1, conf.high = 1,
        n.event = 0L, n.censor = 0L, n.risk = NA_integer_,
        std.error = 0, strata = paste0('source=', src_names))
      overlay_km_df <- bind_rows(ov_origins, overlay_km_df)
      overlay_km_df <- overlay_km_df %>%
        group_by(source) %>%
        group_modify(~ {
          lr <- .x %>% slice_tail(n = 1)
          if (lr$time_x < x_max) {
            ext <- lr; ext$time_x <- x_max; ext$time <- x_max
            ext$n.event <- 0L; ext$n.censor <- 0L
            bind_rows(.x, ext)
          } else .x
        }) %>% ungroup()
    } else {
      overlay_km_df <- overlay_km_df %>% mutate(time_x = time / 12)
      ov_origins <- tibble(
        source = factor(src_names, levels = names(pal)),
        time = 0, time_x = 0, estimate = 1, conf.low = 1, conf.high = 1,
        n.event = 0L, n.censor = 0L, n.risk = NA_integer_,
        std.error = 0, strata = paste0('source=', src_names))
      overlay_km_df <- bind_rows(ov_origins, overlay_km_df)
      overlay_km_df <- overlay_km_df %>%
        group_by(source) %>%
        group_modify(~ {
          lr <- .x %>% slice_tail(n = 1)
          if (lr$time_x < x_max) {
            ext <- lr; ext$time_x <- x_max; ext$time <- x_max * 12
            ext$n.event <- 0L; ext$n.censor <- 0L
            bind_rows(.x, ext)
          } else .x
        }) %>% ungroup()
    }
  }

  # Median survival from full (uncensored) data
  med_tbl <- summary(km_fit_full)$table
  if (is.null(dim(med_tbl)))
    med_tbl <- matrix(med_tbl, nrow = 1,
                      dimnames = list(names(km_fit_full$strata), names(med_tbl)))
  med_df <- tibble(
    source = str_remove(rownames(med_tbl), 'source='),
    median = med_tbl[, 'median'],
    lcl    = med_tbl[, '0.95LCL'],
    ucl    = med_tbl[, '0.95UCL']
  ) %>%
    mutate(
      med_str = if_else(is.na(median), 'NR', sprintf('%.1f', median)),
      ci_str  = sprintf('[%s\u2013%s]',
                        if_else(is.na(lcl), 'NR', sprintf('%.1f', lcl)),
                        if_else(is.na(ucl), 'NR', sprintf('%.1f', ucl))))

  # Annotation box -- two-column: source name left, median [CI] right
  box_left <- x_max * 0.66
  box_right <- x_max * 0.99
  box_pad  <- x_max * 0.01
  box_top  <- 0.99
  line_h   <- 0.036
  box_h    <- line_h * (nrow(med_df) + 1)

  med_annotations <- list(
    annotate('rect',
             xmin = box_left, xmax = box_right,
             ymin = box_top - box_h, ymax = box_top,
             fill = 'white', color = 'grey40', linewidth = 0.35),
    annotate('text', x = (box_left + box_right) / 2,
             y = box_top - line_h * 0.5,
             hjust = 0.5, vjust = 0.5, size = 2.8, fontface = 'bold',
             family = 'Helvetica', color = 'black',
             label = 'Median Survival (months)'))

  for (i in seq_len(nrow(med_df))) {
    src  <- med_df$source[i]
    col  <- pal[src]
    vals <- sprintf('%s  %s', med_df$med_str[i], med_df$ci_str[i])
    y_row <- box_top - line_h * (i + 0.5)
    med_annotations <- c(med_annotations, list(
      annotate('text', x = box_left + box_pad, y = y_row,
               hjust = 0, vjust = 0.5, size = 2.6,
               family = 'Helvetica', color = col, label = src),
      annotate('text', x = box_right - box_pad, y = y_row,
               hjust = 1, vjust = 0.5, size = 2.6,
               family = 'Helvetica', color = col, label = vals)))
  }

  # Risk table (from plot-censored data)
  risk_tbl <- expand_grid(src = src_names, time_x = x_breaks) %>%
    rowwise() %>%
    mutate(n_risk = sum(data_plot$T_x[data_plot$source == src] >= time_x,
                        na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(source = factor(src, levels = names(pal)))

  p_main <- ggplot(km_df,
                   aes(x = time_x, y = estimate,
                       color = source, fill = source)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                alpha = 0.07, color = NA) +
    { if (!is.null(overlay_km_df))
        geom_step(data = overlay_km_df,
                  aes(x = time_x, y = estimate, color = source),
                  linetype = 'dashed', linewidth = 0.5, alpha = 0.6,
                  inherit.aes = FALSE) } +
    geom_step(linewidth = 0.8) +
    geom_segment(data = censor_df,
                 aes(x = time_x, xend = time_x,
                     y = estimate - 0.020, yend = estimate + 0.020),
                 linewidth = 0.6, show.legend = FALSE) +
    geom_hline(yintercept = 0.5, linetype = 'dashed',
               color = 'grey30', linewidth = 0.4) +
    med_annotations +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    scale_x_continuous(breaks = x_breaks,
                       expand = expansion(mult = c(0, 0))) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, 1.01)) +
    scale_y_continuous(labels = scales::percent_format(1)) +
    labs(title = title_str,
         x = x_lab, y = y_label,
         color = NULL, fill = NULL) +
    theme_nature(base_size = 8) +
    theme(legend.position = 'none',
          plot.margin = margin(5.5, 5.5, 0, 5.5))

  p_risk <- ggplot(risk_tbl,
                   aes(x = time_x, y = fct_rev(source),
                       label = n_risk, color = source)) +
    geom_text(size = 2.0) +
    scale_color_manual(values = pal) +
    scale_x_continuous(breaks = x_breaks,
                       expand = expansion(mult = c(0, 0))) +
    coord_cartesian(xlim = c(0, x_max), clip = 'off') +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 4) +
    theme(
      axis.text.y      = if (show_risk_labels)
                           element_text(size = 4.5, hjust = 1, color = 'grey40',
                                        margin = margin(r = 5))
                         else element_blank(),
      axis.title       = element_blank(),
      legend.position  = 'none',
      plot.margin      = margin(0, 5.5, 0, 5.5))

  p_main / p_risk + plot_layout(heights = c(16, 1))
}

plot_confusion <- function(pred, truth, levels, title_str, fill_color = '#2166AC',
                           x_label = 'LLM') {
  valid <- !is.na(pred) & !is.na(truth)
  tbl <- table(
    Oncologist = factor(truth[valid], levels = levels),
    LLM        = factor(pred[valid],  levels = levels)
  ) %>%
    as.data.frame() %>%
    rename(count = Freq) %>%
    mutate(
      Oncologist = factor(as.character(Oncologist), levels = levels),
      LLM        = factor(as.character(LLM),        levels = levels))

  ggplot(tbl, aes(x = LLM, y = fct_rev(Oncologist), fill = count)) +
    geom_tile(color = 'white', linewidth = 0.5) +
    geom_text(aes(label = count), size = 3.5) +
    scale_fill_gradient(low = 'white', high = fill_color, guide = 'none') +
    labs(title = title_str, x = x_label, y = 'Oncologist') +
    theme_nature(base_size = 8) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          axis.line = element_blank(),
          plot.margin = margin(8, 10, 8, 10))
}

plot_confusion_heatmap <- function(var_col, gt_col, levels, title_str,
                                   all_data, gt_df, llm_ids) {
  combo <- map_dfr(llm_ids, function(llm_id) {
    m <- inner_join(all_data[[llm_id]], gt_df, by = 'gsID')
    pred  <- m[[var_col]]
    truth <- m[[gt_col]]
    valid <- !is.na(pred) & !is.na(truth)
    if (sum(valid) == 0) return(tibble())
    table(
      Oncologist = factor(truth[valid], levels = levels),
      Predicted  = factor(pred[valid],  levels = levels)
    ) %>% as.data.frame() %>% rename(count = Freq)
  }) %>%
    mutate(
      Oncologist = factor(as.character(Oncologist), levels = levels),
      Predicted  = factor(as.character(Predicted),  levels = levels)
    ) %>%
    group_by(Oncologist, Predicted) %>%
    summarise(count = sum(count), .groups = 'drop')

  combo <- combo %>%
    group_by(Oncologist) %>%
    mutate(pct = count / sum(count)) %>%
    ungroup() %>%
    mutate(
      label   = if_else(count > 0, sprintf('%d\n(%.0f%%)', count, pct * 100), ''),
      on_diag = as.character(Oncologist) == as.character(Predicted))

  ggplot(combo, aes(x = Predicted, y = fct_rev(Oncologist))) +
    geom_tile(aes(fill = pct), color = 'white', linewidth = 0.8) +
    geom_tile(data = combo %>% filter(on_diag),
              fill = NA, color = 'black', linewidth = 0.9) +
    geom_text(aes(label = label), size = 3, lineheight = 0.85,
              color = if_else(combo$pct > 0.5, 'white', 'black')) +
    scale_fill_gradient(low = 'white', high = '#2166AC',
                        limits = c(0, 1), guide = 'none') +
    labs(title = title_str, x = 'LLM Extraction', y = 'Oncologist (Expert-derived)') +
    theme_nature(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'none', panel.grid = element_blank())
}

build_concordance <- function(llm_id, llm_label, gt_df, all_llm_data, drug_df) {
  d <- all_llm_data[[llm_id]]
  m <- left_join(gt_df, d, by = 'gsID', suffix = c('', '_llm'))

  date_check <- function(llm_col, gt_col) {
    if (!llm_col %in% names(m)) return(rep('missing', nrow(m)))
    llm_v <- tryCatch(as.Date(m[[llm_col]]), error = function(e) rep(as.Date(NA), nrow(m)))
    gt_v  <- tryCatch(as.Date(m[[gt_col]]),  error = function(e) rep(as.Date(NA), nrow(m)))
    case_when(
      is.na(gt_v) & is.na(llm_v) ~ 'correct',
      is.na(gt_v) | is.na(llm_v) ~ 'missing',
      abs(as.numeric(llm_v - gt_v)) <= TOL_PRIMARY ~ 'correct',
      TRUE ~ 'incorrect')
  }

  cat_check <- function(llm_col, gt_col) {
    case_when(
      is.na(m[[gt_col]]) & is.na(m[[llm_col]]) ~ 'correct',
      is.na(m[[gt_col]]) | is.na(m[[llm_col]]) ~ 'missing',
      m[[llm_col]] == m[[gt_col]] ~ 'correct',
      TRUE ~ 'incorrect')
  }

  # Genetics: mutated vs other -- binarize_mutation maps NA -> 'other', so columns
  # are never NA here.  case_when is used (not if_else) to avoid NA propagation
  # if that invariant is ever broken upstream.
  gene_check <- function(llm_col, gt_col) {
    if (!llm_col %in% names(m) || !gt_col %in% names(m))
      return(rep('missing', nrow(m)))
    case_when(
      m[[llm_col]] == m[[gt_col]] ~ 'correct',
      TRUE                        ~ 'incorrect'
    )
  }

  conc <- tibble(gsID = m$gsID)
  conc$init_dx    <- date_check('init_dx', 'gt_init_dx')
  conc$lrr_dx     <- date_check('lrr_dx',  'gt_lrr_dx')
  conc$met_dx     <- date_check('met_dx',  'gt_met_dx')
  conc$death_date <- date_check('death_date', 'gt_death')
  conc$stage      <- if ('stage' %in% names(m)) cat_check('stage', 'gt_stage') else 'missing'
  conc$subtype    <- if ('subtype' %in% names(m)) cat_check('subtype', 'gt_subtype') else 'missing'
  conc$BRCA       <- gene_check('brca', 'gt_brca')
  conc$PIK3CA     <- gene_check('pik3ca', 'gt_pik3ca')
  conc$ESR1       <- gene_check('esr1', 'gt_esr1')

  # Treatment line count
  gt_lc  <- drug_df %>% filter(version == 'v1') %>%
    count(gsID = as.character(gsID), name = 'gt_n')
  llm_lc <- drug_df %>% filter(version == llm_id) %>%
    count(gsID = as.character(gsID), name = 'llm_n')
  lc <- full_join(gt_lc, llm_lc, by = 'gsID')
  lc_map <- setNames(
    case_when(
      is.na(lc$gt_n) & is.na(lc$llm_n) ~ 'correct',
      is.na(lc$gt_n) | is.na(lc$llm_n) ~ 'missing',
      lc$gt_n == lc$llm_n ~ 'correct',
      TRUE ~ 'incorrect'), lc$gsID)
  conc$line_count <- lc_map[conc$gsID]
  conc$line_count[is.na(conc$line_count)] <- 'missing'

  conc$llm <- llm_label
  conc
}

plot_forest <- function(data, strat_order, palette, shapes, title,
                        dodge_width = 0.40, output_width = 8.5,
                        output_height = 7.5, base_size = 9, caption = NULL) {

  Y_GAP <- 3.0

  plot_df <- data %>%
    filter(!is.na(est), strategy %in% strat_order) %>%
    mutate(
      n_suffix  = str_extract(variable, '\\s*\\(n=\\d+\\)'),
      var_clean = str_remove(variable, '\\s*\\(n=\\d+\\)'),
      var_label = if_else(var_clean %in% names(VAR_LABELS),
                          VAR_LABELS[var_clean], var_clean)) %>%
    group_by(var_clean) %>%
    mutate(n_display = max(n, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      var_label = if_else(n_display != 100 & var_label %in% c(
        'Line of Therapy Intent (curative vs palliative)',
        'Treatment Stop Reason (progression vs other)'),
        paste0(var_label, '*'), var_label),
      strategy  = factor(strategy, levels = strat_order),
      lo = pmax(lo, 0), hi = pmin(hi, 1))

  # Numeric y positions: DOB pinned to top, rest by GPT-5 est
  PIN_TOP <- c('Date of Birth')
  gpt5_rank <- plot_df %>%
    filter(strategy == 'GPT-5', !var_label %in% PIN_TOP,
           !str_remove(var_label, '\\*$') %in% PIN_TOP) %>%
    distinct(var_label, est) %>%
    arrange(est)   # ascending: lowest est gets y=1 (bottom), highest gets top
  all_labels <- unique(plot_df$var_label)
  missing <- setdiff(all_labels, c(gpt5_rank$var_label,
                                    all_labels[str_remove(all_labels, '\\*$') %in% PIN_TOP]))
  pinned <- all_labels[str_remove(all_labels, '\\*$') %in% PIN_TOP]
  # Order: bottom → top = worst GPT-5 … best GPT-5 … pinned items
  var_order <- c(missing, gpt5_rank$var_label, rev(pinned))
  y_pos_map <- setNames(seq_along(var_order) * Y_GAP, var_order)

  # Explicit per-model y offsets -- top of strat_order gets highest y (most positive)
  y_offsets <- setNames(
    seq(dodge_width / 2, -dodge_width / 2, length.out = length(strat_order)),
    strat_order)

  plot_df <- plot_df %>%
    mutate(y_base   = y_pos_map[var_label],
           y_dodged = y_base + y_offsets[as.character(strategy)])

  ggplot(plot_df,
         aes(color = strategy, shape = strategy)) +
    geom_vline(xintercept = 1, linetype = 'dashed',
               color = 'grey70', linewidth = 0.25) +
    geom_segment(aes(x = lo, xend = hi, y = y_dodged, yend = y_dodged),
                 linewidth = 1.15) +
    geom_point(aes(x = est, y = y_dodged), size = 2.8) +
    scale_color_manual(values = palette, breaks = strat_order) +
    scale_shape_manual(values = shapes, breaks = strat_order) +
    scale_x_continuous(labels = scales::percent_format(),
                       breaks = seq(0.5, 1.0, 0.1)) +
    coord_cartesian(xlim = c(0.48, 1.02)) +   # zoom, don't drop out-of-range CIs
    scale_y_continuous(breaks = y_pos_map,
                       labels = names(y_pos_map),
                       expand = expansion(add = Y_GAP * 0.6)) +
    labs(title = title,
         x = 'Agreement (95% Bootstrap CI)',
         y = NULL, color = NULL, shape = NULL,
         caption = caption) +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 2.5)),
           shape = guide_legend(nrow = 1)) +
    theme_nature(base_size = base_size) +
    theme(legend.position   = 'bottom',
          panel.grid.major.y = element_blank(),
          axis.text.y        = element_text(size = 8.5, face = 'bold'),
          plot.title         = element_text(hjust = 0.5),
          plot.caption       = element_text(size = 6, color = 'grey40', hjust = 0))
}

plot_jaccard_violin <- function(data, summary_data, title, palette,
                                level_order,
                                subtitle = 'Diamond = median; error bar = IQR',
                                y_label = 'Jaccard Similarity') {
  ggplot(data, aes(x = label, y = jaccard,
                   fill = label, color = label)) +
    geom_violin(alpha = 0.30, linewidth = 0.45, trim = TRUE) +
    geom_jitter(width = 0.10, size = 1.3, alpha = 0.40, shape = 16) +
    geom_point(data = summary_data,
               mapping = aes(x = label, y = m),
               shape = 18, size = 3.5, color = 'black', inherit.aes = FALSE) +
    geom_errorbar(data = summary_data,
                  mapping = aes(x = label, ymin = lo, ymax = hi),
                  color = 'black', width = 0.14, linewidth = 0.55,
                  inherit.aes = FALSE) +
    scale_fill_manual(values  = palette, guide = 'none') +
    scale_color_manual(values = palette, guide = 'none') +
    scale_y_continuous(labels = scales::percent_format(),
                       limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title    = title,
         subtitle = subtitle,
         x = NULL, y = y_label) +
    theme_nature(base_size = 8) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x        = element_text(size = 8.8))
}

### Inputs -- all data reads at top

# Final dfs
non_drug_long <- read_excel(file.path(FINAL_DIR, 'non_drug_long.xlsx'))
drug_df       <- read_excel(file.path(FINAL_DIR, 'final_drug_df_for_comparisonfull.xlsx'))
pat_acc       <- read_excel(file.path(FINAL_DIR, 'patient_level_accuracy_by_version.xlsx'))
vec_rank      <- read_excel(file.path(FINAL_DIR, 'llm_vector_rank.xlsx'))
per_pat_err   <- read_excel(file.path(LAB_MTG_DIR, 'non_drug_per_patient_errors.xlsx'))

# Apply censor date: remove lines started after each source's review cutoff
# (CRC_CENSOR_DATE for CRC rows; patient's last_note_date for v1/v2/LLMs)
drug_df <- drug_df %>%
  mutate(Start       = as.Date(Start),
         censor_date = as.Date(censor_date)) %>%
  filter(is.na(Start) | is.na(censor_date) | Start <= censor_date)

# Cleaned intermediates -- all data reads from intermediate/ or final_dfs/ only
gt_cleaned  <- read_xlsx(GT_NON_DRUG_CLEANED)
crc_cleaned <- read_xlsx(CRC_NON_DRUG_CLEANED)
llm_cleaned <- read_xlsx(LLM_NON_DRUG_CLEANED)

# LLM death dates (from cleaned intermediate)
llm_death <- llm_cleaned %>%
  transmute(gsID = as.character(gsID), version = version,
            death_date = as.Date(date_of_death))

### Transformations

# GT reference with survival data
gt <- non_drug_long %>%
  filter(version == 'gt') %>%
  transmute(
    gsID,
    gt_init_dx  = as.Date(date_initial_diagnosis),
    gt_lrr_dx   = as.Date(date_first_locoregional_recurrence),
    gt_met_dx   = as.Date(date_first_metastatic_diagnosis),
    gt_birth    = as.Date(birth_date),
    gt_stage    = as.character(stage_at_diagnosis),
    gt_subtype  = hr_her2_status,
    gt_hr = case_when(
      str_detect(coalesce(hr_her2_status, ''), 'HR\\+') ~ 'HR+',
      str_detect(coalesce(hr_her2_status, ''), 'HR-')   ~ 'HR-',
      TRUE ~ NA_character_),
    gt_her2 = case_when(
      str_detect(coalesce(hr_her2_status, ''), 'HER2\\+') ~ 'HER2+',
      str_detect(coalesce(hr_her2_status, ''), 'HER2-')   ~ 'HER2-',
      TRUE ~ NA_character_),
    gt_primary_onc = tolower(trimws(as.character(primary_oncologist_shc))),
    gt_brca   = binarize_mutation(brca_status),
    gt_pik3ca = binarize_mutation(pik3ca_status),
    gt_esr1   = binarize_mutation(esr1_status)
  ) %>%
  left_join(
    gt_cleaned %>% transmute(
      gsID         = as.character(gsID),
      gt_death     = parse_gt_date(death_date_clinician),
      gt_last_note = parse_gt_date(last_note_date)),
    by = 'gsID')

# Per-LLM data list
all_llm_data <- map(setNames(ALL_LLMS, ALL_LLMS), function(llm_id) {
  d <- non_drug_long %>%
    filter(version == llm_id) %>%
    transmute(
      gsID,
      init_dx  = as.Date(date_initial_diagnosis),
      lrr_dx   = as.Date(date_first_locoregional_recurrence),
      met_dx   = as.Date(date_first_metastatic_diagnosis),
      birth    = as.Date(birth_date),
      stage    = as.character(stage_at_diagnosis),
      subtype  = hr_her2_status,
      cap_hr = case_when(
        str_detect(coalesce(hr_her2_status, ''), 'HR\\+') ~ 'HR+',
        str_detect(coalesce(hr_her2_status, ''), 'HR-')   ~ 'HR-',
        TRUE ~ NA_character_),
      cap_her2 = case_when(
        str_detect(coalesce(hr_her2_status, ''), 'HER2\\+') ~ 'HER2+',
        str_detect(coalesce(hr_her2_status, ''), 'HER2-')   ~ 'HER2-',
        TRUE ~ NA_character_),
      brca   = binarize_mutation(brca_status),
      pik3ca = binarize_mutation(pik3ca_status),
      esr1   = binarize_mutation(esr1_status),
      primary_onc = tolower(trimws(as.character(primary_oncologist_shc))))

  # Add death dates (from cleaned intermediate)
  death <- llm_death %>% filter(version == llm_id)
  if (nrow(death) > 0) {
    d <- d %>% left_join(death %>% select(gsID, death_date), by = 'gsID')
  } else {
    d$death_date <- as.Date(NA)
  }

  d
})

# CRC data for Fig 2 (from cleaned intermediate)
crc_fig2 <- crc_cleaned %>%
  transmute(
    gsID    = as.character(gsID),
    init_dx = as.Date(date_initial_diagnosis),
    met_dx  = as.Date(date_first_metastatic_diagnosis))

# Drug comparison sets (for Jaccard in Fig 2 and S8)
drug_sets <- drug_df %>%
  filter(!is.na(regimen_bag_of_drugs) & regimen_bag_of_drugs != '') %>%
  mutate(gsID = as.character(gsID)) %>%
  group_by(gsID, version) %>%
  summarise(
    drugs = list(unique(trimws(tolower(unlist(str_split(regimen_bag_of_drugs, ',')))))),
    .groups = 'drop')

# Survival datasets
stage4_pids <- gt %>% filter(gt_stage == '4') %>% pull(gsID)

### Shared palettes and constants

# Palette and shapes
PAL_FIG2 <- c('Second Oncologist' = '#000000',
              PAL_LLM,
              'Research Coordinator' = COL_RC)

SHAPE_FIG2 <- c('Second Oncologist'    = 4,    # X
                'Gemini 2.5 Pro'           = 8,    # asterisk
                'GPT-5'                = 15,   # filled square
                'GPT-4o'               = 17,   # filled triangle
                'DeepSeek-R1'          = 18,   # filled diamond
                'Research Coordinator' = 1)    # open circle

# Clean variable labels -- descriptive enough without domain headers
VAR_LABELS <- c(
  'Date of Birth'                   = 'Date of Birth',
  'Initial Date Dx'                 = 'Initial Diagnosis Date',
  'Patient Had Recurrence'          = 'Recurrence Occurred (Yes/No)',
  'Locoregional Recurrence Date Dx' = 'Locoregional Recurrence Date',
  'Metastatic Date Dx'              = 'Metastatic Recurrence Date',
  'Date of Death'                   = 'Date of Death Documented',
  'Initial Cancer Stage'            = 'Clinical Stage',
  'Hormone Receptor Status'         = 'Hormone Receptor Status',
  'HER2 Receptor Status'            = 'HER2 Status',
  'Primary Oncologist'              = 'Treating Oncologist Identified',
  'Germline BRCA'                            = 'Germline BRCA Pathogenic Variant',
  'Somatic PIK3CA'                           = 'Somatic PIK3CA Mutated',
  'Somatic ESR1'                             = 'Somatic ESR1 Mutated',
  'Curative v Palliative Intent'             = 'Line of Therapy Intent (curative vs palliative)',
  'Progression vs Other Reason for Stopping' = 'Treatment Stop Reason (progression vs other)')

# Shared violin plot labels and palettes for Fig 4 + lab meeting taxane/AI version
S16_LABELS <- c(v2          = 'Second\nOncologist',
                gemini_2_5  = 'Gemini 2.5 Pro',
                gpt5        = 'GPT-5',
                gpt4o       = 'GPT-4o',
                deepseek_r1 = 'DeepSeek-R1',
                crc         = 'Research\nCoordinator')

S16_PAL <- c('Second\nOncologist' = '#FF7F00',   # orange -- distinct from GPT-4o pink
             PAL_LLM,
             'Research\nCoordinator' = COL_RC)

S16_ORDER <- c('Second\nOncologist',
               'Gemini 2.5 Pro', 'GPT-5', 'GPT-4o', 'DeepSeek-R1',
               'Research\nCoordinator')

### Fig 2: Forest plot -- extraction agreement across variables by LLM

# Estimates and CIs from manuscript_tables.R
table_variables_for_figure <- c(
  'Date of birth'                   = 'Date of Birth',
  'Date of initial diagnosis'       = 'Initial Date Dx',
  'Recurrence occurred'             = 'Patient Had Recurrence',
  'Date of locoregional recurrence' = 'Locoregional Recurrence Date Dx',
  'Date of metastatic diagnosis'    = 'Metastatic Date Dx',
  'Exact date of death documented'  = 'Date of Death',
  'Initial cancer stage'            = 'Initial Cancer Stage',
  'HR status'                       = 'Hormone Receptor Status',
  'HER2 status'                     = 'HER2 Receptor Status',
  'Primary oncologist'              = 'Primary Oncologist',
  'Germline BRCA'                   = 'Germline BRCA',
  'Somatic PIK3CA'                  = 'Somatic PIK3CA',
  'Somatic ESR1'                    = 'Somatic ESR1')

variable_rows <- read_xlsx(file.path(TBL_DIR, 'TableS2_full_model_performance.xlsx'),
                           sheet = 'detail_with_counts') %>%
  filter(Variable %in% names(table_variables_for_figure),
         Tolerance %in% c(paste0('\u00b1', TOL_PRIMARY, ' days'), 'Exact match')) %>%
  transmute(variable = unname(table_variables_for_figure[Variable]), domain = Category,
            strategy = recode(LLM, 'Gemini 2.5' = 'Gemini 2.5 Pro'),
            n = `Evaluable N`, est, lo, hi)

line_strategy_labels <- c(v2 = 'Second Oncologist', LLM_LABELS, crc = 'Research Coordinator')
matched_line_rows <- read_xlsx(file.path(TBL_DIR, 'Table3_line_metrics.xlsx'),
                               sheet = 'matched_line_agreement') %>%
  transmute(variable = recode(attribute, 'Treatment intent' = 'Curative v Palliative Intent',
                              'Reason for discontinuation' = 'Progression vs Other Reason for Stopping'),
            domain = 'Systemic Therapy', strategy = unname(line_strategy_labels[model]),
            n = matched_lines, est, lo, hi)

fig2_df <- bind_rows(variable_rows, matched_line_rows)
# -- Console: print n values per variable per strategy (for manuscript footnotes) --
cat('\n=== Fig 2 / FigS1 -- n values by variable and strategy ===\n')
n_tbl <- fig2_df %>%
  filter(!is.na(est)) %>%
  mutate(var_clean = str_remove(variable, '\\s*\\(n=\\d+\\)')) %>%
  select(var_clean, strategy, n) %>%
  arrange(var_clean, strategy)
print(as.data.frame(n_tbl), row.names = FALSE)

cat('\n--- Summary: max n per variable (used for axis label) ---\n')
print(as.data.frame(
  n_tbl %>% group_by(var_clean) %>% summarise(max_n = max(n, na.rm = TRUE), .groups = 'drop')
), row.names = FALSE)

# -- Inline stats: LRR false positive / false negative counts --
# Needed for manuscript text: "0 and 1 false positive locoregional events
# and 3 and 5 false negative (missed events)"
cat('\n=== LRR FP/FN counts per LLM (for manuscript inline stats) ===\n')
lrr_fp_fn <- map_dfr(names(all_llm_data), function(llm_id) {
  d <- all_llm_data[[llm_id]]
  m <- inner_join(d %>% select(gsID, lrr_dx),
                  gt %>% select(gsID, gt_lrr_dx), by = 'gsID')
  tibble(
    llm      = LLM_LABELS[[llm_id]],
    n_gt_lrr = sum(!is.na(m$gt_lrr_dx)),                          # Total GT LRR events
    n_fp_lrr = sum(!is.na(m$lrr_dx) & is.na(m$gt_lrr_dx)),       # LLM says LRR, GT has none
    n_fn_lrr = sum( is.na(m$lrr_dx) & !is.na(m$gt_lrr_dx)),      # GT has LRR, LLM missed it
    n        = nrow(m)
  )
})
print(as.data.frame(lrr_fp_fn), row.names = FALSE)

# Fig 2: Gemini 2.5 Pro and GPT-5 only (two best-performing LLMs).
# All four LLMs + CRC + second oncologist are shown in FigS1 (expanded forest; Supplementary Figure S1 in manuscript).
if (nrow(fig2_df) > 0) {
  STRAT_ORDER <- c('Second Oncologist', 'Gemini 2.5 Pro', 'GPT-5', 'Research Coordinator')

  p_fig2 <- plot_forest(fig2_df, STRAT_ORDER, PAL_FIG2, SHAPE_FIG2,
                         'Agreement with Study Reference Comparator',
                         dodge_width = 1.40, output_width = 8.5,
                         caption = NULL)

  # Print n values for starred (*) variables -- for figure legend
  fig2_plot_tmp <- fig2_df %>%
    filter(!is.na(est), strategy %in% STRAT_ORDER) %>%
    mutate(
      n_suffix  = str_extract(variable, '\\s*\\(n=\\d+\\)'),
      var_clean = str_remove(variable, '\\s*\\(n=\\d+\\)'),
      var_label = if_else(var_clean %in% names(VAR_LABELS),
                          VAR_LABELS[var_clean], var_clean)) %>%
    group_by(var_clean) %>%
    mutate(n_display = max(n, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(var_label = if_else(n_display != 100 & var_label %in% c(
      'Line of Therapy Intent (curative vs palliative)',
      'Treatment Stop Reason (progression vs other)'),
      paste0(var_label, '*'), var_label))

  starred_n_tbl <- fig2_plot_tmp %>%
    filter(n_display != 100) %>%
    distinct(strategy, var_label, n) %>%
    arrange(var_label, strategy)
  cat('\n=== Fig 2 -- n values for starred (*) variables (for figure legend) ===\n')
  print(as.data.frame(starred_n_tbl), row.names = FALSE)

  ggsave(file.path(FIG_DIR, 'Fig2_forest_plot.png'), p_fig2,
         width = 8.5, height = 9.5, dpi = 300, bg = 'white')

  # FigS1: expanded forest -- all four LLMs + second oncologist + CRC
  STRAT_ORDER_ALL <- c('Second Oncologist', 'Gemini 2.5 Pro', 'GPT-5', 'GPT-4o', 'DeepSeek-R1', 'Research Coordinator')

  p_figs1 <- plot_forest(fig2_df, STRAT_ORDER_ALL, PAL_FIG2, SHAPE_FIG2,
                          'Agreement with Study Reference Comparator -- All LLMs',
                          dodge_width = 1.80, output_width = 9.5,
                          output_height = 7.5, caption = NULL)

  ggsave(file.path(FIG_DIR, 'FigS1_expanded_forest.png'), p_figs1,
         width = 9.5, height = 9.5, dpi = 300, bg = 'white')
}

### Fig 3: KM survival curves (OS + RFS) -- lower bound (censored -> dead)
# Patients without a documented death date are treated as having died at their
# last clinical follow-up (worst-case / lower-bound survival).  Dashed overlay
# curves show the standard Kaplan-Meier (censored at last follow-up).

KM_LLMS   <- c('gemini_2_5', 'gpt5')
PAL_KM_F3 <- c('Expert-derived' = '#000000',
                'Gemini 2.5 Pro'  = unname(PAL_LLM['Gemini 2.5 Pro']),
                'GPT-5'       = unname(PAL_LLM['GPT-5']))

# GT KM dataset -- all dates from ground truth
gt_km <- gt %>% transmute(
  gsID,
  init_dx    = gt_init_dx,
  lrr_dx     = gt_lrr_dx,
  met_dx     = gt_met_dx,
  death_date = gt_death,
  last_note  = gt_last_note)

# -- Standard censoring (retained for the overlay) ----------------------------
gt_os  <- build_os(gt_km, 'init_dx', 'Expert-derived')
gt_rfs <- build_rfs(gt_km, 'init_dx', 'lrr_dx', 'met_dx', 'Expert-derived', stage4_pids)

os_list  <- list(gt_os)
rfs_list <- list(gt_rfs)
for (llm_id in KM_LLMS) {
  d   <- all_llm_data[[llm_id]]
  lab <- LLM_LABELS[llm_id]
  # Each LLM uses its own dx, death, and recurrence dates;
  # last_note is a chart-level property (not source-specific)
  d_km <- d %>%
    left_join(gt %>% select(gsID, last_note = gt_last_note), by = 'gsID')
  if ('init_dx' %in% names(d_km))
    os_list <- c(os_list, list(build_os(d_km, 'init_dx', lab)))
  if (all(c('lrr_dx', 'met_dx') %in% names(d_km))) {
    stage4_llm <- d_km$gsID[coalesce(d_km$stage == '4', FALSE)]  # each LLM's own Stage 4
    rfs_list <- c(rfs_list,
                  list(build_rfs(d_km, 'init_dx', 'lrr_dx', 'met_dx',
                                 lab, stage4_llm)))
  }
}

all_os  <- bind_rows(os_list)  %>% mutate(source = factor(source, levels = names(PAL_KM_F3)))
all_rfs <- bind_rows(rfs_list) %>% mutate(source = factor(source, levels = names(PAL_KM_F3)))

# -- Lower-bound builds (censored patients treated as dead) -------------------
build_os_lb <- function(df, dx_col, source_label) {
  df %>%
    mutate(
      dx_date  = .data[[dx_col]],
      event    = 1L,
      end_date = if_else(!is.na(death_date), death_date, last_note),
      T_days   = as.numeric(difftime(end_date, dx_date, units = 'days')),
      T_months = T_days / 30.4375
    ) %>%
    filter(!is.na(dx_date), !is.na(end_date), T_days >= 0) %>%
    mutate(source = source_label) %>%
    select(gsID, source, event, T_months)
}

build_rfs_lb <- function(df, dx_col, lrr_col, met_col,
                         source_label, stage4_pids) {
  df %>%
    filter(!(gsID %in% stage4_pids)) %>%
    rowwise() %>%
    mutate(
      dx_date  = .data[[dx_col]],
      recur    = {
        cands <- na.omit(c(.data[[lrr_col]], .data[[met_col]]))
        if (length(cands) > 0) min(cands) else as.Date(NA)
      },
      event    = 1L,
      end_date = if (!is.na(recur) & !is.na(dx_date) & recur >= dx_date)
                   recur else coalesce(death_date, last_note),
      T_days   = as.numeric(difftime(end_date, dx_date, units = 'days')),
      T_months = T_days / 30.4375
    ) %>%
    ungroup() %>%
    filter(!is.na(dx_date), !is.na(end_date), T_days >= 0) %>%
    mutate(source = source_label) %>%
    select(gsID, source, event, T_months)
}

gt_os_lb  <- build_os_lb(gt_km, 'init_dx', 'Expert-derived')
gt_rfs_lb <- build_rfs_lb(gt_km, 'init_dx', 'lrr_dx', 'met_dx', 'Expert-derived', stage4_pids)

os_lb_list  <- list(gt_os_lb)
rfs_lb_list <- list(gt_rfs_lb)
for (llm_id in KM_LLMS) {
  d   <- all_llm_data[[llm_id]]
  lab <- LLM_LABELS[llm_id]
  d_km <- d %>%
    left_join(gt %>% select(gsID, last_note = gt_last_note), by = 'gsID')
  if ('init_dx' %in% names(d_km))
    os_lb_list <- c(os_lb_list, list(build_os_lb(d_km, 'init_dx', lab)))
  if (all(c('lrr_dx', 'met_dx') %in% names(d_km))) {
    stage4_llm <- d_km$gsID[coalesce(d_km$stage == '4', FALSE)]  # each LLM's own Stage 4
    rfs_lb_list <- c(rfs_lb_list,
                     list(build_rfs_lb(d_km, 'init_dx', 'lrr_dx', 'met_dx',
                                       lab, stage4_llm)))
  }
}

all_os_lb  <- bind_rows(os_lb_list)  %>% mutate(source = factor(source, levels = names(PAL_KM_F3)))
all_rfs_lb <- bind_rows(rfs_lb_list) %>% mutate(source = factor(source, levels = names(PAL_KM_F3)))

# -- Plot lower-bound as Fig3, overlay standard curves as dashed --------------
p_os <- plot_km(all_os_lb, 'Overall Survival',
  PAL_KM_F3, 'OS Probability', max_yr = 12, use_months = FALSE,
  show_risk_labels = TRUE, overlay_data = all_os)

p_rfs <- plot_km(all_rfs_lb, 'Recurrence-Free Survival (Stage 1\u20133 only)',
  PAL_KM_F3, 'RFS Probability', max_yr = 5, use_months = TRUE,
  show_risk_labels = FALSE, overlay_data = all_rfs)

p_fig3 <- (p_os | p_rfs)
ggsave(file.path(FIG_DIR, 'Fig3_KM_curves.png'), p_fig3,
       width = 16, height = 6.5, dpi = 300, bg = 'white')

# -- Inline stats: OS and RFS medians (for manuscript text)
cat('\n=== Fig 3 KM inline stats (lower bound, for manuscript text) ===\n')
for (km_obj in list(list(data = all_os_lb, label = 'Overall Survival'),
                    list(data = all_rfs_lb, label = 'Recurrence-Free Survival'))) {
  km_data <- km_obj$data
  km_fit_full <- survfit(Surv(T_months, event) ~ source, data = km_data)
  med_tbl <- summary(km_fit_full)$table
  if (is.null(dim(med_tbl)))
    med_tbl <- matrix(med_tbl, nrow = 1,
                      dimnames = list(names(km_fit_full$strata), names(med_tbl)))

  cat(sprintf('\n  %s:\n', km_obj$label))
  for (i in seq_len(nrow(med_tbl))) {
    src <- str_remove(rownames(med_tbl)[i], 'source=')
    cat(sprintf('    %-25s median = %.1f months\n', src, med_tbl[i, 'median']))
  }

  # 5-year survival estimates (per source)
  cat('    5-year estimates:\n')
  for (src in levels(km_data$source)) {
    src_data <- km_data %>% filter(source == src)
    km_src <- survfit(Surv(T_months, event) ~ 1, data = src_data)
    s5 <- summary(km_src, times = 60, extend = TRUE)
    cat(sprintf('      %-23s %.1f%% [%.1f\u2013%.1f]\n',
                src, s5$surv*100, s5$lower*100, s5$upper*100))
  }
}

### Fig S2: Dispersion of event dates

date_pairs <- list(
  'Initial Diagnosis'       = c('init_dx', 'gt_init_dx'),
  'Locoregional Recurrence' = c('lrr_dx',  'gt_lrr_dx'),
  'Metastatic Diagnosis'    = c('met_dx',  'gt_met_dx'))

disp_rows <- list()

# All LLMs
for (llm_id in ALL_LLMS) {
  m <- inner_join(all_llm_data[[llm_id]], gt, by = 'gsID')
  for (event_name in names(date_pairs)) {
    cols  <- date_pairs[[event_name]]
    if (!cols[1] %in% names(m)) next
    valid <- !is.na(m[[cols[1]]]) & !is.na(m[[cols[2]]])
    if (sum(valid) == 0) next
    disp_rows <- c(disp_rows, list(tibble(
      event = event_name, strategy = LLM_LABELS[llm_id],
      diff = as.numeric(m[[cols[1]]][valid] - m[[cols[2]]][valid]),
      gsID = m$gsID[valid])))
  }
}

# Research Coordinator
{
  crc_m <- crc_fig2 %>%
    inner_join(gt, by = 'gsID')
  for (pair in list(c('Initial Diagnosis', 'init_dx', 'gt_init_dx'),
                    c('Metastatic Diagnosis', 'met_dx', 'gt_met_dx'))) {
    valid <- !is.na(crc_m[[pair[2]]]) & !is.na(crc_m[[pair[3]]])
    if (sum(valid) == 0) next
    disp_rows <- c(disp_rows, list(tibble(
      event = pair[1], strategy = 'Research Coordinator',
      diff = as.numeric(crc_m[[pair[2]]][valid] - crc_m[[pair[3]]][valid]),
      gsID = crc_m$gsID[valid])))
  }
}

if (length(disp_rows) > 0) {
  pal_s2 <- c(PAL_LLM, 'Research Coordinator' = COL_RC)
  clip_mo <- 60

  disp_df <- bind_rows(disp_rows) %>%
    mutate(
      event    = factor(event, levels = c('Initial Diagnosis',
                                          'Locoregional Recurrence',
                                          'Metastatic Diagnosis')),
      strategy = factor(strategy, levels = names(pal_s2)),
      diff_mo  = pmax(pmin(diff / 30.4375, clip_mo), -clip_mo))

  make_disp_panel <- function(evt_label) {
    dd <- disp_df %>% filter(event == evt_label)
    ggplot(dd, aes(x = strategy, y = diff_mo, color = strategy)) +
      geom_hline(yintercept = 0, linetype = 'solid',
                 color = 'grey50', linewidth = 0.4) +
      geom_hline(yintercept = c(-3, 3), linetype = 'dashed',
                 color = 'grey50', linewidth = 0.4) +
      geom_hline(yintercept = c(-12, 12), linetype = 'dashed',
                 color = '#CC0000', linewidth = 0.4) +
      geom_jitter(width = 0.15, size = 1.3, alpha = 0.60) +
      scale_color_manual(values = pal_s2) +
      scale_y_continuous(breaks = seq(-60, 60, by = 12)) +
      coord_cartesian(ylim = c(-clip_mo, clip_mo)) +
      labs(title = evt_label, x = NULL, y = NULL) +
      theme_nature(base_size = 8) +
      theme(legend.position = 'none',
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  }

  panels <- map(levels(disp_df$event), make_disp_panel)

  p_s2 <- (panels[[1]] | panels[[2]] | panels[[3]]) +
    plot_annotation(
      title    = 'Concordance with Oncologist on Diagnosis Timing',
      theme = theme(
        plot.title    = element_text(face = 'bold', size = 11, hjust = 0.5)))

  p_s2 <- wrap_elements(p_s2) +
    labs(tag = 'Difference from Oncologist (months)') +
    theme(plot.tag = element_text(size = 9, angle = 90),
          plot.tag.position = 'left')

  ggsave(file.path(FIG_DIR, 'FigS2_dispersion.png'), p_s2,
         width = 14, height = 6, dpi = 300, bg = 'white')
}

### Fig S3: Confusion matrices (stage, subtype)

# FigS3: Gemini 2.5 Pro and GPT-5 confusion matrices, color-coded by LLM
CM_LLMS <- c('gemini_2_5', 'gpt5')

stage_plots <- map(CM_LLMS, function(llm_id) {
  m <- inner_join(all_llm_data[[llm_id]], gt, by = 'gsID')
  plot_confusion(m$stage, m$gt_stage, as.character(1:4),
                 '', fill_color = PAL_LLM[LLM_LABELS[llm_id]],
                 x_label = LLM_LABELS[llm_id])
})

subtype_plots <- map(CM_LLMS, function(llm_id) {
  m <- inner_join(all_llm_data[[llm_id]], gt, by = 'gsID')
  plot_confusion(m$subtype, m$gt_subtype,
                 c('HR+/HER2-', 'HR+/HER2+', 'HR-/HER2+', 'HR-/HER2-'), '',
                 fill_color = PAL_LLM[LLM_LABELS[llm_id]],
                 x_label = LLM_LABELS[llm_id])
})

p_s3 <- wrap_plots(c(stage_plots, subtype_plots), ncol = 2) +
  plot_annotation(
    title   = 'Confusion Matrices: Oncologist vs LLMs',
    caption = 'Top row: Cancer Stage (1-4)  |  Bottom row: HR/HER2 Subtype',
    theme = theme(
      plot.title   = element_text(face = 'bold', size = 11, hjust = 0.5,
                                  margin = margin(b = 6)),
      plot.caption = element_text(size = 8, color = 'grey40', hjust = 0.5)))

ggsave(file.path(FIG_DIR, 'FigS3_confusion_matrices.png'), p_s3,
       width = 10, height = 9, dpi = 300, bg = 'white')

### Fig S4: Complex patient treatment timeline -- Oncologist vs LLMs

# Single-patient timeline. No patient ID is hard-coded; supply one through the
# FIGS4_CASE_PID environment variable (the figure is skipped when unset).
# The case-specific configuration below has been removed to protect HPI, so
# this block will not run as published.
CASE_PID <- Sys.getenv('FIGS4_CASE_PID', unset = '')

# -- Regimen display labels --
# Keys are the raw regimen_standardized values from the data; values are the
# display labels drawn on the y-axis. Any regimen not listed displays as-is.
### = ###   # removed to protect HPI
S4_REGIMEN_LABELS <- character(0)

# Sources to overlay -- Oncologist on top, then Gemini, then GPT-5
S4_VERSIONS <- c('v1', 'gemini_2_5', 'gpt5')
S4_LABELS   <- c(v1 = 'Oncologist', gemini_2_5 = 'Gemini 2.5 Pro', gpt5 = 'GPT-5')
S4_COLORS   <- c('Oncologist' = '#000000',
                 'Gemini 2.5 Pro' = '#377EB8',
                 'GPT-5'      = '#E41A1C')
S4_LWIDTH   <- c('Oncologist' = 2.2, 'Gemini 2.5 Pro' = 1.5, 'GPT-5' = 1.5)
# Slight y-dodge: Oncologist on top (negative = up in reversed y-axis), LLMs below
S4_DODGE    <- c('Oncologist' = -0.25, 'Gemini 2.5 Pro' = 0.0, 'GPT-5' = 0.25)

# Helper: build and save a treatment timeline comparison plot for one patient
make_s4_plot <- function(pid, gt_df, gt_cleaned_df, drug_df,
                         versions = S4_VERSIONS, labels = S4_LABELS,
                         colors = S4_COLORS, lwidths = S4_LWIDTH,
                         dodge = S4_DODGE,
                         regimen_labels = S4_REGIMEN_LABELS) {
  cd     <- gt_df %>% filter(gsID == pid)
  anch   <- cd$gt_init_dx
  met_mo <- if (!is.na(cd$gt_met_dx))
               as.numeric(cd$gt_met_dx - anch) / 30.4375 else NA_real_
  lrr_mo <- if (!is.na(cd$gt_lrr_dx))
               as.numeric(cd$gt_lrr_dx - anch) / 30.4375 else NA_real_

  # Death date from cleaned intermediate
  death_date <- gt_cleaned_df %>% filter(as.character(gsID) == pid) %>%
    pull(death_date_clinician) %>% parse_gt_date()
  death_mo <- if (length(death_date) > 0 && !is.na(death_date))
                as.numeric(death_date - anch) / 30.4375 else NA_real_

  df <- drug_df %>%
    filter(gsID == pid, version %in% versions) %>%
    mutate(
      Start      = as.Date(Start),
      Stop       = as.Date(Stop),
      start_mo   = as.numeric(Start - anch) / 30.4375,
      stop_mo    = as.numeric(Stop  - anch) / 30.4375,
      source     = factor(labels[version], levels = names(colors)),
      stop_type  = case_when(
        grepl('progress', tolower(coalesce(as.character(ReasonStop), ''))) ~ 'Progression',
        is.na(Stop) ~ 'Ongoing',
        TRUE        ~ 'Other')
    ) %>%
    filter(!is.na(start_mo)) %>%
    arrange(source, Tx_LineN)

  # -- Align by regimen similarity instead of Tx_LineN --
  # Synonym map: LLM regimen variants -> oncologist regimen for alignment
  ### = ###   # removed to protect HPI
  regimen_synonyms <- character(0)

  # Normalized regimen name for matching purposes
  df <- df %>%
    mutate(regimen_align = coalesce(regimen_synonyms[regimen_standardized],
                                     regimen_standardized))

  # Use oncologist lines as reference y-axis; match LLM lines by regimen name
  onc_ref <- df %>%
    filter(source == 'Oncologist') %>%
    distinct(regimen_align, .keep_all = TRUE) %>%
    arrange(Tx_LineN) %>%
    mutate(y_ref = row_number())

  df <- df %>%
    left_join(onc_ref %>% select(regimen_align, y_ref),
              by = 'regimen_align')

  # LLM-only regimens (not in oncologist) get appended below
  max_ref <- if (nrow(onc_ref) > 0) max(onc_ref$y_ref) else 0L
  extra <- df %>%
    filter(is.na(y_ref)) %>%
    distinct(regimen_align) %>%
    mutate(y_ref_extra = max_ref + row_number())

  if (nrow(extra) > 0) {
    df <- df %>%
      left_join(extra, by = 'regimen_align') %>%
      mutate(y_ref = coalesce(y_ref, y_ref_extra)) %>%
      select(-y_ref_extra)
  }

  df <- df %>%
    mutate(y_pos = y_ref + dodge[as.character(source)]) %>%
    arrange(source, y_ref)

  y_lbl <- df %>%
    group_by(y_ref) %>%
    summarise(
      label = coalesce(
        first(regimen_standardized[source == 'Oncologist']),
        first(regimen_standardized)),
      .groups = 'drop') %>%
    arrange(y_ref) %>%
    mutate(label = coalesce(regimen_labels[label], label))

  ggplot(df) +
    geom_vline(xintercept = 0, linetype = 'dashed',
               color = '#117733', linewidth = 0.55) +
    { if (!is.na(lrr_mo))
        geom_vline(xintercept = lrr_mo, linetype = 'dashed',
                   color = '#CC6600', linewidth = 0.55) } +
    { if (!is.na(met_mo))
        geom_vline(xintercept = met_mo, linetype = 'dashed',
                   color = '#882255', linewidth = 0.55) } +
    { if (!is.na(death_mo))
        geom_vline(xintercept = death_mo, linetype = 'dashed',
                   color = '#CC0000', linewidth = 0.55) } +
    geom_segment(aes(x = start_mo, xend = stop_mo,
                     y = y_pos, yend = y_pos,
                     color = source, linewidth = source),
                 lineend = 'butt') +
    geom_point(aes(x = start_mo, y = y_pos, color = source),
               shape = 16, size = 1.8) +
    geom_point(data = df %>% filter(stop_type == 'Progression'),
               aes(x = stop_mo, y = y_pos, color = source),
               shape = 4, size = 2.2, stroke = 1.1) +
    geom_point(data = df %>% filter(stop_type == 'Other' & !is.na(stop_mo)),
               aes(x = stop_mo, y = y_pos, color = source),
               shape = 1, size = 2.0, stroke = 0.8) +
    scale_color_manual(values = colors, name = NULL) +
    scale_linewidth_manual(values = lwidths, guide = 'none') +
    scale_y_reverse(
      breaks = y_lbl$y_ref,
      labels = y_lbl$label,
      expand = expansion(add = 0.8)) +
    scale_x_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0.01, 0.03))) +
    labs(
      x       = 'Months from Initial Diagnosis',
      y       = NULL,
      title   = 'Treatment Timeline for Example Patient',
      caption = NULL) +
    theme_nature(base_size = 8) +
    theme(
      panel.grid.major.x = element_line(color = 'grey90', linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      axis.text.y        = element_text(size = 6.5, family = 'Helvetica'),
      legend.position    = 'top',
      plot.caption       = element_text(size = 6, color = 'grey40', hjust = 0),
      plot.margin        = margin(8, 15, 5, 5))
}

# FigS4 renders only when a case ID is supplied via FIGS4_CASE_PID.
if (nzchar(CASE_PID)) {
  p_s4 <- make_s4_plot(CASE_PID, gt, gt_cleaned, drug_df)
  ggsave(file.path(FIG_DIR, 'FigS4_treatment_timeline.png'), p_s4,
         width = 9, height = 8, dpi = 300, bg = 'white')
}

### Fig 4 (top panel): All Anti-Cancer Drugs Prescribed vs Oncologist
### Drug bag Jaccard similarity vs Oncologist (v1) -- component of Fig 4

# Normalize only generic AI shorthand (keep specific drug names intact)
norm_drug_fig <- function(x) {
  ai_re <- regex('^(ai|aromatase inhibitor)$', ignore_case = TRUE)
  if_else(str_detect(x, ai_re), 'aromatase_inhibitor', x)
}

# Per-patient drug bags with normalization
bags_s16 <- drug_df %>%
  mutate(gsID = as.character(gsID)) %>%
  separate_rows(regimen_bag_of_drugs, sep = ',') %>%
  mutate(drug = norm_drug_fig(str_trim(regimen_bag_of_drugs))) %>%
  filter(!is.na(drug) & drug != '') %>%
  distinct(gsID, version, drug) %>%
  group_by(gsID, version) %>%
  summarise(drugs = list(sort(unique(drug))), .groups = 'drop')

jac_fn_s16 <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NA_real_)
  i <- length(intersect(a, b)); u <- length(union(a, b))
  if (u == 0) return(NA_real_); i / u
}

wide_s16 <- bags_s16 %>%
  select(gsID, version, drugs) %>%
  pivot_wider(names_from = version, values_from = drugs)

# Versions to show: v2 as Human IRR reference, then all LLMs + CRC
s16_versions <- intersect(c('v2', 'gemini_2_5', 'gpt5', 'gpt4o', 'deepseek_r1', 'crc'),
                           names(wide_s16))

# CRC-censored v1 drug bags for fair CRC comparison in Fig 4
bags_crc_v1_ref <- drug_df %>%
  filter(version == 'v1', is.na(Start) | Start <= CRC_CENSOR_DATE) %>%
  mutate(gsID = as.character(gsID)) %>%
  separate_rows(regimen_bag_of_drugs, sep = ',') %>%
  mutate(drug = norm_drug_fig(str_trim(regimen_bag_of_drugs))) %>%
  filter(!is.na(drug) & drug != '') %>%
  distinct(gsID, drug) %>%
  group_by(gsID) %>%
  summarise(drugs_v1_crc = list(sort(unique(drug))), .groups = 'drop')

s16_long <- map_dfr(s16_versions, function(v) {
  if (!v %in% names(wide_s16)) return(tibble())
  if (v == 'crc') {
    # Use CRC-censored v1 reference for fair comparison
    tibble(gsID = wide_s16$gsID, version = v) %>%
      left_join(bags_crc_v1_ref, by = 'gsID') %>%
      mutate(jaccard = map2_dbl(wide_s16[['crc']][match(gsID, wide_s16$gsID)],
                                drugs_v1_crc, jac_fn_s16)) %>%
      select(gsID, version, jaccard)
  } else {
    tibble(gsID    = wide_s16$gsID,
           version = v,
           jaccard = map2_dbl(wide_s16[[v]], wide_s16[['v1']], jac_fn_s16))
  }
}) %>% filter(!is.na(jaccard))

s16_long <- s16_long %>%
  mutate(label = factor(S16_LABELS[version], levels = S16_ORDER))

# Median + IQR overlay
s16_sum <- s16_long %>%
  group_by(label) %>%
  summarise(m  = median(jaccard,             na.rm = TRUE),
            lo = quantile(jaccard, 0.25,     na.rm = TRUE),
            hi = quantile(jaccard, 0.75,     na.rm = TRUE),
            n  = n(),
            .groups = 'drop')

if (nrow(s16_long) > 0) {
  p_s16 <- plot_jaccard_violin(s16_long, s16_sum,
                                'Similarity in All Anti-Cancer Drugs Prescribed Compared to Oncologist',
                                S16_PAL, S16_ORDER)
}

### Fig 4 (bottom panel): Lines of Therapy Jaccard vs Oncologist
### Same design as Fig 4 top panel (drug bags) but for treatment line sets vs v1

# Panel B: matched-line regimen Jaccard from the line_metrics.R per-patient sheet,
# the same metric as Table 3 and Table S2. A line counts as shared only when the
# chronological match pairs it and the regimen content is identical.
line_metrics_per_patient <- read_xlsx(file.path(TBL_DIR, 'Table3_line_metrics.xlsx'),
                                      sheet = 'per_patient')
s15_versions <- intersect(c('v2', 'gemini_2_5', 'gpt5', 'gpt4o', 'deepseek_r1', 'crc'),
                          unique(line_metrics_per_patient$model))

s15_long <- line_metrics_per_patient %>%
  filter(model %in% s15_versions) %>%
  transmute(gsID = as.character(gsID), version = model, jaccard = jaccard_matched) %>%
  filter(!is.na(jaccard))

if (nrow(s15_long) > 0) {
  s15_long <- s15_long %>%
    mutate(label = factor(S16_LABELS[version], levels = S16_ORDER))

  s15_sum <- s15_long %>%
    group_by(label) %>%
    summarise(m  = median(jaccard,            na.rm = TRUE),
              lo = quantile(jaccard, 0.25,    na.rm = TRUE),
              hi = quantile(jaccard, 0.75,    na.rm = TRUE),
              n  = n(),
              .groups = 'drop')

  p_lines_jac <- plot_jaccard_violin(s15_long, s15_sum,
                                      'Therapy-Line Regimen Agreement vs Oncologist (order-aware matched Jaccard)',
                                      S16_PAL, S16_ORDER)

  # FigS6: all 6 strategies (expanded violin -- mirrors current Fig4)
  ggsave(file.path(FIG_DIR, 'FigS6_drug_lines_jaccard_all.png'),
         (p_s16 / p_lines_jac) + plot_annotation(tag_levels = 'A'),
         width = 9, height = 10, dpi = 300, bg = 'white')

  # Fig4: subset to Gemini 2.5 Pro, GPT-5, Research Coordinator, Second Oncologist
  FIG4_ORDER <- c('Second\nOncologist', 'Gemini 2.5 Pro', 'GPT-5', 'Research\nCoordinator')
  FIG4_PAL   <- S16_PAL[FIG4_ORDER]

  s16_long_f4 <- s16_long %>% filter(label %in% FIG4_ORDER) %>%
    mutate(label = factor(label, levels = FIG4_ORDER))
  s16_sum_f4  <- s16_sum  %>% filter(label %in% FIG4_ORDER) %>%
    mutate(label = factor(label, levels = FIG4_ORDER))
  p_s16_f4 <- plot_jaccard_violin(s16_long_f4, s16_sum_f4,
                                   'Similarity in All Anti-Cancer Drugs Prescribed Compared to Oncologist',
                                   FIG4_PAL, FIG4_ORDER)

  s15_long_f4 <- s15_long %>% filter(label %in% FIG4_ORDER) %>%
    mutate(label = factor(label, levels = FIG4_ORDER))
  s15_sum_f4  <- s15_sum  %>% filter(label %in% FIG4_ORDER) %>%
    mutate(label = factor(label, levels = FIG4_ORDER))
  p_lines_jac_f4 <- plot_jaccard_violin(s15_long_f4, s15_sum_f4,
                                          'Therapy-Line Regimen Agreement vs Oncologist (order-aware matched Jaccard)',
                                          FIG4_PAL, FIG4_ORDER)

  ggsave(file.path(FIG_DIR, 'Fig4_drug_lines_jaccard.png'),
         (p_s16_f4 / p_lines_jac_f4) + plot_annotation(tag_levels = 'A'),
         width = 9, height = 10, dpi = 300, bg = 'white')
}

### FigS5 — Chart Span vs Drug Extraction Accuracy (scatter + per-source trend)

s5_versions <- c('gemini_2_5', 'gpt5', 'gpt4o', 'deepseek_r1', 'crc')
s5_labels   <- c(gemini_2_5 = 'Gemini 2.5 Pro', gpt5 = 'GPT-5',
                 gpt4o = 'GPT-4o', deepseek_r1 = 'DeepSeek-R1',
                 crc = 'Research coordinator')
s5_pal      <- c('Gemini 2.5 Pro' = PAL_LLM[['Gemini 2.5 Pro']],
                 'GPT-5'       = PAL_LLM[['GPT-5']],
                 'GPT-4o'      = PAL_LLM[['GPT-4o']],
                 'DeepSeek-R1' = PAL_LLM[['DeepSeek-R1']],
                 'Research coordinator' = COL_RC)

s5_df <- pat_acc %>%
  filter(version %in% s5_versions,
         !is.na(chart_span_years), !is.na(jaccard_bag_vs_v1)) %>%
  mutate(label = factor(s5_labels[version],
                        levels = c('Gemini 2.5 Pro', 'GPT-5', 'GPT-4o',
                                   'DeepSeek-R1', 'Research coordinator')))

p_s5 <- ggplot(s5_df, aes(x = chart_span_years, y = jaccard_bag_vs_v1,
                           colour = label, fill = label)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = 'lm', se = TRUE, alpha = 0.12, linewidth = 0.8,
              linetype = 'dashed') +
  scale_colour_manual(values = s5_pal) +
  scale_fill_manual(values = s5_pal) +
  labs(x = 'Chart Span (years, first diagnosis to last note)',
       y = 'Patient-level drug Jaccard similarity vs first oncologist',
       title = 'Chart Span vs Drug Extraction Accuracy',
       colour = NULL, fill = NULL) +
  theme_nature(base_size = 8) +
  theme(legend.position = 'bottom')

ggsave(file.path(FIG_DIR, 'FigS5_chart_span_vs_jaccard.png'),
       p_s5, width = 7, height = 5, dpi = 300, bg = 'white')

### Session info

dir.create(LOGS_DIR, showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()),
           file.path(LOGS_DIR, 'session_info_figures.txt'))
