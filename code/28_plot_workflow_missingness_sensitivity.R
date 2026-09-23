suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

options(mc.cores = 1L)

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1) stop("Run this file with Rscript.")
script_path <- normalizePath(
  sub("^--file=", "", gsub("~\\+~", " ", script_argument)),
  mustWork = TRUE
)
study <- dirname(dirname(script_path))

date_tag <- "2026-07-06"
table_dir <- file.path(study, "results", "tables")
figure_dir <- file.path(study, "results", "figures")
derived_dir <- file.path(study, "data", "derived")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

palette <- c(
  blue = "#0F4D92",
  blue_mid = "#5B82B5",
  blue_light = "#D9E6F2",
  teal = "#3E8E8C",
  rose = "#B55D6A",
  rose_light = "#F2DCE0",
  gray_dark = "#4A4A4A",
  gray_mid = "#7A7A7A",
  gray_light = "#E7E7E7",
  gold = "#D69B2D"
)

theme_b1 <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#202020"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#202020"),
      axis.ticks.length = grid::unit(1.4, "mm"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.3, colour = "#202020"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.1, face = "bold"),
      plot.title = element_text(size = base_size + 0.2, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#4A4A4A"),
      plot.tag = element_text(size = 8.2, face = "bold"),
      plot.margin = margin(5, 7, 5, 5),
      panel.grid = element_blank()
    )
}
theme_set(theme_b1())

save_pub_r <- function(
  plot,
  filename,
  width_mm = 183,
  height_mm = 120,
  dpi = 600
) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(
    paste0(filename, ".svg"),
    width = width_in,
    height = height_in
  )
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(filename, ".pdf"),
    width = width_in,
    height = height_in,
    family = "Arial"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(filename, ".tiff"),
    width = width_in,
    height = height_in,
    units = "in",
    res = dpi,
    compression = "lzw"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(filename, ".png"),
    width = width_in,
    height = height_in,
    units = "in",
    res = 300
  )
  print(plot)
  grDevices::dev.off()
}

# Figure 1: study design and calibration workflow --------------------------

figure1_source <- data.frame(
  stage = c(
    rep("national_wave", 4),
    "linked_hcap_frame",
    "hcap_any_indicator",
    "hcap_all_missing",
    "hcap_missing_weight",
    "transport_calibration"
  ),
  year = c(2011, 2013, 2015, 2018, 2018, 2018, 2018, 2018, 2018),
  n = c(7290, 8524, 9845, 10800, 10765, 9755, 1010, 508, 9247),
  definition = c(
    rep("CHARLS participants aged >=60 years", 4),
    "Linked 2018 HCAP analysis frame",
    "At least one of 12 HCAP indicators observed",
    "All 12 HCAP indicators missing",
    "Missing 2018 cross-sectional respondent weight",
    "Calibration sample per routine-cognition imputation"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  figure1_source,
  file.path(table_dir, paste0("b1_figure1_source_data_", date_tag, ".csv")),
  row.names = FALSE
)

box_data <- data.frame(
  xmin = c(0.2, 2.95, 5.25, 7.65, 10.0),
  xmax = c(2.7, 5.0, 7.4, 9.75, 12.5),
  ymin = c(1.25, 1.25, 1.25, 1.25, 1.25),
  ymax = c(5.65, 5.65, 5.65, 5.65, 5.65),
  fill = c(
    palette[["gray_light"]],
    palette[["blue_light"]],
    "#E8E2F0",
    "#DCEBE8",
    palette[["rose_light"]]
  ),
  heading = c(
    "1  National\nwaves",
    "2  2018 HCAP\nanchor",
    "3  Measurement\nmodel",
    "4  Frozen\ntransport",
    "5  Study\noutputs"
  ),
  body = c(
    "Age >=60 years\n2011  n=7,290\n2013  n=8,524\n2015  n=9,845\n2018  n=10,800",
    "Linked frame  n=10,765\n>=1 of 12 indicators\nn=9,755\nAll 12 missing  n=1,010\nWeighted 2018 metric",
    "Education\npartial-invariance\nbifactor model\n\nGeneral cognition anchor\nFIML primary\n50 HCAP item-level\nimputations",
    "Five routine components\nNonlinear age + sex\nThree education groups\nEducation-specific slopes\nMissing 2018 weight n=508\nCalibration n=9,247 per MI",
    "Continuous cognitive\nperformance\nFixed lower-15%\nlow-performance probability\n2011-2018 trends\nSocial inequalities\nResearch-use scoring tool"
  ),
  stringsAsFactors = FALSE
)

figure1 <- ggplot() +
  geom_rect(
    data = box_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    colour = "#A0A0A0",
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  geom_text(
    data = box_data,
    aes(x = (xmin + xmax) / 2, y = 5.12, label = heading),
    family = "Arial",
    fontface = "bold",
    size = 2.55,
    lineheight = 0.95,
    colour = "#202020"
  ) +
  geom_text(
    data = box_data,
    aes(x = (xmin + xmax) / 2, y = 3.25, label = body),
    family = "Arial",
    size = 2.2,
    lineheight = 1.18,
    colour = "#303030"
  ) +
  geom_segment(
    data = data.frame(
      x = c(2.72, 5.02, 7.42, 9.77),
      xend = c(2.91, 5.21, 7.61, 9.96),
      y = 3.45,
      yend = 3.45
    ),
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = grid::arrow(length = grid::unit(2.1, "mm"), type = "closed"),
    linewidth = 0.6,
    colour = palette[["blue"]]
  ) +
  annotate(
    "text",
    x = 6.35,
    y = 0.72,
    label = paste(
      "The weighted 2018 HCAP metric and thresholds were fixed,",
      "then transported backward without wave-specific restandardization."
    ),
    family = "Arial",
    fontface = "italic",
    size = 2.5,
    colour = palette[["gray_dark"]]
  ) +
  coord_cartesian(xlim = c(0, 12.7), ylim = c(0.35, 6), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(5, 5, 5, 5))

save_pub_r(
  figure1,
  file.path(
    figure_dir,
    paste0("b1_figure1_study_design_calibration_", date_tag)
  ),
  width_mm = 183,
  height_mm = 82
)

# Figure 2: corrected HCAP item missingness -------------------------------

hcap_path <- file.path(
  derived_dir,
  "b1_hcap_fiml_input_2026-06-19.csv"
)
hcap <- read.csv(hcap_path, check.names = FALSE)
hcap_items <- c(
  "z_test_short_learning_trial1",
  "z_test_short_delayed",
  "z_test_ten_learning_trial1",
  "z_test_ten_learning_trial2",
  "z_test_ten_learning_trial3",
  "z_test_ten_delayed",
  "z_test_recognition",
  "z_test_serial7",
  "z_test_animal_fluency",
  "z_test_language",
  "z_test_visuospatial",
  "z_test_orientation"
)
item_labels <- c(
  "3-word learning",
  "3-word delayed",
  "10-word learning T1",
  "10-word learning T2",
  "10-word learning T3",
  "10-word delayed",
  "Word recognition",
  "Serial 7",
  "Animal fluency",
  "Language",
  "Visuospatial",
  "Orientation"
)
names(item_labels) <- hcap_items

overall_missing <- do.call(
  rbind,
  lapply(hcap_items, function(item) {
    data.frame(
      source_section = "overall",
      subgroup_domain = "Overall",
      subgroup = "Overall",
      item = item,
      item_label = item_labels[[item]],
      n = nrow(hcap),
      missing_n = sum(is.na(hcap[[item]])),
      missing_pct = 100 * mean(is.na(hcap[[item]])),
      stringsAsFactors = FALSE
    )
  })
)

subgroup_specs <- list(
  Education = list(
    variable = "edu3",
    values = c(1, 2, 3),
    labels = c(
      "No formal schooling",
      "Some schooling below middle",
      "Middle school+"
    )
  ),
  Residence = list(
    variable = "rural",
    values = c(0, 1),
    labels = c("Urban", "Rural")
  ),
  Vision = list(
    variable = "vision_problem",
    values = c(0, 1),
    labels = c("No vision problem", "Vision problem")
  ),
  Hearing = list(
    variable = "hearing_problem",
    values = c(0, 1),
    labels = c("No hearing problem", "Hearing problem")
  )
)

subgroup_missing_rows <- list()
row_index <- 1L
for (domain in names(subgroup_specs)) {
  specification <- subgroup_specs[[domain]]
  for (value_index in seq_along(specification$values)) {
    value <- specification$values[value_index]
    subgroup_label <- specification$labels[value_index]
    target <- !is.na(hcap[[specification$variable]]) &
      hcap[[specification$variable]] == value
    for (item in hcap_items) {
      subgroup_missing_rows[[row_index]] <- data.frame(
        source_section = "subgroup",
        subgroup_domain = domain,
        subgroup = subgroup_label,
        item = item,
        item_label = item_labels[[item]],
        n = sum(target),
        missing_n = sum(is.na(hcap[[item]][target])),
        missing_pct = 100 * mean(is.na(hcap[[item]][target])),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
}
subgroup_missing <- bind_rows(subgroup_missing_rows) %>%
  left_join(
    select(overall_missing, item, overall_missing_pct = missing_pct),
    by = "item"
  ) %>%
  mutate(deviation_pp = missing_pct - overall_missing_pct)

figure2_source <- bind_rows(
  mutate(overall_missing, overall_missing_pct = missing_pct, deviation_pp = 0),
  subgroup_missing
)
write.csv(
  figure2_source,
  file.path(table_dir, paste0("b1_figure2_source_data_", date_tag, ".csv")),
  row.names = FALSE
)

item_order <- overall_missing %>%
  arrange(missing_pct) %>%
  pull(item_label)
overall_missing$item_label <- factor(
  overall_missing$item_label,
  levels = item_order
)

p2a <- ggplot(
  overall_missing,
  aes(x = missing_pct, y = item_label)
) +
  geom_col(width = 0.68, fill = palette[["blue_mid"]]) +
  geom_text(
    aes(label = sprintf("%.1f%%", missing_pct)),
    hjust = -0.12,
    size = 2.15,
    family = "Arial",
    colour = "#303030"
  ) +
  scale_x_continuous(
    limits = c(0, max(overall_missing$missing_pct) * 1.35),
    labels = label_number(suffix = "%", accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Overall item missingness",
    x = "Participants with missing item",
    y = NULL
  ) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

subgroup_order <- c(
  "No formal schooling",
  "Some schooling below middle",
  "Middle school+",
  "Urban",
  "Rural",
  "No vision problem",
  "Vision problem",
  "No hearing problem",
  "Hearing problem"
)
subgroup_missing$subgroup <- factor(
  subgroup_missing$subgroup,
  levels = rev(subgroup_order)
)
subgroup_missing$item_label <- factor(
  subgroup_missing$item_label,
  levels = item_labels[hcap_items]
)
fill_limit <- ceiling(max(abs(subgroup_missing$deviation_pp), na.rm = TRUE))
subgroup_missing$text_colour <- if_else(
  abs(subgroup_missing$deviation_pp) >= 12,
  "white",
  "#252525"
)

p2b <- ggplot(
  subgroup_missing,
  aes(x = item_label, y = subgroup, fill = deviation_pp)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(
    aes(
      label = sprintf("%+.1f", deviation_pp),
      colour = text_colour
    ),
    size = 1.75,
    family = "Arial"
  ) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#5B82B5",
    mid = "white",
    high = "#C76B72",
    midpoint = 0,
    limits = c(-fill_limit, fill_limit),
    name = "Deviation\n(pp)"
  ) +
  labs(
    title = "Subgroup deviation from overall item missingness",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 48,
      hjust = 1,
      vjust = 1,
      size = 5.6
    ),
    axis.text.y = element_text(size = 6.1),
    legend.position = "right",
    legend.key.height = grid::unit(12, "mm")
  )

figure2 <- (p2a | p2b) +
  plot_layout(widths = c(0.78, 2.1)) +
  plot_annotation(
    tag_levels = "a",
    caption = paste(
      "Positive heatmap values indicate more missingness than the item's",
      "overall rate. Education uses the corrected three-group definition.",
      "Analysis frame n=10,765."
    )
  ) &
  theme(
    plot.caption = element_text(
      size = 5.8,
      colour = "#4A4A4A",
      hjust = 0
    )
  )

save_pub_r(
  figure2,
  file.path(
    figure_dir,
    paste0("b1_figure2_hcap_item_missingness_", date_tag)
  ),
  width_mm = 183,
  height_mm = 126
)

# Figure 5: robustness and scoring-model validation ------------------------

formal <- read.csv(
  file.path(
    table_dir,
    "b1_frozen_transport_mi_pooled_estimands_2026-06-28.csv"
  ),
  check.names = FALSE
)
weight_mnar <- read.csv(
  file.path(table_dir, "b1_weight_mnar_sensitivity_2026-06-28.csv"),
  check.names = FALSE
)
item_mi <- read.csv(
  file.path(
    table_dir,
    "b1_hcap_item_mi_vs_fiml_longitudinal_2026-06-28.csv"
  ),
  check.names = FALSE
)
bootstrap <- read.csv(
  file.path(
    table_dir,
    "b1_community_cluster_bootstrap_summary_2026-06-28.csv"
  ),
  check.names = FALSE
)
all_missing <- read.csv(
  file.path(
    table_dir,
    "b1_all_missing_longitudinal_bounds_2026-06-28.csv"
  ),
  check.names = FALSE
)
scoring <- read.csv(
  file.path(
    table_dir,
    "b1_scoring_tool_validation_summary_2026-06-28.csv"
  ),
  check.names = FALSE
)

scenario_labels <- c(
  primary = "Primary MI + survey",
  unweighted = "Unweighted",
  weight_trim_99 = "Weights trimmed at 99th percentile",
  mnar_delta_0.25 = "Routine cognition MNAR -0.25 SD",
  mnar_delta_0.50 = "Routine cognition MNAR -0.50 SD",
  item_mi = "HCAP item-level MI",
  bootstrap = "Community bootstrap",
  complete_predictor = "Complete predictors only"
)
scenario_order <- unname(scenario_labels[c(
  "primary",
  "unweighted",
  "weight_trim_99",
  "mnar_delta_0.25",
  "mnar_delta_0.50",
  "item_mi",
  "bootstrap",
  "complete_predictor"
)])

primary_sensitivity <- formal %>%
  filter(
    model == "primary",
    dimension == "overall_trend",
    group == "2018_minus_2011",
    estimand %in% c(
      "standardized_mean_cognition_change",
      "standardized_low15_probability_change"
    )
  ) %>%
  transmute(
    scenario = "primary",
    estimand,
    estimate,
    ci_lower,
    ci_upper,
    interval_type = "Rubin-pooled 95% CI"
  )

weight_sensitivity <- weight_mnar %>%
  filter(
    model == "primary",
    dimension == "overall_trend",
    group == "2018_minus_2011",
    estimand %in% c(
      "standardized_mean_cognition_change",
      "standardized_low15_probability_change"
    )
  ) %>%
  transmute(
    scenario,
    estimand,
    estimate,
    ci_lower,
    ci_upper,
    interval_type = if_else(
      scenario == "complete_predictor",
      "Model-based 95% CI; selection diagnostic",
      "Rubin-pooled 95% CI"
    )
  )

item_sensitivity <- item_mi %>%
  filter(
    dimension == "overall_trend",
    group == "2018_minus_2011",
    estimand %in% c(
      "standardized_mean_cognition_change",
      "standardized_low15_probability_change"
    )
  ) %>%
  transmute(
    scenario = "item_mi",
    estimand,
    estimate = estimate_item_mi,
    ci_lower = ci_lower_item_mi,
    ci_upper = ci_upper_item_mi,
    interval_type = "Rubin-pooled 95% CI"
  )

bootstrap_sensitivity <- bootstrap %>%
  filter(variable %in% c(
    "mean_cognition_change",
    "low15_probability_change"
  )) %>%
  transmute(
    scenario = "bootstrap",
    estimand = if_else(
      variable == "mean_cognition_change",
      "standardized_mean_cognition_change",
      "standardized_low15_probability_change"
    ),
    estimate = bootstrap_mean,
    ci_lower = percentile_ci_lower,
    ci_upper = percentile_ci_upper,
    interval_type = "Bootstrap percentile 95% interval"
  )

sensitivity <- bind_rows(
  primary_sensitivity,
  weight_sensitivity,
  item_sensitivity,
  bootstrap_sensitivity
) %>%
  mutate(
    scenario_label = unname(scenario_labels[scenario]),
    outcome = if_else(
      estimand == "standardized_mean_cognition_change",
      "Mean cognitive performance",
      "Low cognitive performance"
    ),
    display_estimate = if_else(
      outcome == "Low cognitive performance",
      100 * estimate,
      estimate
    ),
    display_ci_lower = if_else(
      outcome == "Low cognitive performance",
      100 * ci_lower,
      ci_lower
    ),
    display_ci_upper = if_else(
      outcome == "Low cognitive performance",
      100 * ci_upper,
      ci_upper
    ),
    evidence_family = case_when(
      scenario == "primary" ~ "Primary",
      scenario == "complete_predictor" ~ "Selection diagnostic",
      TRUE ~ "Robustness"
    ),
    source_section = "trend_sensitivity"
  )

bound_labels <- c(
  primary = "Primary: observed HCAP indicators",
  MAR = "Include all-missing group: MAR",
  delta_0.25_SD_lower = "All-missing group: -0.25 SD",
  delta_0.50_SD_lower = "All-missing group: -0.50 SD"
)
bound_order <- unname(bound_labels[c(
  "primary",
  "MAR",
  "delta_0.25_SD_lower",
  "delta_0.50_SD_lower"
)])

mean_bounds <- all_missing %>%
  filter(outcome == "mean_cognition") %>%
  transmute(
    scenario,
    outcome = "Mean cognitive performance",
    estimate = bounded_change_2018_minus_2011
  )
low_bounds <- all_missing %>%
  filter(
    outcome == "low_cognitive_performance",
    lower_tail_pct == 15
  ) %>%
  transmute(
    scenario,
    outcome = "Low cognitive performance",
    estimate = 100 * bounded_change_2018_minus_2011
  )
primary_bounds <- data.frame(
  scenario = rep("primary", 2),
  outcome = c("Mean cognitive performance", "Low cognitive performance"),
  estimate = c(
    unique(
      all_missing$primary_change_2018_minus_2011[
        all_missing$outcome == "mean_cognition"
      ]
    )[1],
    100 * unique(
      all_missing$primary_change_2018_minus_2011[
        all_missing$outcome == "low_cognitive_performance" &
          all_missing$lower_tail_pct == 15
      ]
    )[1]
  )
)
selection_bounds <- bind_rows(
  primary_bounds,
  mean_bounds,
  low_bounds
) %>%
  mutate(
    scenario_label = unname(bound_labels[scenario]),
    source_section = "all_hcap_missing_bounds"
  )

scoring_labels <- c(
  overall = "Overall",
  no_formal_schooling = "No formal schooling",
  some_schooling_below_middle = "Some schooling below middle",
  middle_school_plus = "Middle school+"
)
scoring_validation <- scoring %>%
  transmute(
    group,
    group_label = unname(scoring_labels[group]),
    n = n_mean,
    rmse = weighted_rmse_mean,
    r_squared = weighted_r_squared_mean,
    calibration_intercept = calibration_intercept_mean,
    calibration_slope = calibration_slope_mean,
    source_section = "community_separated_validation"
  )

figure5_source <- bind_rows(
  sensitivity %>%
    transmute(
      source_section,
      group = scenario,
      label = scenario_label,
      outcome,
      estimate,
      ci_lower,
      ci_upper,
      display_estimate,
      display_ci_lower,
      display_ci_upper,
      metric_2 = NA_real_,
      n = NA_real_,
      interval_type
    ),
  selection_bounds %>%
    transmute(
      source_section,
      group = scenario,
      label = scenario_label,
      outcome,
      estimate,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      display_estimate = estimate,
      display_ci_lower = NA_real_,
      display_ci_upper = NA_real_,
      metric_2 = NA_real_,
      n = NA_real_,
      interval_type = "Selection point bound; no sampling interval"
    ),
  scoring_validation %>%
    transmute(
      source_section,
      group,
      label = group_label,
      outcome = "Calibration slope vs R-squared",
      estimate = calibration_slope,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      display_estimate = calibration_slope,
      display_ci_lower = NA_real_,
      display_ci_upper = NA_real_,
      metric_2 = r_squared,
      n,
      interval_type = "Mean across 50 imputations"
    )
)
write.csv(
  figure5_source,
  file.path(table_dir, paste0("b1_figure5_source_data_", date_tag, ".csv")),
  row.names = FALSE
)

sensitivity$scenario_label <- factor(
  sensitivity$scenario_label,
  levels = rev(scenario_order)
)
evidence_colors <- c(
  Primary = palette[["blue"]],
  Robustness = palette[["teal"]],
  `Selection diagnostic` = palette[["rose"]]
)
evidence_shapes <- c(
  Primary = 16,
  Robustness = 17,
  `Selection diagnostic` = 15
)

sensitivity_plot <- function(data, title, x_label, accuracy) {
  ggplot(
    data,
    aes(
      x = display_estimate,
      y = scenario_label,
      colour = evidence_family,
      shape = evidence_family
    )
  ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.4,
      linetype = "dashed",
      colour = palette[["gray_mid"]]
    ) +
    geom_errorbar(
      aes(xmin = display_ci_lower, xmax = display_ci_upper),
      orientation = "y",
      width = 0,
      linewidth = 0.5
    ) +
    geom_point(size = 2.35, stroke = 0.3) +
    scale_colour_manual(values = evidence_colors) +
    scale_shape_manual(values = evidence_shapes) +
    scale_x_continuous(labels = label_number(accuracy = accuracy)) +
    labs(title = title, x = x_label, y = NULL) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    )
}

p5a <- sensitivity_plot(
  filter(sensitivity, outcome == "Mean cognitive performance"),
  "Mean cognitive performance trend",
  "Change in mean cognitive performance, 2011-2018 (SD)",
  0.05
)
p5b <- sensitivity_plot(
  filter(sensitivity, outcome == "Low cognitive performance"),
  "Low-performance trend",
  "Change, 2011-2018 (percentage points)",
  1
)

selection_bounds$scenario_label <- factor(
  selection_bounds$scenario_label,
  levels = rev(bound_order)
)
p5c <- ggplot(
  selection_bounds,
  aes(x = estimate, y = scenario_label)
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    colour = palette[["gray_mid"]]
  ) +
  geom_point(
    aes(fill = outcome),
    shape = 21,
    size = 2.4,
    stroke = 0.35,
    colour = "#303030"
  ) +
  facet_wrap(
    ~outcome,
    scales = "free_x",
    labeller = as_labeller(c(
      `Mean cognitive performance` =
        "Mean cognitive\nperformance (SD)\n\u2190 More adverse",
      `Low cognitive performance` =
        "Low performance (pp)\nMore adverse \u2192"
    ))
  ) +
  scale_fill_manual(
    values = c(
      "Mean cognitive performance" = palette[["blue"]],
      "Low cognitive performance" = palette[["rose"]]
    ),
    guide = "none"
  ) +
  labs(
    title = "All-HCAP-indicator-missing selection bounds",
    subtitle = "Point bounds; no sampling confidence intervals",
    x = "Change, 2011-2018",
    y = NULL
  ) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 5.8, lineheight = 0.95)
  )

scoring_validation$group_label <- factor(
  scoring_validation$group_label,
  levels = c(
    "Overall",
    "No formal schooling",
    "Some schooling below middle",
    "Middle school+"
  )
)
validation_colors <- c(
  "Overall" = palette[["gray_dark"]],
  "No formal schooling" = "#4676A9",
  "Some schooling below middle" = "#7298BF",
  "Middle school+" = "#9AB7D2"
)
p5d <- ggplot(
  scoring_validation,
  aes(
    x = r_squared,
    y = calibration_slope,
    colour = group_label,
    size = n
  )
) +
  geom_hline(
    yintercept = 1,
    linewidth = 0.4,
    linetype = "dashed",
    colour = palette[["gray_mid"]]
  ) +
  geom_point(alpha = 0.95) +
  geom_text(
    aes(label = group_label),
    size = 2.05,
    family = "Arial",
    colour = "#252525",
    show.legend = FALSE,
    nudge_y = c(0.006, 0.008, -0.008, -0.009)
  ) +
  scale_colour_manual(values = validation_colors, guide = "none") +
  scale_size_continuous(range = c(2.3, 4.4), guide = "none") +
  scale_x_continuous(
    limits = c(0.5, 0.84),
    labels = label_number(accuracy = 0.1)
  ) +
  scale_y_continuous(
    limits = c(0.95, 1.03),
    labels = label_number(accuracy = 0.01)
  ) +
  labs(
    title = "Community-separated internal validation",
    subtitle = "Point size reflects group n; overall RMSE = 0.455 SD",
    x = expression(R^2),
    y = "Calibration slope"
  )

figure5 <- (
  (p5a | p5b) /
    (p5c | p5d)
) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(
    tag_levels = "a",
    caption = paste(
      "Panels a-b show model-based 95% confidence intervals except the",
      "community bootstrap, which uses percentile intervals.\n",
      "Complete predictors are a selection diagnostic, not an alternative primary analysis."
    )
  ) &
  theme(
    plot.caption = element_text(
      size = 5.8,
      colour = "#4A4A4A",
      hjust = 0
    )
  )

save_pub_r(
  figure5,
  file.path(
    figure_dir,
    paste0("b1_figure5_robustness_validation_", date_tag)
  ),
  width_mm = 183,
  height_mm = 142
)

cat("wrote Figure 1, Figure 2, and Figure 5 publication bundles\n")
