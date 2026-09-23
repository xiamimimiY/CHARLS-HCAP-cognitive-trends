suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
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
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(
  table_dir,
  "b1_frozen_transport_mi_pooled_estimands_2026-06-28.csv"
)
estimands <- read.csv(input_path, check.names = FALSE)

palette <- c(
  "HCAP-calibrated" = "#0F4D92",
  "Social-blind" = "#777777",
  "Uncalibrated" = "#B55D6A"
)
model_labels <- c(
  primary = "HCAP-calibrated",
  social_blind = "Social-blind",
  uncalibrated = "Uncalibrated"
)
model_shapes <- c(
  "HCAP-calibrated" = 16,
  "Social-blind" = 17,
  "Uncalibrated" = 15
)
model_linetypes <- c(
  "HCAP-calibrated" = "solid",
  "Social-blind" = "dotdash",
  "Uncalibrated" = "dashed"
)

theme_b1 <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#202020"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#202020"),
      axis.ticks.length = grid::unit(1.4, "mm"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.3, colour = "#202020"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.4),
      legend.key.width = grid::unit(8, "mm"),
      legend.spacing.x = grid::unit(1.5, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
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
  height_mm = 125,
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

model_factor <- function(x) {
  factor(
    unname(model_labels[x]),
    levels = c("HCAP-calibrated", "Social-blind", "Uncalibrated")
  )
}

# Figure 3: national trends and calibration impact -------------------------

figure3_levels <- estimands %>%
  filter(
    year %in% c(2011, 2013, 2015, 2018),
    dimension == "overall",
    group == "overall",
    model %in% c("primary", "uncalibrated"),
    estimand %in% c(
      "standardized_mean_cognition",
      "standardized_low15_probability"
    )
  ) %>%
  mutate(
    method = model_factor(model),
    outcome = if_else(
      estimand == "standardized_mean_cognition",
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
    source_section = "wave_level"
  )

figure3_changes <- estimands %>%
  filter(
    year == 0,
    dimension == "overall_trend",
    group == "2018_minus_2011",
    model %in% c("primary", "social_blind", "uncalibrated"),
    estimand %in% c(
      "standardized_mean_cognition_change",
      "standardized_low15_probability_change"
    )
  ) %>%
  mutate(
    method = model_factor(model),
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
    source_section = "change_2011_2018"
  )

figure3_source <- bind_rows(figure3_levels, figure3_changes) %>%
  select(
    source_section,
    wave,
    year,
    model,
    method,
    outcome,
    estimand,
    estimate,
    ci_lower,
    ci_upper,
    display_estimate,
    display_ci_lower,
    display_ci_upper
  )
write.csv(
  figure3_source,
  file.path(
    table_dir,
    paste0("b1_figure3_source_data_", date_tag, ".csv")
  ),
  row.names = FALSE
)

p3a_data <- filter(
  figure3_levels,
  outcome == "Mean cognitive performance"
)
p3a <- ggplot(
  p3a_data,
  aes(
    x = year,
    y = display_estimate,
    colour = method,
    shape = method,
    linetype = method,
    group = method
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#B9B9B9") +
  geom_errorbar(
    aes(ymin = display_ci_lower, ymax = display_ci_upper),
    width = 0,
    linewidth = 0.45
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2, stroke = 0.3) +
  scale_x_continuous(breaks = c(2011, 2013, 2015, 2018)) +
  scale_colour_manual(values = palette) +
  scale_shape_manual(values = model_shapes) +
  scale_linetype_manual(values = model_linetypes) +
  labs(
    title = "Standardized mean cognitive performance",
    x = NULL,
    y = "Weighted 2018 SD units"
  ) +
  theme(legend.position = "top")

p3b_data <- filter(
  figure3_levels,
  outcome == "Low cognitive performance"
)
p3b <- ggplot(
  p3b_data,
  aes(
    x = year,
    y = display_estimate,
    colour = method,
    shape = method,
    linetype = method,
    group = method
  )
) +
  geom_errorbar(
    aes(ymin = display_ci_lower, ymax = display_ci_upper),
    width = 0,
    linewidth = 0.45
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2, stroke = 0.3) +
  scale_x_continuous(breaks = c(2011, 2013, 2015, 2018)) +
  scale_y_continuous(labels = label_number(suffix = "%", accuracy = 1)) +
  scale_colour_manual(values = palette) +
  scale_shape_manual(values = model_shapes) +
  scale_linetype_manual(values = model_linetypes) +
  labs(
    title = "Low cognitive performance",
    x = NULL,
    y = "Weighted probability"
  ) +
  theme(legend.position = "top")

change_plot <- function(data, title, x_label, accuracy) {
  ggplot(
    data,
    aes(
      x = display_estimate,
      y = method,
      colour = method,
      shape = method
    )
  ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.4,
      linetype = "dashed",
      colour = "#777777"
    ) +
    geom_errorbar(
      aes(xmin = display_ci_lower, xmax = display_ci_upper),
      orientation = "y",
      width = 0,
      linewidth = 0.55
    ) +
    geom_point(size = 2.5, stroke = 0.3) +
    scale_colour_manual(values = palette, guide = "none") +
    scale_shape_manual(values = model_shapes, guide = "none") +
    scale_y_discrete(
      limits = c("Uncalibrated", "Social-blind", "HCAP-calibrated")
    ) +
    scale_x_continuous(labels = label_number(accuracy = accuracy)) +
    labs(title = title, x = x_label, y = NULL) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
}

p3c <- change_plot(
  filter(figure3_changes, outcome == "Mean cognitive performance"),
  "Change, 2011-2018",
  "Mean cognitive performance change (SD)",
  0.01
)
p3d <- change_plot(
  filter(figure3_changes, outcome == "Low cognitive performance"),
  "Change, 2011-2018",
  "Low-performance change (percentage points)",
  1
)

figure3 <- (
  (p3a | p3b) /
    (p3c | p3d)
) +
  plot_layout(guides = "collect", heights = c(1.25, 1)) +
  plot_annotation(
    tag_levels = "a",
    caption = paste(
      "Points are age- and sex-standardized estimates;",
      "error bars are Rubin-pooled 95% confidence intervals.",
      "Low performance uses the fixed weighted 2018 lower-15% threshold."
    )
  ) &
  theme(
    legend.position = "top",
    plot.caption = element_text(
      size = 5.8,
      colour = "#4A4A4A",
      hjust = 0
    )
  )

figure3_base <- file.path(
  figure_dir,
  paste0("b1_figure3_calibrated_trends_", date_tag)
)
save_pub_r(figure3, figure3_base, width_mm = 183, height_mm = 127)

# Figure 4: inequality change and calibration impact -----------------------

gap_spec <- data.frame(
  dimension = c("education", "residence", "hukou", "consumption_rank"),
  trend_dimension = c(
    "education_gap_trend",
    "residence_gap_trend",
    "hukou_gap_trend",
    "consumption_rank_trend"
  ),
  group = c(
    "no_formal_minus_middle_school_plus",
    "rural_minus_urban",
    "agricultural_minus_non_agricultural",
    "poorest_minus_richest"
  ),
  display_group = c(
    "No formal schooling vs middle school+",
    "Rural vs urban",
    "Agricultural vs non-agricultural hukou",
    "Consumption-rank SII"
  ),
  stringsAsFactors = FALSE
)

mean_gap_estimands <- c(
  "standardized_mean_cognition_absolute_gap",
  "age_sex_adjusted_mean_cognition_sii"
)
low_gap_estimands <- c(
  "standardized_low15_probability_absolute_gap",
  "age_sex_adjusted_low15_probability_sii"
)
mean_change_estimands <- c(
  "standardized_mean_cognition_absolute_gap_change",
  "age_sex_adjusted_mean_cognition_sii_change"
)
low_change_estimands <- c(
  "standardized_low15_probability_absolute_gap_change",
  "age_sex_adjusted_low15_probability_sii_change"
)

figure4_changes <- estimands %>%
  inner_join(
    select(gap_spec, trend_dimension, display_group),
    by = c("dimension" = "trend_dimension")
  ) %>%
  filter(
    year == 0,
    group == "2018_minus_2011",
    model == "primary",
    estimand %in% c(mean_change_estimands, low_change_estimands)
  ) %>%
  mutate(
    outcome = if_else(
      estimand %in% mean_change_estimands,
      "Mean cognitive performance",
      "Low cognitive performance"
    ),
    display_estimate = if_else(
      outcome == "Mean cognitive performance",
      -estimate,
      100 * estimate
    ),
    display_ci_lower = if_else(
      outcome == "Mean cognitive performance",
      -ci_upper,
      100 * ci_lower
    ),
    display_ci_upper = if_else(
      outcome == "Mean cognitive performance",
      -ci_lower,
      100 * ci_upper
    ),
    method = "HCAP-calibrated",
    source_section = "gap_change_2011_2018"
  )

figure4_2018 <- estimands %>%
  inner_join(
    select(gap_spec, dimension, group, display_group),
    by = c("dimension", "group")
  ) %>%
  filter(
    year == 2018,
    model %in% c("primary", "social_blind", "uncalibrated"),
    estimand %in% c(mean_gap_estimands, low_gap_estimands)
  ) %>%
  mutate(
    outcome = if_else(
      estimand %in% mean_gap_estimands,
      "Mean cognitive performance",
      "Low cognitive performance"
    ),
    display_estimate = if_else(
      outcome == "Mean cognitive performance",
      -estimate,
      100 * estimate
    ),
    display_ci_lower = if_else(
      outcome == "Mean cognitive performance",
      -ci_upper,
      100 * ci_lower
    ),
    display_ci_upper = if_else(
      outcome == "Mean cognitive performance",
      -ci_lower,
      100 * ci_upper
    ),
    method = model_factor(model),
    source_section = "gap_level_2018"
  )

figure4_source <- bind_rows(figure4_changes, figure4_2018) %>%
  select(
    source_section,
    wave,
    year,
    model,
    method,
    dimension,
    group,
    display_group,
    outcome,
    estimand,
    estimate,
    ci_lower,
    ci_upper,
    display_estimate,
    display_ci_lower,
    display_ci_upper
  )
write.csv(
  figure4_source,
  file.path(
    table_dir,
    paste0("b1_figure4_source_data_", date_tag, ".csv")
  ),
  row.names = FALSE
)

gap_levels <- rev(gap_spec$display_group)
figure4_changes$display_group <- factor(
  figure4_changes$display_group,
  levels = gap_levels
)
figure4_2018$display_group <- factor(
  figure4_2018$display_group,
  levels = gap_levels
)

single_gap_plot <- function(data, title, subtitle, x_label, accuracy) {
  ggplot(
    data,
    aes(x = display_estimate, y = display_group)
  ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.4,
      linetype = "dashed",
      colour = "#777777"
    ) +
    geom_errorbar(
      aes(xmin = display_ci_lower, xmax = display_ci_upper),
      orientation = "y",
      width = 0,
      linewidth = 0.55,
      colour = palette[["HCAP-calibrated"]]
    ) +
    geom_point(
      size = 2.5,
      shape = 16,
      colour = palette[["HCAP-calibrated"]]
    ) +
    scale_x_continuous(labels = label_number(accuracy = accuracy)) +
    labs(title = title, subtitle = subtitle, x = x_label, y = NULL) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
}

model_gap_plot <- function(data, title, subtitle, x_label, accuracy) {
  ggplot(
    data,
    aes(
      x = display_estimate,
      y = display_group,
      colour = method,
      shape = method
    )
  ) +
    geom_errorbar(
      aes(xmin = display_ci_lower, xmax = display_ci_upper),
      orientation = "y",
      width = 0,
      linewidth = 0.45,
      position = position_dodge(width = 0.58)
    ) +
    geom_point(
      size = 2.2,
      stroke = 0.3,
      position = position_dodge(width = 0.58)
    ) +
    scale_colour_manual(values = palette) +
    scale_shape_manual(values = model_shapes) +
    scale_x_continuous(labels = label_number(accuracy = accuracy)) +
    labs(title = title, subtitle = subtitle, x = x_label, y = NULL) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "top"
    )
}

p4a <- single_gap_plot(
  filter(figure4_changes, outcome == "Mean cognitive performance"),
  "Inequality change in mean cognitive performance",
  "HCAP-calibrated metric only",
  "Gap change, 2011-2018 (SD)",
  0.05
)
p4b <- single_gap_plot(
  filter(figure4_changes, outcome == "Low cognitive performance"),
  "Low-performance inequality change",
  "HCAP-calibrated metric only",
  "Gap change, 2011-2018 (percentage points)",
  1
)
p4c <- model_gap_plot(
  filter(figure4_2018, outcome == "Mean cognitive performance"),
  "Gaps in mean cognitive performance, 2018",
  "Comparison across three metrics",
  "Adverse gap (SD)",
  0.1
)
p4d <- model_gap_plot(
  filter(figure4_2018, outcome == "Low cognitive performance"),
  "Low-performance gaps in 2018",
  "Comparison across three metrics",
  "Adverse gap (percentage points)",
  10
)

figure4_top <- p4a | p4b
figure4_bottom <- (p4c | p4d) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

figure4 <- (
  figure4_top /
    figure4_bottom
) +
  plot_layout(heights = c(1, 1.18)) +
  plot_annotation(
    tag_levels = "a",
    caption = paste(
      "Panels a-b show HCAP-calibrated changes; panels c-d compare three metrics in 2018.\n",
      "Positive values indicate worse cognitive performance or more low performance",
      "in the disadvantaged group.\nConsumption estimates are age- and sex-adjusted",
      "slope indices of inequality.\nError bars are Rubin-pooled 95% confidence intervals."
    )
  ) &
  theme(
    plot.caption = element_text(
      size = 5.8,
      colour = "#4A4A4A",
      hjust = 0
    )
  )

figure4_base <- file.path(
  figure_dir,
  paste0("b1_figure4_social_inequalities_", date_tag)
)
save_pub_r(figure4, figure4_base, width_mm = 183, height_mm = 145)

cat("wrote Figure 3 and Figure 4 publication bundles\n")
