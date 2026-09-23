suppressPackageStartupMessages({
  library(dplyr)
})

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1) stop("Run this file with Rscript.")
script_path <- normalizePath(
  sub("^--file=", "", gsub("~\\+~", " ", script_argument)),
  mustWork = TRUE
)
study <- dirname(dirname(script_path))

date_tag <- "2026-06-29"
table_dir <- file.path(study, "results", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

weighted_mean <- function(x, weight) {
  keep <- is.finite(x) & is.finite(weight) & weight > 0
  if (!any(keep)) {
    return(NA_real_)
  }
  sum(x[keep] * weight[keep]) / sum(weight[keep])
}

weighted_sd <- function(x, weight) {
  keep <- is.finite(x) & is.finite(weight) & weight > 0
  if (!any(keep)) {
    return(NA_real_)
  }
  center <- weighted_mean(x[keep], weight[keep])
  sqrt(sum(weight[keep] * (x[keep] - center)^2) / sum(weight[keep]))
}

kish_effective_n <- function(weight) {
  keep <- is.finite(weight) & weight > 0
  weight <- weight[keep]
  if (length(weight) == 0) {
    return(NA_real_)
  }
  sum(weight)^2 / sum(weight^2)
}

long <- read.csv(
  file.path(
    study,
    "data",
    "derived",
    "b1_frozen_transport_mi_input_2026-06-28.csv"
  ),
  check.names = FALSE,
  colClasses = c(ID = "character")
)

wave_years <- c("1" = 2011, "2" = 2013, "3" = 2015, "4" = 2018)
routine_items <- c("imrc", "dlrc", "ser7", "orient", "draw")

characteristic_specs <- list(
  list(
    key = "unweighted_n",
    label = "Participants, unweighted n",
    statistic = "count"
  ),
  list(
    key = "kish_effective_n",
    label = "Kish effective sample size",
    statistic = "effective_n"
  ),
  list(
    key = "age",
    label = "Age, weighted mean (SD), years",
    statistic = "mean_sd"
  ),
  list(
    key = "female",
    label = "Women, weighted %",
    statistic = "percent"
  ),
  list(
    key = "education_no_formal",
    label = "No formal schooling, weighted %",
    statistic = "percent"
  ),
  list(
    key = "education_some_below_middle",
    label = "Some schooling below middle school, weighted %",
    statistic = "percent"
  ),
  list(
    key = "education_middle_plus",
    label = "Middle school or higher, weighted %",
    statistic = "percent"
  ),
  list(
    key = "rural",
    label = "Rural residence, weighted %",
    statistic = "percent"
  ),
  list(
    key = "agricultural_hukou",
    label = "Agricultural hukou, weighted %",
    statistic = "percent"
  ),
  list(
    key = "any_adl",
    label = "At least one ADL difficulty, weighted %",
    statistic = "percent"
  ),
  list(
    key = "any_iadl",
    label = "At least one IADL difficulty, weighted %",
    statistic = "percent"
  ),
  list(
    key = "any_routine_cognition",
    label = "At least one routine cognition component observed, weighted %",
    statistic = "percent"
  ),
  list(
    key = "all_routine_cognition",
    label = "All five routine cognition components observed, weighted %",
    statistic = "percent"
  )
)

table1_rows <- list()
row_index <- 1L
for (wave_value in 1:4) {
  target <- long[long$wave == wave_value, ]
  weight <- target$weight
  valid_hukou <- target$hukou %in% c(1, 2, 3)
  values <- list(
    unweighted_n = rep(1, nrow(target)),
    kish_effective_n = weight,
    age = target$age,
    female = target$female,
    education_no_formal = ifelse(
      is.na(target$education3),
      NA_real_,
      as.numeric(target$education3 == 1)
    ),
    education_some_below_middle = ifelse(
      is.na(target$education3),
      NA_real_,
      as.numeric(target$education3 == 2)
    ),
    education_middle_plus = ifelse(
      is.na(target$education3),
      NA_real_,
      as.numeric(target$education3 == 3)
    ),
    rural = target$rural,
    agricultural_hukou = ifelse(
      valid_hukou,
      as.numeric(target$hukou == 1),
      NA_real_
    ),
    any_adl = ifelse(
      is.na(target$adl),
      NA_real_,
      as.numeric(target$adl > 0)
    ),
    any_iadl = ifelse(
      is.na(target$iadl),
      NA_real_,
      as.numeric(target$iadl > 0)
    ),
    any_routine_cognition = as.numeric(
      rowSums(!is.na(target[routine_items])) >= 1
    ),
    all_routine_cognition = as.numeric(
      rowSums(!is.na(target[routine_items])) == length(routine_items)
    )
  )

  for (specification in characteristic_specs) {
    key <- specification$key
    statistic <- specification$statistic
    x <- values[[key]]
    if (statistic == "count") {
      estimate <- nrow(target)
      standard_deviation <- NA_real_
      denominator <- nrow(target)
    } else if (statistic == "effective_n") {
      estimate <- kish_effective_n(weight)
      standard_deviation <- NA_real_
      denominator <- sum(is.finite(weight) & weight > 0)
    } else if (statistic == "mean_sd") {
      estimate <- weighted_mean(x, weight)
      standard_deviation <- weighted_sd(x, weight)
      denominator <- sum(
        is.finite(x) & is.finite(weight) & weight > 0
      )
    } else {
      estimate <- 100 * weighted_mean(x, weight)
      standard_deviation <- NA_real_
      denominator <- sum(
        is.finite(x) & is.finite(weight) & weight > 0
      )
    }
    table1_rows[[row_index]] <- data.frame(
      wave = wave_value,
      year = unname(wave_years[as.character(wave_value)]),
      key = key,
      characteristic = specification$label,
      statistic = statistic,
      estimate = estimate,
      standard_deviation = standard_deviation,
      unweighted_denominator = denominator,
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}
table1_source <- bind_rows(table1_rows)

format_table1_value <- function(statistic, estimate, standard_deviation) {
  if (statistic == "count") {
    return(format(round(estimate), big.mark = ",", scientific = FALSE))
  }
  if (statistic == "effective_n") {
    return(format(round(estimate), big.mark = ",", scientific = FALSE))
  }
  if (statistic == "mean_sd") {
    return(sprintf("%.1f (%.1f)", estimate, standard_deviation))
  }
  sprintf("%.1f", estimate)
}

table1_source$display_value <- mapply(
  format_table1_value,
  table1_source$statistic,
  table1_source$estimate,
  table1_source$standard_deviation,
  USE.NAMES = FALSE
)
table1 <- table1_source %>%
  select(key, characteristic, year, display_value) %>%
  tidyr::pivot_wider(
    names_from = year,
    values_from = display_value,
    names_prefix = "year_"
  )

write.csv(
  table1_source,
  file.path(
    table_dir,
    paste0("b1_table1_characteristics_source_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  table1,
  file.path(
    table_dir,
    paste0("b1_table1_characteristics_", date_tag, ".csv")
  ),
  row.names = FALSE
)

model_fit <- read.csv(
  file.path(table_dir, "b1_hcap_fiml_model_fit_2026-06-19.csv"),
  check.names = FALSE
)
bifactor <- read.csv(
  file.path(table_dir, "b1_hcap_fiml_bifactor_indices_2026-06-19.csv"),
  check.names = FALSE
)
partial_metric <- read.csv(
  file.path(table_dir, "b1_hcap_partial_metric_fit_2026-06-19.csv"),
  check.names = FALSE
)
partial_scalar <- read.csv(
  file.path(table_dir, "b1_hcap_partial_scalar_fit_2026-06-19.csv"),
  check.names = FALSE
)
scoring_fit <- read.csv(
  file.path(table_dir, "b1_scoring_tool_fit_summary_2026-06-28.csv"),
  check.names = FALSE
)
scoring_validation <- read.csv(
  file.path(
    table_dir,
    "b1_scoring_tool_validation_summary_2026-06-28.csv"
  ),
  check.names = FALSE
)

anchor <- model_fit[model_fit$model == "bifactor_memory", ]
bifactor_main <- bifactor[bifactor$model == "bifactor_memory", ]
partial_metric_main <- partial_metric[
  partial_metric$model == "partial_metric_core",
]
partial_scalar_main <- partial_scalar[
  partial_scalar$model == "partial_scalar_core",
]
validation_overall <- scoring_validation[
  scoring_validation$group == "overall",
]
validation_education <- scoring_validation[
  scoring_validation$group != "overall",
]

table2 <- data.frame(
  stage = c(
    rep("HCAP bifactor anchor", 8),
    rep("Education partial invariance", 6),
    rep("Frozen transport and validation", 8)
  ),
  metric = c(
    "Participants with at least one HCAP indicator, n",
    "Robust CFI",
    "Robust TLI",
    "Robust RMSEA (90% CI)",
    "SRMR",
    "General-factor explained common variance",
    "Omega hierarchical",
    "General-score determinacy",
    "Full metric delta CFI vs configural",
    "Partial metric delta CFI vs configural",
    "Partial metric delta RMSEA vs configural",
    "Full scalar delta CFI vs selected partial metric",
    "Partial scalar delta CFI vs selected partial metric",
    "Partial scalar delta RMSEA vs selected partial metric",
    "Calibration participants per imputation, n",
    "Apparent weighted R-squared",
    "Weighted residual SD",
    "Validation weighted R-squared",
    "Validation RMSE, SD",
    "Calibration intercept",
    "Calibration slope",
    "Education-group calibration slope range"
  ),
  estimate = c(
    format(anchor$n_used, big.mark = ",", scientific = FALSE),
    sprintf("%.3f", anchor$cfi_robust),
    sprintf("%.3f", anchor$tli_robust),
    sprintf(
      "%.3f (%.3f-%.3f)",
      anchor$rmsea_robust,
      anchor$rmsea_ci_lower_robust,
      anchor$rmsea_ci_upper_robust
    ),
    sprintf("%.3f", anchor$srmr),
    sprintf("%.3f", bifactor_main$ecv_general),
    sprintf("%.3f", bifactor_main$omega_hierarchical),
    sprintf("%.3f", anchor$general_score_determinacy),
    sprintf(
      "%.3f",
      partial_metric$delta_cfi_vs_reference[
        partial_metric$model == "full_metric"
      ]
    ),
    sprintf("%.3f", partial_metric_main$delta_cfi_vs_reference),
    sprintf("%.3f", partial_metric_main$delta_rmsea_vs_reference),
    sprintf(
      "%.3f",
      partial_scalar$delta_cfi_vs_reference[
        partial_scalar$model == "full_scalar"
      ]
    ),
    sprintf("%.3f", partial_scalar_main$delta_cfi_vs_reference),
    sprintf("%.3f", partial_scalar_main$delta_rmsea_vs_reference),
    format(
      scoring_fit$n_calibration,
      big.mark = ",",
      scientific = FALSE
    ),
    sprintf("%.3f", scoring_fit$r_squared),
    sprintf("%.3f", scoring_fit$residual_sd),
    sprintf("%.3f", validation_overall$weighted_r_squared_mean),
    sprintf("%.3f", validation_overall$weighted_rmse_mean),
    sprintf("%.3f", validation_overall$calibration_intercept_mean),
    sprintf("%.3f", validation_overall$calibration_slope_mean),
    sprintf(
      "%.3f-%.3f",
      min(validation_education$calibration_slope_mean),
      max(validation_education$calibration_slope_mean)
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  table2,
  file.path(
    table_dir,
    paste0("b1_table2_measurement_calibration_", date_tag, ".csv")
  ),
  row.names = FALSE
)

markdown_table <- function(data) {
  header <- paste0("| ", paste(names(data), collapse = " | "), " |")
  separator <- paste0(
    "| ",
    paste(rep("---", ncol(data)), collapse = " | "),
    " |"
  )
  rows <- apply(
    data,
    1,
    function(row) paste0("| ", paste(row, collapse = " | "), " |")
  )
  c(header, separator, rows)
}

table1_markdown <- table1 %>%
  select(
    Characteristic = characteristic,
    `2011` = year_2011,
    `2013` = year_2013,
    `2015` = year_2015,
    `2018` = year_2018
  )
table2_markdown <- table2 %>%
  select(Stage = stage, Metric = metric, Estimate = estimate)

table_lines <- c(
  "# Main Analysis Tables",
  "",
  "## Table 1. Characteristics of CHARLS participants aged 60 years or older by survey wave",
  "",
  markdown_table(table1_markdown),
  "",
  paste(
    "Values are weighted percentages unless otherwise indicated.",
    "Weighted statistics use the wave-specific cross-sectional respondent",
    "weight and available observations for each characteristic. Agricultural",
    "hukou percentages exclude missing and unified/no-hukou categories.",
    "ADL indicates activities of daily living; IADL, instrumental activities",
    "of daily living."
  ),
  "",
  "## Table 2. HCAP measurement model, transport model, and internal validation performance",
  "",
  markdown_table(table2_markdown),
  "",
  paste(
    "CFI indicates comparative fit index; TLI, Tucker-Lewis index;",
    "RMSEA, root mean square error of approximation; SRMR, standardized",
    "root mean square residual; RMSE, root mean square error.",
    "Validation used ten community-separated folds across 50 imputations."
  )
)
writeLines(
  table_lines,
  file.path(table_dir, "main_analysis_tables.md"),
  useBytes = TRUE
)

cat("wrote B1 main Tables 1 and 2\n")
