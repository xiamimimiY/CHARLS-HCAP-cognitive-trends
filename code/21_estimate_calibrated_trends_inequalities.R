#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lavaan)
  library(mice)
  library(survey)
})

options(mc.cores = 1L)
options(survey.lonely.psu = "adjust")

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1) stop("Run this file with Rscript.")
script_path <- normalizePath(
  sub("^--file=", "", gsub("~\\+~", " ", script_argument)),
  mustWork = TRUE
)
study <- dirname(dirname(script_path))
long_input <- file.path(
  study,
  "data",
  "derived",
  "b1_frozen_transport_mi_input_2026-06-28.csv"
)
hcap_input <- file.path(
  study,
  "data",
  "derived",
  "b1_hcap_fiml_input_2026-06-19.csv"
)
table_dir <- file.path(study, "results", "tables")
report_path <- file.path(table_dir, "calibrated_trends_report.md")
sensitivity_report_path <- file.path(
  table_dir,
  "weight_missingness_sensitivity_report.md"
)
date_tag <- "2026-06-28"

m <- as.integer(Sys.getenv("B1_M", unset = "50"))
maxit <- as.integer(Sys.getenv("B1_MAXIT", unset = "10"))
run_sensitivity <- identical(
  Sys.getenv("B1_RUN_SENSITIVITY", unset = "1"),
  "1"
)
seed <- 3427L

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

long <- read.csv(long_input, check.names = FALSE, colClasses = c(
  ID = "character",
  pid = "character",
  household_id = "character",
  community_id = "character"
))
hcap <- read.csv(hcap_input, check.names = FALSE, colClasses = c(ID = "character"))

hcap_indicators <- c(
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

memory_indicators <- c(
  "z_test_short_learning_trial1",
  "z_test_short_delayed",
  "z_test_ten_learning_trial1",
  "z_test_ten_learning_trial2",
  "z_test_ten_learning_trial3",
  "z_test_ten_delayed",
  "z_test_recognition"
)

routine_items <- c("imrc", "dlrc", "ser7", "orient", "draw")

hcap <- hcap[
  !is.na(hcap$edu3) & rowSums(!is.na(hcap[hcap_indicators])) >= 1,
]

anchor_syntax <- paste(
  paste0("general =~ ", paste(hcap_indicators, collapse = " + ")),
  paste0("memory_specific =~ ", paste(memory_indicators, collapse = " + ")),
  "general ~~ 0*memory_specific",
  sep = "\n"
)

partial_paths <- c(
  "general=~z_test_orientation",
  "general=~z_test_animal_fluency",
  "memory_specific=~z_test_ten_delayed",
  "general=~z_test_serial7",
  "general=~z_test_recognition",
  "z_test_language~1",
  "z_test_orientation~1",
  "z_test_recognition~1",
  "z_test_animal_fluency~1",
  "z_test_serial7~1"
)

anchor_fit <- cfa(
  anchor_syntax,
  data = hcap,
  group = "edu3",
  group.equal = c("loadings", "intercepts"),
  group.partial = partial_paths,
  estimator = "MLR",
  missing = "fiml",
  std.lv = TRUE,
  meanstructure = TRUE
)

if (!isTRUE(lavInspect(anchor_fit, "converged"))) {
  stop("Frozen HCAP anchor did not converge.")
}

assemble_general_scores <- function(fit, n_rows) {
  scores <- lavPredict(fit, type = "lv")
  case_idx <- lavInspect(fit, "case.idx")
  result <- rep(NA_real_, n_rows)
  for (group_index in seq_along(scores)) {
    result[case_idx[[group_index]]] <- scores[[group_index]][, "general"]
  }
  result
}

weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

weighted_sd <- function(x, w) {
  center <- weighted_mean(x, w)
  sqrt(weighted_mean((x - center)^2, w))
}

weighted_quantile <- function(x, w, probability) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  order_index <- order(x)
  x <- x[order_index]
  w <- w[order_index]
  cumulative <- cumsum(w) / sum(w)
  x[which(cumulative >= probability)[1]]
}

weighted_fractional_rank <- function(x, w) {
  result <- rep(NA_real_, length(x))
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(result)
  values <- x[keep]
  weights <- w[keep]
  group_weight <- tapply(weights, values, sum)
  ordered_values <- as.numeric(names(group_weight))
  order_index <- order(ordered_values)
  ordered_values <- ordered_values[order_index]
  group_weight <- as.numeric(group_weight[order_index])
  midpoint <- (cumsum(group_weight) - 0.5 * group_weight) / sum(group_weight)
  rank_map <- setNames(midpoint, as.character(ordered_values))
  result[keep] <- unname(rank_map[as.character(values)])
  result
}

raw_anchor <- assemble_general_scores(anchor_fit, nrow(hcap))
anchor_mean <- weighted_mean(raw_anchor, hcap$r4wtrespb)
anchor_sd <- weighted_sd(raw_anchor, hcap$r4wtrespb)
hcap$anchor_z <- (raw_anchor - anchor_mean) / anchor_sd

thresholds <- sapply(
  c(0.10, 0.15, 0.20),
  function(probability) weighted_quantile(
    hcap$anchor_z,
    hcap$r4wtrespb,
    probability
  )
)
names(thresholds) <- c("10", "15", "20")

scaling <- do.call(rbind, lapply(routine_items, function(item) {
  hcap_name <- paste0("r4", item)
  data.frame(
    variable = item,
    mean_2018 = mean(hcap[[hcap_name]], na.rm = TRUE),
    sd_2018 = stats::sd(hcap[[hcap_name]], na.rm = TRUE)
  )
}))

anchor_map <- hcap[, c("ID", "anchor_z")]
long <- merge(long, anchor_map, by = "ID", all.x = TRUE)
long$anchor_z[long$wave != 4] <- NA_real_

age_sex_levels <- c(
  "60-69_male", "60-69_female",
  "70-79_male", "70-79_female",
  "80+_male", "80+_female"
)

long$age_band <- cut(
  long$age,
  breaks = c(59, 69, 79, Inf),
  labels = c("60-69", "70-79", "80+")
)
long$age_sex_cell <- factor(
  paste0(
    long$age_band,
    "_",
    ifelse(long$female == 1, "female", "male")
  ),
  levels = age_sex_levels
)

reference <- long[
  long$wave == 1 & is.finite(long$weight) & long$weight > 0,
]
reference_totals <- tapply(
  reference$weight,
  reference$age_sex_cell,
  sum,
  na.rm = TRUE
)
reference_totals <- reference_totals[age_sex_levels]
reference_totals[is.na(reference_totals)] <- 0

mi_variables <- c(
  routine_items,
  "age", "female", "education3", "rural", "hukou",
  "adl", "iadl", "log_consumption"
)

imputations <- list()
missing_rows <- list()

for (wave_value in sort(unique(long$wave))) {
  wave_data <- long[long$wave == wave_value, ]
  cat(
    "starting MICE for wave ",
    wave_value,
    " (",
    unique(wave_data$year),
    "), m=",
    m,
    ", maxit=",
    maxit,
    "\n",
    sep = ""
  )
  missing_rows[[length(missing_rows) + 1]] <- data.frame(
    wave = wave_value,
    year = unique(wave_data$year),
    variable = mi_variables,
    n = nrow(wave_data),
    n_missing = sapply(wave_data[mi_variables], function(x) sum(is.na(x))),
    missing_pct = sapply(wave_data[mi_variables], function(x) mean(is.na(x)) * 100)
  )

  mice_data <- wave_data[mi_variables]
  mice_data$education3 <- factor(
    mice_data$education3,
    levels = c(1, 2, 3)
  )
  mice_data$rural <- factor(mice_data$rural, levels = c(0, 1))
  mice_data$hukou <- factor(mice_data$hukou, levels = c(1, 2, 3, 4))

  method <- rep("", ncol(mice_data))
  names(method) <- names(mice_data)
  method[routine_items] <- "pmm"
  method["education3"] <- "polyreg"
  method["rural"] <- "logreg"
  method["hukou"] <- "polyreg"
  method[c("adl", "iadl", "log_consumption")] <- "pmm"

  predictor <- make.predictorMatrix(mice_data)
  predictor[,] <- 0
  imputed_targets <- names(method)[method != ""]
  fixed_predictors <- c("age", "female")
  for (target in imputed_targets) {
    predictor[target, setdiff(c(imputed_targets, fixed_predictors), target)] <- 1
  }

  imputations[[as.character(wave_value)]] <- mice(
    mice_data,
    m = m,
    maxit = maxit,
    method = method,
    predictorMatrix = predictor,
    seed = seed + wave_value,
    printFlag = FALSE
  )
  cat("completed MICE for wave ", wave_value, "\n", sep = "")
}

missing_table <- do.call(rbind, missing_rows)

apply_frozen_scaling <- function(data) {
  result <- data
  for (item in routine_items) {
    row <- scaling[scaling$variable == item, ]
    result[[paste0("z_", item)]] <- (
      result[[item]] - row$mean_2018
    ) / row$sd_2018
  }
  result$routine_composite <- rowMeans(
    result[paste0("z_", routine_items)]
  )
  result$age_decade <- (
    result$age - mean(hcap$r4agey, na.rm = TRUE)
  ) / 10
  result$education3 <- as.integer(as.character(result$education3))
  result$edu_some_schooling_below_middle <- as.numeric(result$education3 == 2)
  result$edu_middle_school_plus <- as.numeric(result$education3 == 3)
  result$education_label <- factor(
    result$education3,
    levels = c(1, 2, 3),
    labels = c(
      "no_formal",
      "some_schooling_below_middle",
      "middle_school_plus"
    )
  )
  result$residence_label <- factor(
    as.integer(as.character(result$rural)),
    levels = c(0, 1),
    labels = c("urban", "rural")
  )
  hukou_numeric <- as.integer(as.character(result$hukou))
  result$hukou_label <- factor(
    ifelse(
      hukou_numeric == 1,
      "agricultural",
      ifelse(hukou_numeric %in% c(2, 3), "non_agricultural", NA)
    ),
    levels = c("non_agricultural", "agricultural")
  )
  result$age_band <- cut(
    result$age,
    breaks = c(59, 69, 79, Inf),
    labels = c("60-69", "70-79", "80+")
  )
  result$age_sex_cell <- factor(
    paste0(
      result$age_band,
      "_",
      ifelse(result$female == 1, "female", "male")
    ),
    levels = age_sex_levels
  )
  result$wave_factor <- factor(
    result$year,
    levels = c(2011, 2013, 2015, 2018)
  )
  result$year_centered <- result$year - 2011
  result
}

primary_formula <- anchor_z ~ z_imrc + z_dlrc + z_ser7 +
  z_orient + z_draw + age_decade + I(age_decade^2) + female +
  edu_some_schooling_below_middle + edu_middle_school_plus +
  routine_composite:edu_some_schooling_below_middle +
  routine_composite:edu_middle_school_plus

blind_formula <- anchor_z ~ z_imrc + z_dlrc + z_ser7 +
  z_orient + z_draw + age_decade + I(age_decade^2) + female

extract_mean <- function(design, variable, estimand, model, imputation, wave, year) {
  estimate <- svymean(
    as.formula(paste0("~", variable)),
    design,
    na.rm = TRUE
  )
  data.frame(
    imputation = imputation,
    wave = wave,
    year = year,
    model = model,
    dimension = "overall",
    group = "overall",
    estimand = estimand,
    estimate = unname(coef(estimate)[1]),
    variance = unname(vcov(estimate)[1, 1])
  )
}

extract_group <- function(
  design,
  variable,
  group_variable,
  dimension,
  contrast_levels,
  estimand,
  model,
  imputation,
  wave,
  year
) {
  group_fit <- svyglm(
    as.formula(paste0(variable, " ~ 0 + ", group_variable)),
    design = design,
    family = gaussian(),
    na.action = na.omit
  )
  estimates <- coef(group_fit)
  covariance <- vcov(group_fit)
  labels <- sub(
    paste0("^", group_variable),
    "",
    names(estimates)
  )
  rows <- lapply(seq_along(estimates), function(i) {
    data.frame(
      imputation = imputation,
      wave = wave,
      year = year,
      model = model,
      dimension = dimension,
      group = labels[i],
      estimand = estimand,
      estimate = unname(estimates[i]),
      variance = unname(covariance[i, i])
    )
  })
  if (all(contrast_levels %in% labels)) {
    first <- which(labels == contrast_levels[1])[1]
    second <- which(labels == contrast_levels[2])[1]
    contrast_variance <- (
      covariance[first, first] +
        covariance[second, second] -
        2 * covariance[first, second]
    )
    rows[[length(rows) + 1]] <- data.frame(
      imputation = imputation,
      wave = wave,
      year = year,
      model = model,
      dimension = dimension,
      group = paste0(
        contrast_levels[1],
        "_minus_",
        contrast_levels[2]
      ),
      estimand = paste0(estimand, "_absolute_gap"),
      estimate = unname(estimates[first] - estimates[second]),
      variance = unname(contrast_variance)
    )
  }
  do.call(rbind, rows)
}

contrast_row <- function(
  estimates,
  covariance,
  contrast,
  imputation,
  model,
  dimension,
  group,
  estimand
) {
  contrast <- as.numeric(contrast)
  data.frame(
    imputation = imputation,
    wave = 0,
    year = 0,
    model = model,
    dimension = dimension,
    group = group,
    estimand = estimand,
    estimate = unname(sum(contrast * estimates)),
    variance = unname(
      as.numeric(t(contrast) %*% covariance %*% contrast)
    )
  )
}

extract_wave_contrasts <- function(
  design,
  variable,
  estimand,
  model,
  imputation
) {
  wave_fit <- svyglm(
    as.formula(paste0(variable, " ~ 0 + wave_factor")),
    design = design,
    family = gaussian(),
    na.action = na.omit
  )
  estimates <- coef(wave_fit)
  covariance <- vcov(wave_fit)
  labels <- sub("^wave_factor", "", names(estimates))
  contrast_pairs <- list(
    `2013_minus_2011` = c("2013", "2011"),
    `2015_minus_2013` = c("2015", "2013"),
    `2018_minus_2015` = c("2018", "2015"),
    `2018_minus_2011` = c("2018", "2011")
  )
  rows <- lapply(names(contrast_pairs), function(label) {
    pair <- contrast_pairs[[label]]
    contrast <- rep(0, length(estimates))
    contrast[labels == pair[1]] <- 1
    contrast[labels == pair[2]] <- -1
    contrast_row(
      estimates,
      covariance,
      contrast,
      imputation,
      model,
      "overall_trend",
      label,
      paste0(estimand, "_change")
    )
  })

  trend_fit <- svyglm(
    as.formula(paste0(variable, " ~ year_centered")),
    design = design,
    family = gaussian(),
    na.action = na.omit
  )
  rows[[length(rows) + 1]] <- data.frame(
    imputation = imputation,
    wave = 0,
    year = 0,
    model = model,
    dimension = "overall_trend",
    group = "annual_linear",
    estimand = paste0(estimand, "_annual_linear_change"),
    estimate = unname(coef(trend_fit)["year_centered"]),
    variance = unname(vcov(trend_fit)["year_centered", "year_centered"])
  )
  do.call(rbind, rows)
}

extract_group_gap_changes <- function(
  data,
  variable,
  group_variable,
  dimension,
  contrast_levels,
  estimand,
  model,
  imputation,
  population = reference_totals
) {
  target <- data[
    !is.na(data[[group_variable]]) &
      is.finite(data$weight) &
      data$weight > 0 &
      !is.na(data$community_id),
  ]
  target$wave_group <- interaction(
    target$wave_factor,
    target[[group_variable]],
    sep = "__",
    drop = TRUE
  )
  base_design <- svydesign(
    ids = ~community_id,
    weights = ~weight,
    data = target,
    nest = TRUE
  )
  standardized <- svystandardize(
    base_design,
    by = ~age_sex_cell,
    over = ~wave_group,
    population = population
  )
  group_fit <- svyglm(
    as.formula(paste0(variable, " ~ 0 + wave_group")),
    design = standardized,
    family = gaussian(),
    na.action = na.omit
  )
  estimates <- coef(group_fit)
  covariance <- vcov(group_fit)
  labels <- sub("^wave_group", "", names(estimates))

  gap_contrast <- function(year) {
    contrast <- rep(0, length(estimates))
    contrast[labels == paste0(year, "__", contrast_levels[1])] <- 1
    contrast[labels == paste0(year, "__", contrast_levels[2])] <- -1
    contrast
  }
  gaps <- lapply(c("2011", "2013", "2015", "2018"), gap_contrast)
  names(gaps) <- c("2011", "2013", "2015", "2018")
  change_pairs <- list(
    `2013_minus_2011` = c("2013", "2011"),
    `2015_minus_2013` = c("2015", "2013"),
    `2018_minus_2015` = c("2018", "2015"),
    `2018_minus_2011` = c("2018", "2011")
  )
  rows <- lapply(names(change_pairs), function(label) {
    pair <- change_pairs[[label]]
    contrast_row(
      estimates,
      covariance,
      gaps[[pair[1]]] - gaps[[pair[2]]],
      imputation,
      model,
      paste0(dimension, "_gap_trend"),
      label,
      paste0(estimand, "_absolute_gap_change")
    )
  })
  do.call(rbind, rows)
}

extract_sii <- function(
  design,
  variable,
  estimand,
  model,
  imputation,
  wave,
  year
) {
  fit <- svyglm(
    as.formula(paste0(variable, " ~ consumption_rank + age_sex_cell")),
    design = design,
    family = gaussian(),
    na.action = na.omit
  )
  data.frame(
    imputation = imputation,
    wave = wave,
    year = year,
    model = model,
    dimension = "consumption_rank",
    group = "poorest_minus_richest",
    estimand = paste0(estimand, "_sii"),
    estimate = -unname(coef(fit)["consumption_rank"]),
    variance = unname(vcov(fit)["consumption_rank", "consumption_rank"])
  )
}

extract_sii_change <- function(
  design,
  variable,
  estimand,
  model,
  imputation
) {
  fit <- svyglm(
    as.formula(
      paste0(
        variable,
        " ~ consumption_rank * wave_factor + age_sex_cell"
      )
    ),
    design = design,
    family = gaussian(),
    na.action = na.omit
  )
  coefficient_name <- grep(
    "consumption_rank:wave_factor2018|wave_factor2018:consumption_rank",
    names(coef(fit)),
    value = TRUE
  )
  if (length(coefficient_name) != 1) {
    stop("Could not identify the 2018 consumption-rank interaction.")
  }
  data.frame(
    imputation = imputation,
    wave = 0,
    year = 0,
    model = model,
    dimension = "consumption_rank_trend",
    group = "2018_minus_2011",
    estimand = paste0(estimand, "_sii_change"),
    estimate = -unname(coef(fit)[coefficient_name]),
    variance = unname(vcov(fit)[coefficient_name, coefficient_name])
  )
}

estimate_completed_data <- function(data, imputation) {
  wave4 <- data[data$wave == 4 & !is.na(data$anchor_z), ]
  primary_fit <- lm(
    primary_formula,
    data = wave4,
    weights = weight
  )
  blind_fit <- lm(
    blind_formula,
    data = wave4,
    weights = weight
  )
  model_fits <- list(primary = primary_fit, social_blind = blind_fit)
  result_rows <- list()
  fit_rows <- list()

  for (model_name in names(model_fits)) {
    fit <- model_fits[[model_name]]
    residual_sd <- sqrt(weighted_mean(
      residuals(fit)^2,
      model.weights(model.frame(fit))
    ))
    fit_rows[[length(fit_rows) + 1]] <- data.frame(
      imputation = imputation,
      model = model_name,
      n_calibration = nobs(fit),
      residual_sd = residual_sd,
      r_squared = summary(fit)$r.squared
    )
    data[[paste0("pred_", model_name)]] <- predict(fit, newdata = data)
    for (threshold_name in names(thresholds)) {
      data[[paste0("low", threshold_name, "_", model_name)]] <- pnorm(
        (
          thresholds[[threshold_name]] -
            data[[paste0("pred_", model_name)]]
        ) / residual_sd
      )
    }
  }

  wave4_population <- data[
    data$wave == 4 &
      is.finite(data$weight) &
      data$weight > 0,
  ]
  uncalibrated_mean <- weighted_mean(
    wave4_population$routine_composite,
    wave4_population$weight
  )
  uncalibrated_sd <- weighted_sd(
    wave4_population$routine_composite,
    wave4_population$weight
  )
  data$pred_uncalibrated <- (
    data$routine_composite - uncalibrated_mean
  ) / uncalibrated_sd
  uncalibrated_thresholds <- sapply(
    c(0.10, 0.15, 0.20),
    function(probability) weighted_quantile(
      data$pred_uncalibrated[data$wave == 4],
      data$weight[data$wave == 4],
      probability
    )
  )
  names(uncalibrated_thresholds) <- c("10", "15", "20")
  for (threshold_name in names(uncalibrated_thresholds)) {
    data[[paste0("low", threshold_name, "_uncalibrated")]] <- as.numeric(
      data$pred_uncalibrated < uncalibrated_thresholds[[threshold_name]]
    )
  }

  comparison_specs <- list(
    primary_minus_uncalibrated = c("primary", "uncalibrated"),
    primary_minus_social_blind = c("primary", "social_blind")
  )
  for (comparison_name in names(comparison_specs)) {
    pair <- comparison_specs[[comparison_name]]
    data[[paste0("pred_", comparison_name)]] <- (
      data[[paste0("pred_", pair[1])]] -
        data[[paste0("pred_", pair[2])]]
    )
    for (threshold_name in names(thresholds)) {
      data[[paste0("low", threshold_name, "_", comparison_name)]] <- (
        data[[paste0("low", threshold_name, "_", pair[1])]] -
          data[[paste0("low", threshold_name, "_", pair[2])]]
      )
    }
  }

  analysis_models <- c(
    "primary",
    "social_blind",
    "uncalibrated",
    names(comparison_specs)
  )
  group_specs <- list(
    education = list(
      variable = "education_label",
      contrast = c("no_formal", "middle_school_plus")
    ),
    residence = list(
      variable = "residence_label",
      contrast = c("rural", "urban")
    ),
    hukou = list(
      variable = "hukou_label",
      contrast = c("agricultural", "non_agricultural")
    )
  )

  for (wave_value in sort(unique(data$wave))) {
    target <- data[
      data$wave == wave_value &
        is.finite(data$weight) &
        data$weight > 0 &
        !is.na(data$community_id),
    ]
    year_value <- unique(target$year)
    design <- svydesign(
      ids = ~community_id,
      weights = ~weight,
      data = target,
      nest = TRUE
    )
    standardized <- svystandardize(
      design,
      by = ~age_sex_cell,
      over = ~1,
      population = reference_totals
    )

    target$consumption_rank <- weighted_fractional_rank(
      target$log_consumption,
      target$weight
    )
    consumption_design <- svydesign(
      ids = ~community_id,
      weights = ~weight,
      data = target[is.finite(target$consumption_rank), ],
      nest = TRUE
    )

    for (model_name in analysis_models) {
      pred_variable <- paste0("pred_", model_name)
      result_rows[[length(result_rows) + 1]] <- extract_mean(
        standardized,
        pred_variable,
        "standardized_mean_cognition",
        model_name,
        imputation,
        wave_value,
        year_value
      )
      for (threshold_name in names(thresholds)) {
        low_variable <- paste0("low", threshold_name, "_", model_name)
        result_rows[[length(result_rows) + 1]] <- extract_mean(
          standardized,
          low_variable,
          paste0("standardized_low", threshold_name, "_probability"),
          model_name,
          imputation,
          wave_value,
          year_value
        )
      }

      result_rows[[length(result_rows) + 1]] <- extract_sii(
        consumption_design,
        pred_variable,
        "age_sex_adjusted_mean_cognition",
        model_name,
        imputation,
        wave_value,
        year_value
      )
      result_rows[[length(result_rows) + 1]] <- extract_sii(
        consumption_design,
        paste0("low15_", model_name),
        "age_sex_adjusted_low15_probability",
        model_name,
        imputation,
        wave_value,
        year_value
      )

      for (dimension_name in names(group_specs)) {
        spec <- group_specs[[dimension_name]]
        group_target <- target[!is.na(target[[spec$variable]]), ]
        group_base_design <- svydesign(
          ids = ~community_id,
          weights = ~weight,
          data = group_target,
          nest = TRUE
        )
        group_design <- svystandardize(
          group_base_design,
          by = ~age_sex_cell,
          over = as.formula(paste0("~", spec$variable)),
          population = reference_totals
        )
        result_rows[[length(result_rows) + 1]] <- extract_group(
          group_design,
          pred_variable,
          spec$variable,
          dimension_name,
          spec$contrast,
          "standardized_mean_cognition",
          model_name,
          imputation,
          wave_value,
          year_value
        )
        result_rows[[length(result_rows) + 1]] <- extract_group(
          group_design,
          paste0("low15_", model_name),
          spec$variable,
          dimension_name,
          spec$contrast,
          "standardized_low15_probability",
          model_name,
          imputation,
          wave_value,
          year_value
        )
      }
    }
  }

  joint_target <- data[
    is.finite(data$weight) &
      data$weight > 0 &
      !is.na(data$community_id),
  ]
  joint_design <- svydesign(
    ids = ~community_id,
    weights = ~weight,
    data = joint_target,
    nest = TRUE
  )
  joint_standardized <- svystandardize(
    joint_design,
    by = ~age_sex_cell,
    over = ~wave_factor,
    population = reference_totals
  )

  joint_target$consumption_rank <- ave(
    seq_len(nrow(joint_target)),
    joint_target$wave,
    FUN = function(index) weighted_fractional_rank(
      joint_target$log_consumption[index],
      joint_target$weight[index]
    )
  )
  joint_consumption_design <- svydesign(
    ids = ~community_id,
    weights = ~weight,
    data = joint_target[is.finite(joint_target$consumption_rank), ],
    nest = TRUE
  )

  for (model_name in analysis_models) {
    pred_variable <- paste0("pred_", model_name)
    result_rows[[length(result_rows) + 1]] <- extract_wave_contrasts(
      joint_standardized,
      pred_variable,
      "standardized_mean_cognition",
      model_name,
      imputation
    )
    result_rows[[length(result_rows) + 1]] <- extract_wave_contrasts(
      joint_standardized,
      paste0("low15_", model_name),
      "standardized_low15_probability",
      model_name,
      imputation
    )
    result_rows[[length(result_rows) + 1]] <- extract_sii_change(
      joint_consumption_design,
      pred_variable,
      "age_sex_adjusted_mean_cognition",
      model_name,
      imputation
    )
    result_rows[[length(result_rows) + 1]] <- extract_sii_change(
      joint_consumption_design,
      paste0("low15_", model_name),
      "age_sex_adjusted_low15_probability",
      model_name,
      imputation
    )

    for (dimension_name in names(group_specs)) {
      spec <- group_specs[[dimension_name]]
      result_rows[[length(result_rows) + 1]] <- extract_group_gap_changes(
        joint_target,
        pred_variable,
        spec$variable,
        dimension_name,
        spec$contrast,
        "standardized_mean_cognition",
        model_name,
        imputation
      )
      result_rows[[length(result_rows) + 1]] <- extract_group_gap_changes(
        joint_target,
        paste0("low15_", model_name),
        spec$variable,
        dimension_name,
        spec$contrast,
        "standardized_low15_probability",
        model_name,
        imputation
      )
    }
  }

  list(
    estimates = do.call(rbind, result_rows),
    fits = do.call(rbind, fit_rows)
  )
}

apply_weight_scenario <- function(data, scenario) {
  result <- data
  if (scenario == "unweighted") {
    result$weight <- 1
  } else if (scenario == "weight_trim_99") {
    for (wave_value in sort(unique(result$wave))) {
      index <- result$wave == wave_value &
        is.finite(result$weight) &
        result$weight > 0
      cap <- stats::quantile(
        result$weight[index],
        probs = 0.99,
        na.rm = TRUE,
        names = FALSE
      )
      result$weight[index] <- pmin(result$weight[index], cap)
    }
  }
  result
}

sensitivity_reference <- function(data) {
  reference <- data[
    data$wave == 1 &
      is.finite(data$weight) &
      data$weight > 0,
  ]
  totals <- tapply(
    reference$weight,
    reference$age_sex_cell,
    sum,
    na.rm = TRUE
  )
  totals <- totals[age_sex_levels]
  totals[is.na(totals)] <- 0
  totals
}

estimate_sensitivity_completed_data <- function(
  data,
  imputation,
  scenario,
  prediction_delta = 0
) {
  data <- apply_weight_scenario(data, scenario)
  scenario_reference <- sensitivity_reference(data)
  wave4 <- data[
    data$wave == 4 &
      !is.na(data$anchor_z) &
      is.finite(data$weight) &
      data$weight > 0,
  ]
  fit <- lm(
    primary_formula,
    data = wave4,
    weights = weight
  )
  residual_sd <- sqrt(weighted_mean(
    residuals(fit)^2,
    model.weights(model.frame(fit))
  ))
  data$pred_primary <- predict(fit, newdata = data)
  if (prediction_delta > 0) {
    high_missing <- data$original_cognition_missing_count >= 3
    data$pred_primary[high_missing] <- (
      data$pred_primary[high_missing] - prediction_delta
    )
  }
  data$low15_primary <- pnorm(
    (thresholds[["15"]] - data$pred_primary) / residual_sd
  )

  rows <- list()
  for (wave_value in sort(unique(data$wave))) {
    target <- data[
      data$wave == wave_value &
        is.finite(data$weight) &
        data$weight > 0 &
        !is.na(data$community_id),
    ]
    year_value <- unique(target$year)
    design <- svydesign(
      ids = ~community_id,
      weights = ~weight,
      data = target,
      nest = TRUE
    )
    standardized <- svystandardize(
      design,
      by = ~age_sex_cell,
      over = ~1,
      population = scenario_reference
    )
    rows[[length(rows) + 1]] <- extract_mean(
      standardized,
      "pred_primary",
      "standardized_mean_cognition",
      "primary",
      imputation,
      wave_value,
      year_value
    )
    rows[[length(rows) + 1]] <- extract_mean(
      standardized,
      "low15_primary",
      "standardized_low15_probability",
      "primary",
      imputation,
      wave_value,
      year_value
    )
  }

  joint_target <- data[
    is.finite(data$weight) &
      data$weight > 0 &
      !is.na(data$community_id),
  ]
  joint_design <- svydesign(
    ids = ~community_id,
    weights = ~weight,
    data = joint_target,
    nest = TRUE
  )
  joint_standardized <- svystandardize(
    joint_design,
    by = ~age_sex_cell,
    over = ~wave_factor,
    population = scenario_reference
  )
  rows[[length(rows) + 1]] <- extract_wave_contrasts(
    joint_standardized,
    "pred_primary",
    "standardized_mean_cognition",
    "primary",
    imputation
  )
  rows[[length(rows) + 1]] <- extract_wave_contrasts(
    joint_standardized,
    "low15_primary",
    "standardized_low15_probability",
    "primary",
    imputation
  )

  group_specs <- list(
    education = list(
      variable = "education_label",
      contrast = c("no_formal", "middle_school_plus")
    ),
    residence = list(
      variable = "residence_label",
      contrast = c("rural", "urban")
    ),
    hukou = list(
      variable = "hukou_label",
      contrast = c("agricultural", "non_agricultural")
    )
  )
  for (dimension_name in names(group_specs)) {
    spec <- group_specs[[dimension_name]]
    rows[[length(rows) + 1]] <- extract_group_gap_changes(
      joint_target,
      "pred_primary",
      spec$variable,
      dimension_name,
      spec$contrast,
      "standardized_mean_cognition",
      "primary",
      imputation,
      population = scenario_reference
    )
    rows[[length(rows) + 1]] <- extract_group_gap_changes(
      joint_target,
      "low15_primary",
      spec$variable,
      dimension_name,
      spec$contrast,
      "standardized_low15_probability",
      "primary",
      imputation,
      population = scenario_reference
    )
  }

  result <- do.call(rbind, rows)
  result$scenario <- scenario
  result$n_calibration <- nobs(fit)
  result
}

all_estimates <- list()
all_fits <- list()
all_sensitivity_estimates <- list()

for (imputation_index in seq_len(m)) {
  completed_waves <- lapply(sort(unique(long$wave)), function(wave_value) {
    base <- long[long$wave == wave_value, ]
    base$original_cognition_missing_count <- rowSums(
      is.na(base[routine_items])
    )
    for (item in routine_items) {
      base[[paste0("original_missing_", item)]] <- is.na(base[[item]])
    }
    completed <- complete(
      imputations[[as.character(wave_value)]],
      imputation_index
    )
    base[mi_variables] <- completed[mi_variables]
    apply_frozen_scaling(base)
  })
  completed_long <- do.call(rbind, completed_waves)
  result <- estimate_completed_data(completed_long, imputation_index)
  all_estimates[[imputation_index]] <- result$estimates
  all_fits[[imputation_index]] <- result$fits
  if (run_sensitivity) {
    sensitivity_rows <- list()
    for (scenario in c("unweighted", "weight_trim_99")) {
      sensitivity_rows[[length(sensitivity_rows) + 1]] <-
        estimate_sensitivity_completed_data(
          completed_long,
          imputation_index,
          scenario
        )
    }
    for (delta in c(0.25, 0.50)) {
      scenario <- paste0("mnar_delta_", format(delta, nsmall = 2))
      sensitivity_rows[[length(sensitivity_rows) + 1]] <-
        estimate_sensitivity_completed_data(
          completed_long,
          imputation_index,
          scenario,
          prediction_delta = delta
        )
    }
    all_sensitivity_estimates[[imputation_index]] <- do.call(
      rbind,
      sensitivity_rows
    )
  }
  cat("completed imputation ", imputation_index, " of ", m, "\n", sep = "")
}

estimate_long <- do.call(rbind, all_estimates)
fit_long <- do.call(rbind, all_fits)

rubin_pool <- function(data) {
  qbar <- mean(data$estimate)
  ubar <- mean(data$variance)
  imputations_n <- nrow(data)
  between <- if (imputations_n > 1) stats::var(data$estimate) else 0
  added_variance <- (1 + 1 / imputations_n) * between
  total <- ubar + added_variance
  relative_increase <- ifelse(
    ubar > 0,
    added_variance / ubar,
    Inf
  )
  degrees_freedom <- ifelse(
    between <= 0,
    Inf,
    (imputations_n - 1) * (1 + 1 / relative_increase)^2
  )
  critical_value <- ifelse(
    is.finite(degrees_freedom),
    stats::qt(0.975, degrees_freedom),
    stats::qnorm(0.975)
  )
  statistic <- qbar / sqrt(total)
  data.frame(
    estimate = qbar,
    se = sqrt(total),
    degrees_freedom = degrees_freedom,
    statistic = statistic,
    p_value = 2 * stats::pt(
      -abs(statistic),
      df = degrees_freedom
    ),
    ci_lower = qbar - critical_value * sqrt(total),
    ci_upper = qbar + critical_value * sqrt(total),
    within_variance = ubar,
    between_variance = between,
    fraction_missing_information = ifelse(
      total > 0,
      added_variance / total,
      0
    )
  )
}

if (run_sensitivity) {
  sensitivity_long <- do.call(rbind, all_sensitivity_estimates)

  complete_data <- long[
    stats::complete.cases(
      long[, c(
        routine_items,
        "age",
        "female",
        "education3",
        "weight",
        "community_id"
      )]
    ) &
      long$weight > 0,
  ]
  complete_data$original_cognition_missing_count <- 0
  for (item in routine_items) {
    complete_data[[paste0("original_missing_", item)]] <- FALSE
  }
  complete_data <- apply_frozen_scaling(complete_data)
  complete_result <- estimate_sensitivity_completed_data(
    complete_data,
    1,
    "complete_predictor"
  )
  sensitivity_long <- rbind(sensitivity_long, complete_result)

  sensitivity_pool_keys <- c(
    "scenario",
    "wave",
    "year",
    "model",
    "dimension",
    "group",
    "estimand"
  )
  sensitivity_pooled <- do.call(rbind, lapply(
    split(
      sensitivity_long,
      interaction(
        sensitivity_long[sensitivity_pool_keys],
        drop = TRUE,
        lex.order = TRUE
      )
    ),
    function(target) {
      pooled_stats <- rubin_pool(target)
      data.frame(
        target[1, sensitivity_pool_keys],
        m = nrow(target),
        n_calibration_min = min(target$n_calibration),
        n_calibration_max = max(target$n_calibration),
        pooled_stats
      )
    }
  ))
  row.names(sensitivity_pooled) <- NULL
}

pool_keys <- c(
  "wave", "year", "model", "dimension", "group", "estimand"
)
pooled <- do.call(rbind, lapply(
  split(
    estimate_long,
    interaction(estimate_long[pool_keys], drop = TRUE, lex.order = TRUE)
  ),
  function(target) {
    pooled_stats <- rubin_pool(target)
    data.frame(target[1, pool_keys], m = nrow(target), pooled_stats)
  }
))
row.names(pooled) <- NULL

fit_summary <- aggregate(
  fit_long[, c("n_calibration", "residual_sd", "r_squared")],
  by = fit_long[, c("model"), drop = FALSE],
  FUN = function(x) c(
    mean = mean(x),
    between_sd = stats::sd(x),
    min = min(x),
    max = max(x)
  )
)

formal_summary <- pooled[
  (
    pooled$dimension == "overall" &
      pooled$model %in% c("primary", "uncalibrated")
  ) |
    (
      pooled$group == "2018_minus_2011" &
        pooled$model %in% c(
          "primary",
          "uncalibrated",
          "primary_minus_uncalibrated",
          "primary_minus_social_blind"
        )
    ) |
    (
      pooled$year == 2018 &
        pooled$dimension == "education" &
        pooled$group == "no_formal_minus_middle_school_plus" &
        pooled$model %in% c(
          "primary",
          "social_blind",
          "uncalibrated",
          "primary_minus_uncalibrated",
          "primary_minus_social_blind"
        )
    ),
]

write.csv(
  missing_table,
  file.path(
    table_dir,
    paste0("b1_frozen_transport_mi_missingness_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  fit_long,
  file.path(
    table_dir,
    paste0("b1_frozen_transport_mi_fit_long_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  fit_summary,
  file.path(
    table_dir,
    paste0("b1_frozen_transport_mi_fit_summary_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  pooled,
  file.path(
    table_dir,
    paste0("b1_frozen_transport_mi_pooled_estimands_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  formal_summary,
  file.path(
    table_dir,
    paste0("b1_formal_trend_calibration_impact_", date_tag, ".csv")
  ),
  row.names = FALSE
)
if (run_sensitivity) {
  write.csv(
    sensitivity_pooled,
    file.path(
      table_dir,
      paste0("b1_weight_mnar_sensitivity_", date_tag, ".csv")
    ),
    row.names = FALSE
  )
}

display <- pooled[
  pooled$dimension == "overall" &
    pooled$model == "primary" &
    pooled$estimand %in% c(
      "standardized_mean_cognition",
      "standardized_low15_probability"
    ),
]
display <- display[order(display$year, display$estimand), ]

fmt <- function(x, digits = 3) {
  ifelse(
    is.na(x) | !is.finite(x),
    "NA",
    formatC(x, digits = digits, format = "f")
  )
}

display_lines <- c(
  "| Year | Estimand | Estimate | 95% CI | FMI |",
  "| ---: | --- | ---: | --- | ---: |"
)
for (i in seq_len(nrow(display))) {
  r <- display[i, ]
  multiplier <- if (grepl("probability", r$estimand)) 100 else 1
  suffix <- if (multiplier == 100) "%" else " SD"
  display_lines <- c(
    display_lines,
    sprintf(
      "| %s | %s | %s%s | %s to %s%s | %s |",
      r$year,
      r$estimand,
      fmt(r$estimate * multiplier, 2),
      suffix,
      fmt(r$ci_lower * multiplier, 2),
      fmt(r$ci_upper * multiplier, 2),
      suffix,
      fmt(r$fraction_missing_information, 3)
    )
  )
}

result_lines <- function(data, labels) {
  lines <- c(
    "| Result | Estimate | 95% CI | P | FMI |",
    "| --- | ---: | --- | ---: | ---: |"
  )
  for (i in seq_len(nrow(data))) {
    row <- data[i, ]
    multiplier <- if (grepl("probability", row$estimand)) 100 else 1
    suffix <- if (multiplier == 100) " pp" else " SD"
    lines <- c(
      lines,
      sprintf(
        "| %s | %s%s | %s to %s%s | %s | %s |",
        labels[i],
        fmt(row$estimate * multiplier, 2),
        suffix,
        fmt(row$ci_lower * multiplier, 2),
        fmt(row$ci_upper * multiplier, 2),
        suffix,
        ifelse(
          row$p_value < 0.001,
          "<0.001",
          fmt(row$p_value, 3)
        ),
        fmt(row$fraction_missing_information, 3)
      )
    )
  }
  lines
}

trend_display <- pooled[
  pooled$model == "primary" &
    pooled$dimension == "overall_trend" &
    pooled$group == "2018_minus_2011" &
    pooled$estimand %in% c(
      "standardized_mean_cognition_change",
      "standardized_low15_probability_change"
    ),
]
trend_display <- trend_display[
  match(
    c(
      "standardized_mean_cognition_change",
      "standardized_low15_probability_change"
    ),
    trend_display$estimand
  ),
]
trend_lines <- result_lines(
  trend_display,
  c(
    "Mean cognition: 2018 minus 2011",
    "Low-performance probability: 2018 minus 2011"
  )
)

calibration_display <- pooled[
  pooled$model == "primary_minus_uncalibrated" &
    pooled$dimension == "overall_trend" &
    pooled$group == "2018_minus_2011" &
    pooled$estimand %in% c(
      "standardized_mean_cognition_change",
      "standardized_low15_probability_change"
    ),
]
calibration_display <- calibration_display[
  match(
    c(
      "standardized_mean_cognition_change",
      "standardized_low15_probability_change"
    ),
    calibration_display$estimand
  ),
]
calibration_lines <- result_lines(
  calibration_display,
  c(
    "Calibration impact on mean-cognition change",
    "Calibration impact on low-performance change"
  )
)

inequality_display <- pooled[
  pooled$model == "primary" &
    pooled$group == "2018_minus_2011" &
    pooled$dimension %in% c(
      "education_gap_trend",
      "residence_gap_trend",
      "hukou_gap_trend",
      "consumption_rank_trend"
    ) &
    grepl(
      "mean_cognition_absolute_gap_change|mean_cognition_sii_change",
      pooled$estimand
    ),
]
inequality_display <- inequality_display[
  match(
    c(
      "education_gap_trend",
      "residence_gap_trend",
      "hukou_gap_trend",
      "consumption_rank_trend"
    ),
    inequality_display$dimension
  ),
]
inequality_lines <- result_lines(
  inequality_display,
  c(
    "Education mean gap change",
    "Rural-urban mean gap change",
    "Agricultural/non-agricultural hukou mean gap change",
    "Consumption-rank SII change"
  )
)

if (run_sensitivity) {
  scenario_labels <- c(
    main_mi = "Primary MI and survey weights",
    unweighted = "Unweighted",
    weight_trim_99 = "Weights trimmed at wave-specific 99th percentile",
    mnar_delta_0.25 = "MNAR delta -0.25 SD",
    mnar_delta_0.50 = "MNAR delta -0.50 SD",
    complete_predictor = "Complete predictor"
  )
  sensitivity_table <- c(
    "| Scenario | Mean cognition change | Low-performance change | Education mean-gap change |",
    "| --- | ---: | ---: | ---: |"
  )
  for (scenario in names(scenario_labels)) {
    if (scenario == "main_mi") {
      source <- pooled[
        pooled$model == "primary" &
          pooled$group == "2018_minus_2011",
      ]
    } else {
      source <- sensitivity_pooled[
        sensitivity_pooled$scenario == scenario &
          sensitivity_pooled$group == "2018_minus_2011",
      ]
    }
    mean_row <- source[
      source$dimension == "overall_trend" &
        source$estimand == "standardized_mean_cognition_change",
    ]
    low_row <- source[
      source$dimension == "overall_trend" &
        source$estimand == "standardized_low15_probability_change",
    ]
    education_row <- source[
      source$dimension == "education_gap_trend" &
        source$estimand ==
          "standardized_mean_cognition_absolute_gap_change",
    ]
    sensitivity_table <- c(
      sensitivity_table,
      sprintf(
        "| %s | %s (%s to %s) SD | %s (%s to %s) pp | %s (%s to %s) SD |",
        scenario_labels[[scenario]],
        fmt(mean_row$estimate, 3),
        fmt(mean_row$ci_lower, 3),
        fmt(mean_row$ci_upper, 3),
        fmt(low_row$estimate * 100, 2),
        fmt(low_row$ci_lower * 100, 2),
        fmt(low_row$ci_upper * 100, 2),
        fmt(education_row$estimate, 3),
        fmt(education_row$ci_lower, 3),
        fmt(education_row$ci_upper, 3)
      )
    )
  }

  original_missing_count <- rowSums(is.na(long[routine_items]))
  mnar_count_lines <- c(
    "| Year | N with >=3/5 cognition items missing | Percent |",
    "| ---: | ---: | ---: |"
  )
  for (wave_value in sort(unique(long$wave))) {
    index <- long$wave == wave_value
    count <- sum(original_missing_count[index] >= 3)
    mnar_count_lines <- c(
      mnar_count_lines,
      sprintf(
        "| %s | %s | %s%% |",
        unique(long$year[index]),
        count,
        fmt(100 * count / sum(index), 2)
      )
    )
  }

  sensitivity_report <- c(
    "# B1 Weight and Routine-Cognition MNAR Sensitivity",
    "",
    paste0("- Date: ", date_tag),
    paste0("- Multiple imputations for MI scenarios: ", m),
    paste0("- MICE iterations: ", maxit),
    "- Status: formal sensitivity analysis",
    "",
    "## Definitions",
    "",
    "- Unweighted: all analysis and 2018 calibration weights set to one.",
    "- Weight trimmed: weights capped at the wave-specific 99th percentile; the 2011 standardization distribution was recalculated with trimmed weights.",
    "- MNAR delta: among respondents with at least three of five routine cognition items originally missing, the HCAP-calibrated predicted cognition was shifted downward by 0.25 or 0.50 SD after MAR imputation.",
    "- Complete predictor: all five routine cognition items and primary calibration predictors observed; retained only as a selection-prone comparator.",
    "",
    "## MNAR Target Group",
    "",
    mnar_count_lines,
    "",
    "## 2011-2018 Robustness",
    "",
    sensitivity_table,
    "",
    "Negative mean-cognition and education-gap changes indicate deterioration and widening disadvantage, respectively. Positive low-performance changes indicate increasing burden.",
    "",
    "## Interpretation Guardrail",
    "",
    "The complete-predictor analysis changes the target population by selecting respondents able to complete all routine cognition items. A divergent complete-predictor result is evidence of selection sensitivity, not evidence that the MI result is invalid. MNAR scenarios are stress tests, not estimated missing-data mechanisms.",
    "",
    "## Remaining Sensitivities",
    "",
    "1. HCAP item-level MI and all-indicator-missing anchor bounds integrated with longitudinal estimation.",
    "2. Community-cluster bootstrap of anchor, transport, and trend estimation.",
    "3. Fixed/rolling cohort IPCW and mortality-bounded analyses."
  )
  writeLines(
    sensitivity_report,
    sensitivity_report_path,
    useBytes = TRUE
  )
  cat("wrote ", sensitivity_report_path, "\n", sep = "")
}

report <- c(
  "# B1 Frozen-Transport MI and Survey Inference",
  "",
  paste0("- Date: ", date_tag),
  paste0("- Multiple imputations: ", m),
  paste0("- MICE iterations: ", maxit),
  "- Survey variance: community-clustered, wave-specific respondent weights",
  "- Standardization: fixed 2011 weighted age/sex distribution",
  "- Status: executable main pipeline with separate sensitivity analyses",
  "",
  "## Main Estimates",
  "",
  display_lines,
  "",
  "Low-performance burden is estimated as the survey mean of the model-based probability that latent cognition is below the fixed 2018 threshold. It is not obtained by hard-thresholding shrunken predicted means.",
  "",
  "## Model Implementation",
  "",
  "The primary education-calibrated transport model and the required social-blind comparator were refit within every imputation using 2018 participants with a frozen FIML anchor. The same 2018 item scaling and transport coefficients were then applied to all waves within that imputation.",
  "",
  "## Formal 2011-2018 Change",
  "",
  trend_lines,
  "",
  "Wave contrasts were estimated inside a stacked survey design, retaining covariance across waves induced by repeated communities and participants. Confidence intervals use Rubin pooling with imputation degrees of freedom.",
  "",
  "## Calibration Impact on Change",
  "",
  calibration_lines,
  "",
  "The uncalibrated comparator is the five-item routine-cognition composite placed on a fixed weighted 2018 mean/SD metric, with fixed weighted 2018 lower-tail thresholds. Positive calibration-impact estimates indicate that the HCAP-calibrated change is more positive, or less negative, than the uncalibrated change.",
  "",
  "## Change in Absolute Inequality",
  "",
  inequality_lines,
  "",
  "Negative mean-gap changes indicate widening disadvantage for the first-named group. The consumption SII is coded poorest minus richest; its change is 2018 minus 2011.",
  "",
  "## Interpretation Boundary",
  "",
  "These are HCAP-calibrated cognitive-performance surveillance estimates, not clinical MCI or dementia prevalence. Education enters the primary model only as a measurement-calibration modifier. Education inequality conclusions must be accompanied by the social-blind calibration-impact sensitivity.",
  "",
  "## Additional Sensitivities",
  "",
  "1. HCAP item-level MI and all-indicator-missing bounds.",
  "2. Community-cluster bootstrap of the full calibration and trend pipeline.",
  "3. Fixed/rolling cohort IPCW and mortality-bounded analyses.",
  "4. 2020 extension only after harmonization audit."
)

writeLines(report, report_path, useBytes = TRUE)
cat("wrote ", report_path, "\n", sep = "")
