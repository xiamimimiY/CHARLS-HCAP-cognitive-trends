#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lavaan)
  library(mice)
})

options(mc.cores = 1L)

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1) stop("Run this file with Rscript.")
script_path <- normalizePath(
  sub("^--file=", "", gsub("~\\+~", " ", script_argument)),
  mustWork = TRUE
)
study <- dirname(dirname(script_path))
hcap_path <- file.path(
  study,
  "data",
  "derived",
  "b1_hcap_fiml_input_2026-06-19.csv"
)
long_path <- file.path(
  study,
  "data",
  "derived",
  "b1_frozen_transport_mi_input_2026-06-28.csv"
)
table_dir <- file.path(study, "results", "tables")
tool_dir <- file.path(study, "scoring_tool")
report_path <- file.path(table_dir, "scoring_model_validation_report.md")
date_tag <- "2026-06-28"

m <- as.integer(Sys.getenv("B1_M", unset = "50"))
maxit <- as.integer(Sys.getenv("B1_MAXIT", unset = "10"))
folds <- as.integer(Sys.getenv("B1_FOLDS", unset = "10"))
seed <- 3427L
if (
  any(!is.finite(c(m, maxit, folds))) ||
    any(c(m, maxit, folds) < 2)
) {
  stop("B1_M, B1_MAXIT, and B1_FOLDS must be integers of at least 2.")
}

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tool_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

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

anchor_syntax <- paste(
  paste0("general =~ ", paste(hcap_indicators, collapse = " + ")),
  paste0(
    "memory_specific =~ ",
    paste(memory_indicators, collapse = " + ")
  ),
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
assemble_general_scores <- function(fit, n_rows) {
  scores <- lavPredict(fit, type = "lv")
  case_idx <- lavInspect(fit, "case.idx")
  result <- rep(NA_real_, n_rows)
  for (group_index in seq_along(scores)) {
    result[case_idx[[group_index]]] <- scores[[group_index]][, "general"]
  }
  result
}

hcap <- read.csv(
  hcap_path,
  check.names = FALSE,
  colClasses = c(ID = "character")
)
hcap <- hcap[
  !is.na(hcap$edu3) &
    rowSums(!is.na(hcap[hcap_indicators])) >= 1,
]
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
  stop("FIML HCAP anchor did not converge.")
}
anchor_raw <- assemble_general_scores(anchor_fit, nrow(hcap))
anchor_mean <- weighted_mean(anchor_raw, hcap$r4wtrespb)
anchor_sd <- weighted_sd(anchor_raw, hcap$r4wtrespb)
hcap$anchor_z <- (anchor_raw - anchor_mean) / anchor_sd
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
age_center <- mean(hcap$r4agey, na.rm = TRUE)

long <- read.csv(
  long_path,
  check.names = FALSE,
  colClasses = c(ID = "character", community_id = "character")
)
wave4 <- long[long$wave == 4, ]
wave4 <- merge(
  wave4,
  hcap[, c("ID", "anchor_z")],
  by = "ID",
  all.x = TRUE
)

mi_variables <- c(
  routine_items,
  "age",
  "female",
  "education3",
  "rural",
  "hukou",
  "adl",
  "iadl",
  "log_consumption"
)
mice_data <- wave4[mi_variables]
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
targets <- names(method)[method != ""]
fixed <- c("age", "female")
for (target in targets) {
  predictor[target, setdiff(c(targets, fixed), target)] <- 1
}

cat(
  "starting scoring-tool MICE: m=",
  m,
  ", maxit=",
  maxit,
  "\n",
  sep = ""
)
imputation <- mice(
  mice_data,
  m = m,
  maxit = maxit,
  method = method,
  predictorMatrix = predictor,
  seed = seed + 4L,
  printFlag = FALSE
)
cat("completed scoring-tool MICE\n")

prepare_data <- function(data) {
  result <- data
  for (item in routine_items) {
    scaling_row <- scaling[scaling$variable == item, ]
    result[[paste0("z_", item)]] <- (
      result[[item]] - scaling_row$mean_2018
    ) / scaling_row$sd_2018
  }
  result$routine_composite <- rowMeans(
    result[paste0("z_", routine_items)]
  )
  result$age_decade <- (result$age - age_center) / 10
  result$education3 <- as.integer(as.character(result$education3))
  result$edu_some_schooling_below_middle <- as.numeric(
    result$education3 == 2
  )
  result$edu_middle_school_plus <- as.numeric(
    result$education3 == 3
  )
  result
}

primary_formula <- anchor_z ~ z_imrc + z_dlrc + z_ser7 +
  z_orient + z_draw + age_decade + I(age_decade^2) + female +
  edu_some_schooling_below_middle + edu_middle_school_plus +
  routine_composite:edu_some_schooling_below_middle +
  routine_composite:edu_middle_school_plus

communities <- sort(unique(wave4$community_id[!is.na(wave4$community_id)]))
set.seed(seed + 26L)
communities <- sample(communities)
community_fold <- setNames(
  rep(seq_len(folds), length.out = length(communities)),
  communities
)

coefficient_rows <- list()
fit_rows <- list()
validation_rows <- list()

validation_metrics <- function(observed, predicted, weight, group, imputation_index) {
  keep <- is.finite(observed) &
    is.finite(predicted) &
    is.finite(weight) &
    weight > 0
  observed <- observed[keep]
  predicted <- predicted[keep]
  weight <- weight[keep]
  error <- predicted - observed
  calibration <- lm(
    observed ~ predicted,
    weights = weight
  )
  data.frame(
    imputation = imputation_index,
    group = group,
    n = length(observed),
    weighted_mean_error = weighted_mean(error, weight),
    weighted_mae = weighted_mean(abs(error), weight),
    weighted_rmse = sqrt(weighted_mean(error^2, weight)),
    weighted_r_squared = 1 - (
      weighted_mean(error^2, weight) /
        weighted_mean(
          (observed - weighted_mean(observed, weight))^2,
          weight
        )
    ),
    calibration_intercept = unname(coef(calibration)["(Intercept)"]),
    calibration_slope = unname(coef(calibration)["predicted"])
  )
}

for (imputation_index in seq_len(m)) {
  completed <- wave4
  completed[mi_variables] <- complete(imputation, imputation_index)
  completed <- prepare_data(completed)
  calibration <- completed[
    !is.na(completed$anchor_z) &
      is.finite(completed$weight) &
      completed$weight > 0 &
      !is.na(completed$community_id),
  ]
  fit <- lm(primary_formula, data = calibration, weights = weight)
  coefficient_rows[[imputation_index]] <- data.frame(
    imputation = imputation_index,
    term = names(coef(fit)),
    estimate = unname(coef(fit)),
    within_variance = diag(vcov(fit))
  )
  fit_rows[[imputation_index]] <- data.frame(
    imputation = imputation_index,
    n_calibration = nobs(fit),
    residual_sd = sqrt(weighted_mean(
      residuals(fit)^2,
      model.weights(model.frame(fit))
    )),
    r_squared = summary(fit)$r.squared
  )

  fold_prediction <- rep(NA_real_, nrow(calibration))
  fold_id <- unname(community_fold[calibration$community_id])
  for (fold in seq_len(folds)) {
    train <- calibration[fold_id != fold, ]
    test <- calibration[fold_id == fold, ]
    fold_fit <- lm(primary_formula, data = train, weights = weight)
    fold_prediction[fold_id == fold] <- predict(fold_fit, newdata = test)
  }
  validation_rows[[length(validation_rows) + 1]] <-
    validation_metrics(
      calibration$anchor_z,
      fold_prediction,
      calibration$weight,
      "overall",
      imputation_index
    )
  education_group_labels <- c(
    "no_formal_schooling",
    "some_schooling_below_middle",
    "middle_school_plus"
  )
  for (education_group in 1:3) {
    target <- calibration$education3 == education_group
    validation_rows[[length(validation_rows) + 1]] <-
      validation_metrics(
        calibration$anchor_z[target],
        fold_prediction[target],
        calibration$weight[target],
        education_group_labels[education_group],
        imputation_index
      )
  }
  cat("completed scoring imputation ", imputation_index, " of ", m, "\n", sep = "")
}

coefficient_long <- do.call(rbind, coefficient_rows)
fit_long <- do.call(rbind, fit_rows)
validation_long <- do.call(rbind, validation_rows)

pool_coefficient <- function(data) {
  qbar <- mean(data$estimate)
  ubar <- mean(data$within_variance)
  between <- stats::var(data$estimate)
  total <- ubar + (1 + 1 / nrow(data)) * between
  data.frame(
    estimate = qbar,
    se = sqrt(total),
    ci_lower = qbar - stats::qnorm(0.975) * sqrt(total),
    ci_upper = qbar + stats::qnorm(0.975) * sqrt(total),
    within_variance = ubar,
    between_variance = between
  )
}
coefficient_pooled <- do.call(rbind, lapply(
  split(coefficient_long, coefficient_long$term),
  function(target) {
    data.frame(
      term = target$term[1],
      pool_coefficient(target)
    )
  }
))
row.names(coefficient_pooled) <- NULL

validation_summary <- do.call(rbind, lapply(
  split(validation_long, validation_long$group),
  function(target) {
    values <- setdiff(
      names(target),
      c("imputation", "group")
    )
    out <- data.frame(group = target$group[1])
    for (value in values) {
      out[[paste0(value, "_mean")]] <- mean(target[[value]])
      out[[paste0(value, "_between_sd")]] <- stats::sd(target[[value]])
    }
    out
  }
))
row.names(validation_summary) <- NULL

fit_summary <- data.frame(
  m = m,
  n_calibration = round(mean(fit_long$n_calibration)),
  residual_sd = mean(fit_long$residual_sd),
  residual_sd_between = stats::sd(fit_long$residual_sd),
  r_squared = mean(fit_long$r_squared),
  r_squared_between = stats::sd(fit_long$r_squared)
)

constants <- data.frame(
  key = c(
    paste0("mean_", routine_items),
    paste0("sd_", routine_items),
    "age_center_years",
    "residual_sd",
    "threshold_low10",
    "threshold_low15",
    "threshold_low20",
    "development_min_age",
    "development_max_age"
  ),
  value = c(
    scaling$mean_2018,
    scaling$sd_2018,
    age_center,
    fit_summary$residual_sd,
    thresholds[c("10", "15", "20")],
    60,
    max(wave4$age, na.rm = TRUE)
  ),
  description = c(
    paste("2018 mean for", routine_items),
    paste("2018 SD for", routine_items),
    "Mean age in the HCAP anchor sample",
    "Mean weighted residual SD across 50 imputations",
    "Fixed weighted 2018 lower 10% HCAP threshold",
    "Fixed weighted 2018 lower 15% HCAP threshold",
    "Fixed weighted 2018 lower 20% HCAP threshold",
    "Minimum development age",
    "Maximum observed development age"
  )
)

write.csv(
  coefficient_pooled,
  file.path(tool_dir, "b1_hcap_scoring_coefficients_v1.csv"),
  row.names = FALSE
)
write.csv(
  constants,
  file.path(tool_dir, "b1_hcap_scoring_constants_v1.csv"),
  row.names = FALSE
)
write.csv(
  coefficient_long,
  file.path(
    table_dir,
    paste0("b1_scoring_tool_coefficients_long_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  fit_summary,
  file.path(
    table_dir,
    paste0("b1_scoring_tool_fit_summary_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  validation_long,
  file.path(
    table_dir,
    paste0("b1_scoring_tool_validation_long_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  validation_summary,
  file.path(
    table_dir,
    paste0("b1_scoring_tool_validation_summary_", date_tag, ".csv")
  ),
  row.names = FALSE
)

source(file.path(tool_dir, "b1_hcap_score.R"))
test_input <- data.frame(
  imrc = c(4, NA, 4),
  dlrc = c(3, 3, 3),
  ser7 = c(4, 4, 4),
  orient = c(4, 4, 4),
  draw = c(1, 1, 1),
  age = c(68, 68, 58),
  female = c(1, 1, 1),
  education3 = c(2, 2, 2)
)
test_scored <- b1_hcap_score(test_input, strict = FALSE)
if (
  !is.finite(test_scored$hcap_calibrated_cognition_z[1]) ||
    !is.na(test_scored$hcap_calibrated_cognition_z[2]) ||
    isTRUE(test_scored$within_frozen_scope[3])
) {
  stop("Scoring-tool unit checks failed.")
}

fmt <- function(x, digits = 3) {
  ifelse(
    is.na(x) | !is.finite(x),
    "NA",
    formatC(x, digits = digits, format = "f")
  )
}
validation_lines <- c(
  "| Group | N | Mean error | RMSE | R2 | Calibration intercept | Calibration slope |",
  "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"
)
validation_group_display <- c(
  no_formal_schooling = "no formal schooling",
  some_schooling_below_middle = "some schooling below middle school",
  middle_school_plus = "middle school or higher",
  overall = "overall"
)
for (index in seq_len(nrow(validation_summary))) {
  row <- validation_summary[index, ]
  validation_lines <- c(
    validation_lines,
    sprintf(
      "| %s | %s | %s | %s | %s | %s | %s |",
      validation_group_display[[row$group]],
      fmt(row$n_mean, 0),
      fmt(row$weighted_mean_error_mean),
      fmt(row$weighted_rmse_mean),
      fmt(row$weighted_r_squared_mean),
      fmt(row$calibration_intercept_mean),
      fmt(row$calibration_slope_mean)
    )
  )
}

report <- c(
  "# HCAP-Anchored Scoring Model",
  "",
  paste0("- Date: ", date_tag),
  "- Version: 1.0 research-use calibration model",
  paste0("- Routine-cognition imputations: ", m),
  paste0("- Community-separated validation folds: ", folds),
  "- Status: deterministic coefficients, constants, thresholds, R function, and internal validation",
  "",
  "## Frozen Tool",
  "",
  paste(
    "- Input: five CHARLS routine cognition components, age, sex, and three-level",
    "education (`1` no formal schooling, `2` some schooling below middle school,",
    "`3` middle school or higher)."
  ),
  "- Primary output: continuous HCAP-calibrated cognition in weighted 2018 SD units.",
  "- Secondary outputs: model-based probabilities below fixed 2018 lower 10%, 15%, and 20% thresholds.",
  paste0("- Calibration N: ", fit_summary$n_calibration),
  paste0("- Mean residual SD: ", fmt(fit_summary$residual_sd)),
  paste0("- Mean apparent R2 across imputations: ", fmt(fit_summary$r_squared)),
  "",
  "## Community-Separated Validation",
  "",
  validation_lines,
  "",
  "## Deployment Decision",
  "",
  paste(
    "The published research-use tool uses Rubin-averaged transport",
    paste0(
      "coefficients from the ",
      m,
      "-imputation 2018 calibration pipeline. Standardizing"
    ),
    "constants, age centering, residual SD, and thresholds are fixed and versioned.",
    "The deterministic score requires complete inputs; cohort studies with",
    "missing items must use multiple imputation."
  ),
  "",
  "## Interpretation Boundary",
  "",
  paste(
    "This is a population-research calibration tool, not a clinical diagnostic",
    "instrument. Threshold outputs are low-cognitive-performance probabilities,",
    "not MCI, dementia, or Alzheimer disease probabilities. Independent external",
    "validation is required before use outside CHARLS or closely harmonized",
    "research settings."
  )
)
writeLines(report, report_path, useBytes = TRUE)
cat("wrote ", report_path, "\n", sep = "")
