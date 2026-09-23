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
report_path <- file.path(
  table_dir,
  "hcap_item_mi_sensitivity_report.md"
)
date_tag <- "2026-06-28"

m <- as.integer(Sys.getenv("B1_M", unset = "50"))
maxit <- as.integer(Sys.getenv("B1_MAXIT", unset = "10"))
seed <- 3427L

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

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

fit_anchor <- function(data, missing = "listwise") {
  cfa(
    anchor_syntax,
    data = data,
    group = "edu3",
    group.equal = c("loadings", "intercepts"),
    group.partial = partial_paths,
    estimator = "MLR",
    missing = missing,
    std.lv = TRUE,
    meanstructure = TRUE
  )
}

hcap_all <- read.csv(
  hcap_path,
  check.names = FALSE,
  colClasses = c(ID = "character")
)
hcap <- hcap_all[
  !is.na(hcap_all$edu3) &
    rowSums(!is.na(hcap_all[hcap_indicators])) >= 1,
]
long <- read.csv(
  long_path,
  check.names = FALSE,
  colClasses = c(
    ID = "character",
    pid = "character",
    household_id = "character",
    community_id = "character"
  )
)

fiml_fit <- fit_anchor(hcap, missing = "fiml")
if (!isTRUE(lavInspect(fiml_fit, "converged"))) {
  stop("FIML reference anchor did not converge.")
}
fiml_raw <- assemble_general_scores(fiml_fit, nrow(hcap))
fiml_mean <- weighted_mean(fiml_raw, hcap$r4wtrespb)
fiml_sd <- weighted_sd(fiml_raw, hcap$r4wtrespb)
fiml_scores <- (fiml_raw - fiml_mean) / fiml_sd
fiml_cut15 <- weighted_quantile(fiml_scores, hcap$r4wtrespb, 0.15)

scaling <- do.call(rbind, lapply(routine_items, function(item) {
  hcap_name <- paste0("r4", item)
  data.frame(
    variable = item,
    mean_2018 = mean(hcap[[hcap_name]], na.rm = TRUE),
    sd_2018 = stats::sd(hcap[[hcap_name]], na.rm = TRUE)
  )
}))

age_sex_levels <- c(
  "60-69_male", "60-69_female",
  "70-79_male", "70-79_female",
  "80+_male", "80+_female"
)

prepare_long <- function(data) {
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
  result$age_decade <- (
    result$age - mean(hcap$r4agey, na.rm = TRUE)
  ) / 10
  result$education3 <- as.integer(as.character(result$education3))
  result$edu_some_schooling_below_middle <- as.numeric(
    result$education3 == 2
  )
  result$edu_middle_school_plus <- as.numeric(
    result$education3 == 3
  )
  result$education_label <- factor(
    result$education3,
    levels = c(1, 2, 3),
    labels = c(
      "no_formal",
      "some_schooling_below_middle",
      "middle_school_plus"
    )
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
  result
}

long <- prepare_long(long)
reference <- long[
  long$wave == 1 &
    is.finite(long$weight) &
    long$weight > 0,
]
reference_totals <- tapply(
  reference$weight,
  reference$age_sex_cell,
  sum,
  na.rm = TRUE
)
reference_totals <- reference_totals[age_sex_levels]
reference_totals[is.na(reference_totals)] <- 0

hcap_mi_columns <- c(
  hcap_indicators,
  "r4agey",
  "female",
  "edu3",
  "iqcode_mean",
  "blessed_score"
)
hcap_mi_data <- hcap[hcap_mi_columns]
hcap_mi_data$edu3 <- factor(hcap_mi_data$edu3, levels = c(1, 2, 3))
hcap_method <- make.method(hcap_mi_data)
hcap_method[] <- ""
hcap_imputed_targets <- c(
  hcap_indicators,
  "iqcode_mean",
  "blessed_score"
)
hcap_method[hcap_imputed_targets] <- "pmm"
hcap_predictor <- make.predictorMatrix(hcap_mi_data)
hcap_predictor[,] <- 0
hcap_fixed <- c(
  "r4agey",
  "female",
  "edu3"
)
for (target in hcap_imputed_targets) {
  hcap_predictor[
    target,
    setdiff(c(hcap_imputed_targets, hcap_fixed), target)
  ] <- 1
}

cat(
  "starting HCAP item-level MICE: m=",
  m,
  ", maxit=",
  maxit,
  "\n",
  sep = ""
)
hcap_imp <- mice(
  hcap_mi_data,
  m = m,
  maxit = maxit,
  method = hcap_method,
  predictorMatrix = hcap_predictor,
  seed = seed + 100L,
  printFlag = FALSE
)
cat("completed HCAP item-level MICE\n")
if (!is.null(hcap_imp$loggedEvents) && nrow(hcap_imp$loggedEvents) > 0) {
  write.csv(
    hcap_imp$loggedEvents,
    file.path(
      table_dir,
      paste0("b1_hcap_item_mi_longitudinal_logged_events_", date_tag, ".csv")
    ),
    row.names = FALSE
  )
}

routine_mi_variables <- c(
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
routine_imputations <- list()
for (wave_value in sort(unique(long$wave))) {
  wave_data <- long[long$wave == wave_value, ]
  mice_data <- wave_data[routine_mi_variables]
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
    predictor[
      target,
      setdiff(c(imputed_targets, fixed_predictors), target)
    ] <- 1
  }

  cat(
    "starting routine MICE for ",
    unique(wave_data$year),
    "\n",
    sep = ""
  )
  routine_imputations[[as.character(wave_value)]] <- mice(
    mice_data,
    m = m,
    maxit = maxit,
    method = method,
    predictorMatrix = predictor,
    seed = seed + wave_value,
    printFlag = FALSE
  )
  cat("completed routine MICE for ", unique(wave_data$year), "\n", sep = "")
}

primary_formula <- anchor_z ~ z_imrc + z_dlrc + z_ser7 +
  z_orient + z_draw + age_decade + I(age_decade^2) + female +
  edu_some_schooling_below_middle + edu_middle_school_plus +
  routine_composite:edu_some_schooling_below_middle +
  routine_composite:edu_middle_school_plus

contrast_row <- function(
  estimates,
  covariance,
  contrast,
  imputation,
  dimension,
  group,
  estimand
) {
  data.frame(
    imputation = imputation,
    dimension = dimension,
    group = group,
    estimand = estimand,
    estimate = unname(sum(contrast * estimates)),
    variance = unname(
      as.numeric(t(contrast) %*% covariance %*% contrast)
    )
  )
}

extract_wave_change <- function(
  design,
  variable,
  estimand,
  imputation
) {
  fit <- svyglm(
    as.formula(paste0(variable, " ~ 0 + wave_factor")),
    design = design,
    family = gaussian(),
    na.action = na.omit
  )
  estimates <- coef(fit)
  covariance <- vcov(fit)
  labels <- sub("^wave_factor", "", names(estimates))
  contrast <- rep(0, length(estimates))
  contrast[labels == "2018"] <- 1
  contrast[labels == "2011"] <- -1
  contrast_row(
    estimates,
    covariance,
    contrast,
    imputation,
    "overall_trend",
    "2018_minus_2011",
    paste0(estimand, "_change")
  )
}

extract_2018_level <- function(
  design,
  variable,
  estimand,
  imputation
) {
  fit <- svymean(
    as.formula(paste0("~", variable)),
    design,
    na.rm = TRUE
  )
  data.frame(
    imputation = imputation,
    dimension = "overall",
    group = "2018",
    estimand = estimand,
    estimate = unname(coef(fit)[1]),
    variance = unname(vcov(fit)[1, 1])
  )
}

extract_education_gap_change <- function(
  data,
  variable,
  estimand,
  imputation
) {
  target <- data[
    !is.na(data$education_label) &
      is.finite(data$weight) &
      data$weight > 0 &
      !is.na(data$community_id),
  ]
  target$wave_group <- interaction(
    target$wave_factor,
    target$education_label,
    sep = "__",
    drop = TRUE
  )
  design <- svydesign(
    ids = ~community_id,
    weights = ~weight,
    data = target,
    nest = TRUE
  )
  standardized <- svystandardize(
    design,
    by = ~age_sex_cell,
    over = ~wave_group,
    population = reference_totals
  )
  fit <- svyglm(
    as.formula(paste0(variable, " ~ 0 + wave_group")),
    design = standardized,
    family = gaussian(),
    na.action = na.omit
  )
  estimates <- coef(fit)
  covariance <- vcov(fit)
  labels <- sub("^wave_group", "", names(estimates))
  gap <- function(year) {
    contrast <- rep(0, length(estimates))
    contrast[
      labels == paste0(year, "__no_formal")
    ] <- 1
    contrast[
      labels == paste0(year, "__middle_school_plus")
    ] <- -1
    contrast
  }
  contrast_row(
    estimates,
    covariance,
    gap("2018") - gap("2011"),
    imputation,
    "education_gap_trend",
    "2018_minus_2011",
    paste0(estimand, "_absolute_gap_change")
  )
}

rubin_pool <- function(data) {
  qbar <- mean(data$estimate)
  ubar <- mean(data$variance)
  m_local <- nrow(data)
  between <- if (m_local > 1) stats::var(data$estimate) else 0
  added <- (1 + 1 / m_local) * between
  total <- ubar + added
  relative <- ifelse(ubar > 0, added / ubar, Inf)
  degrees_freedom <- ifelse(
    between <= 0,
    Inf,
    (m_local - 1) * (1 + 1 / relative)^2
  )
  critical <- ifelse(
    is.finite(degrees_freedom),
    stats::qt(0.975, degrees_freedom),
    stats::qnorm(0.975)
  )
  statistic <- qbar / sqrt(total)
  data.frame(
    m = m_local,
    estimate = qbar,
    se = sqrt(total),
    degrees_freedom = degrees_freedom,
    p_value = 2 * stats::pt(-abs(statistic), df = degrees_freedom),
    ci_lower = qbar - critical * sqrt(total),
    ci_upper = qbar + critical * sqrt(total),
    within_variance = ubar,
    between_variance = between,
    fraction_missing_information = ifelse(total > 0, added / total, 0)
  )
}

estimate_rows <- list()
diagnostic_rows <- list()

for (imputation_index in seq_len(m)) {
  hcap_completed <- complete(hcap_imp, imputation_index)
  hcap_completed$edu3 <- as.integer(
    as.character(hcap_completed$edu3)
  )
  if (anyNA(hcap_completed[hcap_indicators])) {
    stop(
      "HCAP indicators remained missing after imputation ",
      imputation_index
    )
  }
  anchor_fit <- fit_anchor(hcap_completed)
  if (!isTRUE(lavInspect(anchor_fit, "converged"))) {
    stop(
      "HCAP item-MI anchor failed at imputation ",
      imputation_index
    )
  }
  raw_scores <- assemble_general_scores(
    anchor_fit,
    nrow(hcap_completed)
  )
  anchor_mean <- weighted_mean(raw_scores, hcap$r4wtrespb)
  anchor_sd <- weighted_sd(raw_scores, hcap$r4wtrespb)
  anchor_scores <- (raw_scores - anchor_mean) / anchor_sd
  cut15 <- weighted_quantile(
    anchor_scores,
    hcap$r4wtrespb,
    0.15
  )

  fiml_flag <- fiml_scores <= fiml_cut15
  mi_flag <- anchor_scores <= cut15
  diagnostic_rows[[imputation_index]] <- data.frame(
    imputation = imputation_index,
    converged = lavInspect(anchor_fit, "converged"),
    cfi_robust = unname(fitMeasures(anchor_fit, "cfi.robust")),
    rmsea_robust = unname(fitMeasures(anchor_fit, "rmsea.robust")),
    srmr = unname(fitMeasures(anchor_fit, "srmr")),
    score_correlation = stats::cor(
      anchor_scores,
      fiml_scores,
      use = "complete.obs"
    ),
    threshold_agreement = weighted_mean(
      mi_flag == fiml_flag,
      hcap$r4wtrespb
    )
  )

  completed_waves <- lapply(
    sort(unique(long$wave)),
    function(wave_value) {
      base <- long[long$wave == wave_value, ]
      completed <- complete(
        routine_imputations[[as.character(wave_value)]],
        imputation_index
      )
      base[routine_mi_variables] <- completed[routine_mi_variables]
      prepare_long(base)
    }
  )
  completed_long <- do.call(rbind, completed_waves)
  anchor_map <- data.frame(
    ID = hcap$ID,
    anchor_z = anchor_scores
  )
  completed_long <- merge(
    completed_long,
    anchor_map,
    by = "ID",
    all.x = TRUE
  )
  completed_long$anchor_z[completed_long$wave != 4] <- NA_real_

  wave4 <- completed_long[
    completed_long$wave == 4 &
      !is.na(completed_long$anchor_z),
  ]
  transport_fit <- lm(
    primary_formula,
    data = wave4,
    weights = weight
  )
  residual_sd <- sqrt(weighted_mean(
    residuals(transport_fit)^2,
    model.weights(model.frame(transport_fit))
  ))
  completed_long$pred_primary <- predict(
    transport_fit,
    newdata = completed_long
  )
  completed_long$low15_primary <- pnorm(
    (
      cut15 - completed_long$pred_primary
    ) / residual_sd
  )

  target <- completed_long[
    is.finite(completed_long$weight) &
      completed_long$weight > 0 &
      !is.na(completed_long$community_id),
  ]
  design <- svydesign(
    ids = ~community_id,
    weights = ~weight,
    data = target,
    nest = TRUE
  )
  standardized <- svystandardize(
    design,
    by = ~age_sex_cell,
    over = ~wave_factor,
    population = reference_totals
  )
  estimate_rows[[length(estimate_rows) + 1]] <- extract_wave_change(
    standardized,
    "pred_primary",
    "standardized_mean_cognition",
    imputation_index
  )
  estimate_rows[[length(estimate_rows) + 1]] <- extract_wave_change(
    standardized,
    "low15_primary",
    "standardized_low15_probability",
    imputation_index
  )
  estimate_rows[[length(estimate_rows) + 1]] <-
    extract_education_gap_change(
      target,
      "pred_primary",
      "standardized_mean_cognition",
      imputation_index
    )
  estimate_rows[[length(estimate_rows) + 1]] <-
    extract_education_gap_change(
      target,
      "low15_primary",
      "standardized_low15_probability",
      imputation_index
    )

  wave4_design <- subset(standardized, wave_factor == "2018")
  estimate_rows[[length(estimate_rows) + 1]] <- extract_2018_level(
    wave4_design,
    "pred_primary",
    "standardized_mean_cognition",
    imputation_index
  )
  estimate_rows[[length(estimate_rows) + 1]] <- extract_2018_level(
    wave4_design,
    "low15_primary",
    "standardized_low15_probability",
    imputation_index
  )

  diagnostic_rows[[imputation_index]]$n_calibration <- nobs(
    transport_fit
  )
  diagnostic_rows[[imputation_index]]$transport_r_squared <- summary(
    transport_fit
  )$r.squared
  diagnostic_rows[[imputation_index]]$residual_sd <- residual_sd
  cat(
    "completed paired HCAP/routine imputation ",
    imputation_index,
    " of ",
    m,
    "\n",
    sep = ""
  )
}

estimate_long <- do.call(rbind, estimate_rows)
diagnostics <- do.call(rbind, diagnostic_rows)

pool_keys <- c("dimension", "group", "estimand")
pooled <- do.call(rbind, lapply(
  split(
    estimate_long,
    interaction(
      estimate_long[pool_keys],
      drop = TRUE,
      lex.order = TRUE
    )
  ),
  function(target) {
    data.frame(target[1, pool_keys], rubin_pool(target))
  }
))
row.names(pooled) <- NULL

primary_path <- file.path(
  table_dir,
  "b1_frozen_transport_mi_pooled_estimands_2026-06-28.csv"
)
primary_results <- read.csv(primary_path, check.names = FALSE)
primary_core <- primary_results[
  primary_results$model == "primary" &
    (
      (
        primary_results$group == "2018_minus_2011" &
          primary_results$dimension %in% c(
            "overall_trend",
            "education_gap_trend"
          ) &
          primary_results$estimand %in% pooled$estimand
      ) |
        (
          primary_results$year == 2018 &
            primary_results$dimension == "overall" &
            primary_results$estimand %in% c(
              "standardized_mean_cognition",
              "standardized_low15_probability"
            )
        )
    ),
  c(
    "dimension",
    "group",
    "estimand",
    "estimate",
    "ci_lower",
    "ci_upper"
  )
]
primary_core$group[
  primary_core$dimension == "overall" &
    primary_core$group == "overall"
] <- "2018"

comparison <- merge(
  primary_core,
  pooled,
  by = c("dimension", "group", "estimand"),
  suffixes = c("_fiml", "_item_mi"),
  all.y = TRUE
)
comparison$difference_item_mi_minus_fiml <- (
  comparison$estimate_item_mi - comparison$estimate_fiml
)

write.csv(
  diagnostics,
  file.path(
    table_dir,
    paste0("b1_hcap_item_mi_longitudinal_diagnostics_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  estimate_long,
  file.path(
    table_dir,
    paste0("b1_hcap_item_mi_longitudinal_long_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  pooled,
  file.path(
    table_dir,
    paste0("b1_hcap_item_mi_longitudinal_pooled_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  comparison,
  file.path(
    table_dir,
    paste0("b1_hcap_item_mi_vs_fiml_longitudinal_", date_tag, ".csv")
  ),
  row.names = FALSE
)

fmt <- function(x, digits = 3) {
  ifelse(
    is.na(x) | !is.finite(x),
    "NA",
    formatC(x, digits = digits, format = "f")
  )
}

comparison_lines <- c(
  "| Estimand | FIML anchor | Item-MI anchor | Difference |",
  "| --- | ---: | ---: | ---: |"
)
for (i in seq_len(nrow(comparison))) {
  row <- comparison[i, ]
  multiplier <- if (grepl("probability", row$estimand)) 100 else 1
  suffix <- if (multiplier == 100) " pp" else " SD"
  comparison_lines <- c(
    comparison_lines,
    sprintf(
      "| %s / %s | %s%s | %s (%s to %s)%s | %s%s |",
      row$dimension,
      row$estimand,
      fmt(row$estimate_fiml * multiplier),
      suffix,
      fmt(row$estimate_item_mi * multiplier),
      fmt(row$ci_lower_item_mi * multiplier),
      fmt(row$ci_upper_item_mi * multiplier),
      suffix,
      fmt(row$difference_item_mi_minus_fiml * multiplier),
      suffix
    )
  )
}

report <- c(
  "# B1 HCAP Item-MI Longitudinal Propagation",
  "",
  paste0("- Date: ", date_tag),
  paste0("- HCAP item-level imputations: ", m),
  paste0("- Routine-cognition imputations: ", m),
  paste0("- MICE iterations: ", maxit),
  "- Combination: independently generated HCAP and routine-cognition imputations paired by imputation index",
  "- Status: formal anchor-missingness sensitivity",
  "",
  "## Method",
  "",
  "Each HCAP item-level imputation was fit with the frozen education partial-invariance bifactor model. The resulting weighted 2018 anchor and 15% threshold were paired with one wave-specific routine-cognition imputation. The frozen education-calibrated transport model was then refit and transported across 2011-2018. Rubin pooling therefore propagates both HCAP item and routine-cognition missing-data uncertainty.",
  "",
  "## Diagnostics",
  "",
  paste0(
    "- All anchor models converged: ",
    ifelse(all(diagnostics$converged), "yes", "no")
  ),
  paste0(
    "- Mean FIML/item-MI score correlation: ",
    fmt(mean(diagnostics$score_correlation))
  ),
  paste0(
    "- Mean weighted 15% threshold agreement: ",
    fmt(100 * mean(diagnostics$threshold_agreement), 1),
    "%"
  ),
  paste0(
    "- Mean transport R2: ",
    fmt(mean(diagnostics$transport_r_squared))
  ),
  "",
  "## Longitudinal Comparison",
  "",
  comparison_lines,
  "",
  "## Interpretation",
  "",
  "The item-MI anchor sensitivity tests whether the longitudinal conclusions depend on FIML treatment of partially missing HCAP tests. It does not identify cognition for participants with all HCAP indicators missing; that group remains a separately labelled selection-bound analysis.",
  "",
  "## Decision Rule",
  "",
  "The FIML anchor remains primary if item-MI preserves the direction of the national trend and education-gap change, with no practically material shift in the principal estimates."
)

writeLines(report, report_path, useBytes = TRUE)
cat("wrote ", report_path, "\n", sep = "")
