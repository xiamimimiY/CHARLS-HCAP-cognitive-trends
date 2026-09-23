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
report_path <- file.path(table_dir, "cluster_bootstrap_report.md")
date_tag <- "2026-06-28"

bootstrap_reps <- as.integer(Sys.getenv("B1_BOOT_REPS", unset = "200"))
m <- as.integer(Sys.getenv("B1_M", unset = "10"))
maxit <- as.integer(Sys.getenv("B1_MAXIT", unset = "10"))
seed <- 3427L

if (
  any(!is.finite(c(bootstrap_reps, m, maxit))) ||
    any(c(bootstrap_reps, m, maxit) < 1)
) {
  stop("B1_BOOT_REPS, B1_M, and B1_MAXIT must be positive integers.")
}

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
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

fit_anchor <- function(data) {
  cfa(
    anchor_syntax,
    data = data,
    group = "edu3",
    group.equal = c("loadings", "intercepts"),
    group.partial = partial_paths,
    estimator = "MLR",
    missing = "fiml",
    std.lv = TRUE,
    meanstructure = TRUE
  )
}

age_sex_levels <- c(
  "60-69_male", "60-69_female",
  "70-79_male", "70-79_female",
  "80+_male", "80+_female"
)

prepare_long <- function(data, scaling, age_center) {
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
  result
}

standardized_mean <- function(data, variable, reference_proportions) {
  cell_means <- tapply(
    data[[variable]] * data$weight,
    data$age_sex_cell,
    sum,
    na.rm = TRUE
  ) / tapply(
    data$weight,
    data$age_sex_cell,
    sum,
    na.rm = TRUE
  )
  cell_means <- cell_means[age_sex_levels]
  if (any(!is.finite(cell_means[reference_proportions > 0]))) {
    return(NA_real_)
  }
  sum(cell_means * reference_proportions, na.rm = TRUE)
}

standardized_group_gap <- function(
  data,
  variable,
  reference_proportions,
  first = "no_formal",
  second = "middle_school_plus"
) {
  group_mean <- function(label) {
    target <- data[data$education_label == label, ]
    standardized_mean(target, variable, reference_proportions)
  }
  group_mean(first) - group_mean(second)
}

primary_formula <- anchor_z ~ z_imrc + z_dlrc + z_ser7 +
  z_orient + z_draw + age_decade + I(age_decade^2) + female +
  edu_some_schooling_below_middle + edu_middle_school_plus +
  routine_composite:edu_some_schooling_below_middle +
  routine_composite:edu_middle_school_plus

hcap <- read.csv(
  hcap_path,
  check.names = FALSE,
  colClasses = c(ID = "character")
)
hcap <- hcap[
  !is.na(hcap$edu3) &
    rowSums(!is.na(hcap[hcap_indicators])) >= 1,
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

wave4_community <- unique(
  long[long$wave == 4, c("ID", "community_id")]
)
hcap <- merge(hcap, wave4_community, by = "ID", all.x = TRUE)
hcap <- hcap[!is.na(hcap$community_id), ]

original_anchor_fit <- fit_anchor(hcap)
original_anchor_raw <- assemble_general_scores(
  original_anchor_fit,
  nrow(hcap)
)
original_anchor_mean <- weighted_mean(
  original_anchor_raw,
  hcap$r4wtrespb
)
original_anchor_sd <- weighted_sd(
  original_anchor_raw,
  hcap$r4wtrespb
)
original_anchor_z <- (
  original_anchor_raw - original_anchor_mean
) / original_anchor_sd

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

imputations <- list()
for (wave_value in sort(unique(long$wave))) {
  wave_data <- long[long$wave == wave_value, ]
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
  targets <- names(method)[method != ""]
  fixed <- c("age", "female")
  for (target in targets) {
    predictor[target, setdiff(c(targets, fixed), target)] <- 1
  }

  cat("starting MICE for ", unique(wave_data$year), "\n", sep = "")
  imputations[[as.character(wave_value)]] <- mice(
    mice_data,
    m = m,
    maxit = maxit,
    method = method,
    predictorMatrix = predictor,
    seed = seed + wave_value,
    printFlag = FALSE
  )
  cat("completed MICE for ", unique(wave_data$year), "\n", sep = "")
}

completed_sets <- lapply(seq_len(m), function(imputation_index) {
  completed_waves <- lapply(
    sort(unique(long$wave)),
    function(wave_value) {
      base <- long[long$wave == wave_value, ]
      completed <- complete(
        imputations[[as.character(wave_value)]],
        imputation_index
      )
      base[mi_variables] <- completed[mi_variables]
      base
    }
  )
  do.call(rbind, completed_waves)
})

communities <- sort(unique(long$community_id))
communities <- communities[
  !is.na(communities) &
    communities != "" &
    communities != "nan"
]

bootstrap_rows <- list()
failure_rows <- list()

for (bootstrap_index in seq_len(bootstrap_reps)) {
  selected <- sample(
    communities,
    length(communities),
    replace = TRUE
  )
  imputation_index <- ((bootstrap_index - 1) %% m) + 1
  completed <- completed_sets[[imputation_index]]

  long_parts <- vector("list", length(selected))
  hcap_parts <- vector("list", length(selected))
  for (draw_index in seq_along(selected)) {
    community <- selected[draw_index]
    long_part <- completed[completed$community_id == community, ]
    long_part$bootstrap_cluster <- rep.int(draw_index, nrow(long_part))
    long_parts[[draw_index]] <- long_part
    hcap_part <- hcap[hcap$community_id == community, ]
    hcap_part$bootstrap_cluster <- rep.int(draw_index, nrow(hcap_part))
    hcap_parts[[draw_index]] <- hcap_part
  }
  boot_long <- do.call(rbind, long_parts)
  boot_hcap <- do.call(rbind, hcap_parts)

  result <- tryCatch({
    anchor_fit <- fit_anchor(boot_hcap)
    if (!isTRUE(lavInspect(anchor_fit, "converged"))) {
      stop("anchor_nonconvergence")
    }
    anchor_raw <- assemble_general_scores(
      anchor_fit,
      nrow(boot_hcap)
    )
    original_match <- original_anchor_z[
      match(boot_hcap$ID, hcap$ID)
    ]
    orientation <- stats::cor(
      anchor_raw,
      original_match,
      use = "complete.obs"
    )
    if (is.finite(orientation) && orientation < 0) {
      anchor_raw <- -anchor_raw
    }
    anchor_mean <- weighted_mean(anchor_raw, boot_hcap$r4wtrespb)
    anchor_sd <- weighted_sd(anchor_raw, boot_hcap$r4wtrespb)
    boot_hcap$anchor_z <- (
      anchor_raw - anchor_mean
    ) / anchor_sd
    threshold15 <- weighted_quantile(
      boot_hcap$anchor_z,
      boot_hcap$r4wtrespb,
      0.15
    )

    scaling <- do.call(rbind, lapply(routine_items, function(item) {
      hcap_name <- paste0("r4", item)
      data.frame(
        variable = item,
        mean_2018 = mean(boot_hcap[[hcap_name]], na.rm = TRUE),
        sd_2018 = stats::sd(boot_hcap[[hcap_name]], na.rm = TRUE)
      )
    }))
    age_center <- mean(boot_hcap$r4agey, na.rm = TRUE)
    boot_long <- prepare_long(boot_long, scaling, age_center)

    anchor_map <- boot_hcap[
      , c("ID", "bootstrap_cluster", "anchor_z")
    ]
    boot_long <- merge(
      boot_long,
      anchor_map,
      by = c("ID", "bootstrap_cluster"),
      all.x = TRUE
    )
    boot_long$anchor_z[boot_long$wave != 4] <- NA_real_
    calibration <- boot_long[
      boot_long$wave == 4 &
        !is.na(boot_long$anchor_z) &
        is.finite(boot_long$weight) &
        boot_long$weight > 0,
    ]
    transport_fit <- lm(
      primary_formula,
      data = calibration,
      weights = weight
    )
    residual_sd <- sqrt(weighted_mean(
      residuals(transport_fit)^2,
      model.weights(model.frame(transport_fit))
    ))
    boot_long$pred_primary <- predict(
      transport_fit,
      newdata = boot_long
    )
    boot_long$low15_primary <- pnorm(
      (
        threshold15 - boot_long$pred_primary
      ) / residual_sd
    )

    wave1 <- boot_long[
      boot_long$wave == 1 &
        is.finite(boot_long$weight) &
        boot_long$weight > 0,
    ]
    reference_weights <- tapply(
      wave1$weight,
      wave1$age_sex_cell,
      sum,
      na.rm = TRUE
    )
    reference_weights <- reference_weights[age_sex_levels]
    reference_weights[is.na(reference_weights)] <- 0
    reference_proportions <- reference_weights / sum(reference_weights)

    wave_estimates <- lapply(c(1, 4), function(wave_value) {
      target <- boot_long[
        boot_long$wave == wave_value &
          is.finite(boot_long$weight) &
          boot_long$weight > 0,
      ]
      data.frame(
        wave = wave_value,
        mean_cognition = standardized_mean(
          target,
          "pred_primary",
          reference_proportions
        ),
        low15_probability = standardized_mean(
          target,
          "low15_primary",
          reference_proportions
        ),
        education_mean_gap = standardized_group_gap(
          target,
          "pred_primary",
          reference_proportions
        ),
        education_low15_gap = standardized_group_gap(
          target,
          "low15_primary",
          reference_proportions
        )
      )
    })
    wave_estimates <- do.call(rbind, wave_estimates)
    wave2011 <- wave_estimates[wave_estimates$wave == 1, ]
    wave2018 <- wave_estimates[wave_estimates$wave == 4, ]

    data.frame(
      bootstrap = bootstrap_index,
      imputation = imputation_index,
      converged = TRUE,
      n_communities = length(unique(selected)),
      n_calibration = nobs(transport_fit),
      anchor_cfi_robust = unname(
        fitMeasures(anchor_fit, "cfi.robust")
      ),
      transport_r_squared = summary(transport_fit)$r.squared,
      mean_cognition_change = (
        wave2018$mean_cognition - wave2011$mean_cognition
      ),
      low15_probability_change = (
        wave2018$low15_probability - wave2011$low15_probability
      ),
      education_mean_gap_change = (
        wave2018$education_mean_gap - wave2011$education_mean_gap
      ),
      education_low15_gap_change = (
        wave2018$education_low15_gap - wave2011$education_low15_gap
      )
    )
  }, error = function(error) {
    failure_rows[[length(failure_rows) + 1]] <<- data.frame(
      bootstrap = bootstrap_index,
      imputation = imputation_index,
      message = conditionMessage(error)
    )
    NULL
  })

  if (!is.null(result)) {
    bootstrap_rows[[length(bootstrap_rows) + 1]] <- result
  }
  if (
    bootstrap_index %% 10 == 0 ||
      bootstrap_index == bootstrap_reps
  ) {
    cat(
      "completed bootstrap ",
      bootstrap_index,
      " of ",
      bootstrap_reps,
      "\n",
      sep = ""
    )
  }
}

bootstrap_long <- do.call(rbind, bootstrap_rows)
if (is.null(bootstrap_long) || nrow(bootstrap_long) == 0) {
  stop("All bootstrap replicates failed.")
}

main_results <- read.csv(
  file.path(
    table_dir,
    "b1_frozen_transport_mi_pooled_estimands_2026-06-28.csv"
  ),
  check.names = FALSE
)
main_value <- function(dimension, estimand) {
  target <- main_results[
    main_results$model == "primary" &
      main_results$dimension == dimension &
      main_results$group == "2018_minus_2011" &
      main_results$estimand == estimand,
  ]
  target$estimate[1]
}

estimand_map <- data.frame(
  variable = c(
    "mean_cognition_change",
    "low15_probability_change",
    "education_mean_gap_change",
    "education_low15_gap_change"
  ),
  dimension = c(
    "overall_trend",
    "overall_trend",
    "education_gap_trend",
    "education_gap_trend"
  ),
  estimand = c(
    "standardized_mean_cognition_change",
    "standardized_low15_probability_change",
    "standardized_mean_cognition_absolute_gap_change",
    "standardized_low15_probability_absolute_gap_change"
  )
)

summary_rows <- lapply(seq_len(nrow(estimand_map)), function(index) {
  spec <- estimand_map[index, ]
  values <- bootstrap_long[[spec$variable]]
  main_estimate <- main_value(spec$dimension, spec$estimand)
  data.frame(
    variable = spec$variable,
    main_estimate = main_estimate,
    bootstrap_reps_requested = bootstrap_reps,
    bootstrap_reps_successful = sum(is.finite(values)),
    bootstrap_success_rate = mean(is.finite(values)),
    bootstrap_mean = mean(values, na.rm = TRUE),
    bootstrap_bias = mean(values, na.rm = TRUE) - main_estimate,
    bootstrap_se = stats::sd(values, na.rm = TRUE),
    percentile_ci_lower = unname(
      stats::quantile(values, 0.025, na.rm = TRUE)
    ),
    percentile_ci_upper = unname(
      stats::quantile(values, 0.975, na.rm = TRUE)
    )
  )
})
bootstrap_summary <- do.call(rbind, summary_rows)

half_cut <- floor(bootstrap_reps / 2)
stability_rows <- lapply(seq_len(nrow(estimand_map)), function(index) {
  spec <- estimand_map[index, ]
  do.call(rbind, lapply(
    list(
      first_half = seq_len(half_cut),
      second_half = seq.int(half_cut + 1, bootstrap_reps)
    ),
    function(replicates) {
      values <- bootstrap_long[
        bootstrap_long$bootstrap %in% replicates,
        spec$variable
      ]
      values <- values[is.finite(values)]
      data.frame(
        variable = spec$variable,
        half = if (max(replicates) <= half_cut) {
          "first_half"
        } else {
          "second_half"
        },
        n = length(values),
        mean = mean(values),
        percentile_ci_lower = unname(stats::quantile(values, 0.025)),
        percentile_ci_upper = unname(stats::quantile(values, 0.975))
      )
    }
  ))
})
bootstrap_stability <- do.call(rbind, stability_rows)

write.csv(
  bootstrap_long,
  file.path(
    table_dir,
    paste0("b1_community_cluster_bootstrap_long_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  bootstrap_stability,
  file.path(
    table_dir,
    paste0("b1_community_cluster_bootstrap_stability_", date_tag, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  bootstrap_summary,
  file.path(
    table_dir,
    paste0("b1_community_cluster_bootstrap_summary_", date_tag, ".csv")
  ),
  row.names = FALSE
)
failure_table <- if (length(failure_rows) > 0) {
  do.call(rbind, failure_rows)
} else {
  data.frame(
    bootstrap = integer(),
    imputation = integer(),
    message = character()
  )
}
write.csv(
  failure_table,
  file.path(
    table_dir,
    paste0("b1_community_cluster_bootstrap_failures_", date_tag, ".csv")
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

summary_lines <- c(
  "| Estimand | Valid replicates | Main estimate | Bootstrap mean | Bias | Bootstrap SE | Percentile 95% CI |",
  "| --- | ---: | ---: | ---: | ---: | ---: | --- |"
)
for (index in seq_len(nrow(bootstrap_summary))) {
  row <- bootstrap_summary[index, ]
  multiplier <- if (grepl("probability|low15", row$variable)) 100 else 1
  suffix <- if (multiplier == 100) " pp" else " SD"
  summary_lines <- c(
    summary_lines,
    sprintf(
      "| %s | %d/%d | %s%s | %s%s | %s%s | %s%s | %s to %s%s |",
      row$variable,
      row$bootstrap_reps_successful,
      row$bootstrap_reps_requested,
      fmt(row$main_estimate * multiplier),
      suffix,
      fmt(row$bootstrap_mean * multiplier),
      suffix,
      fmt(row$bootstrap_bias * multiplier),
      suffix,
      fmt(row$bootstrap_se * multiplier),
      suffix,
      fmt(row$percentile_ci_lower * multiplier),
      fmt(row$percentile_ci_upper * multiplier),
      suffix
    )
  )
}

report <- c(
  "# B1 Community-Cluster Bootstrap",
  "",
  paste0("- Date: ", date_tag),
  paste0("- Requested bootstrap replicates: ", bootstrap_reps),
  paste0("- Anchor/transport model successful replicates: ", nrow(bootstrap_long)),
  paste0("- Failed model replicates: ", nrow(failure_table)),
  paste0("- Routine-cognition imputations cycled across replicates: ", m),
  paste0("- MICE iterations: ", maxit),
  "- Resampling unit: community, with all waves and participants retained within each selected community draw",
  "- Refit per replicate: FIML partial-invariance HCAP anchor, 2018 threshold, routine-item scaling, education-calibrated transport model, and age/sex standardization",
  "- Status: full-pipeline cluster-bootstrap robustness analysis",
  "",
  "## Results",
  "",
  summary_lines,
  "",
  "## Diagnostics",
  "",
  paste0(
    "- Anchor convergence rate: ",
    fmt(100 * nrow(bootstrap_long) / bootstrap_reps, 1),
    "%"
  ),
  paste0(
    "- Median unique communities per replicate: ",
    fmt(stats::median(bootstrap_long$n_communities), 0)
  ),
  paste0(
    "- Median calibration N: ",
    fmt(stats::median(bootstrap_long$n_calibration), 0)
  ),
  paste0(
    "- Mean anchor robust CFI: ",
    fmt(mean(bootstrap_long$anchor_cfi_robust))
  ),
  "- First-half and second-half estimates are stored in the stability table to assess Monte Carlo drift.",
  "",
  "## Interpretation",
  "",
  "The percentile intervals quantify cluster-sampling and calibration-pipeline instability after resampling entire communities. Routine-cognition MI is represented by cycling pre-generated completed datasets across bootstrap replicates. Sparse education-by-age/sex cells can make an education-gap estimand unavailable in an otherwise converged replicate; valid counts are therefore reported separately for every estimand. This is a robustness analysis rather than a replacement for the primary Rubin-pooled survey inference."
)

writeLines(report, report_path, useBytes = TRUE)
cat("wrote ", report_path, "\n", sep = "")
