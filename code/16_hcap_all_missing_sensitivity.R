#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lavaan)
  library(mice)
})

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1) stop("Run this file with Rscript.")
script_path <- normalizePath(
  sub("^--file=", "", gsub("~\\+~", " ", script_argument)),
  mustWork = TRUE
)
study <- dirname(dirname(script_path))
input_path <- file.path(
  study,
  "data",
  "derived",
  "b1_hcap_fiml_input_2026-06-19.csv"
)
table_dir <- file.path(study, "results", "tables")
report_path <- file.path(table_dir, "hcap_all_missing_report.md")
date_tag <- "2026-06-19"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(3427)

indicators <- c(
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

model_syntax <- paste(
  paste0("general =~ ", paste(indicators, collapse = " + ")),
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

raw <- read.csv(input_path, check.names = FALSE)
raw$observed_indicator_count <- rowSums(!is.na(raw[indicators]))
raw$all_indicators_missing <- raw$observed_indicator_count == 0
analysis <- raw[!is.na(raw$edu3), ]

weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
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

profile_rows <- list()
for (missing_value in c(FALSE, TRUE)) {
  target <- analysis[analysis$all_indicators_missing == missing_value, ]
  profile_rows[[length(profile_rows) + 1]] <- data.frame(
    group = ifelse(missing_value, "all_indicators_missing", "any_indicator_observed"),
    n = nrow(target),
    weighted_age_mean = weighted_mean(target$r4agey, target$r4wtrespb),
    weighted_female_pct = weighted_mean(target$female, target$r4wtrespb) * 100,
    weighted_edu1_pct = weighted_mean(target$edu3 == 1, target$r4wtrespb) * 100,
    weighted_edu2_pct = weighted_mean(target$edu3 == 2, target$r4wtrespb) * 100,
    weighted_edu3_pct = weighted_mean(target$edu3 == 3, target$r4wtrespb) * 100,
    iqcode_available_pct = mean(!is.na(target$iqcode_mean)) * 100,
    weighted_iqcode_mean = weighted_mean(target$iqcode_mean, target$r4wtrespb),
    blessed_available_pct = mean(!is.na(target$blessed_score)) * 100,
    weighted_blessed_mean = weighted_mean(target$blessed_score, target$r4wtrespb)
  )
}
profile_table <- do.call(rbind, profile_rows)

mi_columns <- c(
  indicators,
  "r4agey",
  "female",
  "edu3",
  "iqcode_mean",
  "blessed_score"
)
mi_data <- analysis[mi_columns]
mi_data$edu3 <- factor(mi_data$edu3, levels = c(1, 2, 3))

method <- make.method(mi_data)
method[] <- ""
method[indicators] <- "pmm"
method[c("iqcode_mean", "blessed_score")] <- "pmm"

predictor <- make.predictorMatrix(mi_data)
predictor[,] <- 0
imputed_targets <- c(indicators, "iqcode_mean", "blessed_score")
fixed_predictors <- c("r4agey", "female", "edu3")
for (target in imputed_targets) {
  predictor[target, setdiff(c(imputed_targets, fixed_predictors), target)] <- 1
}

m <- 10
maxit <- 10
imp <- mice(
  mi_data,
  m = m,
  maxit = maxit,
  method = method,
  predictorMatrix = predictor,
  printFlag = FALSE,
  seed = 3428
)

fit_partial <- function(data) {
  cfa(
    model_syntax,
    data = data,
    group = "edu3",
    group.equal = c("loadings", "intercepts"),
    group.partial = partial_paths,
    estimator = "MLR",
    missing = "listwise",
    std.lv = TRUE,
    meanstructure = TRUE
  )
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

result_rows <- list()
education_rows <- list()
score_rows <- list()
fit_rows <- list()

for (i in seq_len(m)) {
  completed <- complete(imp, i)
  completed$edu3 <- as.integer(as.character(completed$edu3))
  fit <- fit_partial(completed)
  measures <- fitMeasures(fit, c("cfi.robust", "rmsea.robust", "srmr"))
  fit_rows[[i]] <- data.frame(
    imputation = i,
    converged = lavInspect(fit, "converged"),
    cfi_robust = unname(measures["cfi.robust"]),
    rmsea_robust = unname(measures["rmsea.robust"]),
    srmr = unname(measures["srmr"])
  )
  scores <- assemble_general_scores(fit, nrow(analysis))
  any_observed <- !analysis$all_indicators_missing
  all_missing <- analysis$all_indicators_missing
  reference_mean <- weighted_mean(
    scores[any_observed],
    analysis$r4wtrespb[any_observed]
  )
  reference_sd <- sqrt(weighted_mean(
    (scores[any_observed] - reference_mean)^2,
    analysis$r4wtrespb[any_observed]
  ))
  aligned_scores <- (scores - reference_mean) / reference_sd

  score_rows[[i]] <- data.frame(
    imputation = i,
    any_observed_mean = weighted_mean(
      aligned_scores[any_observed],
      analysis$r4wtrespb[any_observed]
    ),
    all_missing_mar_mean = weighted_mean(
      aligned_scores[all_missing],
      analysis$r4wtrespb[all_missing]
    )
  )

  for (delta in c(0, -0.25, -0.50)) {
    scenario_scores <- aligned_scores
    scenario_scores[all_missing] <- scenario_scores[all_missing] + delta
    scenario_label <- if (delta == 0) {
      "MAR"
    } else {
      paste0("delta_", format(abs(delta), nsmall = 2), "_SD_lower")
    }

    for (probability in c(0.10, 0.15, 0.20)) {
      cut <- weighted_quantile(
        aligned_scores[any_observed],
        analysis$r4wtrespb[any_observed],
        probability
      )
      impaired <- scenario_scores <= cut
      result_rows[[length(result_rows) + 1]] <- data.frame(
        imputation = i,
        scenario = scenario_label,
        lower_tail_pct = probability * 100,
        reference_cut = cut,
        any_observed_prevalence_pct = weighted_mean(
          impaired[any_observed],
          analysis$r4wtrespb[any_observed]
        ) * 100,
        all_missing_prevalence_pct = weighted_mean(
          impaired[all_missing],
          analysis$r4wtrespb[all_missing]
        ) * 100,
        expanded_sample_prevalence_pct = weighted_mean(
          impaired,
          analysis$r4wtrespb
        ) * 100
      )

      for (education_group in 1:3) {
        target <- analysis$edu3 == education_group
        education_rows[[length(education_rows) + 1]] <- data.frame(
          imputation = i,
          scenario = scenario_label,
          lower_tail_pct = probability * 100,
          education_group = education_group,
          prevalence_pct = weighted_mean(
            impaired[target],
            analysis$r4wtrespb[target]
          ) * 100
        )
      }
    }
  }
}

fit_table <- do.call(rbind, fit_rows)
score_table <- do.call(rbind, score_rows)
prevalence_long <- do.call(rbind, result_rows)
education_long <- do.call(rbind, education_rows)

summarize_columns <- function(data, group_columns, value_columns) {
  groups <- interaction(data[group_columns], drop = TRUE, lex.order = TRUE)
  split_data <- split(data, groups)
  do.call(rbind, lapply(split_data, function(target) {
    out <- target[1, group_columns, drop = FALSE]
    for (value in value_columns) {
      out[[paste0(value, "_mean")]] <- mean(target[[value]])
      out[[paste0(value, "_between_sd")]] <- sd(target[[value]])
      out[[paste0(value, "_min")]] <- min(target[[value]])
      out[[paste0(value, "_max")]] <- max(target[[value]])
    }
    out
  }))
}

prevalence_summary <- summarize_columns(
  prevalence_long,
  c("scenario", "lower_tail_pct"),
  c(
    "any_observed_prevalence_pct",
    "all_missing_prevalence_pct",
    "expanded_sample_prevalence_pct"
  )
)
education_summary <- summarize_columns(
  education_long,
  c("scenario", "lower_tail_pct", "education_group"),
  "prevalence_pct"
)

write.csv(
  profile_table,
  file.path(table_dir, paste0("b1_hcap_all_missing_profile_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  fit_table,
  file.path(table_dir, paste0("b1_hcap_all_missing_mi_fit_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  score_table,
  file.path(table_dir, paste0("b1_hcap_all_missing_scores_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  prevalence_summary,
  file.path(table_dir, paste0("b1_hcap_all_missing_mnar_prevalence_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  education_summary,
  file.path(table_dir, paste0("b1_hcap_all_missing_mnar_education_", date_tag, ".csv")),
  row.names = FALSE
)

fmt <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

profile_lines <- c(
  "| Group | N | Age | Female | No formal education | IQCODE | Informant score |",
  "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"
)
for (i in seq_len(nrow(profile_table))) {
  r <- profile_table[i, ]
  profile_lines <- c(
    profile_lines,
    sprintf(
      "| %s | %s | %s | %s%% | %s%% | %s | %s |",
      r$group,
      r$n,
      fmt(r$weighted_age_mean, 1),
      fmt(r$weighted_female_pct, 1),
      fmt(r$weighted_edu1_pct, 1),
      fmt(r$weighted_iqcode_mean, 2),
      fmt(r$weighted_blessed_mean, 2)
    )
  )
}

prevalence_lines <- c(
  "| Scenario | Threshold | Any-item sample | All-missing group | Expanded sample |",
  "| --- | ---: | ---: | ---: | ---: |"
)
for (i in seq_len(nrow(prevalence_summary))) {
  r <- prevalence_summary[i, ]
  prevalence_lines <- c(
    prevalence_lines,
    sprintf(
      "| %s | %s%% | %s%% | %s%% | %s%% |",
      r$scenario,
      fmt(r$lower_tail_pct, 0),
      fmt(r$any_observed_prevalence_pct_mean, 1),
      fmt(r$all_missing_prevalence_pct_mean, 1),
      fmt(r$expanded_sample_prevalence_pct_mean, 1)
    )
  )
}

mar_15 <- prevalence_summary[
  prevalence_summary$scenario == "MAR" &
    prevalence_summary$lower_tail_pct == 15,
]
delta50_15 <- prevalence_summary[
  prevalence_summary$scenario == "delta_0.50_SD_lower" &
    prevalence_summary$lower_tail_pct == 15,
]

report <- c(
  "# B1 HCAP All-Indicator-Missing and MNAR Sensitivity",
  "",
  paste0("- Date: ", date_tag),
  paste0("- Method: expanded-sample MICE-PMM, ", m, " imputations, ", maxit, " iterations"),
  "- Status: expanded-sample/MNAR sensitivity; not a clinical prevalence estimate",
  "",
  "## Why This Group Requires Separate Treatment",
  "",
  profile_lines,
  "",
  paste(
    "Participants with all 12 indicators missing are older, less educated, and",
    "have substantially worse informant-reported cognitive and functional",
    "change. Treating them as ordinary complete-case exclusions is therefore",
    "not defensible."
  ),
  "",
  "## Expanded-Sample Burden",
  "",
  paste(
    "Thresholds are defined within the any-item-observed sample for each",
    "imputation, then transported unchanged to the all-missing group and full",
    "expanded sample. Delta scenarios lower only the all-missing general-factor",
    "scores after MAR imputation."
  ),
  "",
  prevalence_lines,
  "",
  "## Key 15% Threshold Result",
  "",
  sprintf(
    "Under MAR, the all-missing group had an estimated low-performance proportion of %s%% and inclusion raised the expanded-sample burden to %s%%.",
    fmt(mar_15$all_missing_prevalence_pct_mean, 1),
    fmt(mar_15$expanded_sample_prevalence_pct_mean, 1)
  ),
  sprintf(
    "With a -0.50 SD delta, the all-missing low-performance proportion was %s%% and expanded-sample burden was %s%%.",
    fmt(delta50_15$all_missing_prevalence_pct_mean, 1),
    fmt(delta50_15$expanded_sample_prevalence_pct_mean, 1)
  ),
  "",
  "## Decision",
  "",
  paste(
    "The primary HCAP anchor remains the any-item-observed FIML sample because",
    "the all-missing respondents have no direct item information. However, the",
    "expanded MAR and delta estimates must be reported prominently as selection",
    "bounds. Excluding the all-missing group without such bounds would likely",
    "underestimate population low-cognitive-performance burden."
  ),
  "",
  "## Interpretation Boundary",
  "",
  paste(
    "The delta values are sensitivity parameters, not empirically identified",
    "corrections. Expanded estimates depend strongly on auxiliary-variable",
    "transport and must not be presented as validated dementia prevalence."
  )
)

writeLines(report, report_path, useBytes = TRUE)
cat("wrote", report_path, "\n")
