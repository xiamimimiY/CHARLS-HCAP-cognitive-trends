#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1) stop("Run this file with Rscript.")
script_path <- normalizePath(
  sub("^--file=", "", gsub("~\\+~", " ", script_argument)),
  mustWork = TRUE
)
study <- dirname(dirname(script_path))
table_dir <- file.path(study, "results", "tables")
date_tag <- "2026-06-28"

hcap_path <- file.path(
  study,
  "data",
  "derived",
  "b1_hcap_fiml_input_2026-06-19.csv"
)
score_path <- file.path(
  table_dir,
  "b1_hcap_all_missing_scores_2026-06-19.csv"
)
prevalence_path <- file.path(
  table_dir,
  "b1_hcap_all_missing_mnar_prevalence_2026-06-19.csv"
)
primary_path <- file.path(
  table_dir,
  "b1_frozen_transport_mi_pooled_estimands_2026-06-28.csv"
)
output_path <- file.path(
  table_dir,
  paste0("b1_all_missing_longitudinal_bounds_", date_tag, ".csv")
)
report_path <- file.path(
  table_dir,
  "all_missing_selection_bounds_report.md"
)

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

hcap <- read.csv(hcap_path, check.names = FALSE)
hcap <- hcap[!is.na(hcap$edu3), ]
hcap$all_indicators_missing <- rowSums(!is.na(hcap[indicators])) == 0
valid_weight <- is.finite(hcap$r4wtrespb) & hcap$r4wtrespb > 0
all_missing_weight_share <- sum(
  hcap$r4wtrespb[valid_weight & hcap$all_indicators_missing]
) / sum(hcap$r4wtrespb[valid_weight])

score_table <- read.csv(score_path, check.names = FALSE)
prevalence <- read.csv(prevalence_path, check.names = FALSE)
primary <- read.csv(primary_path, check.names = FALSE)

primary_value <- function(year, estimand) {
  target <- primary[
    primary$model == "primary" &
      primary$dimension == "overall" &
      primary$year == year &
      primary$estimand == estimand,
  ]
  if (nrow(target) != 1) {
    stop("Primary endpoint lookup did not return exactly one row.")
  }
  target$estimate
}

scenario_delta <- c(
  MAR = 0,
  delta_0.25_SD_lower = 0.25,
  delta_0.50_SD_lower = 0.50
)

rows <- list()
for (scenario in names(scenario_delta)) {
  shifted_all_missing_mean <- (
    score_table$all_missing_mar_mean - scenario_delta[[scenario]]
  )
  mean_shift <- all_missing_weight_share * shifted_all_missing_mean
  primary_2011 <- primary_value(2011, "standardized_mean_cognition")
  primary_2018 <- primary_value(2018, "standardized_mean_cognition")
  rows[[length(rows) + 1]] <- data.frame(
    outcome = "mean_cognition",
    scenario = scenario,
    lower_tail_pct = NA_real_,
    all_missing_weight_share = all_missing_weight_share,
    primary_2011 = primary_2011,
    primary_2018 = primary_2018,
    selection_shift_2018 = mean(mean_shift),
    selection_shift_between_sd = stats::sd(mean_shift),
    bounded_2018 = primary_2018 + mean(mean_shift),
    primary_change_2018_minus_2011 = primary_2018 - primary_2011,
    bounded_change_2018_minus_2011 = (
      primary_2018 - primary_2011 + mean(mean_shift)
    )
  )

  for (threshold in c(10, 15, 20)) {
    target <- prevalence[
      prevalence$scenario == scenario &
        prevalence$lower_tail_pct == threshold,
    ]
    if (nrow(target) != 1) {
      stop("All-missing prevalence lookup did not return exactly one row.")
    }
    estimand <- paste0("standardized_low", threshold, "_probability")
    primary_2011 <- primary_value(2011, estimand)
    primary_2018 <- primary_value(2018, estimand)
    selection_shift <- (
      target$expanded_sample_prevalence_pct_mean -
        target$any_observed_prevalence_pct_mean
    ) / 100
    rows[[length(rows) + 1]] <- data.frame(
      outcome = "low_cognitive_performance",
      scenario = scenario,
      lower_tail_pct = threshold,
      all_missing_weight_share = all_missing_weight_share,
      primary_2011 = primary_2011,
      primary_2018 = primary_2018,
      selection_shift_2018 = selection_shift,
      selection_shift_between_sd = NA_real_,
      bounded_2018 = primary_2018 + selection_shift,
      primary_change_2018_minus_2011 = primary_2018 - primary_2011,
      bounded_change_2018_minus_2011 = (
        primary_2018 - primary_2011 + selection_shift
      )
    )
  }
}
bound_table <- do.call(rbind, rows)

write.csv(bound_table, output_path, row.names = FALSE)

fmt <- function(x, digits = 2) {
  ifelse(
    is.na(x) | !is.finite(x),
    "NA",
    formatC(x, digits = digits, format = "f")
  )
}

main_rows <- bound_table[
  bound_table$outcome == "mean_cognition" |
    (
      bound_table$outcome == "low_cognitive_performance" &
        bound_table$lower_tail_pct == 15
    ),
]
report_lines <- c(
  "| Outcome | Scenario | 2018 selection shift | Bounded 2018 level | Bounded 2011-2018 change |",
  "| --- | --- | ---: | ---: | ---: |"
)
for (index in seq_len(nrow(main_rows))) {
  row <- main_rows[index, ]
  is_low <- row$outcome == "low_cognitive_performance"
  multiplier <- if (is_low) 100 else 1
  suffix <- if (is_low) " pp" else " SD"
  report_lines <- c(
    report_lines,
    sprintf(
      "| %s | %s | %s%s | %s%s | %s%s |",
      ifelse(is_low, "Low cognitive performance (15%)", "Mean cognition"),
      row$scenario,
      fmt(row$selection_shift_2018 * multiplier, 2),
      suffix,
      fmt(row$bounded_2018 * multiplier, 2),
      suffix,
      fmt(row$bounded_change_2018_minus_2011 * multiplier, 2),
      suffix
    )
  )
}

report <- c(
  "# B1 All-HCAP-Indicator-Missing Longitudinal Bounds",
  "",
  paste0("- Date: ", date_tag),
  "- Status: selection-bound integration; not a replacement for primary Rubin-pooled inference",
  paste0(
    "- Weighted 2018 share with all 12 HCAP indicators missing: ",
    fmt(100 * all_missing_weight_share, 1),
    "%"
  ),
  "",
  "## Method",
  "",
  paste(
    "The primary 2011 estimate and the formal transport model are held fixed.",
    "For 2018, the expanded-sample minus any-item-observed HCAP difference is",
    "added to the primary transported endpoint. This isolates the plausible",
    "selection effect of the group with no directly observed HCAP indicator",
    "without redefining the primary longitudinal metric."
  ),
  "",
  "For continuous cognition, the shift equals the weighted all-missing share",
  "multiplied by its MAR-imputed mean score, with additional -0.25 or -0.50 SD",
  "scenario shifts. For low cognitive performance, the shift is the difference",
  "between expanded-sample and any-item-observed burden at each fixed threshold.",
  "",
  "## Main 15% Results",
  "",
  report_lines,
  "",
  "## Interpretation",
  "",
  paste(
    "Every all-missing scenario strengthens rather than removes the estimated",
    "2011-2018 deterioration. Under MAR, the bounded mean-cognition decline is",
    "more negative and the low-performance increase is larger than in the",
    "primary analysis. Pessimistic delta scenarios further strengthen both",
    "patterns."
  ),
  "",
  paste(
    "These are additive selection bounds, not jointly pooled confidence",
    "intervals and not clinical impairment or dementia prevalence estimates.",
    "Their purpose is to show the direction and plausible magnitude of bias",
    "from excluding respondents with no observed HCAP test."
  )
)

writeLines(report, report_path, useBytes = TRUE)
cat("wrote ", report_path, "\n", sep = "")
