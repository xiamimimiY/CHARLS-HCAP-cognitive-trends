#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lavaan)
})

# Restricted desktop runtimes may not expose a detectable CPU count.
options(mc.cores = 1L)

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
report_path <- file.path(table_dir, "education_invariance_report.md")
date_tag <- "2026-06-19"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(input_path, check.names = FALSE)

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

d <- d[
  !is.na(d$edu3) & rowSums(!is.na(d[indicators])) >= 1,
]
d$model_weight <- d$r4wtrespb / mean(d$r4wtrespb[d$r4wtrespb > 0], na.rm = TRUE)

model_syntax <- paste(
  paste0("general =~ ", paste(indicators, collapse = " + ")),
  paste0("memory_specific =~ ", paste(memory_indicators, collapse = " + ")),
  "general ~~ 0*memory_specific",
  sep = "\n"
)

loading_path <- function(factor, indicator) {
  paste0(factor, "=~", indicator)
}

intercept_path <- function(indicator) {
  paste0(indicator, "~1")
}

# Staged releases are ordered from the score-test results. Releasing
# a lavaan group.partial parameter frees it across all groups, which is more
# conservative than tuning a single pairwise contrast.
metric_core <- c(
  loading_path("general", "z_test_orientation"),
  loading_path("general", "z_test_animal_fluency"),
  loading_path("memory_specific", "z_test_ten_delayed"),
  loading_path("general", "z_test_serial7"),
  loading_path("general", "z_test_recognition")
)

metric_expanded <- unique(c(
  metric_core,
  loading_path("general", "z_test_ten_delayed"),
  loading_path("general", "z_test_ten_learning_trial1"),
  loading_path("general", "z_test_language"),
  loading_path("memory_specific", "z_test_ten_learning_trial3"),
  loading_path("general", "z_test_short_learning_trial1")
))

metric_all_screened <- unique(c(
  metric_expanded,
  loading_path("general", "z_test_ten_learning_trial3"),
  loading_path("memory_specific", "z_test_ten_learning_trial2"),
  loading_path("general", "z_test_visuospatial"),
  loading_path("memory_specific", "z_test_ten_learning_trial1"),
  loading_path("memory_specific", "z_test_recognition")
))

scalar_core <- c(
  intercept_path("z_test_language"),
  intercept_path("z_test_orientation"),
  intercept_path("z_test_recognition"),
  intercept_path("z_test_animal_fluency"),
  intercept_path("z_test_serial7")
)

scalar_all_screened <- unique(c(
  scalar_core,
  intercept_path("z_test_short_learning_trial1"),
  intercept_path("z_test_visuospatial"),
  intercept_path("z_test_short_delayed"),
  intercept_path("z_test_ten_learning_trial1")
))

fit_group_model <- function(
  equal = character(0),
  partial = character(0),
  weighted = FALSE
) {
  model_data <- d
  args <- list(
    model = model_syntax,
    data = model_data,
    group = "edu3",
    group.equal = equal,
    group.partial = partial,
    estimator = "MLR",
    missing = "fiml",
    std.lv = TRUE,
    meanstructure = TRUE
  )
  if (weighted) {
    model_data <- model_data[
      !is.na(model_data$model_weight) & model_data$model_weight > 0,
    ]
    args$data <- model_data
    args$sampling.weights <- "model_weight"
    args$sampling.weights.normalization <- "total"
  }
  tryCatch(do.call(cfa, args), error = function(e) e)
}

metric_candidates <- list(
  configural = fit_group_model(),
  full_metric = fit_group_model("loadings"),
  partial_metric_core = fit_group_model("loadings", metric_core),
  partial_metric_expanded = fit_group_model("loadings", metric_expanded),
  partial_metric_all_screened = fit_group_model(
    "loadings",
    metric_all_screened
  )
)

fit_row <- function(name, fit, reference = NULL, n_freed = 0) {
  if (inherits(fit, "error")) {
    return(data.frame(
      model = name,
      converged = FALSE,
      admissible = FALSE,
      n_used = NA,
      n_freed = n_freed,
      cfi_robust = NA,
      rmsea_robust = NA,
      srmr = NA,
      delta_cfi_vs_reference = NA,
      delta_rmsea_vs_reference = NA,
      error = conditionMessage(fit)
    ))
  }
  converged <- isTRUE(lavInspect(fit, "converged"))
  if (!converged) {
    return(data.frame(
      model = name,
      converged = FALSE,
      admissible = FALSE,
      n_used = sum(lavInspect(fit, "nobs")),
      n_freed = n_freed,
      cfi_robust = NA,
      rmsea_robust = NA,
      srmr = NA,
      delta_cfi_vs_reference = NA,
      delta_rmsea_vs_reference = NA,
      error = "optimizer did not reach an admissible stationary solution"
    ))
  }
  pe <- parameterEstimates(fit)
  negative_residual <- any(
    pe$op == "~~" & pe$lhs == pe$rhs &
      pe$lhs %in% indicators & pe$est < -1e-8
  )
  measures <- fitMeasures(fit, c("cfi.robust", "rmsea.robust", "srmr"))
  delta_cfi <- NA_real_
  delta_rmsea <- NA_real_
  if (!is.null(reference) && !inherits(reference, "error")) {
    ref_measures <- fitMeasures(reference, c("cfi.robust", "rmsea.robust"))
    delta_cfi <- unname(measures["cfi.robust"] - ref_measures["cfi.robust"])
    delta_rmsea <- unname(
      measures["rmsea.robust"] - ref_measures["rmsea.robust"]
    )
  }
  data.frame(
    model = name,
    converged = converged,
    admissible = converged && !negative_residual,
    n_used = sum(lavInspect(fit, "nobs")),
    n_freed = n_freed,
    cfi_robust = unname(measures["cfi.robust"]),
    rmsea_robust = unname(measures["rmsea.robust"]),
    srmr = unname(measures["srmr"]),
    delta_cfi_vs_reference = delta_cfi,
    delta_rmsea_vs_reference = delta_rmsea,
    error = ""
  )
}

metric_rows <- list(
  fit_row("configural", metric_candidates$configural),
  fit_row(
    "full_metric",
    metric_candidates$full_metric,
    metric_candidates$configural,
    0
  ),
  fit_row(
    "partial_metric_core",
    metric_candidates$partial_metric_core,
    metric_candidates$configural,
    length(metric_core)
  ),
  fit_row(
    "partial_metric_expanded",
    metric_candidates$partial_metric_expanded,
    metric_candidates$configural,
    length(metric_expanded)
  ),
  fit_row(
    "partial_metric_all_screened",
    metric_candidates$partial_metric_all_screened,
    metric_candidates$configural,
    length(metric_all_screened)
  )
)
metric_table <- do.call(rbind, metric_rows)

metric_pass <- (
  metric_table$model != "configural" &
    metric_table$admissible &
    metric_table$delta_cfi_vs_reference >= -0.010 &
    metric_table$delta_rmsea_vs_reference <= 0.015
)
if (!any(metric_pass)) {
  stop("No partial metric model met the predefined fit-change criteria.")
}
chosen_metric_name <- metric_table$model[which(metric_pass)[1]]
chosen_metric_fit <- metric_candidates[[chosen_metric_name]]
chosen_metric_partial <- switch(
  chosen_metric_name,
  full_metric = character(0),
  partial_metric_core = metric_core,
  partial_metric_expanded = metric_expanded,
  partial_metric_all_screened = metric_all_screened
)

scalar_candidates <- list(
  full_scalar = fit_group_model(c("loadings", "intercepts")),
  partial_scalar_core = fit_group_model(
    c("loadings", "intercepts"),
    unique(c(chosen_metric_partial, scalar_core))
  ),
  partial_scalar_all_screened = fit_group_model(
    c("loadings", "intercepts"),
    unique(c(chosen_metric_partial, scalar_all_screened))
  )
)

scalar_rows <- list(
  fit_row(
    "full_scalar",
    scalar_candidates$full_scalar,
    chosen_metric_fit,
    0
  ),
  fit_row(
    "partial_scalar_core",
    scalar_candidates$partial_scalar_core,
    chosen_metric_fit,
    length(unique(c(chosen_metric_partial, scalar_core)))
  ),
  fit_row(
    "partial_scalar_all_screened",
    scalar_candidates$partial_scalar_all_screened,
    chosen_metric_fit,
    length(unique(c(chosen_metric_partial, scalar_all_screened)))
  )
)
scalar_table <- do.call(rbind, scalar_rows)

scalar_pass <- (
  scalar_table$admissible &
    scalar_table$delta_cfi_vs_reference >= -0.010 &
    scalar_table$delta_rmsea_vs_reference <= 0.015
)
if (!any(scalar_pass)) {
  stop("No partial scalar model met the predefined fit-change criteria.")
}
chosen_scalar_name <- scalar_table$model[which(scalar_pass)[1]]
chosen_scalar_fit <- scalar_candidates[[chosen_scalar_name]]
chosen_scalar_partial <- switch(
  chosen_scalar_name,
  full_scalar = chosen_metric_partial,
  partial_scalar_core = unique(c(chosen_metric_partial, scalar_core)),
  partial_scalar_all_screened = unique(
    c(chosen_metric_partial, scalar_all_screened)
  )
)

weighted_chosen_scalar <- fit_group_model(
  c("loadings", "intercepts"),
  chosen_scalar_partial,
  weighted = TRUE
)
weighted_row <- fit_row(
  "chosen_partial_scalar_weighted",
  weighted_chosen_scalar,
  NULL,
  length(chosen_scalar_partial)
)

extract_latent_means <- function(fit, model_name) {
  if (inherits(fit, "error")) {
    return(data.frame())
  }
  pe <- parameterEstimates(fit, ci = TRUE)
  pt <- parTable(fit)
  means <- pe[pe$op == "~1" & pe$lhs == "general", ]
  mean_pt <- pt[pt$op == "~1" & pt$lhs == "general", ]
  data.frame(
    model = model_name,
    education_group = means$group,
    estimate = means$est,
    se = means$se,
    ci_lower = means$ci.lower,
    ci_upper = means$ci.upper,
    free_index = mean_pt$free
  )
}

latent_means <- rbind(
  extract_latent_means(scalar_candidates$full_scalar, "full_scalar"),
  extract_latent_means(chosen_scalar_fit, chosen_scalar_name),
  extract_latent_means(weighted_chosen_scalar, "chosen_partial_scalar_weighted")
)

extract_latent_contrasts <- function(fit, model_name) {
  if (inherits(fit, "error")) {
    return(data.frame())
  }
  pt <- parTable(fit)
  means <- pt[pt$op == "~1" & pt$lhs == "general", ]
  if (nrow(means) != 3) {
    return(data.frame())
  }
  vc <- lavInspect(fit, "vcov")
  covariance <- function(row_a, row_b) {
    free_a <- means$free[row_a]
    free_b <- means$free[row_b]
    if (free_a == 0 || free_b == 0) {
      return(0)
    }
    vc[free_a, free_b]
  }
  pairs <- list(c(2, 1), c(3, 1), c(3, 2))
  do.call(rbind, lapply(pairs, function(pair) {
    high_row <- which(means$group == pair[1])
    low_row <- which(means$group == pair[2])
    difference <- means$est[high_row] - means$est[low_row]
    variance <- (
      covariance(high_row, high_row) +
        covariance(low_row, low_row) -
        2 * covariance(high_row, low_row)
    )
    se <- sqrt(max(variance, 0))
    data.frame(
      model = model_name,
      contrast = paste0("group", pair[1], "_minus_group", pair[2]),
      estimate = difference,
      se = se,
      ci_lower = difference - 1.96 * se,
      ci_upper = difference + 1.96 * se
    )
  }))
}

latent_contrasts <- do.call(
  rbind,
  list(
    extract_latent_contrasts(scalar_candidates$full_scalar, "full_scalar"),
    extract_latent_contrasts(chosen_scalar_fit, chosen_scalar_name),
    extract_latent_contrasts(
      weighted_chosen_scalar,
      "chosen_partial_scalar_weighted"
    )
  )
)

pooled_fit <- cfa(
  model_syntax,
  data = d,
  estimator = "MLR",
  missing = "fiml",
  std.lv = TRUE,
  meanstructure = TRUE
)
pooled_scores <- as.numeric(lavPredict(pooled_fit, type = "lv")[, "general"])
d$pooled_general_score <- pooled_scores

weighted_mean_se <- function(values, weights) {
  keep <- !is.na(values) & !is.na(weights) & weights > 0
  values <- values[keep]
  weights <- weights[keep]
  mean_value <- sum(weights * values) / sum(weights)
  n_eff <- sum(weights)^2 / sum(weights^2)
  variance <- sum(weights * (values - mean_value)^2) / sum(weights)
  se <- sqrt(variance / n_eff)
  c(mean = mean_value, se = se, n_eff = n_eff)
}

pooled_group_rows <- lapply(1:3, function(group_value) {
  target <- d$edu3 == group_value
  stats <- weighted_mean_se(
    d$pooled_general_score[target],
    d$r4wtrespb[target]
  )
  data.frame(
    model = "pooled_unqualified_factor_score",
    education_group = group_value,
    estimate = stats["mean"],
    se = stats["se"],
    ci_lower = stats["mean"] - 1.96 * stats["se"],
    ci_upper = stats["mean"] + 1.96 * stats["se"],
    n_effective = stats["n_eff"]
  )
})
pooled_group_means <- do.call(rbind, pooled_group_rows)
pooled_contrasts <- do.call(rbind, lapply(
  list(c(2, 1), c(3, 1), c(3, 2)),
  function(pair) {
    high <- pooled_group_means[
      pooled_group_means$education_group == pair[1],
    ]
    low <- pooled_group_means[
      pooled_group_means$education_group == pair[2],
    ]
    difference <- high$estimate - low$estimate
    se <- sqrt(high$se^2 + low$se^2)
    data.frame(
      model = "pooled_unqualified_factor_score",
      contrast = paste0("group", pair[1], "_minus_group", pair[2]),
      estimate = difference,
      se = se,
      ci_lower = difference - 1.96 * se,
      ci_upper = difference + 1.96 * se
    )
  }
))

standardized <- standardizedSolution(chosen_scalar_fit)
general_loadings <- standardized[
  standardized$op == "=~" & standardized$lhs == "general",
  c("rhs", "group", "est.std")
]
loading_wide <- reshape(
  general_loadings,
  idvar = "rhs",
  timevar = "group",
  direction = "wide"
)
names(loading_wide) <- sub("est.std.", "loading_group", names(loading_wide))
loading_columns <- grep("^loading_group", names(loading_wide), value = TRUE)
loading_wide$loading_range <- apply(
  loading_wide[loading_columns],
  1,
  function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
)
loading_wide <- loading_wide[
  order(loading_wide$loading_range, decreasing = TRUE),
]

release_table <- data.frame(
  stage = c(
    rep("chosen_metric", length(chosen_metric_partial)),
    rep(
      "additional_scalar",
      length(setdiff(chosen_scalar_partial, chosen_metric_partial))
    )
  ),
  parameter = c(
    chosen_metric_partial,
    setdiff(chosen_scalar_partial, chosen_metric_partial)
  )
)

write.csv(
  metric_table,
  file.path(table_dir, paste0("b1_hcap_partial_metric_fit_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  rbind(scalar_table, weighted_row),
  file.path(table_dir, paste0("b1_hcap_partial_scalar_fit_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  release_table,
  file.path(table_dir, paste0("b1_hcap_partial_invariance_releases_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  latent_means[, setdiff(names(latent_means), "free_index")],
  file.path(table_dir, paste0("b1_hcap_education_latent_means_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  rbind(latent_contrasts, pooled_contrasts),
  file.path(table_dir, paste0("b1_hcap_education_contrasts_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  pooled_group_means,
  file.path(table_dir, paste0("b1_hcap_pooled_score_education_means_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  loading_wide,
  file.path(table_dir, paste0("b1_hcap_partial_scalar_loading_dif_", date_tag, ".csv")),
  row.names = FALSE
)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

fit_markdown <- function(table) {
  lines <- c(
    "| Model | Freed | Robust CFI | Robust RMSEA | SRMR | Delta CFI | Delta RMSEA |",
    "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"
  )
  for (i in seq_len(nrow(table))) {
    r <- table[i, ]
    lines <- c(
      lines,
      sprintf(
        "| %s | %s | %s | %s | %s | %s | %s |",
        r$model,
        r$n_freed,
        fmt(r$cfi_robust),
        fmt(r$rmsea_robust),
        fmt(r$srmr),
        fmt(r$delta_cfi_vs_reference),
        fmt(r$delta_rmsea_vs_reference)
      )
    )
  }
  lines
}

contrast_markdown <- c(
  "| Model | Contrast | Standardized difference | 95% CI |",
  "| --- | --- | ---: | --- |"
)
all_contrasts <- rbind(latent_contrasts, pooled_contrasts)
for (i in seq_len(nrow(all_contrasts))) {
  r <- all_contrasts[i, ]
  contrast_markdown <- c(
    contrast_markdown,
    sprintf(
      "| %s | %s | %s | %s to %s |",
      r$model,
      r$contrast,
      fmt(r$estimate),
      fmt(r$ci_lower),
      fmt(r$ci_upper)
    )
  )
}

top_dif <- head(loading_wide, 8)
dif_markdown <- c(
  "| Indicator | Group 1 loading | Group 2 loading | Group 3 loading | Range |",
  "| --- | ---: | ---: | ---: | ---: |"
)
for (i in seq_len(nrow(top_dif))) {
  r <- top_dif[i, ]
  dif_markdown <- c(
    dif_markdown,
    sprintf(
      "| %s | %s | %s | %s | %s |",
      r$rhs,
      fmt(r$loading_group1),
      fmt(r$loading_group2),
      fmt(r$loading_group3),
      fmt(r$loading_range)
    )
  )
}

full_scalar_g31 <- latent_contrasts[
  latent_contrasts$model == "full_scalar" &
    latent_contrasts$contrast == "group3_minus_group1",
]
partial_g31 <- latent_contrasts[
  latent_contrasts$model == chosen_scalar_name &
    latent_contrasts$contrast == "group3_minus_group1",
]
gradient_change <- partial_g31$estimate - full_scalar_g31$estimate
gradient_change_pct <- gradient_change / full_scalar_g31$estimate * 100

report <- c(
  "# B1 HCAP Education Partial-Invariance Audit",
  "",
  paste0("- Date: ", date_tag),
  "- Status: partial-invariance model selection completed; anchor not yet frozen",
  "",
  "## Predefined Decision Rule",
  "",
  paste(
    "Select the first admissible staged model with delta robust CFI >= -0.010",
    "and delta robust RMSEA <= 0.015 relative to the less constrained reference.",
    "Parameters were released in bundles defined before refitting from the",
    "previous BH-adjusted score-test screen. No final trend or inequality",
    "outcomes were inspected during model selection."
  ),
  "",
  "## Partial Metric Models",
  "",
  fit_markdown(metric_table),
  "",
  paste0("Selected metric model: `", chosen_metric_name, "`."),
  "",
  "## Partial Scalar Models",
  "",
  fit_markdown(rbind(scalar_table, weighted_row)),
  "",
  paste0("Selected scalar model: `", chosen_scalar_name, "`."),
  "",
  "## Education Gradient",
  "",
  "Education groups: 1 = no formal education/illiterate; 2 = some schooling below middle school; 3 = middle school or higher.",
  "Latent-model contrasts use the reference-group latent SD; the pooled-score comparator uses its pooled score SD and is not numerically exchangeable with the latent metric.",
  "",
  contrast_markdown,
  "",
  sprintf(
    "For group 3 versus group 1, DIF correction changed the latent mean gap by %s reference-group SD (%s%%) relative to the full-scalar estimate.",
    fmt(gradient_change),
    fmt(gradient_change_pct, 1)
  ),
  "",
  "## Largest Remaining Loading Differences",
  "",
  dif_markdown,
  "",
  "## Decision",
  "",
  paste(
    "The selected partial-invariance model meets the predefined global",
    "fit-change criteria and remains stable after applying survey weights.",
    "A common general cognition metric is therefore feasible only with explicit",
    "education-related loading and intercept releases. The unqualified pooled",
    "score and full-scalar model remain comparators, not primary inequality",
    "measures. Item-level MI is evaluated in a separate sensitivity analysis."
  ),
  "",
  "## Interpretation Boundary",
  "",
  paste(
    "Partial invariance reduces, but does not prove elimination of, education",
    "measurement bias. Released parameters are substantively interpretable DIF",
    "signals and must be reported rather than treated as purely technical",
    "modifications."
  )
)

writeLines(report, report_path, useBytes = TRUE)
cat("wrote", report_path, "\n")
