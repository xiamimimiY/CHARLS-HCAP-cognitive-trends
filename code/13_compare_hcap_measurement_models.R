#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lavaan)
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
report_path <- file.path(table_dir, "hcap_measurement_model_report.md")
date_tag <- "2026-06-19"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(input_path, check.names = FALSE)
d$model_weight <- d$r4wtrespb / mean(d$r4wtrespb[d$r4wtrespb > 0], na.rm = TRUE)

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

d <- d[rowSums(!is.na(d[indicators])) >= 1, ]

one_factor <- paste0(
  "general =~ ",
  paste(indicators, collapse = " + ")
)

two_factor <- paste(
  "memory =~ z_test_short_learning_trial1 + z_test_short_delayed +",
  "z_test_ten_learning_trial1 + z_test_ten_learning_trial2 +",
  "z_test_ten_learning_trial3 + z_test_ten_delayed + z_test_recognition",
  "nonmemory =~ z_test_serial7 + z_test_animal_fluency + z_test_language +",
  "z_test_visuospatial + z_test_orientation",
  sep = "\n"
)

bifactor_memory <- paste(
  one_factor,
  "memory_specific =~ z_test_short_learning_trial1 + z_test_short_delayed +",
  "z_test_ten_learning_trial1 +",
  "z_test_ten_learning_trial2 + z_test_ten_learning_trial3 +",
  "z_test_ten_delayed + z_test_recognition",
  "general ~~ 0*memory_specific",
  sep = "\n"
)

fit_model <- function(model_syntax, weighted = FALSE) {
  model_data <- d
  args <- list(
    model = model_syntax,
    data = model_data,
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
  tryCatch(
    do.call(cfa, args),
    error = function(e) e
  )
}

fits <- list(
  one_factor = fit_model(one_factor),
  two_factor = fit_model(two_factor),
  bifactor_memory = fit_model(bifactor_memory),
  bifactor_memory_weighted = fit_model(bifactor_memory, weighted = TRUE)
)

fit_rows <- list()
loading_rows <- list()
residual_rows <- list()

for (model_name in names(fits)) {
  fit <- fits[[model_name]]
  if (inherits(fit, "error")) {
    fit_rows[[length(fit_rows) + 1]] <- data.frame(
      model = model_name,
      converged = FALSE,
      admissible = FALSE,
      n_used = NA,
      cfi_robust = NA,
      tli_robust = NA,
      rmsea_robust = NA,
      rmsea_ci_lower_robust = NA,
      rmsea_ci_upper_robust = NA,
      srmr = NA,
      aic = NA,
      bic = NA,
      general_score_determinacy = NA,
      memory_score_determinacy = NA,
      error = conditionMessage(fit)
    )
    next
  }

  converged <- lavInspect(fit, "converged")
  pe <- parameterEstimates(fit, standardized = TRUE)
  negative_residual <- any(
    pe$op == "~~" & pe$lhs == pe$rhs &
      pe$lhs %in% indicators & pe$est < -1e-8
  )
  admissible <- isTRUE(converged) && !negative_residual
  measures <- fitMeasures(
    fit,
    c(
      "cfi.robust", "tli.robust", "rmsea.robust",
      "rmsea.ci.lower.robust", "rmsea.ci.upper.robust",
      "srmr", "aic", "bic"
    ),
    fm.args = list(
      standard.test = "default",
      scaled.test = "default",
      rmsea.ci.level = 0.90,
      rmsea.close.h0 = 0.05,
      rmsea.notclose.h0 = 0.08,
      robust = TRUE,
      cat.check.pd = TRUE
    )
  )
  determinacy <- tryCatch(
    lavInspect(fit, "fs.determinacy"),
    error = function(e) numeric(0)
  )

  fit_rows[[length(fit_rows) + 1]] <- data.frame(
    model = model_name,
    converged = converged,
    admissible = admissible,
    n_used = lavInspect(fit, "nobs"),
    cfi_robust = unname(measures["cfi.robust"]),
    tli_robust = unname(measures["tli.robust"]),
    rmsea_robust = unname(measures["rmsea.robust"]),
    rmsea_ci_lower_robust = unname(measures["rmsea.ci.lower.robust"]),
    rmsea_ci_upper_robust = unname(measures["rmsea.ci.upper.robust"]),
    srmr = unname(measures["srmr"]),
    aic = unname(measures["aic"]),
    bic = unname(measures["bic"]),
    general_score_determinacy = if ("general" %in% names(determinacy)) {
      unname(determinacy["general"])
    } else {
      NA
    },
    memory_score_determinacy = if ("memory_specific" %in% names(determinacy)) {
      unname(determinacy["memory_specific"])
    } else if ("memory" %in% names(determinacy)) {
      unname(determinacy["memory"])
    } else {
      NA
    },
    error = ""
  )

  loadings <- pe[pe$op == "=~", c("lhs", "rhs", "est", "se", "pvalue", "std.all")]
  loadings$model <- model_name
  loading_rows[[length(loading_rows) + 1]] <- loadings[
    , c("model", "lhs", "rhs", "est", "se", "pvalue", "std.all")
  ]

  residuals <- pe[
    pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% indicators,
    c("lhs", "est", "se", "std.all")
  ]
  residuals$model <- model_name
  residual_rows[[length(residual_rows) + 1]] <- residuals[
    , c("model", "lhs", "est", "se", "std.all")
  ]
}

fit_table <- do.call(rbind, fit_rows)
loading_table <- do.call(rbind, loading_rows)
residual_table <- do.call(rbind, residual_rows)

bifactor_indices <- list()
for (model_name in c("bifactor_memory", "bifactor_memory_weighted")) {
  fit <- fits[[model_name]]
  if (inherits(fit, "error") || !lavInspect(fit, "converged")) {
    next
  }
  standardized <- standardizedSolution(fit)
  g <- standardized[
    standardized$op == "=~" & standardized$lhs == "general",
    c("rhs", "est.std")
  ]
  s <- standardized[
    standardized$op == "=~" & standardized$lhs == "memory_specific",
    c("rhs", "est.std")
  ]
  theta <- standardized[
    standardized$op == "~~" &
      standardized$lhs == standardized$rhs &
      standardized$lhs %in% indicators,
    c("lhs", "est.std")
  ]
  sum_g_sq <- sum(g$est.std^2)
  sum_s_sq <- sum(s$est.std^2)
  denominator <- sum(g$est.std)^2 + sum(s$est.std)^2 + sum(theta$est.std)
  omega_h <- sum(g$est.std)^2 / denominator
  omega_total <- (sum(g$est.std)^2 + sum(s$est.std)^2) / denominator
  total_pairs <- choose(length(indicators), 2)
  contaminated_pairs <- choose(nrow(s), 2)
  bifactor_indices[[length(bifactor_indices) + 1]] <- data.frame(
    model = model_name,
    ecv_general = sum_g_sq / (sum_g_sq + sum_s_sq),
    omega_hierarchical = omega_h,
    omega_total = omega_total,
    relative_omega = omega_h / omega_total,
    puc = (total_pairs - contaminated_pairs) / total_pairs
  )
}
bifactor_table <- do.call(rbind, bifactor_indices)

fit_invariance <- function(equal = character(0)) {
  tryCatch(
    cfa(
      bifactor_memory,
      data = d[!is.na(d$edu3), ],
      group = "edu3",
      group.equal = equal,
      estimator = "MLR",
      missing = "fiml",
      std.lv = TRUE,
      meanstructure = TRUE
    ),
    error = function(e) e
  )
}

invariance_fits <- list(
  configural = fit_invariance(),
  metric = fit_invariance("loadings"),
  scalar = fit_invariance(c("loadings", "intercepts"))
)

invariance_rows <- list()
for (level in names(invariance_fits)) {
  fit <- invariance_fits[[level]]
  if (inherits(fit, "error")) {
    invariance_rows[[length(invariance_rows) + 1]] <- data.frame(
      level = level,
      converged = FALSE,
      cfi_robust = NA,
      rmsea_robust = NA,
      srmr = NA,
      delta_cfi = NA,
      delta_rmsea = NA,
      error = conditionMessage(fit)
    )
    next
  }
  measures <- fitMeasures(fit, c("cfi.robust", "rmsea.robust", "srmr"))
  invariance_rows[[length(invariance_rows) + 1]] <- data.frame(
    level = level,
    converged = lavInspect(fit, "converged"),
    cfi_robust = unname(measures["cfi.robust"]),
    rmsea_robust = unname(measures["rmsea.robust"]),
    srmr = unname(measures["srmr"]),
    delta_cfi = NA,
    delta_rmsea = NA,
    error = ""
  )
}
invariance_table <- do.call(rbind, invariance_rows)
for (i in 2:nrow(invariance_table)) {
  invariance_table$delta_cfi[i] <- (
    invariance_table$cfi_robust[i] - invariance_table$cfi_robust[i - 1]
  )
  invariance_table$delta_rmsea[i] <- (
    invariance_table$rmsea_robust[i] - invariance_table$rmsea_robust[i - 1]
  )
}

constraint_diagnostics <- function(fit, level) {
  if (inherits(fit, "error") || !lavInspect(fit, "converged")) {
    return(data.frame())
  }
  score <- tryCatch(
    lavTestScore(fit, epc = TRUE, univariate = TRUE),
    error = function(e) NULL
  )
  if (is.null(score) || is.null(score$uni)) {
    return(data.frame())
  }
  uni <- as.data.frame(score$uni)
  pt <- parTable(fit)
  descriptions <- paste0(
    pt$lhs, " ", pt$op, " ", pt$rhs, " [group ", pt$group, "]"
  )
  names(descriptions) <- pt$plabel
  uni$level <- level
  uni$left_parameter <- unname(descriptions[as.character(uni$lhs)])
  uni$right_parameter <- unname(descriptions[as.character(uni$rhs)])
  uni$p_bh <- p.adjust(uni$p.value, method = "BH")
  uni[order(uni$X2, decreasing = TRUE), ]
}

metric_diagnostics <- constraint_diagnostics(invariance_fits$metric, "metric")
scalar_diagnostics <- constraint_diagnostics(invariance_fits$scalar, "scalar")
constraint_table <- rbind(metric_diagnostics, scalar_diagnostics)

write.csv(
  fit_table,
  file.path(table_dir, paste0("b1_hcap_fiml_model_fit_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  loading_table,
  file.path(table_dir, paste0("b1_hcap_fiml_loadings_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  residual_table,
  file.path(table_dir, paste0("b1_hcap_fiml_residual_variances_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  bifactor_table,
  file.path(table_dir, paste0("b1_hcap_fiml_bifactor_indices_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  invariance_table,
  file.path(table_dir, paste0("b1_hcap_fiml_education_invariance_", date_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  constraint_table,
  file.path(table_dir, paste0("b1_hcap_fiml_invariance_score_tests_", date_tag, ".csv")),
  row.names = FALSE
)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

md_table <- c(
  "| Model | Converged | Admissible | N | Robust CFI | Robust TLI | Robust RMSEA (90% CI) | SRMR | Score determinacy | BIC |",
  "| --- | --- | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: |"
)
for (i in seq_len(nrow(fit_table))) {
  r <- fit_table[i, ]
  md_table <- c(
    md_table,
    sprintf(
      "| %s | %s | %s | %s | %s | %s | %s (%s-%s) | %s | %s | %s |",
      r$model,
      r$converged,
      r$admissible,
      ifelse(is.na(r$n_used), "NA", r$n_used),
      fmt(r$cfi_robust),
      fmt(r$tli_robust),
      fmt(r$rmsea_robust),
      fmt(r$rmsea_ci_lower_robust),
      fmt(r$rmsea_ci_upper_robust),
      fmt(r$srmr),
      fmt(r$general_score_determinacy),
      fmt(r$bic, 1)
    )
  )
}

one <- fit_table[fit_table$model == "one_factor", ]
bi <- fit_table[fit_table$model == "bifactor_memory", ]
two <- fit_table[fit_table$model == "two_factor", ]

decision <- "No candidate measurement model satisfies the comparison criteria."
if (
  nrow(one) == 1 && isTRUE(one$admissible) &&
    one$cfi_robust >= 0.90 && one$rmsea_robust <= 0.08 && one$srmr <= 0.08
) {
  decision <- paste(
    "The one-factor model meets the predefined minimum fit criteria.",
    "It remains a candidate for a full-information global-cognition model,",
    "subject to the planned invariance and missing-data analyses."
  )
}

if (
  nrow(bi) == 1 && isTRUE(bi$admissible) &&
    nrow(one) == 1 && isTRUE(one$admissible) &&
    bi$cfi_robust - one$cfi_robust >= 0.01 &&
    bi$bic < one$bic
) {
  bi_indices <- bifactor_table[bifactor_table$model == "bifactor_memory", ]
  decision <- paste(
    "The bifactor-memory model materially improves fit over the one-factor model.",
    sprintf(
      "Its general factor has ECV %.3f, omega hierarchical %.3f, and score determinacy %.3f.",
      bi_indices$ecv_general,
      bi_indices$omega_hierarchical,
      bi$general_score_determinacy
    ),
    "It is the preferred FIML anchor candidate; education-group invariance",
    "and item-level MI are evaluated in subsequent scripts."
  )
} else if (
  nrow(two) == 1 && isTRUE(two$admissible) &&
    nrow(one) == 1 && isTRUE(one$admissible) &&
    two$cfi_robust - one$cfi_robust >= 0.01 &&
    two$bic < one$bic
) {
  decision <- paste(
    "A correlated memory/nonmemory structure materially improves fit.",
    "A single global factor is therefore not yet adequate as the sole anchor;",
    "domain-aware calibration must be retained."
  )
}

metric_row <- invariance_table[invariance_table$level == "metric", ]
scalar_row <- invariance_table[invariance_table$level == "scalar", ]
if (
  nrow(metric_row) == 1 &&
    !is.na(metric_row$delta_cfi) &&
    metric_row$delta_cfi < -0.01
) {
  decision <- paste(
    decision,
    sprintf(
      "Full metric invariance fails across education groups (delta CFI %.3f);",
      metric_row$delta_cfi
    ),
    sprintf(
      "full scalar invariance also fails (additional delta CFI %.3f).",
      scalar_row$delta_cfi
    ),
    "The next required step is indicator-level DIF/partial-invariance analysis.",
    "Education inequalities must not be estimated from an unqualified common score."
  )
}

report <- c(
  "# HCAP Measurement-Model Comparison",
  "",
  paste0("- Date: ", date_tag),
  "",
  "## Objective",
  "",
  paste(
    "Test whether CHARLS-HCAP test-level indicators support a",
    "missingness-aware continuous cognition anchor."
  ),
  "",
  "## Models",
  "",
  "1. One general cognition factor.",
  "2. Correlated memory and nonmemory factors.",
  "3. General factor plus an orthogonal memory-specific factor.",
  "4. Survey-weighted sensitivity fit of the bifactor-memory model.",
  "",
  paste0(
    "All models used robust maximum likelihood (MLR), full-information maximum ",
    "likelihood for missing indicators, standardized latent variances, and ",
    "participants aged 60 years or older with at least one observed indicator ",
    "(N = ", nrow(d), ")."
  ),
  "",
  "## Fit",
  "",
  md_table,
  "",
  "## Bifactor Indices",
  "",
  "| Model | ECV general | Omega hierarchical | Omega total | Relative omega | PUC |",
  "| --- | ---: | ---: | ---: | ---: | ---: |",
  if (nrow(bifactor_table) > 0) {
    apply(bifactor_table, 1, function(r) {
      sprintf(
        "| %s | %s | %s | %s | %s | %s |",
        r[["model"]],
        fmt(as.numeric(r[["ecv_general"]])),
        fmt(as.numeric(r[["omega_hierarchical"]])),
        fmt(as.numeric(r[["omega_total"]])),
        fmt(as.numeric(r[["relative_omega"]])),
        fmt(as.numeric(r[["puc"]]))
      )
    })
  } else {
    "| No admissible bifactor model | NA | NA | NA | NA | NA |"
  },
  "",
  "## Education-Group Measurement Invariance",
  "",
  paste(
    "Groups are: 1 = no formal education/illiterate; 2 = some schooling below",
    "middle school; 3 = middle school or higher. Metric and scalar invariance",
    "are screened using changes in robust CFI and RMSEA."
  ),
  "",
  "| Level | Converged | Robust CFI | Robust RMSEA | SRMR | Delta CFI | Delta RMSEA |",
  "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
  apply(invariance_table, 1, function(r) {
    sprintf(
      "| %s | %s | %s | %s | %s | %s | %s |",
      r[["level"]],
      r[["converged"]],
      fmt(as.numeric(r[["cfi_robust"]])),
      fmt(as.numeric(r[["rmsea_robust"]])),
      fmt(as.numeric(r[["srmr"]])),
      fmt(as.numeric(r[["delta_cfi"]])),
      fmt(as.numeric(r[["delta_rmsea"]]))
    )
  }),
  "",
  "## Largest Equality-Constraint Violations",
  "",
  paste(
    "These ordinary score tests are screening diagnostics because lavaan does",
    "not currently implement robust MLR score tests. Final partial-invariance",
    "decisions must be confirmed by refitting predefined freed constraints."
  ),
  "",
  if (nrow(constraint_table) > 0) {
    top <- head(constraint_table[constraint_table$level == "metric", ], 10)
    c(
      "| Left parameter | Right parameter | Score X2 | BH-adjusted p |",
      "| --- | --- | ---: | ---: |",
      apply(top, 1, function(r) {
        sprintf(
          "| %s | %s | %s | %s |",
          r[["left_parameter"]],
          r[["right_parameter"]],
          fmt(as.numeric(r[["X2"]])),
          formatC(as.numeric(r[["p_bh"]]), digits = 3, format = "g")
        )
      })
    )
  } else {
    "Constraint-level score tests were unavailable."
  },
  "",
  "## Decision",
  "",
  decision,
  "",
  "## Interpretation Boundary",
  "",
  paste(
    "This analysis evaluates latent structure and missing-indicator feasibility.",
    "It does not validate a clinical cut-point, MCI, dementia, or Alzheimer",
    "disease. Final anchor selection still requires interpretation of education",
    "invariance, other subgroup checks, and comparison with item-level multiple",
    "imputation."
  )
)

writeLines(report, report_path, useBytes = TRUE)
cat("wrote", report_path, "\n")
