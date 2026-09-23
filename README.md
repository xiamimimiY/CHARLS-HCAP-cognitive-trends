# Measurement calibration and selective non-completion in estimates of cognitive trends and social inequalities: a CHARLS-HCAP analysis

## Overview

This repository contains the principal analysis code for a study of cognitive
measurement, population trends, and social inequalities using repeated China
Health and Retirement Longitudinal Study (CHARLS) surveys and the 2018
CHARLS Harmonized Cognitive Assessment Protocol (HCAP). The workflow constructs
12 HCAP indicators, compares latent measurement models, evaluates education
partial invariance, links routine CHARLS cognition to the HCAP metric, and
estimates survey-weighted trends and inequalities. It also implements
wave-specific multiple imputation, item-missingness and selection analyses,
community-cluster bootstrap inference, attrition and mortality analyses, and
internal validation of a research-use scoring model. The repository contains
code only and does not include participant-level data or analysis outputs.

## Data Availability

CHARLS data are not included in this repository because they are subject to the
data-use requirements of the original data provider. Researchers wishing to
run the analyses must obtain the corresponding CHARLS data independently.
The repository contains analysis code but does not redistribute
participant-level CHARLS data.

By default, scripts expect authorized inputs under:

```text
data/
├── raw/
│   ├── 2018 CHARLS Wave 4/
│   │   ├── Cognition.dta
│   │   ├── Insider.dta
│   │   └── Health_Status_and_Functioning.dta
│   └── Harmonized CHARLS/
│       └── H_CHARLS_D_Data.zip
└── derived/                 # created locally; excluded from version control
```

The `data/raw/` location can be overridden without editing code:

```bash
export CHARLS_DATA_ROOT="/path/to/authorized/CHARLS/data"
```

All raw files are read only. Participant-level model inputs are written to
`data/derived/`, and numerical and graphical outputs are written to `results/`.
Both directories are excluded from version control.

## Code Structure

| Script | Purpose |
|---|---|
| `code/project_paths.py` | Defines repository-relative raw, derived, result, and temporary paths. |
| `code/11_construct_hcap_indicators.py` | Constructs the 12 scored HCAP indicators and final `raeduc_c` education groups. |
| `code/12_build_hcap_model_input.py` | Standardizes HCAP indicators and builds the local measurement-model input. |
| `code/13_compare_hcap_measurement_models.R` | Compares one-factor, two-factor, and bifactor HCAP models using robust FIML. |
| `code/14_test_education_measurement_invariance.R` | Tests metric and scalar invariance and estimates the education partial-invariance model. |
| `code/16_hcap_all_missing_sensitivity.R` | Evaluates expanded-sample and lower-score scenarios for participants with all HCAP indicators missing. |
| `code/20_prepare_longitudinal_cognition.py` | Builds the repeated cross-sectional routine-cognition input for 2011-2018. |
| `code/21_estimate_calibrated_trends_inequalities.R` | Runs wave-specific imputation, HCAP-to-routine calibration, survey-weighted trends, inequalities, and principal sensitivity analyses. |
| `code/22_hcap_item_mi_sensitivity.R` | Propagates HCAP item-level imputation through calibration and longitudinal estimation. |
| `code/23_full_pipeline_cluster_bootstrap.R` | Refits the measurement-to-estimation pipeline in a community-cluster bootstrap. |
| `code/24_all_missing_selection_bounds.R` | Translates all-HCAP-missing scenarios into longitudinal selection bounds. |
| `code/25_attrition_mortality_sensitivity.R` | Evaluates fixed/rolling cohorts, inverse-probability weighting, and mortality-composite outcomes. |
| `code/26_validate_scoring_model.R` | Estimates and internally validates the research-use HCAP-linked scoring model. |
| `code/27_plot_trends_inequalities.R` | Produces the principal trend and social-inequality displays from aggregate estimates. |
| `code/28_plot_workflow_missingness_sensitivity.R` | Produces workflow, HCAP missingness, and robustness displays. |
| `code/29_build_main_tables.R` | Produces participant-characteristic and measurement/calibration summary tables. |

The separate sensitivity scripts are retained because each reproduces a
distinct reported analysis; exploratory and superseded scripts are not
included.

## Running the Analysis

Create the local output directories and install the Python dependencies:

```bash
mkdir -p data/derived results/tables results/figures
python -m pip install -r requirements.txt
```

Run the scripts from the repository root in this order:

```bash
python code/12_build_hcap_model_input.py
Rscript code/13_compare_hcap_measurement_models.R
Rscript code/14_test_education_measurement_invariance.R
Rscript code/16_hcap_all_missing_sensitivity.R
python code/20_prepare_longitudinal_cognition.py
Rscript code/21_estimate_calibrated_trends_inequalities.R
Rscript code/22_hcap_item_mi_sensitivity.R
Rscript code/23_full_pipeline_cluster_bootstrap.R
Rscript code/24_all_missing_selection_bounds.R
Rscript code/25_attrition_mortality_sensitivity.R
Rscript code/26_validate_scoring_model.R
Rscript code/27_plot_trends_inequalities.R
Rscript code/28_plot_workflow_missingness_sensitivity.R
Rscript code/29_build_main_tables.R
```

`11_construct_hcap_indicators.py` defines the indicator-construction functions;
`12_build_hcap_model_input.py` imports and executes them, so the first command
performs both indicator construction and model-input preparation. The primary analysis uses 50
imputations by default, and the full-pipeline bootstrap uses 200 replicates.
Environment variables such as `B1_M`, `B1_MAXIT`, and `B1_BOOT_REPS` are
available in the relevant scripts for local computational checks; the defaults
correspond to the final analysis.

## Software

The analysis was run with Python 3.8.16 and R 4.1.3. Python dependencies are
listed in `requirements.txt`.

Principal R packages:

- `lavaan` 0.6-21
- `mice` 3.15.0
- `survey` 4.2.1
- `haven` 2.5.4
- `dplyr` 1.1.4
- `tidyr` 1.3.0
- `ggplot2` 4.0.3
- `patchwork` 1.3.2
- `scales` 1.4.0
- `ragg` 1.2.7
- `svglite` 2.1.0

The HCAP-linked score is a population research metric. It is not a clinical
diagnostic instrument or an estimate of individual dementia probability.
