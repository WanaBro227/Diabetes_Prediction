# PROJECT: DIABETES BIOINFORMATICS AND COMPUTATIONAL BIOLOGY
# SCRIPT: 02_data_cleaning.R
# PURPOSE: Clinical Data Cleaning and Preprocessing
# STUDENT REGISTRATION: 2021/ICT/34
# UNIVERSITY: UNIVERSITY OF VAVUNIYA

# 1. INITIAL CONFIGURATION

cat("\n")
cat("============================================================\n")
cat("       CLINICAL DATA CLEANING AND PREPROCESSING\n")
cat("============================================================\n")

# Confirm working directory
cat("\nProject directory:", getwd(), "\n")

# Check input file
input_file <- "data/01_raw_diabetes.rds"

if (!file.exists(input_file)) {
  
  stop(
    paste(
      "Raw dataset not found.",
      "Execute code/01_data_import.R first."
    )
  )
  
}

# 2. REQUIRED PACKAGES

required_packages <- c(
  "tidyverse",
  "caret"
)

for (package in required_packages) {
  
  if (!requireNamespace(package, quietly = TRUE)) {
    
    cat("\nInstalling package:", package, "\n")
    
    install.packages(
      package,
      repos = "https://cloud.r-project.org"
    )
    
  }
  
}

library(tidyverse)
library(caret)

# 3. LOAD RAW DATASET

cat("\n[STEP 1] Loading raw clinical dataset...\n")

diabetes_raw <- readRDS(input_file)

cat("\nDataset loaded successfully.\n")

cat("Observations:", nrow(diabetes_raw), "\n")
cat("Variables:", ncol(diabetes_raw), "\n")

# 4. VERIFY DATA STRUCTURE

expected_columns <- c(
  "Pregnancies",
  "Glucose",
  "BloodPressure",
  "SkinThickness",
  "Insulin",
  "BMI",
  "DiabetesPedigreeFunction",
  "Age",
  "Outcome"
)

if (!all(expected_columns %in% names(diabetes_raw))) {
  
  stop("The imported dataset has an unexpected column structure.")
  
}

# Convert historical package outcome labels when necessary.
# Positive outcome = 1; negative outcome = 0.

if (
  is.character(diabetes_raw$Outcome) ||
  is.factor(diabetes_raw$Outcome)
) {
  
  outcome_labels <- tolower(
    as.character(diabetes_raw$Outcome)
  )
  
  diabetes_raw$Outcome <- ifelse(
    outcome_labels %in% c("pos", "positive", "1"),
    1,
    ifelse(
      outcome_labels %in% c("neg", "negative", "0"),
      0,
      NA_real_
    )
  )
  
}

diabetes_raw$Outcome <- as.integer(
  diabetes_raw$Outcome
)

# 5. DEFINE CLINICAL VARIABLES

clinical_variables <- c(
  "Pregnancies",
  "Glucose",
  "BloodPressure",
  "SkinThickness",
  "Insulin",
  "BMI",
  "DiabetesPedigreeFunction",
  "Age"
)

biologically_impossible_zero <- c(
  "Glucose",
  "BloodPressure",
  "SkinThickness",
  "Insulin",
  "BMI"
)

target_variable <- "Outcome"

# 6. INITIAL DATA QUALITY ASSESSMENT

cat("\n[STEP 2] Initial data quality assessment\n")

initial_missing <- colSums(
  is.na(diabetes_raw)
)

initial_zero_count <- sapply(
  biologically_impossible_zero,
  function(variable) {
    sum(
      diabetes_raw[[variable]] == 0,
      na.rm = TRUE
    )
  }
)

initial_duplicates <- sum(
  duplicated(diabetes_raw)
)

initial_quality <- data.frame(
  Variable = names(diabetes_raw),
  Missing_Values = initial_missing,
  Zero_Values = sapply(
    diabetes_raw,
    function(x) sum(x == 0, na.rm = TRUE)
  )
)

cat("\nInitial missing-value summary:\n")
print(initial_quality)

cat("\nBiologically impossible zero measurements:\n")
print(initial_zero_count)

cat("\nDuplicate observations:", initial_duplicates, "\n")

# 7. HANDLE BIOLOGICALLY IMPOSSIBLE ZERO VALUES

cat("\n[STEP 3] Replacing impossible zero values with NA...\n")

diabetes_clean <- diabetes_raw

for (variable in biologically_impossible_zero) {
  
  zero_count <- sum(
    diabetes_clean[[variable]] == 0,
    na.rm = TRUE
  )
  
  diabetes_clean[[variable]][
    diabetes_clean[[variable]] == 0
  ] <- NA
  
  cat(
    variable,
    ":",
    zero_count,
    "zero values converted to NA\n"
  )
  
}

# 8. MISSING VALUE ASSESSMENT

cat("\n[STEP 4] Missing-value assessment after zero correction\n")

missing_after_zero <- colSums(
  is.na(diabetes_clean)
)

print(
  data.frame(
    Variable = names(diabetes_clean),
    Missing_Values = missing_after_zero
  )
)

# Verify that no entire clinical variable is missing
for (variable in clinical_variables) {
  
  if (all(is.na(diabetes_clean[[variable]]))) {
    
    stop(
      paste("All values are missing in:", variable)
    )
    
  }
  
}

# The outcome must be available for supervised classification.
if (any(is.na(diabetes_clean$Outcome))) {
  
  stop(
    "Missing outcome values detected. Resolve them before analysis."
  )
  
}

# 9. MEDIAN IMPUTATION

cat("\n[STEP 5] Performing median imputation...\n")

median_imputed <- diabetes_clean

imputation_summary <- data.frame(
  Variable = character(),
  Missing_Before = integer(),
  Median_Value = numeric(),
  Values_Imputed = integer()
)

for (variable in clinical_variables) {
  
  missing_count <- sum(
    is.na(median_imputed[[variable]])
  )
  
  median_value <- median(
    median_imputed[[variable]],
    na.rm = TRUE
  )
  
  if (!is.finite(median_value)) {
    stop(paste("Invalid median for:", variable))
  }
  
  median_imputed[[variable]][
    is.na(median_imputed[[variable]])
  ] <- median_value
  
  imputation_summary <- rbind(
    imputation_summary,
    data.frame(
      Variable = variable,
      Missing_Before = missing_count,
      Median_Value = median_value,
      Values_Imputed = missing_count
    )
  )
  
}

cat("\nMedian imputation summary:\n")
print(imputation_summary)

# 10. K-NEAREST-NEIGHBOUR IMPUTATION

cat("\n[STEP 6] Performing k-NN imputation...\n")

# Retain the zero-corrected data as input so that the
# k-NN procedure uses observed clinical measurements.

knn_input <- diabetes_clean %>%
  select(all_of(clinical_variables))

# k-NN preprocessing
# k = 5 nearest observations
# Centering and scaling are applied to the predictors.

knn_preprocessor <- caret::preProcess(
  knn_input,
  method = c(
    "center",
    "scale",
    "knnImpute"
  ),
  k = 5
)

knn_imputed_scaled <- predict(
  knn_preprocessor,
  newdata = knn_input
)

# Verify imputation results
if (anyNA(knn_imputed_scaled)) {
  
  stop(
    "k-NN preprocessing produced missing values."
  )
  
}

cat("\nk-NN imputation completed.\n")
cat("Number of neighbours: 5\n")
cat("Variables processed:", ncol(knn_imputed_scaled), "\n")

# 11. STANDARDIZATION

cat("\n[STEP 7] Generating standardized variables...\n")

# The k-NN preprocessor returns centered and scaled
# predictor measurements.

standardized_data <- as.data.frame(
  knn_imputed_scaled
)

names(standardized_data) <- paste0(
  names(standardized_data),
  "_scaled"
)

# 12. MIN-MAX NORMALIZATION

cat("\n[STEP 8] Generating normalized variables...\n")

normalized_data <- median_imputed %>%
  select(all_of(clinical_variables))

normalized_data <- as.data.frame(
  lapply(
    normalized_data,
    function(x) {
      
      minimum <- min(x, na.rm = TRUE)
      maximum <- max(x, na.rm = TRUE)
      
      if (maximum == minimum) {
        return(rep(0, length(x)))
      }
      
      (x - minimum) / (maximum - minimum)
      
    }
  )
)

names(normalized_data) <- paste0(
  names(normalized_data),
  "_normalized"
)

# 13. COMBINE PROCESSED DATA

cat("\n[STEP 9] Combining processed clinical variables...\n")

# Original clinical values are retained in their
# median-imputed form.

diabetes_clean <- median_imputed %>%
  bind_cols(standardized_data) %>%
  bind_cols(normalized_data)

# 14. FINAL DATA QUALITY ASSESSMENT

cat("\n[STEP 10] Final data quality assessment\n")

final_missing <- colSums(
  is.na(diabetes_clean)
)

final_duplicates <- sum(
  duplicated(
    diabetes_clean %>%
      select(all_of(expected_columns))
  )
)

quality_comparison <- data.frame(
  Measurement = c(
    "Total observations",
    "Total variables",
    "Missing values",
    "Impossible zero measurements",
    "Duplicate observations"
  ),
  
  Before_Cleaning = c(
    nrow(diabetes_raw),
    ncol(diabetes_raw),
    sum(initial_missing),
    sum(initial_zero_count),
    initial_duplicates
  ),
  
  After_Cleaning = c(
    nrow(diabetes_clean),
    ncol(diabetes_clean),
    sum(final_missing),
    0,
    final_duplicates
  )
)

cat("\n============================================\n")
cat("        DATA QUALITY COMPARISON\n")
cat("============================================\n")

print(quality_comparison)

cat("\nRemaining missing values by variable:\n")

print(
  data.frame(
    Variable = names(diabetes_clean),
    Missing_Values = final_missing
  )
)

if (any(final_missing > 0)) {
  
  stop("Missing values remain in the processed dataset.")
  
}

# ------------------------------------------------------------
# 15. PROCESSED DATASET SUMMARY
# ------------------------------------------------------------

cat("\n[STEP 11] Processed dataset summary\n")

cat("\nOriginal clinical variables:\n")

print(
  summary(
    diabetes_clean %>%
      select(all_of(expected_columns))
  )
)

cat("\nOutcome distribution:\n")

print(
  table(
    diabetes_clean$Outcome
  )
)

cat("\nTotal standardized variables:",
    length(clinical_variables), "\n")

cat("Total normalized variables:",
    length(clinical_variables), "\n")

# 16. SAVE PROCESSED DATASET

cat("\n[STEP 12] Saving processed dataset...\n")

saveRDS(
  diabetes_clean,
  file = "data/02_clean_diabetes.rds"
)

# Save preprocessing information for reproducibility
saveRDS(
  knn_preprocessor,
  file = "data/knn_preprocessor.rds"
)

saveRDS(
  imputation_summary,
  file = "data/imputation_summary.rds"
)

# 17. VERIFY OUTPUT FILES

stopifnot(
  file.exists("data/02_clean_diabetes.rds"),
  file.exists("data/knn_preprocessor.rds"),
  file.exists("data/imputation_summary.rds")
)

verification_data <- readRDS(
  "data/02_clean_diabetes.rds"
)

stopifnot(
  nrow(verification_data) == 768,
  !anyNA(verification_data)
)

# 18. FINAL SUMMARY

cat("\n")
cat("============================================================\n")
cat("       DATA CLEANING COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nFinal observations:", nrow(diabetes_clean), "\n")
cat("Final variables:", ncol(diabetes_clean), "\n")

cat("\nGenerated files:\n")

cat("1. data/02_clean_diabetes.rds\n")
cat("2. data/knn_preprocessor.rds\n")
cat("3. data/imputation_summary.rds\n")

cat("\nPreprocessing operations completed:\n")

cat("- Biological zero-value correction\n")
cat("- Median imputation\n")
cat("- k-NN imputation (k = 5)\n")
cat("- Z-score standardization\n")
cat("- Min-Max normalization\n")
cat("- Data quality verification\n")

cat("\nNext stage:\n")
cat("code/03_exploratory_analysis.R\n")

cat("\n============================================================\n")