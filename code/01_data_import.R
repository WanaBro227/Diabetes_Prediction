
# PROJECT: DIABETES BIOINFORMATICS AND COMPUTATIONAL BIOLOGY
# SCRIPT: 01_data_import.R
# PURPOSE: Clinical Dataset Acquisition and Initial Inspection
# STUDENT REGISTRATION: 2021/ICT/34
# UNIVERSITY: UNIVERSITY OF VAVUNIYA

# 1. INITIAL CONFIGURATION

cat("\n")
cat("============================================================\n")
cat("       DIABETES CLINICAL DATA ACQUISITION\n")
cat("============================================================\n")

# Confirm project root
cat("\nProject working directory:\n")
print(getwd())

# Create required directories if they do not exist
required_dirs <- c("data", "code", "figures")

for (folder in required_dirs) {
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
}

# 2. REQUIRED PACKAGES

required_packages <- c("tidyverse", "mlbench")

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

# 3. DATASET SOURCE CONFIGURATION

# Preferred public dataset source
source_urls <- c(
  "https://raw.githubusercontent.com/selva86/datasets/master/PimaIndiansDiabetes.csv",
  "https://raw.githubusercontent.com/plotly/datasets/master/diabetes.csv"
)

dataset_source <- NULL
diabetes_raw <- NULL

# Expected clinical variables
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

# 4. DOWNLOAD AND VALIDATE DATASET

cat("\n[STEP 1] Acquiring the clinical dataset...\n")

for (url in source_urls) {
  
  cat("\nAttempting source:\n")
  cat(url, "\n")
  
  temporary_file <- tempfile(fileext = ".csv")
  
  download_result <- tryCatch({
    
    download.file(
      url = url,
      destfile = temporary_file,
      mode = "wb",
      quiet = TRUE
    )
    
    TRUE
    
  }, error = function(e) {
    FALSE
  })
  
  if (!download_result || !file.exists(temporary_file)) {
    cat("Download unsuccessful. Trying next source.\n")
    next
  }
  
  candidate <- tryCatch({
    
    read.csv(
      temporary_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
  }, error = function(e) {
    NULL
  })
  
  if (is.null(candidate)) {
    cat("Could not read this dataset. Trying next source.\n")
    next
  }
  
  # Normalize column names for validation
  normalized_names <- tolower(
    gsub("[^A-Za-z]", "", names(candidate))
  )
  
  expected_names <- tolower(
    gsub("[^A-Za-z]", "", expected_columns)
  )
  
  if (
    nrow(candidate) == 768 &&
    all(expected_names %in% normalized_names)
  ) {
    
    diabetes_raw <- candidate
    dataset_source <- url
    
    cat("\nDataset successfully downloaded and validated.\n")
    
    break
    
  } else {
    
    cat("Dataset integrity check failed for this source.\n")
    
  }
  
  unlink(temporary_file)
}

# 5. PACKAGE FALLBACK

# Use the mlbench package as a final alternative.
# This prevents the workflow from depending on a single
# external download location.

if (is.null(diabetes_raw)) {
  
  cat("\n[STEP 2] Activating package-based fallback...\n")
  
  if (!requireNamespace("mlbench", quietly = TRUE)) {
    
    install.packages(
      "mlbench",
      repos = "https://cloud.r-project.org"
    )
    
  }
  
  data(
    "PimaIndiansDiabetes",
    package = "mlbench",
    envir = environment()
  )
  
  diabetes_raw <- as.data.frame(
    PimaIndiansDiabetes
  )
  
  dataset_source <- paste(
    "R package: mlbench, dataset: PimaIndiansDiabetes"
  )
  
  # Convert package outcome into numeric binary values
  diabetes_raw$diabetes <- as.character(
    diabetes_raw$diabetes
  )
  
  diabetes_raw$diabetes <- ifelse(
    diabetes_raw$diabetes == "pos",
    1,
    0
  )
  
  cat("Package fallback loaded successfully.\n")
  
}

# 6. STANDARDIZE VARIABLE NAMES

cat("\n[STEP 3] Standardizing clinical variable names...\n")

# Display original downloaded column names
cat("\nOriginal downloaded column names:\n")
print(names(diabetes_raw))

# Normalize names by removing spaces and special characters
normalized_names <- tolower(
  gsub("[^A-Za-z]", "", names(diabetes_raw))
)

# Support different naming conventions used by public mirrors
column_aliases <- list(
  
  Pregnancies = c("pregnancies", "pregnant"),
  
  Glucose = c("glucose"),
  
  BloodPressure = c("bloodpressure", "pressure"),
  
  SkinThickness = c("skinthickness", "triceps"),
  
  Insulin = c("insulin"),
  
  BMI = c("bmi", "mass"),
  
  DiabetesPedigreeFunction = c(
    "diabetespedigreefunction",
    "pedigree"
  ),
  
  Age = c("age"),
  
  Outcome = c(
    "outcome",
    "diabetes",
    "class"
  )
  
)

# Identify the correct position of every required variable
column_positions <- vapply(
  column_aliases,
  function(aliases) {
    
    matched_position <- which(
      normalized_names %in% aliases
    )
    
    if (length(matched_position) == 0) {
      return(NA_integer_)
    }
    
    matched_position[1]
    
  },
  integer(1)
)

# Stop if any required variable is missing
if (anyNA(column_positions)) {
  
  missing_columns <- names(column_positions)[
    is.na(column_positions)
  ]
  
  stop(
    paste(
      "Required clinical variables not found:",
      paste(missing_columns, collapse = ", "),
      "\nAvailable columns:",
      paste(names(diabetes_raw), collapse = ", ")
    )
  )
  
}

# Reorder and rename the dataset
diabetes_raw <- diabetes_raw[
  ,
  column_positions,
  drop = FALSE
]

names(diabetes_raw) <- names(column_aliases)

# Convert all clinical measurements to numeric
diabetes_raw <- diabetes_raw %>%
  mutate(
    across(
      everything(),
      ~ suppressWarnings(
        as.numeric(as.character(.x))
      )
    )
  )

# Convert binary outcome into numeric representation
if (anyNA(diabetes_raw$Outcome)) {
  
  stop(
    "The outcome variable contains missing or invalid values."
  )
  
}

if (!all(diabetes_raw$Outcome %in% c(0, 1))) {
  
  stop(
    "The outcome variable must contain only 0 and 1."
  )
  
}

diabetes_raw$Outcome <- as.integer(
  diabetes_raw$Outcome
)

cat("\nStandardized column names:\n")
print(names(diabetes_raw))

cat("\nVariable standardization completed successfully.\n")


# 7. DATA INTEGRITY VERIFICATION

cat("\n[STEP 4] Performing data integrity verification...\n")

# Check dimensions
dataset_dimensions <- dim(diabetes_raw)

# Check expected structure
stopifnot(
  nrow(diabetes_raw) == 768,
  ncol(diabetes_raw) == 9
)

# Check duplicate observations
duplicate_count <- sum(
  duplicated(diabetes_raw)
)

# Check outcome values
stopifnot(
  all(na.omit(diabetes_raw$Outcome) %in% c(0, 1))
)

# Check numeric columns
stopifnot(
  all(sapply(diabetes_raw, is.numeric))
)

# Check missing values
missing_values <- colSums(
  is.na(diabetes_raw)
)

# 8. DATASET INFORMATION

cat("\n[STEP 5] Dataset information\n")

cat("\nDataset dimensions:\n")
cat("Observations:", nrow(diabetes_raw), "\n")
cat("Variables:", ncol(diabetes_raw), "\n")

cat("\nClinical variables:\n")

print(
  data.frame(
    Variable = names(diabetes_raw),
    Data_Type = sapply(
      diabetes_raw,
      function(x) class(x)[1]
    ),
    Missing_Values = missing_values
  )
)

cat("\nOutcome distribution:\n")

print(
  table(
    diabetes_raw$Outcome,
    useNA = "ifany"
  )
)

cat("\nDuplicate observations:", duplicate_count, "\n")

cat("\nDataset structure:\n")
str(diabetes_raw)

cat("\nStatistical summary:\n")
print(summary(diabetes_raw))

# 9. DATASET SOURCE AND CITATION

cat("\n[STEP 6] Dataset provenance\n")

cat("\nDataset: Pima Indians Diabetes Clinical Dataset\n")

cat("\nPrimary repository:\n")
cat("UCI Repository\n")

cat("\nDataset information page:\n")
cat("https://archive.ics.uci.edu/dataset/34/diabetes\n")

cat("\nDownloaded source:\n")
cat(dataset_source, "\n")

cat("\nDataset citation:\n")

cat(
  "Smith, J. W., Everhart, J. E., Dickson, W. C.,",
  "Knowler, W. C., & Johannes, R. S. (1988).",
  "Using the ADAP learning algorithm to forecast the onset",
  "of diabetes mellitus. Proceedings of the Annual Symposium",
  "on Computer Application in Medical Care, 261-265.\n"
)

cat("\nLicensing note:\n")
cat(
  "Consult the UCI dataset page for its current licensing",
  "and reuse information. Public availability alone does",
  "not establish an open license.\n"
)

# 10. SAVE RAW DATASET

cat("\n[STEP 7] Saving imported dataset...\n")

# Save CSV
write.csv(
  diabetes_raw,
  file = "data/diabetes_raw.csv",
  row.names = FALSE
)

# Save RDS
saveRDS(
  diabetes_raw,
  file = "data/01_raw_diabetes.rds"
)

# 11. VERIFY SAVED FILES

stopifnot(
  file.exists("data/diabetes_raw.csv"),
  file.exists("data/01_raw_diabetes.rds")
)

# Reload saved RDS to verify readability
verification_data <- readRDS(
  "data/01_raw_diabetes.rds"
)

stopifnot(
  identical(
    dim(verification_data),
    dim(diabetes_raw)
  )
)

# 12. FINAL SUMMARY

cat("\n")
cat("============================================================\n")
cat("          DATA IMPORT COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nTotal observations:", nrow(diabetes_raw), "\n")
cat("Total variables:", ncol(diabetes_raw), "\n")

cat("\nGenerated files:\n")
cat("1. data/diabetes_raw.csv\n")
cat("2. data/01_raw_diabetes.rds\n")

cat("\nNext stage:\n")
cat("code/02_data_cleaning.R\n")

cat("\n============================================================\n")
