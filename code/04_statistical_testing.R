
# ============================================================
# IT4133 - Bioinformatics and Computational Biology
# Project: Diabetes Prediction
# Script: 04_statistical_testing.R
# Purpose: Statistical Testing and Group Comparisons
# ============================================================

# ------------------------------------------------------------
# 1. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "tidyverse"
)

for (package_name in required_packages) {
  
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name)
  }
  
  library(
    package_name,
    character.only = TRUE
  )
}

# ------------------------------------------------------------
# 2. DEFINE PROJECT DIRECTORIES
# ------------------------------------------------------------

data_dir <- "data"
figures_dir <- "figures"

if (!dir.exists(data_dir)) {
  stop("The data directory does not exist.")
}

if (!dir.exists(figures_dir)) {
  dir.create(
    figures_dir,
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# 3. LOAD CLEANED DATASET
# ------------------------------------------------------------

input_file <- file.path(
  data_dir,
  "02_clean_diabetes.rds"
)

if (!file.exists(input_file)) {
  stop(
    "Cleaned dataset not found. ",
    "Run 02_data_cleaning.R first."
  )
}

diabetes_data <- readRDS(input_file)

cat("\n============================================\n")
cat("STATISTICAL TESTING - DIABETES DATASET\n")
cat("============================================\n")

cat(
  "\nDataset dimensions:",
  nrow(diabetes_data),
  "rows and",
  ncol(diabetes_data),
  "columns\n"
)

# ------------------------------------------------------------
# 4. DEFINE CLINICAL VARIABLES
# ------------------------------------------------------------

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

outcome_variable <- "Outcome"

required_columns <- c(
  clinical_variables,
  outcome_variable
)

missing_columns <- setdiff(
  required_columns,
  names(diabetes_data)
)

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Retain only the variables needed for this analysis.
analysis_data <- diabetes_data %>%
  select(all_of(required_columns)) %>%
  mutate(
    Outcome = factor(
      Outcome,
      levels = c(0, 1),
      labels = c("No Diabetes", "Diabetes")
    )
  )

# Check missing values.
if (anyNA(analysis_data)) {
  stop(
    "Missing values were detected in the analysis dataset. ",
    "Check 02_data_cleaning.R."
  )
}

if (nlevels(droplevels(analysis_data$Outcome)) != 2) {
  stop("Both outcome groups must be present.")
}

cat("\nOutcome group distribution:\n")

print(
  table(analysis_data$Outcome)
)

cat("\nClinical variables included:\n")

print(clinical_variables)

# ------------------------------------------------------------
# 5. HELPER FUNCTIONS
# ------------------------------------------------------------

# Safely perform statistical tests.
safe_test <- function(test_function) {
  
  tryCatch(
    test_function(),
    error = function(e) {
      warning(
        "A statistical test could not be completed: ",
        conditionMessage(e)
      )
      NULL
    }
  )
}

# Safely calculate Cohen's d.
calculate_cohens_d <- function(data, variable) {
  
  group_no <- data %>%
    filter(Outcome == "No Diabetes") %>%
    pull(all_of(variable))
  
  group_yes <- data %>%
    filter(Outcome == "Diabetes") %>%
    pull(all_of(variable))
  
  n_no <- length(group_no)
  n_yes <- length(group_yes)
  
  if (n_no < 2 || n_yes < 2) {
    return(NA_real_)
  }
  
  variance_no <- var(group_no)
  variance_yes <- var(group_yes)
  
  pooled_variance <- (
    (n_no - 1) * variance_no +
      (n_yes - 1) * variance_yes
  ) / (n_no + n_yes - 2)
  
  if (is.na(pooled_variance) || pooled_variance <= 0) {
    return(NA_real_)
  }
  
  pooled_sd <- sqrt(pooled_variance)
  
  mean_no <- mean(group_no)
  mean_yes <- mean(group_yes)
  
  # Positive values indicate a higher mean in the
  # Diabetes group.
  (mean_yes - mean_no) / pooled_sd
}

# Round numeric columns without affecting text columns.
round_numeric_columns <- function(data, digits = 4) {
  
  data %>%
    mutate(
      across(
        where(is.numeric),
        ~ round(.x, digits)
      )
    )
}

# ------------------------------------------------------------
# 6. SHAPIRO-WILK NORMALITY TESTS
# ------------------------------------------------------------

cat("\n--------------------------------------------\n")
cat("6. SHAPIRO-WILK NORMALITY TESTS\n")
cat("--------------------------------------------\n")

normality_results <- map_dfr(
  clinical_variables,
  function(variable) {
    
    values <- analysis_data[[variable]]
    
    test_result <- safe_test(
      function() {
        shapiro.test(values)
      }
    )
    
    if (is.null(test_result)) {
      
      return(
        tibble(
          Variable = variable,
          W_Statistic = NA_real_,
          P_Value = NA_real_
        )
      )
    }
    
    tibble(
      Variable = variable,
      W_Statistic = unname(test_result$statistic),
      P_Value = test_result$p.value
    )
  }
) %>%
  mutate(
    Adjusted_P_Value = p.adjust(
      P_Value,
      method = "BH"
    ),
    Normality_At_0_05 = case_when(
      is.na(Adjusted_P_Value) ~ "Not evaluated",
      Adjusted_P_Value < 0.05 ~ "Evidence against normality",
      TRUE ~ "No significant evidence against normality"
    )
  )

normality_results <- round_numeric_columns(
  normality_results
)

print(normality_results)

# ------------------------------------------------------------
# 7. WELCH TWO-SAMPLE T-TESTS
# ------------------------------------------------------------

cat("\n--------------------------------------------\n")
cat("7. WELCH TWO-SAMPLE T-TESTS\n")
cat("--------------------------------------------\n")

welch_results <- map_dfr(
  clinical_variables,
  function(variable) {
    
    test_result <- safe_test(
      function() {
        t.test(
          formula = as.formula(
            paste(variable, "~ Outcome")
          ),
          data = analysis_data,
          var.equal = FALSE
        )
      }
    )
    
    if (is.null(test_result)) {
      
      return(
        tibble(
          Variable = variable,
          Mean_No_Diabetes = NA_real_,
          Mean_Diabetes = NA_real_,
          Mean_Difference = NA_real_,
          T_Statistic = NA_real_,
          Degrees_of_Freedom = NA_real_,
          P_Value = NA_real_,
          Cohen_d = NA_real_
        )
      )
    }
    
    group_means <- analysis_data %>%
      group_by(Outcome) %>%
      summarise(
        Mean = mean(.data[[variable]]),
        .groups = "drop"
      )
    
    mean_no <- group_means %>%
      filter(Outcome == "No Diabetes") %>%
      pull(Mean)
    
    mean_yes <- group_means %>%
      filter(Outcome == "Diabetes") %>%
      pull(Mean)
    
    tibble(
      Variable = variable,
      Mean_No_Diabetes = mean_no,
      Mean_Diabetes = mean_yes,
      Mean_Difference = mean_yes - mean_no,
      T_Statistic = unname(test_result$statistic),
      Degrees_of_Freedom = unname(test_result$parameter),
      P_Value = test_result$p.value,
      Cohen_d = calculate_cohens_d(
        analysis_data,
        variable
      )
    )
  }
) %>%
  mutate(
    Adjusted_P_Value = p.adjust(
      P_Value,
      method = "BH"
    ),
    Significant = case_when(
      is.na(Adjusted_P_Value) ~ NA,
      Adjusted_P_Value < 0.05 ~ TRUE,
      TRUE ~ FALSE
    )
  )

welch_results <- round_numeric_columns(
  welch_results
)

print(welch_results)

# ------------------------------------------------------------
# 8. WILCOXON RANK-SUM TESTS
# ------------------------------------------------------------

cat("\n--------------------------------------------\n")
cat("8. WILCOXON RANK-SUM TESTS\n")
cat("--------------------------------------------\n")

wilcoxon_results <- map_dfr(
  clinical_variables,
  function(variable) {
    
    test_result <- safe_test(
      function() {
        wilcox.test(
          formula = as.formula(
            paste(variable, "~ Outcome")
          ),
          data = analysis_data,
          exact = FALSE,
          conf.int = FALSE
        )
      }
    )
    
    if (is.null(test_result)) {
      
      return(
        tibble(
          Variable = variable,
          W_Statistic = NA_real_,
          P_Value = NA_real_
        )
      )
    }
    
    tibble(
      Variable = variable,
      W_Statistic = unname(test_result$statistic),
      P_Value = test_result$p.value
    )
  }
) %>%
  mutate(
    Adjusted_P_Value = p.adjust(
      P_Value,
      method = "BH"
    ),
    Significant = case_when(
      is.na(Adjusted_P_Value) ~ NA,
      Adjusted_P_Value < 0.05 ~ TRUE,
      TRUE ~ FALSE
    )
  )

wilcoxon_results <- round_numeric_columns(
  wilcoxon_results
)

print(wilcoxon_results)

# ------------------------------------------------------------
# 9. FLIGNER-KILLEEN VARIANCE TESTS
# ------------------------------------------------------------

cat("\n--------------------------------------------\n")
cat("9. FLIGNER-KILLEEN VARIANCE TESTS\n")
cat("--------------------------------------------\n")

variance_results <- map_dfr(
  clinical_variables,
  function(variable) {
    
    test_result <- safe_test(
      function() {
        fligner.test(
          formula = as.formula(
            paste(variable, "~ Outcome")
          ),
          data = analysis_data
        )
      }
    )
    
    if (is.null(test_result)) {
      
      return(
        tibble(
          Variable = variable,
          Chi_Square = NA_real_,
          Degrees_of_Freedom = NA_real_,
          P_Value = NA_real_
        )
      )
    }
    
    tibble(
      Variable = variable,
      Chi_Square = unname(test_result$statistic),
      Degrees_of_Freedom = unname(test_result$parameter),
      P_Value = test_result$p.value
    )
  }
) %>%
  mutate(
    Adjusted_P_Value = p.adjust(
      P_Value,
      method = "BH"
    ),
    Variance_Difference = case_when(
      is.na(Adjusted_P_Value) ~ "Not evaluated",
      Adjusted_P_Value < 0.05 ~ "Evidence of unequal variances",
      TRUE ~ "No significant evidence of unequal variances"
    )
  )

variance_results <- round_numeric_columns(
  variance_results
)

print(variance_results)

# ------------------------------------------------------------
# 10. PEARSON AND SPEARMAN CORRELATION TESTS
# ------------------------------------------------------------

cat("\n--------------------------------------------\n")
cat("10. PEARSON AND SPEARMAN CORRELATIONS\n")
cat("--------------------------------------------\n")

# Convert the binary outcome into numeric values:
# No Diabetes = 0
# Diabetes = 1

numeric_outcome <- as.numeric(
  analysis_data$Outcome
) - 1

correlation_results <- map_dfr(
  clinical_variables,
  function(variable) {
    
    values <- analysis_data[[variable]]
    
    pearson_result <- safe_test(
      function() {
        cor.test(
          values,
          numeric_outcome,
          method = "pearson"
        )
      }
    )
    
    spearman_result <- safe_test(
      function() {
        cor.test(
          values,
          numeric_outcome,
          method = "spearman",
          exact = FALSE
        )
      }
    )
    
    tibble(
      Variable = variable,
      
      Pearson_Correlation = if (
        is.null(pearson_result)
      ) {
        NA_real_
      } else {
        unname(pearson_result$estimate)
      },
      
      Pearson_P_Value = if (
        is.null(pearson_result)
      ) {
        NA_real_
      } else {
        pearson_result$p.value
      },
      
      Spearman_Correlation = if (
        is.null(spearman_result)
      ) {
        NA_real_
      } else {
        unname(spearman_result$estimate)
      },
      
      Spearman_P_Value = if (
        is.null(spearman_result)
      ) {
        NA_real_
      } else {
        spearman_result$p.value
      }
    )
  }
) %>%
  mutate(
    Pearson_Adjusted_P = p.adjust(
      Pearson_P_Value,
      method = "BH"
    ),
    
    Spearman_Adjusted_P = p.adjust(
      Spearman_P_Value,
      method = "BH"
    )
  )

correlation_results <- round_numeric_columns(
  correlation_results
)

print(correlation_results)

# ------------------------------------------------------------
# 11. PREPARE GROUP MEANS FOR VISUALIZATION
# ------------------------------------------------------------

group_means <- analysis_data %>%
  pivot_longer(
    cols = all_of(clinical_variables),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  group_by(
    Variable,
    Outcome
  ) %>%
  summarise(
    Mean = mean(Value),
    Standard_Deviation = sd(Value),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 12. FIGURE 10 - GROUP MEAN COMPARISON
# ------------------------------------------------------------

cat("\n--------------------------------------------\n")
cat("12. GENERATING STATISTICAL FIGURES\n")
cat("--------------------------------------------\n")

figure_10 <- ggplot(
  group_means,
  aes(
    x = Outcome,
    y = Mean,
    fill = Outcome
  )
) +
  
  geom_col(
    width = 0.65,
    position = "dodge"
  ) +
  
  facet_wrap(
    ~ Variable,
    scales = "free_y",
    ncol = 2
  ) +
  
  labs(
    title = "Clinical Variable Mean Comparisons",
    subtitle = "Comparison between diabetes outcome groups",
    x = "Diabetes Outcome",
    y = "Mean Value",
    fill = "Outcome"
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    strip.text = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 15,
      hjust = 1
    ),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    figures_dir,
    "10_statistical_comparison.png"
  ),
  plot = figure_10,
  width = 12,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------
# 13. FIGURE 11 - ADJUSTED P-VALUE COMPARISON
# ------------------------------------------------------------

plot_welch <- welch_results %>%
  mutate(
    Plot_P_Value = pmax(
      Adjusted_P_Value,
      .Machine$double.xmin
    ),
    Minus_Log10_P = -log10(Plot_P_Value),
    Variable = factor(
      Variable,
      levels = rev(clinical_variables)
    )
  )

figure_11 <- ggplot(
  plot_welch,
  aes(
    x = Minus_Log10_P,
    y = Variable
  )
) +
  
  geom_col(
    fill = "steelblue",
    width = 0.7
  ) +
  
  geom_vline(
    xintercept = -log10(0.05),
    linetype = "dashed",
    color = "red",
    linewidth = 0.8
  ) +
  
  labs(
    title = "Welch's T-Test: Adjusted P-Values",
    subtitle = "Dashed line represents the 0.05 significance threshold",
    x = expression(-log[10]("Adjusted P-value")),
    y = "Clinical Variable"
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

ggsave(
  filename = file.path(
    figures_dir,
    "11_pvalue_comparison.png"
  ),
  plot = figure_11,
  width = 10,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 14. FIGURE 12 - CORRELATION COMPARISON
# ------------------------------------------------------------

correlation_plot_data <- correlation_results %>%
  select(
    Variable,
    Pearson_Correlation,
    Spearman_Correlation
  ) %>%
  pivot_longer(
    cols = c(
      Pearson_Correlation,
      Spearman_Correlation
    ),
    names_to = "Correlation_Method",
    values_to = "Correlation"
  ) %>%
  mutate(
    Correlation_Method = recode(
      Correlation_Method,
      Pearson_Correlation = "Pearson",
      Spearman_Correlation = "Spearman"
    ),
    Variable = factor(
      Variable,
      levels = clinical_variables
    )
  )

figure_12 <- ggplot(
  correlation_plot_data,
  aes(
    x = Variable,
    y = Correlation,
    fill = Correlation_Method
  )
) +
  
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.5
  ) +
  
  labs(
    title = "Clinical Variables and Diabetes Outcome",
    subtitle = "Pearson and Spearman correlation coefficients",
    x = "Clinical Variable",
    y = "Correlation Coefficient",
    fill = "Method"
  ) +
  
  coord_cartesian(
    ylim = c(-1, 1)
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    figures_dir,
    "12_correlation_comparison.png"
  ),
  plot = figure_12,
  width = 11,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 15. SAVE STATISTICAL RESULTS AS CSV FILES
# ------------------------------------------------------------

write_csv(
  normality_results,
  file.path(
    data_dir,
    "04_normality_tests.csv"
  )
)

write_csv(
  welch_results,
  file.path(
    data_dir,
    "04_welch_t_tests.csv"
  )
)

write_csv(
  wilcoxon_results,
  file.path(
    data_dir,
    "04_wilcoxon_tests.csv"
  )
)

write_csv(
  variance_results,
  file.path(
    data_dir,
    "04_variance_tests.csv"
  )
)

write_csv(
  correlation_results,
  file.path(
    data_dir,
    "04_correlation_tests.csv"
  )
)

# ------------------------------------------------------------
# 16. SAVE ALL RESULTS AS A SINGLE RDS OBJECT
# ------------------------------------------------------------

statistical_results <- list(
  
  normality_tests = normality_results,
  
  welch_t_tests = welch_results,
  
  wilcoxon_tests = wilcoxon_results,
  
  variance_tests = variance_results,
  
  correlation_tests = correlation_results,
  
  group_means = group_means,
  
  analysis_variables = clinical_variables,
  
  outcome_variable = outcome_variable
  
)

saveRDS(
  statistical_results,
  file.path(
    data_dir,
    "04_statistical_results.rds"
  )
)

# ------------------------------------------------------------
# 17. VERIFY OUTPUT FILES
# ------------------------------------------------------------

expected_outputs <- c(
  
  "data/04_normality_tests.csv",
  
  "data/04_welch_t_tests.csv",
  
  "data/04_wilcoxon_tests.csv",
  
  "data/04_variance_tests.csv",
  
  "data/04_correlation_tests.csv",
  
  "data/04_statistical_results.rds",
  
  "figures/10_statistical_comparison.png",
  
  "figures/11_pvalue_comparison.png",
  
  "figures/12_correlation_comparison.png"
  
)

output_status <- tibble(
  File = expected_outputs,
  Exists = file.exists(expected_outputs)
)

cat("\n============================================\n")
cat("OUTPUT FILE VERIFICATION\n")
cat("============================================\n")

print(output_status)

if (!all(output_status$Exists)) {
  stop(
    "One or more output files were not generated. ",
    "Check the errors above."
  )
}

# ------------------------------------------------------------
# 18. FINAL SUMMARY
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STATISTICAL TESTING COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat("\nStatistical methods completed:\n")

cat("1. Shapiro-Wilk normality tests\n")
cat("2. Welch two-sample t-tests\n")
cat("3. Wilcoxon rank-sum tests\n")
cat("4. Fligner-Killeen variance tests\n")
cat("5. Pearson correlation tests\n")
cat("6. Spearman correlation tests\n")
cat("7. Benjamini-Hochberg adjustments\n")
cat("8. Cohen's d effect size calculations\n")

cat("\nGenerated figures:\n")

cat("figures/10_statistical_comparison.png\n")
cat("figures/11_pvalue_comparison.png\n")
cat("figures/12_correlation_comparison.png\n")

cat("\nGenerated statistical result files:\n")

cat("data/04_normality_tests.csv\n")
cat("data/04_welch_t_tests.csv\n")
cat("data/04_wilcoxon_tests.csv\n")
cat("data/04_variance_tests.csv\n")
cat("data/04_correlation_tests.csv\n")
cat("data/04_statistical_results.rds\n")

cat("\nNext step:\n")
cat("Run code/05_predictive_modelling.R\n")

cat("\n============================================\n")
