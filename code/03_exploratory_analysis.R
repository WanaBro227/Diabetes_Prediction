
# ============================================================
# PROJECT: DIABETES BIOINFORMATICS AND COMPUTATIONAL BIOLOGY
# SCRIPT: 03_exploratory_analysis.R
# PURPOSE: Exploratory Data Analysis of Clinical Measurements
# DATASET: Pima Indians Diabetes Dataset
# ============================================================


# ------------------------------------------------------------
# 1. PROJECT CONFIGURATION
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("       EXPLORATORY DATA ANALYSIS\n")
cat("============================================================\n")

cat("\nProject directory:", getwd(), "\n")

input_file <- "data/02_clean_diabetes.rds"

if (!file.exists(input_file)) {
  
  stop(
    paste(
      "Cleaned dataset not found.",
      "Execute code/02_data_cleaning.R first."
    )
  )
  
}

# Create figures directory if it does not exist.

if (!dir.exists("figures")) {
  
  dir.create("figures", recursive = TRUE)
  
}


# ------------------------------------------------------------
# 2. REQUIRED PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "ggplot2",
  "corrplot"
)

for (package in required_packages) {
  
  if (!requireNamespace(package, quietly = TRUE)) {
    
    cat("\nInstalling package:", package, "\n")
    
    install.packages(
      package,
      dependencies = TRUE
    )
    
  }
  
}

library(tidyverse)
library(ggplot2)
library(corrplot)


# ------------------------------------------------------------
# 3. LOAD CLEANED DATASET
# ------------------------------------------------------------

cat("\n[STEP 1] Loading cleaned clinical dataset...\n")

diabetes_data <- readRDS(input_file)

cat("\nDataset loaded successfully.\n")

cat("Observations:", nrow(diabetes_data), "\n")
cat("Variables:", ncol(diabetes_data), "\n")

if (anyNA(diabetes_data)) {
  
  stop(
    "Missing values detected in the cleaned dataset."
  )
  
}


# ------------------------------------------------------------
# 4. DEFINE VARIABLES
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

target_variable <- "Outcome"

if (!all(clinical_variables %in% names(diabetes_data))) {
  
  stop(
    "One or more required clinical variables are missing."
  )
  
}

if (!target_variable %in% names(diabetes_data)) {
  
  stop("Outcome variable is missing.")
  
}

# Use the median-imputed clinical measurements for
# descriptive statistics and exploratory visualizations.

eda_data <- diabetes_data %>%
  select(
    all_of(clinical_variables),
    Outcome
  ) %>%
  mutate(
    Outcome = factor(
      Outcome,
      levels = c(0, 1),
      labels = c(
        "No Diabetes",
        "Diabetes"
      )
    )
  )

cat("\nClinical variables selected:", length(clinical_variables), "\n")


# ------------------------------------------------------------
# 5. DESCRIPTIVE STATISTICAL ANALYSIS
# ------------------------------------------------------------

cat("\n[STEP 2] Calculating descriptive statistics...\n")

descriptive_statistics <- data.frame(
  Variable = clinical_variables,
  
  Mean = sapply(
    eda_data[clinical_variables],
    mean
  ),
  
  Median = sapply(
    eda_data[clinical_variables],
    median
  ),
  
  Standard_Deviation = sapply(
    eda_data[clinical_variables],
    sd
  ),
  
  Minimum = sapply(
    eda_data[clinical_variables],
    min
  ),
  
  Maximum = sapply(
    eda_data[clinical_variables],
    max
  ),
  
  First_Quartile = sapply(
    eda_data[clinical_variables],
    function(x) quantile(x, 0.25)
  ),
  
  Third_Quartile = sapply(
    eda_data[clinical_variables],
    function(x) quantile(x, 0.75)
  )
)

cat("\nDescriptive statistics:\n")

descriptive_statistics_display <- descriptive_statistics

numeric_columns <- sapply(
  descriptive_statistics_display,
  is.numeric
)

descriptive_statistics_display[numeric_columns] <-
  lapply(
    descriptive_statistics_display[numeric_columns],
    round,
    digits = 3
  )

print(descriptive_statistics_display)

write.csv(
  descriptive_statistics,
  "data/descriptive_statistics.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. OUTCOME DISTRIBUTION ANALYSIS
# ------------------------------------------------------------

cat("\n[STEP 3] Analysing diabetes outcome distribution...\n")

outcome_distribution <- eda_data %>%
  count(Outcome) %>%
  mutate(
    Percentage = (n / sum(n)) * 100
  )

cat("\nDiabetes outcome distribution:\n")

print(outcome_distribution)

write.csv(
  outcome_distribution,
  "data/outcome_distribution.csv",
  row.names = FALSE
)


# Figure 1: Outcome distribution

outcome_plot <- ggplot(
  outcome_distribution,
  aes(
    x = Outcome,
    y = n,
    fill = Outcome
  )
) +
  
  geom_col(
    width = 0.65,
    color = "black"
  ) +
  
  geom_text(
    aes(
      label = paste0(
        n,
        "\n",
        round(Percentage, 1),
        "%"
      )
    ),
    vjust = -0.4,
    size = 4
  ) +
  
  scale_fill_manual(
    values = c(
      "No Diabetes" = "#3498DB",
      "Diabetes" = "#E74C3C"
    )
  ) +
  
  labs(
    title = "Distribution of Diabetes Outcomes",
    subtitle = "Number and percentage of observations",
    x = "Diabetes Outcome",
    y = "Number of Observations"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  ) +
  
  expand_limits(
    y = max(outcome_distribution$n) * 1.15
  )

ggsave(
  "figures/01_outcome_distribution.png",
  plot = outcome_plot,
  width = 8,
  height = 6,
  dpi = 300
)

print(outcome_plot)


# ------------------------------------------------------------
# 7. CLINICAL VARIABLE DISTRIBUTIONS
# ------------------------------------------------------------

cat("\n[STEP 4] Examining clinical variable distributions...\n")

# Convert the dataset into long format.

clinical_long <- eda_data %>%
  pivot_longer(
    cols = all_of(clinical_variables),
    names_to = "Clinical_Variable",
    values_to = "Measurement"
  )

# Figure 2: Histograms

histogram_plot <- ggplot(
  clinical_long,
  aes(
    x = Measurement
  )
) +
  
  geom_histogram(
    bins = 30,
    fill = "#3498DB",
    color = "white",
    alpha = 0.85
  ) +
  
  facet_wrap(
    ~ Clinical_Variable,
    scales = "free",
    ncol = 2
  ) +
  
  labs(
    title = "Distribution of Clinical Measurements",
    subtitle = "Histograms of the eight clinical variables",
    x = "Measurement Value",
    y = "Frequency"
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

ggsave(
  "figures/02_clinical_histograms.png",
  plot = histogram_plot,
  width = 12,
  height = 10,
  dpi = 300
)

print(histogram_plot)


# Figure 3: Density plots

density_plot <- ggplot(
  clinical_long,
  aes(
    x = Measurement,
    fill = Clinical_Variable
  )
) +
  
  geom_density(
    alpha = 0.65,
    color = "black"
  ) +
  
  facet_wrap(
    ~ Clinical_Variable,
    scales = "free",
    ncol = 2
  ) +
  
  labs(
    title = "Density Distributions of Clinical Variables",
    subtitle = "Smoothed distributions of clinical measurements",
    x = "Measurement Value",
    y = "Density"
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

ggsave(
  "figures/03_clinical_density_plots.png",
  plot = density_plot,
  width = 12,
  height = 10,
  dpi = 300
)

print(density_plot)


# ------------------------------------------------------------
# 8. COMPARISON OF CLINICAL VARIABLES BY OUTCOME
# ------------------------------------------------------------

cat("\n[STEP 5] Comparing clinical measurements by outcome...\n")

# Calculate group-specific descriptive statistics.

group_statistics <- eda_data %>%
  
  pivot_longer(
    cols = all_of(clinical_variables),
    names_to = "Clinical_Variable",
    values_to = "Measurement"
  ) %>%
  
  group_by(
    Clinical_Variable,
    Outcome
  ) %>%
  
  summarise(
    Count = n(),
    
    Mean = mean(Measurement),
    
    Median = median(Measurement),
    
    Standard_Deviation = sd(Measurement),
    
    .groups = "drop"
  )

cat("\nClinical measurements grouped by diabetes outcome:\n")

group_statistics_display <- group_statistics

numeric_columns <- sapply(
  group_statistics_display,
  is.numeric
)

group_statistics_display[numeric_columns] <-
  lapply(
    group_statistics_display[numeric_columns],
    round,
    digits = 3
  )

print(group_statistics_display)

write.csv(
  group_statistics,
  "data/group_statistics.csv",
  row.names = FALSE
)


# Figure 4: Boxplots

boxplot_data <- clinical_long

boxplot_plot <- ggplot(
  boxplot_data,
  aes(
    x = Outcome,
    y = Measurement,
    fill = Outcome
  )
) +
  
  geom_boxplot(
    alpha = 0.75,
    outlier.color = "black",
    outlier.alpha = 0.5
  ) +
  
  facet_wrap(
    ~ Clinical_Variable,
    scales = "free",
    ncol = 2
  ) +
  
  scale_fill_manual(
    values = c(
      "No Diabetes" = "#3498DB",
      "Diabetes" = "#E74C3C"
    )
  ) +
  
  labs(
    title = "Clinical Measurements by Diabetes Outcome",
    subtitle = "Comparison of clinical variables between outcome groups",
    x = "Diabetes Outcome",
    y = "Measurement Value"
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

ggsave(
  "figures/04_clinical_boxplots_by_outcome.png",
  plot = boxplot_plot,
  width = 12,
  height = 10,
  dpi = 300
)

print(boxplot_plot)


# ------------------------------------------------------------
# 9. CORRELATION ANALYSIS
# ------------------------------------------------------------

cat("\n[STEP 6] Calculating correlations between clinical variables...\n")

correlation_matrix <- cor(
  eda_data %>%
    select(all_of(clinical_variables)),
  method = "pearson"
)

cat("\nPearson correlation matrix:\n")

print(
  round(
    correlation_matrix,
    3
  )
)

write.csv(
  correlation_matrix,
  "data/correlation_matrix.csv"
)


# Figure 5: Correlation heatmap

png(
  filename = "figures/05_correlation_heatmap.png",
  width = 3000,
  height = 2500,
  res = 300
)

corrplot(
  correlation_matrix,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",
  number.cex = 0.7,
  tl.col = "black",
  tl.srt = 45,
  col = colorRampPalette(
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )(200),
  title = "Correlation Between Clinical Variables",
  mar = c(0, 0, 2, 0)
)

dev.off()


# ------------------------------------------------------------
# 10. CORRELATION WITH DIABETES OUTCOME
# ------------------------------------------------------------

cat("\n[STEP 7] Examining associations with diabetes outcome...\n")

# Calculate Pearson correlations between each clinical
# measurement and the binary outcome (0 or 1).

outcome_correlations <- data.frame(
  Clinical_Variable = clinical_variables,
  
  Correlation = sapply(
    clinical_variables,
    function(variable) {
      
      cor(
        eda_data[[variable]],
        as.numeric(eda_data$Outcome) - 1,
        method = "pearson"
      )
      
    }
  )
) %>%
  
  arrange(
    desc(
      abs(Correlation)
    )
  )

cat("\nCorrelation with diabetes outcome:\n")

outcome_correlations_display <- outcome_correlations

numeric_columns <- sapply(
  outcome_correlations_display,
  is.numeric
)

outcome_correlations_display[numeric_columns] <-
  lapply(
    outcome_correlations_display[numeric_columns],
    round,
    digits = 3
  )

print(outcome_correlations_display)

write.csv(
  outcome_correlations,
  "data/outcome_correlations.csv",
  row.names = FALSE
)


# Figure 6: Outcome correlation chart

outcome_correlation_plot <- ggplot(
  outcome_correlations,
  aes(
    x = reorder(
      Clinical_Variable,
      Correlation
    ),
    y = Correlation,
    fill = Correlation
  )
) +
  
  geom_col(
    width = 0.7
  ) +
  
  coord_flip() +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  
  labs(
    title = "Correlation of Clinical Variables with Diabetes Outcome",
    subtitle = "Pearson correlation coefficients",
    x = "Clinical Variable",
    y = "Correlation Coefficient"
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

ggsave(
  "figures/06_outcome_correlation.png",
  plot = outcome_correlation_plot,
  width = 10,
  height = 7,
  dpi = 300
)

print(outcome_correlation_plot)


# ------------------------------------------------------------
# 11. SCATTER PLOT ANALYSIS
# ------------------------------------------------------------

cat("\n[STEP 8] Examining relationships between selected variables...\n")

# Glucose and BMI are examined as two important
# clinical measurements.

scatter_plot <- ggplot(
  eda_data,
  aes(
    x = Glucose,
    y = BMI,
    color = Outcome
  )
) +
  
  geom_point(
    alpha = 0.65,
    size = 2
  ) +
  
  scale_color_manual(
    values = c(
      "No Diabetes" = "#3498DB",
      "Diabetes" = "#E74C3C"
    )
  ) +
  
  labs(
    title = "Relationship Between Glucose and BMI",
    subtitle = "Clinical observations grouped by diabetes outcome",
    x = "Glucose Concentration",
    y = "Body Mass Index (BMI)",
    color = "Outcome"
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

ggsave(
  "figures/07_glucose_bmi_scatterplot.png",
  plot = scatter_plot,
  width = 9,
  height = 7,
  dpi = 300
)

print(scatter_plot)


# ------------------------------------------------------------
# 12. AGE AND GLUCOSE DISTRIBUTIONS BY OUTCOME
# ------------------------------------------------------------

cat("\n[STEP 9] Analysing age and glucose patterns...\n")

# Figure 8: Age distribution by outcome

age_plot <- ggplot(
  eda_data,
  aes(
    x = Age,
    fill = Outcome
  )
) +
  
  geom_density(
    alpha = 0.5
  ) +
  
  scale_fill_manual(
    values = c(
      "No Diabetes" = "#3498DB",
      "Diabetes" = "#E74C3C"
    )
  ) +
  
  labs(
    title = "Age Distribution by Diabetes Outcome",
    x = "Age (Years)",
    y = "Density",
    fill = "Outcome"
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

ggsave(
  "figures/08_age_distribution_by_outcome.png",
  plot = age_plot,
  width = 9,
  height = 6,
  dpi = 300
)

print(age_plot)


# Figure 9: Glucose distribution by outcome

glucose_plot <- ggplot(
  eda_data,
  aes(
    x = Glucose,
    fill = Outcome
  )
) +
  
  geom_density(
    alpha = 0.5
  ) +
  
  scale_fill_manual(
    values = c(
      "No Diabetes" = "#3498DB",
      "Diabetes" = "#E74C3C"
    )
  ) +
  
  labs(
    title = "Glucose Distribution by Diabetes Outcome",
    x = "Glucose Concentration",
    y = "Density",
    fill = "Outcome"
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

ggsave(
  "figures/09_glucose_distribution_by_outcome.png",
  plot = glucose_plot,
  width = 9,
  height = 6,
  dpi = 300
)

print(glucose_plot)


# ------------------------------------------------------------
# 13. SAVE EDA RESULTS
# ------------------------------------------------------------

cat("\n[STEP 10] Saving exploratory analysis results...\n")

eda_results <- list(
  
  descriptive_statistics = descriptive_statistics,
  
  outcome_distribution = outcome_distribution,
  
  group_statistics = group_statistics,
  
  correlation_matrix = correlation_matrix,
  
  outcome_correlations = outcome_correlations
  
)

saveRDS(
  eda_results,
  "data/03_eda_results.rds"
)


# ------------------------------------------------------------
# 14. VERIFY GENERATED FILES
# ------------------------------------------------------------

cat("\n[STEP 11] Verifying generated output files...\n")

required_outputs <- c(
  
  "data/descriptive_statistics.csv",
  
  "data/outcome_distribution.csv",
  
  "data/group_statistics.csv",
  
  "data/correlation_matrix.csv",
  
  "data/outcome_correlations.csv",
  
  "data/03_eda_results.rds",
  
  "figures/01_outcome_distribution.png",
  
  "figures/02_clinical_histograms.png",
  
  "figures/03_clinical_density_plots.png",
  
  "figures/04_clinical_boxplots_by_outcome.png",
  
  "figures/05_correlation_heatmap.png",
  
  "figures/06_outcome_correlation.png",
  
  "figures/07_glucose_bmi_scatterplot.png",
  
  "figures/08_age_distribution_by_outcome.png",
  
  "figures/09_glucose_distribution_by_outcome.png"
  
)

output_status <- data.frame(
  
  File = required_outputs,
  
  Exists = file.exists(required_outputs)
  
)

print(output_status)

if (!all(output_status$Exists)) {
  
  stop(
    "One or more EDA output files were not generated."
  )
  
}


# ------------------------------------------------------------
# 15. FINAL SUMMARY
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("       EXPLORATORY ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nDataset observations:", nrow(eda_data), "\n")

cat("Clinical variables:", length(clinical_variables), "\n")

cat(
  "No Diabetes observations:",
  sum(eda_data$Outcome == "No Diabetes"),
  "\n"
)

cat(
  "Diabetes observations:",
  sum(eda_data$Outcome == "Diabetes"),
  "\n"
)

cat("\nGenerated statistical files:\n")

cat("- data/descriptive_statistics.csv\n")
cat("- data/outcome_distribution.csv\n")
cat("- data/group_statistics.csv\n")
cat("- data/correlation_matrix.csv\n")
cat("- data/outcome_correlations.csv\n")
cat("- data/03_eda_results.rds\n")

cat("\nGenerated figures:\n")

for (i in 1:9) {
  
  cat(
    paste0(
      "- figures/",
      sprintf("%02d", i),
      "_"
    ),
    "\n"
  )
  
}

cat("\nNext stage:\n")
cat("code/04_statistical_testing.R\n")

cat("\n============================================================\n")