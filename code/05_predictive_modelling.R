# IT4133 - Bioinformatics and Computational Biology
# Project: Diabetes Prediction
# Script: 05_predictive_modelling.R
# Purpose: Predictive Modelling and Model Evaluation


# 1. LOAD REQUIRED PACKAGES

required_packages <- c(
  "tidyverse",
  "e1071",
  "rpart",
  "pROC"
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


# 2. SET RANDOM SEED

set.seed(123)


# 3. DEFINE PROJECT DIRECTORIES

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


# 4. LOAD CLEANED DATASET

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
cat("PREDICTIVE MODELLING - DIABETES DATASET\n")
cat("============================================\n")

cat(
  "\nDataset dimensions:",
  nrow(diabetes_data),
  "rows and",
  ncol(diabetes_data),
  "columns\n"
)


# 5. SELECT CLINICAL VARIABLES

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

required_columns <- c(
  clinical_variables,
  "Outcome"
)

missing_columns <- setdiff(
  required_columns,
  names(diabetes_data)
)

if (length(missing_columns) > 0) {
  stop(
    "Required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

analysis_data <- diabetes_data %>%
  select(all_of(required_columns)) %>%
  mutate(
    Outcome = factor(
      as.character(Outcome),
      levels = c("0", "1"),
      labels = c("No Diabetes", "Diabetes")
    )
  )

if (anyNA(analysis_data)) {
  stop(
    "Missing values detected. ",
    "Check the data cleaning script and Outcome encoding."
  )
}

if (!all(is.finite(as.matrix(
  analysis_data %>%
  select(all_of(clinical_variables))
)))) {
  stop("Non-finite values detected in clinical variables.")
}

cat("\nOutcome distribution:\n")
print(table(analysis_data$Outcome))


# 6. CREATE TRAINING AND TESTING DATASETS

cat("\n--------------------------------------------\n")
cat("6. TRAINING AND TESTING DATA PARTITION\n")
cat("--------------------------------------------\n")

# Stratified partitioning preserves the approximate
# outcome distribution in both datasets.

training_indices <- unlist(
  lapply(
    levels(analysis_data$Outcome),
    function(group) {
      
      group_indices <- which(
        analysis_data$Outcome == group
      )
      
      sample(
        group_indices,
        size = floor(
          0.80 * length(group_indices)
        )
      )
    }
  ),
  use.names = FALSE
)

training_indices <- sort(training_indices)

training_data <- analysis_data[
  training_indices,
  ,
  drop = FALSE
]

testing_data <- analysis_data[
  -training_indices,
  ,
  drop = FALSE
]

training_data$Outcome <- droplevels(
  training_data$Outcome
)

testing_data$Outcome <- factor(
  testing_data$Outcome,
  levels = levels(training_data$Outcome)
)

if (nrow(training_data) == 0 || nrow(testing_data) == 0) {
  stop("Training or testing dataset is empty.")
}

cat(
  "\nTraining observations:",
  nrow(training_data),
  "\n"
)

cat(
  "Testing observations:",
  nrow(testing_data),
  "\n"
)

cat("\nTraining outcome distribution:\n")
print(table(training_data$Outcome))

cat("\nTesting outcome distribution:\n")
print(table(testing_data$Outcome))


# 7. STANDARDIZE USING TRAINING DATA ONLY

cat("\n--------------------------------------------\n")
cat("7. DATA STANDARDIZATION\n")
cat("--------------------------------------------\n")

# Calculate means and standard deviations from training data.
# Testing data uses the same training parameters.

training_means <- sapply(
  training_data[clinical_variables],
  mean
)

training_sds <- sapply(
  training_data[clinical_variables],
  sd
)

# Protect against variables with zero or invalid standard deviation.

training_sds[
  training_sds == 0 |
    !is.finite(training_sds)
] <- 1

standardize_data <- function(data) {
  
  standardized <- as.data.frame(
    scale(
      data[, clinical_variables, drop = FALSE],
      center = training_means,
      scale = training_sds
    )
  )
  
  standardized$Outcome <- data$Outcome
  
  standardized
}

training_scaled <- standardize_data(
  training_data
)

testing_scaled <- standardize_data(
  testing_data
)

cat("\nStandardization completed.\n")


# 8. PREPARE MODEL INPUTS

x_train <- training_scaled %>%
  select(all_of(clinical_variables))

x_test <- testing_scaled %>%
  select(all_of(clinical_variables))

y_train <- training_scaled$Outcome
y_test <- testing_scaled$Outcome

cat("\nPredictor variables prepared.\n")


# 9. DEFINE HELPER FUNCTIONS

# Convert probability values into class predictions.

probability_to_class <- function(
    probabilities,
    threshold = 0.5) {
  
  factor(
    ifelse(
      probabilities >= threshold,
      "Diabetes",
      "No Diabetes"
    ),
    levels = c(
      "No Diabetes",
      "Diabetes"
    )
  )
}


# Safely calculate a binary classification model.

safe_model_prediction <- function(
    model_name,
    prediction_function) {
  
  tryCatch(
    
    prediction_function(),
    
    error = function(e) {
      
      warning(
        model_name,
        " prediction failed: ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
}


# 10. LOGISTIC REGRESSION

cat("\n--------------------------------------------\n")
cat("10. LOGISTIC REGRESSION\n")
cat("--------------------------------------------\n")

logistic_model <- glm(
  Outcome ~ .,
  data = training_scaled,
  family = binomial(
    link = "logit"
  )
)

logistic_probabilities <- predict(
  logistic_model,
  newdata = testing_scaled,
  type = "response"
)

logistic_probabilities <- as.numeric(
  logistic_probabilities
)

logistic_predictions <- probability_to_class(
  logistic_probabilities
)

cat("\nLogistic Regression completed.\n")

cat("\nModel summary:\n")
print(summary(logistic_model))


# 11. NAIVE BAYES

cat("\n--------------------------------------------\n")
cat("11. NAIVE BAYES\n")
cat("--------------------------------------------\n")

naive_bayes_model <- e1071::naiveBayes(
  Outcome ~ .,
  data = training_scaled,
  laplace = 1
)

naive_bayes_probability_matrix <- predict(
  naive_bayes_model,
  newdata = testing_scaled,
  type = "raw"
)

naive_bayes_probabilities <- as.numeric(
  naive_bayes_probability_matrix[
    ,
    "Diabetes"
  ]
)

naive_bayes_predictions <- probability_to_class(
  naive_bayes_probabilities
)

cat("\nNaive Bayes completed.\n")


# 12. DECISION TREE

cat("\n--------------------------------------------\n")
cat("12. DECISION TREE\n")
cat("--------------------------------------------\n")

decision_tree_model <- rpart(
  Outcome ~ .,
  data = training_scaled,
  method = "class",
  control = rpart.control(
    cp = 0.01,
    maxdepth = 5,
    minsplit = 20
  )
)

decision_tree_probability_matrix <- predict(
  decision_tree_model,
  newdata = testing_scaled,
  type = "prob"
)

decision_tree_probabilities <- as.numeric(
  decision_tree_probability_matrix[
    ,
    "Diabetes"
  ]
)

decision_tree_predictions <- probability_to_class(
  decision_tree_probabilities
)

cat("\nDecision Tree completed.\n")

cat("\nDecision Tree structure:\n")
print(decision_tree_model)


# 13. k-NEAREST NEIGHBOURS (k-NN)

cat("\n--------------------------------------------\n")
cat("13. k-NEAREST NEIGHBOURS\n")
cat("--------------------------------------------\n")

knn_k <- 5

knn_predictions <- class::knn(
  train = as.matrix(x_train),
  test = as.matrix(x_test),
  cl = y_train,
  k = knn_k,
  prob = TRUE
)

# class::knn returns the probability of the winning class.
# Convert it into the estimated probability of Diabetes.

knn_winning_probability <- as.numeric(
  attr(
    knn_predictions,
    "prob"
  )
)

knn_probabilities <- ifelse(
  knn_predictions == "Diabetes",
  knn_winning_probability,
  1 - knn_winning_probability
)

knn_predictions <- factor(
  knn_predictions,
  levels = c(
    "No Diabetes",
    "Diabetes"
  )
)

cat(
  "\nk-NN completed with k =",
  knn_k,
  "\n"
)


# 14. SUPPORT VECTOR MACHINE (SVM)

cat("\n--------------------------------------------\n")
cat("14. SUPPORT VECTOR MACHINE\n")
cat("--------------------------------------------\n")

svm_model <- e1071::svm(
  Outcome ~ .,
  data = training_scaled,
  kernel = "radial",
  cost = 1,
  gamma = 1 / length(clinical_variables),
  probability = TRUE,
  scale = FALSE
)

svm_predictions_raw <- predict(
  svm_model,
  newdata = testing_scaled,
  probability = TRUE
)

svm_probability_matrix <- attr(
  svm_predictions_raw,
  "probabilities"
)

if (
  is.null(svm_probability_matrix) ||
  !"Diabetes" %in% colnames(svm_probability_matrix)
) {
  stop(
    "SVM probability estimates could not be generated."
  )
}

svm_probabilities <- as.numeric(
  svm_probability_matrix[
    ,
    "Diabetes"
  ]
)

svm_predictions <- factor(
  as.character(svm_predictions_raw),
  levels = c(
    "No Diabetes",
    "Diabetes"
  )
)

cat("\nSupport Vector Machine completed.\n")


# 15. K-MEANS CLUSTERING

cat("\n--------------------------------------------\n")
cat("15. K-MEANS CLUSTERING\n")
cat("--------------------------------------------\n")

# K-means is an unsupervised method.
# The training observations are grouped into two clusters.
# Each cluster is assigned the majority outcome observed
# among its training members.
#
# Testing observations are assigned to their nearest
# training cluster centre.

set.seed(123)

kmeans_model <- kmeans(
  x = as.matrix(x_train),
  centers = 2,
  nstart = 25,
  iter.max = 100
)

# Determine the majority outcome in each training cluster.

cluster_training_data <- tibble(
  Cluster = kmeans_model$cluster,
  Outcome = y_train
)

cluster_mapping <- cluster_training_data %>%
  count(
    Cluster,
    Outcome
  ) %>%
  group_by(Cluster) %>%
  slice_max(
    order_by = n,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

# Ensure every cluster has a mapping.

cluster_labels <- setNames(
  rep(
    "No Diabetes",
    2
  ),
  as.character(1:2)
)

for (i in seq_len(nrow(cluster_mapping))) {
  
  cluster_labels[
    as.character(cluster_mapping$Cluster[i])
  ] <- as.character(cluster_mapping$Outcome[i])
}

# Calculate squared Euclidean distances between each
# testing observation and the training cluster centres.

kmeans_centres <- kmeans_model$centers

test_distances <- sapply(
  seq_len(nrow(kmeans_centres)),
  function(cluster_number) {
    
    rowSums(
      sweep(
        as.matrix(x_test),
        2,
        kmeans_centres[cluster_number, ],
        "-"
      )^2
    )
  }
)

test_distances <- as.matrix(
  test_distances
)

testing_clusters <- max.col(
  -test_distances
)

kmeans_predictions <- factor(
  cluster_labels[
    as.character(testing_clusters)
  ],
  levels = c(
    "No Diabetes",
    "Diabetes"
  )
)

# K-means does not produce calibrated class probabilities.
# Use normalized inverse-distance scores as ranking scores
# for ROC-AUC calculation.

inverse_distances <- 1 / (
  test_distances + 1e-10
)

cluster_diabetes_indicator <- as.numeric(
  cluster_labels[as.character(1:2)] == "Diabetes"
)

kmeans_probabilities <- as.vector(
  (
    inverse_distances %*%
      cluster_diabetes_indicator
  ) /
    rowSums(inverse_distances)
)

cat("\nK-means clustering completed.\n")


# 16. COLLECT MODEL PREDICTIONS

cat("\n--------------------------------------------\n")
cat("16. COLLECTING MODEL PREDICTIONS\n")
cat("--------------------------------------------\n")

model_predictions <- list(
  
  "Logistic Regression" = list(
    predicted = logistic_predictions,
    probability = logistic_probabilities
  ),
  
  "Naive Bayes" = list(
    predicted = naive_bayes_predictions,
    probability = naive_bayes_probabilities
  ),
  
  "Decision Tree" = list(
    predicted = decision_tree_predictions,
    probability = decision_tree_probabilities
  ),
  
  "k-Nearest Neighbours" = list(
    predicted = knn_predictions,
    probability = knn_probabilities
  ),
  
  "Support Vector Machine" = list(
    predicted = svm_predictions,
    probability = svm_probabilities
  ),
  
  "K-means Clustering" = list(
    predicted = kmeans_predictions,
    probability = kmeans_probabilities
  )
)


# Verify that all methods generated valid predictions.

for (model_name in names(model_predictions)) {
  
  current_predictions <- model_predictions[[model_name]]$predicted
  
  current_probabilities <- model_predictions[[model_name]]$probability
  
  if (
    length(current_predictions) != length(y_test) ||
    length(current_probabilities) != length(y_test)
  ) {
    stop(
      "Prediction length mismatch for ",
      model_name
    )
  }
  
  if (
    anyNA(current_predictions) ||
    anyNA(current_probabilities)
  ) {
    stop(
      "Missing predictions detected for ",
      model_name
    )
  }
  
  if (
    any(!is.finite(current_probabilities)) ||
    any(current_probabilities < 0 | current_probabilities > 1)
  ) {
    stop(
      "Invalid probability values detected for ",
      model_name
    )
  }
}

cat("\nAll six methods generated testing predictions.\n")


# 17. CALCULATE CLASSIFICATION METRICS

cat("\n--------------------------------------------\n")
cat("17. MODEL PERFORMANCE EVALUATION\n")
cat("--------------------------------------------\n")

calculate_metrics <- function(
    actual,
    predicted,
    probabilities,
    model_name) {
  
  actual <- factor(
    actual,
    levels = c(
      "No Diabetes",
      "Diabetes"
    )
  )
  
  predicted <- factor(
    predicted,
    levels = c(
      "No Diabetes",
      "Diabetes"
    )
  )
  
  confusion <- table(
    Actual = actual,
    Predicted = predicted
  )
  
  tn <- as.numeric(
    confusion["No Diabetes", "No Diabetes"]
  )
  
  fp <- as.numeric(
    confusion["No Diabetes", "Diabetes"]
  )
  
  fn <- as.numeric(
    confusion["Diabetes", "No Diabetes"]
  )
  
  tp <- as.numeric(
    confusion["Diabetes", "Diabetes"]
  )
  
  accuracy <- (tp + tn) / (
    tp + tn + fp + fn
  )
  
  sensitivity <- if (
    tp + fn == 0
  ) {
    NA_real_
  } else {
    tp / (tp + fn)
  }
  
  specificity <- if (
    tn + fp == 0
  ) {
    NA_real_
  } else {
    tn / (tn + fp)
  }
  
  precision <- if (
    tp + fp == 0
  ) {
    NA_real_
  } else {
    tp / (tp + fp)
  }
  
  f1_score <- if (
    is.na(precision) ||
    is.na(sensitivity) ||
    precision + sensitivity == 0
  ) {
    NA_real_
  } else {
    2 * precision * sensitivity /
      (precision + sensitivity)
  }
  
  roc_result <- tryCatch(
    
    pROC::roc(
      response = actual,
      predictor = as.numeric(probabilities),
      levels = c(
        "No Diabetes",
        "Diabetes"
      ),
      direction = "<",
      quiet = TRUE
    ),
    
    error = function(e) {
      NULL
    }
  )
  
  auc_value <- if (is.null(roc_result)) {
    NA_real_
  } else {
    as.numeric(
      pROC::auc(roc_result)
    )
  }
  
  metrics <- tibble(
    Model = model_name,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1_Score = f1_score,
    ROC_AUC = auc_value
  )
  
  list(
    metrics = metrics,
    confusion_matrix = confusion,
    roc = roc_result
  )
}


# 18. EVALUATE ALL SIX METHODS

evaluation_results <- list()

for (model_name in names(model_predictions)) {
  
  evaluation_results[[model_name]] <- calculate_metrics(
    actual = y_test,
    predicted = model_predictions[[model_name]]$predicted,
    probabilities = model_predictions[[model_name]]$probability,
    model_name = model_name
  )
}

performance_results <- purrr::map_dfr(
  evaluation_results,
  "metrics"
) %>%
  arrange(
    desc(ROC_AUC),
    desc(F1_Score)
  )

cat("\nModel performance comparison:\n")

print(
  performance_results %>%
    mutate(
      across(
        where(is.numeric),
        ~ round(.x, 4)
      )
    )
)


# 19. DISPLAY CONFUSION MATRICES

cat("\n--------------------------------------------\n")
cat("19. CONFUSION MATRICES\n")
cat("--------------------------------------------\n")

for (model_name in names(evaluation_results)) {
  
  cat("\nModel:", model_name, "\n")
  
  print(
    evaluation_results[[model_name]]$confusion_matrix
  )
}


# 20. PREPARE ROC CURVE DATA

roc_data <- purrr::map_dfr(
  names(evaluation_results),
  function(model_name) {
    
    roc_result <- evaluation_results[[model_name]]$roc
    
    if (is.null(roc_result)) {
      return(tibble())
    }
    
    tibble(
      Model = model_name,
      Specificity = roc_result$specificities,
      Sensitivity = roc_result$sensitivities,
      False_Positive_Rate = 1 - roc_result$specificities
    )
  }
)

if (nrow(roc_data) == 0) {
  stop("No valid ROC curves could be generated.")
}


# 21. FIGURE 13 - MODEL PERFORMANCE COMPARISON

figure_13 <- performance_results %>%
  select(
    Model,
    Accuracy,
    Sensitivity,
    Specificity,
    Precision,
    F1_Score
  ) %>%
  pivot_longer(
    cols = -Model,
    names_to = "Metric",
    values_to = "Score"
  ) %>%
  ggplot(
    aes(
      x = Model,
      y = Score,
      fill = Metric
    )
  ) +
  
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  
  labs(
    title = "Classification Performance Comparison",
    subtitle = "Performance metrics evaluated on the testing dataset",
    x = "Computational Method",
    y = "Performance Score",
    fill = "Metric"
  ) +
  
  coord_cartesian(
    ylim = c(0, 1)
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
      angle = 30,
      hjust = 1
    ),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    figures_dir,
    "13_model_performance_comparison.png"
  ),
  plot = figure_13,
  width = 12,
  height = 7,
  dpi = 300
)


# 22. FIGURE 14 - ROC CURVE COMPARISON

figure_14 <- ggplot(
  roc_data,
  aes(
    x = False_Positive_Rate,
    y = Sensitivity,
    color = Model
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "gray"
  ) +
  
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  
  labs(
    title = "Receiver Operating Characteristic Curves",
    subtitle = "Comparison of computational methods",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Method"
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    figures_dir,
    "14_roc_curve_comparison.png"
  ),
  plot = figure_14,
  width = 10,
  height = 8,
  dpi = 300
)


# 23. FIGURE 15 - ROC-AUC COMPARISON

figure_15 <- performance_results %>%
  filter(
    !is.na(ROC_AUC)
  ) %>%
  ggplot(
    aes(
      x = reorder(Model, ROC_AUC),
      y = ROC_AUC
    )
  ) +
  
  geom_col(
    fill = "steelblue",
    width = 0.7
  ) +
  
  geom_text(
    aes(
      label = sprintf("%.3f", ROC_AUC)
    ),
    hjust = -0.15,
    size = 3.5
  ) +
  
  # Set the axis limit before flipping the coordinates.
  coord_flip(
    ylim = c(
      0,
      min(
        1.15,
        max(
          performance_results$ROC_AUC,
          na.rm = TRUE
        ) + 0.15
      )
    )
  ) +
  
  labs(
    title = "ROC-AUC Comparison",
    subtitle = "Area under the ROC curve for each method",
    x = "Computational Method",
    y = "ROC-AUC"
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
    "15_roc_auc_comparison.png"
  ),
  plot = figure_15,
  width = 10,
  height = 6,
  dpi = 300
)


# 24. FIGURE 16 - CONFUSION MATRIX COMPARISON

confusion_data <- purrr::imap_dfr(
  evaluation_results,
  function(result, model_name) {
    
    as.data.frame(
      result$confusion_matrix
    ) %>%
      as_tibble() %>%
      rename(
        Count = Freq
      ) %>%
      mutate(
        Model = model_name
      )
  }
)

figure_16 <- ggplot(
  confusion_data,
  aes(
    x = Predicted,
    y = Actual,
    fill = Count
  )
) +
  
  geom_tile(
    color = "white"
  ) +
  
  geom_text(
    aes(
      label = Count
    ),
    size = 4
  ) +
  
  facet_wrap(
    ~ Model,
    ncol = 3
  ) +
  
  labs(
    title = "Confusion Matrix Comparison",
    subtitle = "Actual versus predicted diabetes outcomes",
    x = "Predicted Outcome",
    y = "Actual Outcome",
    fill = "Observations"
  ) +
  
  theme_minimal(
    base_size = 11
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
      angle = 20,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    figures_dir,
    "16_confusion_matrix_comparison.png"
  ),
  plot = figure_16,
  width = 13,
  height = 9,
  dpi = 300
)


# 25. SAVE MODEL PERFORMANCE RESULTS

write_csv(
  performance_results,
  file.path(
    data_dir,
    "05_model_performance.csv"
  )
)


# 26. SAVE INDIVIDUAL CONFUSION MATRICES

confusion_matrices <- purrr::imap(
  evaluation_results,
  function(result, model_name) {
    result$confusion_matrix
  }
)

confusion_export <- purrr::imap_dfr(
  confusion_matrices,
  function(matrix, model_name) {
    
    as.data.frame(matrix) %>%
      as_tibble() %>%
      rename(
        Count = Freq
      ) %>%
      mutate(
        Model = model_name,
        .before = 1
      )
  }
)

write_csv(
  confusion_export,
  file.path(
    data_dir,
    "05_confusion_matrices.csv"
  )
)


# 27. SAVE TESTING PREDICTIONS

prediction_results <- tibble(
  Actual_Outcome = y_test
)

for (model_name in names(model_predictions)) {
  
  safe_column_name <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    model_name
  )
  
  prediction_results[[
    paste0(
      safe_column_name,
      "_Prediction"
    )
  ]] <- model_predictions[[model_name]]$predicted
  
  prediction_results[[
    paste0(
      safe_column_name,
      "_Probability"
    )
  ]] <- model_predictions[[model_name]]$probability
}

write_csv(
  prediction_results,
  file.path(
    data_dir,
    "05_testing_predictions.csv"
  )
)


# 28. SAVE TRAINED MODELS AND RESULTS

predictive_results <- list(
  
  training_data = training_data,
  
  testing_data = testing_data,
  
  training_indices = training_indices,
  
  training_means = training_means,
  
  training_sds = training_sds,
  
  logistic_regression = logistic_model,
  
  naive_bayes = naive_bayes_model,
  
  decision_tree = decision_tree_model,
  
  knn_predictions = knn_predictions,
  
  knn_k = knn_k,
  
  svm = svm_model,
  
  kmeans = kmeans_model,
  
  kmeans_cluster_mapping = cluster_labels,
  
  predictions = model_predictions,
  
  evaluation = evaluation_results,
  
  performance = performance_results,
  
  confusion_matrices = confusion_matrices,
  
  roc_data = roc_data,
  
  clinical_variables = clinical_variables
)

saveRDS(
  predictive_results,
  file.path(
    data_dir,
    "05_predictive_results.rds"
  )
)


# 29. IDENTIFY THE TOP PERFORMING METHOD

valid_performance <- performance_results %>%
  filter(
    !is.na(ROC_AUC),
    !is.na(F1_Score)
  )

if (nrow(valid_performance) > 0) {
  
  best_model <- valid_performance %>%
    arrange(
      desc(ROC_AUC),
      desc(F1_Score)
    ) %>%
    slice(1)
  
  cat("\n--------------------------------------------\n")
  cat("29. TOP PERFORMING METHOD\n")
  cat("--------------------------------------------\n")
  
  cat(
    "\nTop method based on ROC-AUC:",
    as.character(best_model$Model),
    "\n"
  )
  
  cat(
    "ROC-AUC:",
    round(best_model$ROC_AUC, 4),
    "\n"
  )
  
  cat(
    "Accuracy:",
    round(best_model$Accuracy, 4),
    "\n"
  )
  
  cat(
    "F1-score:",
    round(best_model$F1_Score, 4),
    "\n"
  )
}


# 30. VERIFY GENERATED OUTPUT FILES

expected_outputs <- c(
  
  "data/05_model_performance.csv",
  
  "data/05_confusion_matrices.csv",
  
  "data/05_testing_predictions.csv",
  
  "data/05_predictive_results.rds",
  
  "figures/13_model_performance_comparison.png",
  
  "figures/14_roc_curve_comparison.png",
  
  "figures/15_roc_auc_comparison.png",
  
  "figures/16_confusion_matrix_comparison.png"
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


# 31. FINAL SUMMARY

cat("\n============================================\n")
cat("PREDICTIVE MODELLING COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat("\nMethods evaluated:\n")

cat("1. Logistic Regression\n")
cat("2. Naive Bayes\n")
cat("3. Decision Tree\n")
cat("4. k-Nearest Neighbours\n")
cat("5. Support Vector Machine\n")
cat("6. K-means Clustering\n")

cat("\nEvaluation metrics:\n")

cat("1. Accuracy\n")
cat("2. Sensitivity\n")
cat("3. Specificity\n")
cat("4. Precision\n")
cat("5. F1-score\n")
cat("6. ROC-AUC\n")

cat("\nGenerated figures:\n")

cat("figures/13_model_performance_comparison.png\n")
cat("figures/14_roc_curve_comparison.png\n")
cat("figures/15_roc_auc_comparison.png\n")
cat("figures/16_confusion_matrix_comparison.png\n")

cat("\nGenerated data files:\n")

cat("data/05_model_performance.csv\n")
cat("data/05_confusion_matrices.csv\n")
cat("data/05_testing_predictions.csv\n")
cat("data/05_predictive_results.rds\n")

cat("\nAll six methods have been evaluated.\n")

cat("\n============================================\n")
