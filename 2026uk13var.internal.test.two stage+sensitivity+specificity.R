library(caret)
library(pROC)
library(dplyr)
library(tidyr)

#set.seed(191)
# Load the dataset
# Define predictor variables (X) and response variable (y)
#X <- stage_clean[,optVariables]  # Top selected features from previous feature selection

# Convert y to a factor (for logistic regression)
#y <- as.factor(y)
#df <- cbind(X,y)
# Select the new set of columns from SOF_clinical
df <- read.csv("/Users/mac/Desktop/2026twostage/2026newy/UK_male_newtscore.14var.dataset.csv", header = TRUE, stringsAsFactors = FALSE)

# Columns to convert
cols_to_factor <- c("smoke", "hip", "SLDFX", "WRSTFX", "IADLdiffnew")

# Convert to factor
df[cols_to_factor] <- lapply(df[cols_to_factor], factor)

# Check structure
str(df)

# Define Clinical and DXA Features
clinical_features <- c(  "height","walk","smoke","weight","age","SLDFX","WRSTFX",
                         "IADLdiffnew","grip")


dxa_features <- c("neckbmd_Tscore","totalbmd_Tscore","spinebmd_Tscore")
target_col <- "hip"




# Drop rows with missing values in selected features
# Create a vector of columns to check for NA values

# Remove rows with NA values in any of the specified columns
df_clean <- na.omit(df)
y <- as.factor(df_clean$hip)
# Merge levels 1 and 2 into a single level 1
#df_clean$smoke <- ifelse(df_clean$smoke %in% c(1, 2), 1, 0)

# Convert to factor with levels 0 and 1
#df_clean$smoke <- factor(df_clean$smoke, levels = c(0, 1))

# Check the new counts
#table(df_clean$smoke)

# Merge levels 3 and 4 into a single level 3
#df_clean$IADLdiffnew <- ifelse(df_clean$IADLdiffnew %in% c(3, 4), 3, df_clean$IADLdiffnew)

# Convert to factor (optional, useful for modeling)
#df_clean$IADLdiffnew <- factor(df_clean$IADLdiffnew, levels = sort(unique(df_clean$IADLdiffnew)))

# Check the new counts
#table(df_clean$IADLdiffnew)

# Extract X (Clinical + DXA) and y (Target) using base R
X_clinical <- df_clean[, clinical_features, drop = FALSE]  # Keep as data frame
X_dxa <- df_clean[, dxa_features, drop = FALSE]  # Keep as data frame

# Split Data (80% Train, 20% Test)
#set.seed(42)
n_runs <- 100
auc_values_staged <- numeric(n_runs)
final_accuracy_scores <- numeric(n_runs)
final_sensitivity_scores <- numeric(n_runs)
final_specificity_scores <- numeric(n_runs)
threshold <- numeric(n_runs)
accuracy_stage1 <- numeric(n_runs)
percent_stage1 <- numeric(n_runs)
accuracy_stage2 <- numeric(n_runs)
percent_stage2 <- numeric(n_runs)

ssindx <- vector("list", n_runs)  # Creates an empty list with n_runs elements

for (run in 1:n_runs) {
  trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)
  X_train_c <- X_clinical[trainIndex, ]
  X_test_c <- X_clinical[-trainIndex, ]
  X_train_d <- X_dxa[trainIndex, ]
  X_test_d <- X_dxa[-trainIndex, ]
  X_train_d <- data.frame(X_train_d = X_train_d)
  X_test_d  <- data.frame(X_train_d = X_test_d)  # note: name must match training
  y_train <- y[trainIndex]
  y_test <- y[-trainIndex]
  
  # Train Ensemble 1 (Bootstrapped Logistic Regression on Clinical Features)
  #set.seed(42)
  # Split Train Data into Training and Validation (80% of train for training, 20% for validation)
  k_folds <- 5
  folds <- createFolds(y_train, k = k_folds, list = TRUE,returnTrain = FALSE)
  
  # Store AUCs and thresholds
  auc_values <- numeric(k_folds)
  thresholds <- numeric(k_folds)
  
  # Perform k-fold cross-validation
  for (fold in seq_along(folds)) {
    # Split into training and validation fold
    train_idx <- unlist(folds[-fold])  # Use all other folds as training
    valid_idx <- folds[[fold]]  # Current fold as validation
    
    X_train_cv <- X_train_c[train_idx, ]
    y_train_cv <- y_train[train_idx]
    X_valid_cv <- X_train_c[valid_idx, ]
    y_valid_cv <- y_train[valid_idx]
    # Check if validation fold has at least one case and one control
    if (length(unique(y_valid_cv)) < 2) {
      stop("Validation fold does not have both classes! Reduce k or use stratified sampling.")
    }
    

    
    # Train Bootstrapped Ensemble Logistic Regression on Training Fold
    n_models <- 50
    proba_valid_cv <- matrix(NA, nrow = nrow(X_valid_cv), ncol = n_models)
    
    for (i in 1:n_models) {
      boot_idx <- sample(seq_len(nrow(X_train_cv)), replace = TRUE)
      model <- glm(y_train_cv[boot_idx] ~ ., data = X_train_cv[boot_idx, ], family = binomial)
      proba_valid_cv[, i] <- predict(model, newdata = X_valid_cv, type = "response")
    }
    
    # Compute mean probability for the validation fold
    mean_proba_valid_cv <- rowMeans(proba_valid_cv)
    
    # Compute ROC curve
    roc_obj_cv <- roc(y_valid_cv, mean_proba_valid_cv)
    auc_values[fold] <- auc(roc_obj_cv)  # Store AUC for this fold
    
    # Compute Youden's index to determine the optimal threshold
    youden_index_cv <- roc_obj_cv$sensitivities + roc_obj_cv$specificities - 1
    valid_thresholds_cv <- roc_obj_cv$thresholds[is.finite(roc_obj_cv$thresholds)]
    thresholds[fold] <- valid_thresholds_cv[which.max(youden_index_cv)]
  }
  
  # Select the best threshold corresponding to the best AUC
  best_fold <- which.max(auc_values)
  best_threshold <- thresholds[best_fold]
  
  # Print results
  print(paste("Best AUC:", round(auc_values[best_fold], 4)))
  print(paste("Optimal Threshold from Best Validation Fold:", round(best_threshold, 4)))
  
  # Train Final Ensemble on Full Training Set (X_train_c) and Evaluate on Test Set
  n_models <- 50
  # Train Ensemble on Full Training Data (X_train_c) and Evaluate on Test Set
  proba_test <- matrix(NA, nrow = nrow(X_test_c), ncol = n_models)
  
  for (i in 1:n_models) {
    boot_idx <- sample(seq_len(nrow(X_train_c)), replace = TRUE)
    model <- glm(y_train[boot_idx] ~ ., data = X_train_c[boot_idx, ], family = binomial)
    proba_test[, i] <- predict(model, newdata = X_test_c, type = "response")
  }
  
  # Compute mean probability for test set
  mean_proba_test <- rowMeans(proba_test)
  
  # Compute standard deviation for each row of proba_test
  sd_proba_test <- apply(proba_test, 1, sd)
  
  # Print first few values
  head(mean_proba_test)
  head(sd_proba_test)
  
  # Compute z-score for uncertainty
  uncertainty_ensemble1 <- abs((mean_proba_test-best_threshold) / sd_proba_test)
  
  # Define uncertain cases based on z-score threshold (e.g., 1 standard deviation)
  z_threshold <- 2  # Adjust as needed
  uncertain_cases <- uncertainty_ensemble1 < z_threshold
  sum(uncertain_cases=="TRUE")
  
  # Prepare Data for Stage 2 (Clinical + DXA Features)
  X_train_combined <- cbind(X_train_c, X_train_d)
  X_test_combined <- cbind(X_test_c, X_test_d)
  y_train_aligned <- y_train
  
  # Train Ensemble 2 (Clinical + DXA)
  #set.seed(42)
  n_models <- 50
  bootstrap_models2 <- list()
  proba_ensemble2 <- matrix(NA, nrow = nrow(X_test_combined), ncol = n_models)
  
  
  for (i in 1:n_models) {
    boot_idx <- sample(seq_len(nrow(X_train_combined)), replace = TRUE)
    model <- glm(y_train_aligned[boot_idx] ~ ., data = X_train_combined[boot_idx, ], family = binomial)
    proba_ensemble2[, i] <- predict(model, newdata = X_test_combined, type = "response")
  }
  mean_proba_ensemble2 <- rowMeans(proba_ensemble2)
  
  # Merge Staged Predictions
  final_predictions <- mean_proba_test
  
  # Ensure mean_proba_ensemble2 has the same length as uncertain cases
  if (sum(uncertain_cases) > 0) {
    final_predictions[uncertain_cases] <- mean_proba_ensemble2[uncertain_cases]
  }
  
  ### for group staying in stage 1
  
  # Get the indices of certain cases (uncertain_cases == "FALSE")
  certain_indices <- which(uncertain_cases == "FALSE")
  y_test[certain_indices]
  mean_proba_test[certain_indices]
  
  
  library(pROC)
  
  # Extract relevant data
  y_true <- y_test[certain_indices]  # Actual labels
  y_proba <- mean_proba_test[certain_indices]  # Predicted probabilities
  
  
  percentage_stage1 <- length(y_true)/length(y_test)
  
  ## enter stage 2 person and percentage
  
  y_true2 <- y_test[-certain_indices]  # Actual labels
  y_proba2 <- mean_proba_ensemble2[-certain_indices]  # Predicted probabilities
  
  percentage_stage2 <- length(y_true2)/length(y_test)
  
  
  # 3️⃣ Accuracy Using Optimal Threshold
  predicted_labels_opt <- ifelse(y_proba >= best_threshold, 1, 0)
  accuracy_opt <- mean(predicted_labels_opt == y_true)
  
  ## stage 2 accutacy
  predicted_labels_opt2 <- ifelse(y_proba2 >= best_threshold, 1, 0)
  accuracy_opt2 <- mean(predicted_labels_opt2 == y_true2)
  
  print(paste("Accuracy with optimal threshold:", round(accuracy_opt, 4)))
  sum(y_proba > best_threshold & y_proba < 0.5)
  
  testIndex_c <- setdiff(1:nrow(X_clinical), trainIndex)
  
  stage1_sample_index <- testIndex_c[certain_indices]
  
  auc_values_staged[run] <- auc(roc(y_test, final_predictions))
  
  roc_obj <- roc(y_test, final_predictions)
  # Get the optimal threshold using Youden's index
  best_threshold <- best_threshold
  #best_threshold <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")
  # Ensure y_pred_prob is a numeric vector
  final_predictions <- as.numeric(final_predictions)  
  
  # Convert predicted probabilities into binary class labels using the best threshold
  predicted_labels <- rep(0, length(final_predictions))  # Initialize with 0s
  predicted_labels[final_predictions >= best_threshold] <- 1  # Assign 1 where condition is met
  # Compute confusion matrix
  cm <- confusionMatrix(factor(predicted_labels), factor(y_test), positive = "1")
  
  # Extract accuracy, sensitivity, and specificity
  final_accuracy_scores[run] <- cm$overall["Accuracy"]
  final_sensitivity_scores[run]<- cm$byClass["Sensitivity"]
  final_specificity_scores[run] <- cm$byClass["Specificity"]
  
  threshold[run] <- best_threshold
  accuracy_stage1[run] <- accuracy_opt
  accuracy_stage2[run] <- accuracy_opt2
  
  percent_stage1[run] <- percentage_stage1
  percent_stage2[run] <- percentage_stage2
  
  ssindx[[run]] <- stage1_sample_index
}

# Compute the average AUC over 100 runs
# Compute the average AUC over 100 runs
average_auc <- mean(auc_values_staged)
cat("\u2705 Average AUC over 100 runs =", round(average_auc, 3), "confirming robust performance.\n")
mean_auc <- round(mean(auc_values_staged), 3)
sd_auc <- round(sd(auc_values_staged), 3)

mean_accuracy <- round(mean(final_accuracy_scores), 3)
sd_accuracy <- round(sd(final_accuracy_scores), 3)

mean_sensitivity <- round(mean(final_sensitivity_scores), 3)
sd_sensitivity <- round(sd(final_sensitivity_scores), 3)

mean_specificity <- round(mean(final_specificity_scores), 3)
sd_specificity <- round(sd(final_specificity_scores), 3)

# Create a data frame with the results
performance_summary <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity"),
  Mean = c(mean_auc, mean_accuracy, mean_sensitivity, mean_specificity),
  SD = c(sd_auc, sd_accuracy, sd_sensitivity, sd_specificity)
)

# Print the summary table
print(performance_summary)

average_accuracy <- mean(accuracy_stage1)
average_accuracy
average_percent <- mean(percent_stage1)
average_percent
average_accuracy <- mean(accuracy_stage2)
average_accuracy
average_percent <- mean(percent_stage2)
average_percent

