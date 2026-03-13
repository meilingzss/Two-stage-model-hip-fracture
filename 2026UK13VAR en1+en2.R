##project

{
  #install.packages(c("caret", "corrplot", "e1071", "lattice", "AppliedPredictiveModeling"))
  #install.packages("mlbench")
  #install.packages( "earth")
  #install.packages( "kernlab")
  #install.packages( "nnet")
  #install.packages( "caret")
  #install.packages("klaR")
  #install.packages("MASS")
  #install.packages("xgboost")
  #install.packages("pls")
  library(AppliedPredictiveModeling)
  library(corrplot)
  library(e1071)
  library(caret) 
  library(VIM) 
  library(dplyr)
  library(mda)
  library(MASS)
  library(klaR)
  library(nnet)
  library(kernlab)
  library(mlbench)
  library(earth)
  library(nnet)   
  library(fastDummies)
  library(glmnet)
  library(ROCR)
  library(pROC)
  library(pls)
  library(xgboost)
}
#mrosfrax_data <- read.csv("/Users/mac/Desktop/dissertation/cooperation/FOS/frax/mros+sof/frxa.mros(male).tscore.csv", header = TRUE, stringsAsFactors = FALSE)
#MrOS_clinical <- mrosfrax_data
#soffrax_data <- read.csv("/Users/mac/Desktop/dissertation/cooperation/FOS/frax/mros+sof/frxa.sof(female).tscore.csv", header = TRUE, stringsAsFactors = FALSE)
#sof_clinical <- soffrax_data

#df <- sof_clinical[, c("HA_AGE", "HA_WEIGHT", "HA_HEIGHT", "HA_FND",
  #                     "HA_SMOKE")]
#y <- sof_clinical[,c("HA_HIPFX")]
DATA <- read.csv("/Users/mac/Desktop/2026twostage/2026newy/UK_male_newtscore.14var.dataset.csv", header = TRUE, stringsAsFactors = FALSE)
# Columns to convert
cols_to_factor <- c("smoke", "hip", "SLDFX", "WRSTFX", "IADLdiffnew")

# Convert to factor
DATA[cols_to_factor] <- lapply(DATA[cols_to_factor], factor)

# Check structure
str(DATA)

# Define Clinical and DXA Features
en1_features <- c(  "height","walk","smoke","weight","age",
                         "IADLdiffnew","grip")


en2_features <- c("height","walk","smoke","weight","age","SLDFX","WRSTFX",
                  "IADLdiffnew","grip","neckbmd_Tscore","totalbmd_Tscore","spinebmd_Tscore")

df_clean <- na.omit(DATA)

df <- df_clean[, en1_features, drop = FALSE]  # Keep as data frame
y <- df_clean[,c("hip")]


## feature selection

# Load required libraries
library(caret) 
library(pROC)
# Set up parameters
n_runs <- 100  # Run the entire process 100 times
final_auc_scores <- numeric(n_runs)  # Store final AUC scores for each run
final_accuracy_scores <- numeric(n_runs)
final_sensitivity_scores <- numeric(n_runs)
final_specificity_scores <- numeric(n_runs)

for (run in 1:n_runs) {
  
  # Define predictor variables (X) and response variable (y)
  X <- df
  
  # Convert y to a factor (for logistic regression)
  y <- as.factor(y)
  
  # Set up parameters
  n_models <- 100  # Number of bootstrap samples (ensemble size)
  n_folds <- 5    # Number of cross-validation folds
  
  # Store AUC scores for each model
  auc_scores <- numeric(n_models)
  accuracy_scores <- numeric(n_models)
  sensitivity_scores <- numeric(n_models)
  specificity_scores <- numeric(n_models)
  
  # Perform Bootstrap Sampling and Logistic Regression Training
  #set.seed(42)  # For reproducibility
  for (i in 1:n_models) {
    
    # Bootstrap sample with replacement
    boot_idx <- sample(1:nrow(X), replace = TRUE)
    X_boot <- X[boot_idx, ]
    y_boot <- y[boot_idx]
    
    # Perform stratified cross-validation
    cv_folds <- createFolds(y_boot, k = n_folds, list = TRUE)
    auc_values <- numeric(n_folds)
    accuracy_values <- numeric(n_folds)
    sensitivity_values <- numeric(n_folds)
    specificity_values <- numeric(n_folds)
    
    for (j in 1:n_folds) {
      # Split data into training and validation sets
      train_idx <- unlist(cv_folds[-j])
      test_idx <- unlist(cv_folds[j])
      
      X_train <- X_boot[train_idx, ]
      y_train <- y_boot[train_idx]
      X_test <- X_boot[test_idx, ]
      y_test <- y_boot[test_idx]
      
      # Train Logistic Regression Model
      model <- glm(y_train ~ ., data = data.frame(y_train, X_train), family = "binomial")
      
      # Get predictions (probabilities)
      y_pred_prob <- predict(model, newdata = data.frame(X_test), type = "response")
      
      # Compute AUC
      roc_obj <- roc(y_test, y_pred_prob, quiet = TRUE)
      auc_values[j] <- roc_obj$auc
      # Youden's index = sensitivity + specificity - 1
      youden_index <- roc_obj$sensitivities + roc_obj$specificities - 1
      
      # Find the index of the maximum Youden's index
      best_idx <- which.max(youden_index)
      
      # Get the threshold corresponding to that index
      best_threshold <- roc_obj$thresholds[best_idx]
      
      
      # Ensure y_pred_prob is a numeric vector
      y_pred_prob <- as.numeric(y_pred_prob)  
      
      # Convert predicted probabilities into binary class labels using the best threshold
      predicted_labels <- rep(0, length(y_pred_prob))  # Initialize with 0s
      predicted_labels[y_pred_prob >= best_threshold] <- 1  # Assign 1 where condition is met
      # Compute confusion matrix
      cm <- confusionMatrix(factor(predicted_labels), factor(y_test), positive = "1")
      
      # Extract accuracy, sensitivity, and specificity
      accuracy_values[j] <- cm$overall["Accuracy"]
      sensitivity_values[j] <- cm$byClass["Sensitivity"]
      specificity_values[j] <- cm$byClass["Specificity"]
    }
    
    # Store the mean AUC across folds for this bootstrap sample
    auc_scores[i] <- mean(auc_values)
    accuracy_scores[i] <- mean(accuracy_values)
    sensitivity_scores[i] <- mean(sensitivity_values)
    specificity_scores[i] <- mean(specificity_values)
  }
  
  # Compute the final average AUC for this run
  final_auc_scores[run] <- mean(auc_scores)
  final_accuracy_scores[run] <- mean(accuracy_scores)
  final_sensitivity_scores[run] <- mean(sensitivity_scores)
  final_specificity_scores[run] <- mean(specificity_scores)
  
}

# Compute overall statistics for AUC, Accuracy, Sensitivity, and Specificity
mean_auc <- round(mean(final_auc_scores), 3)
sd_auc <- round(sd(final_auc_scores), 3)

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


