# Load necessary libraries
library(OmicsPLS)
library(caret)
library(readxl)
library(writexl)

# Load your datasets
df1 <- read_excel("PhD/OPLS-DA_training/Pareto_scaled_data_log10_scaled.xlsx") #metabolome
df2 <- read_excel("PhD/OPLS-DA_training/pareto_scaled_data_prot_knn_cv_scaled.xlsx") #proteome

# Extract features, assuming the first two columns are sample names and class labels
features_df1 <- as.matrix(df1[, -c(1, 2)])
features_df2 <- as.matrix(df2[, -c(1, 2)])
class_labels <- as.factor(df1[[2]])  # Extract the second column correctly as a vector

# Define function for MCCV
perform_mccv <- function(features_df1, features_df2, class_labels, n, nx, ny, num_iterations = 100, train_ratio = 0.8) {
  set.seed(123)
  results <- data.frame(n = integer(), nx = integer(), ny = integer(), R2X = double(), R2Y = double())
  
  for (i in 1:num_iterations) {
    # Create balanced train/test split
    train_indices <- createDataPartition(class_labels, p = train_ratio, list = FALSE)
    X_train <- features_df1[train_indices, ]
    Y_train <- features_df2[train_indices, ]
    X_test <- features_df1[-train_indices, ]
    Y_test <- features_df2[-train_indices, ]
    
    # Ensure number of components is less than sample size
    sample_size <- nrow(X_train)
    cat("Iteration", i, "Sample size:", sample_size, "n =", n, "nx =", nx, "ny =", ny, "\n")
    
    if (n + max(nx, ny) >= sample_size) {
      cat("Skipping iteration: n =", n, "nx =", nx, "ny =", ny, "sample_size =", sample_size, "\n")
      next
    }
    
    # Fit O2PLS model
    cat("Fitting O2PLS model: n =", n, "nx =", nx, "ny =", ny, "sample_size =", sample_size, "\n")
    o2pls_model <- o2m(X_train, Y_train, n = n, nx = nx, ny = ny)
    
    # Calculate R2 metrics for test data
    R2X <- sum(o2pls_model$explvar$Xjoint) / sum(o2pls_model$explvar$Xtotal)
    R2Y <- sum(o2pls_model$explvar$Yjoint) / sum(o2pls_model$explvar$Ytotal)
    
    # Store results
    results <- rbind(results, data.frame(n = n, nx = nx, ny = ny, R2X = R2X, R2Y = R2Y))
    cat("Iteration", i, "completed: n =", n, "nx =", nx, "ny =", ny, "R2X =", R2X, "R2Y =", R2Y, "\n")
  }
  
  return(results)
}

# Define grid of parameters to test
n_values <- c(1, 2, 3, 4)
nx_values <- c(0, 1)
ny_values <- c(0, 1)
results_list <- list()

# Perform MCCV for each combination of parameters
for (n in n_values) {
  for (nx in nx_values) {
    for (ny in ny_values) {
      cat("Testing combination: n =", n, "nx =", nx, "ny =", ny, "\n")
      mccv_results <- perform_mccv(features_df1, features_df2, class_labels, n, nx, ny, num_iterations = 100, train_ratio = 0.8)
      if (nrow(mccv_results) > 0) {
        results_list[[paste(n, nx, ny, sep = "_")]] <- mccv_results
      } else {
        cat("No results for combination: n =", n, "nx =", nx, "ny =", ny, "\n")
      }
    }
  }
}

# Combine and analyze results
combined_results <- do.call(rbind, results_list)

# Check if combined_results is empty
if (nrow(combined_results) == 0) {
  stop("No results to aggregate. Please check the parameter grid and sample size constraints.")
}

summary_results <- aggregate(cbind(R2X, R2Y) ~ n + nx + ny, data = combined_results, mean)

# Save results to Excel
write_xlsx(summary_results, "PhD/OPLS-DA_training/mccv_results.xlsx")

# Print the summary of the best model based on R2X and R2Y
print(summary_results)
