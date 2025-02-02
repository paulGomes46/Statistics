# Load necessary libraries
library(OmicsPLS)
library(openxlsx)
library(writexl)
library(readxl)
library(ggplot2)

# Load your datasets
folder <- "Li_Oi-DAP_DAM"
proteome <- sprintf('O2PLS_DA/%s/data/proteome/proteome_centered.xlsx', folder)
metabolome <- sprintf('O2PLS_DA/%s/data/metabolome/metabolome_centered.xlsx', folder)

df1 <- read_excel(proteome)
df2 <- read_excel(metabolome)

# Extract features, assuming the first two columns are sample names and class labels
features_df1 <- as.matrix(df1[, -c(1, 2)])
features_df2 <- as.matrix(df2[, -c(1, 2)])

# Define permutation function
perform_permutation <- function(X, Y, n, nx, ny, n_permutations) {
  original_model <- o2m(X = X, Y = Y, n = n, nx = nx, ny = ny)
  original_R2X <- original_model$R2X
  original_R2Y <- original_model$R2Y
  
  permuted_R2X <- numeric(n_permutations)
  permuted_R2Y <- numeric(n_permutations)
  
  for (i in 1:n_permutations) {
    permuted_Y <- Y[sample(nrow(Y)), ]
    permuted_model <- o2m(X = X, Y = permuted_Y, n = n, nx = nx, ny = ny)
    permuted_R2X[i] <- permuted_model$R2X
    permuted_R2Y[i] <- permuted_model$R2Y
  }
  
  p_value_R2X <- mean(permuted_R2X >= original_R2X)
  p_value_R2Y <- mean(permuted_R2Y >= original_R2Y)
  
  return(list(p_value_R2X = p_value_R2X, p_value_R2Y = p_value_R2Y))
}

# Perform O2PLS with the selected number of components
joint <- 1
ortho_x <- 1
ortho_y <- 1
o2pls_model <- o2m(X = features_df1, Y = features_df2, n = joint, nx = ortho_x, ny = ortho_y)
summary(o2pls_model)

# Perform permutation testing
n_permutations <- 1000
perm_test_results <- perform_permutation(features_df1, features_df2, joint, ortho_x, ortho_y, n_permutations)

# Print permutation test results
print(perm_test_results)

# Extract loadings and scores
loadings_list_X <- loadings(o2pls_model, "Xjoint")
loadings_list_Y <- loadings(o2pls_model, "Yjoint")

# Convert row names to a column
loadings_list_X_df <- as.data.frame(loadings_list_X)
loadings_list_Y_df <- as.data.frame(loadings_list_Y)
loadings_list_X_df$Variable <- rownames(loadings_list_X_df)
loadings_list_Y_df$Variable <- rownames(loadings_list_Y_df)

# Reorder columns to have Variable names first
loadings_list_X_df <- loadings_list_X_df[, c("Variable", setdiff(colnames(loadings_list_X_df), "Variable"))]
loadings_list_Y_df <- loadings_list_Y_df[, c("Variable", setdiff(colnames(loadings_list_Y_df), "Variable"))]

# Write the data frame to an Excel file
write_xlsx(loadings_list_X_df, "PhD/OPLS-DA_training/loadings_list_X_df.xlsx")
write_xlsx(loadings_list_Y_df, "PhD/OPLS-DA_training/loadings_list_Y_df.xlsx")

# Extract scores
Score_X <- scores(o2pls_model, which_part = "Xjoint")
Score_Y <- scores(o2pls_model, which_part = "Yjoint")

# Convert matrix to df
score_X_df <- as.data.frame(Score_X)
score_Y_df <- as.data.frame(Score_Y)

# Export df as excel file
write_xlsx(score_X_df, "PhD/OPLS-DA_training/score_X_df.xlsx")
write_xlsx(score_Y_df, "PhD/OPLS-DA_training/score_Y_df.xlsx")
