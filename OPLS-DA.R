library(openxlsx)
library(ropls)
library(readxl)
library(gridExtra)
library(ggplot2)

# Read the data

###### define the folders
folder <- "Li_Oi-DAP_DAM"
#folder <- "Li-Oi-metbolome_réduit"
#folder <- "Li_Lh-total_metabolome"
#folder <- "6 samples"

#define the omic
omic <- "metabolome"



################################################################  Loop-OPLS-DA ################################################################

#define the datasets to use
#metabolome_dataset <- list("log2", "log2_PS", "log2_UV", "log2_centered", "log10", "log10_PS", "log10_UV", "log10_centered") # log2+log10
metabolome_preprocess <- list("log2_PS", "log2_UV", "log2_centered") #uniquement log2
proteome_preprocess <- list("KNN_CV", "KNN_CV_PS", "KNN_CV_UV", "KNN_CV_centered") 

# Initialize an empty data frame to collect all metrics
all_metrics_df <- data.frame()

for (i in metabolome_preprocess) {m
  print(i)
  file_path <- sprintf("O2PLS_DA/%s/data/%s/%s_%s.xlsx", folder, omic, omic, i)
  data <- read_excel(file_path)
  #print(head(data))
  
  # Prepare the data
  # The first column is sample names, the second column is class labels, the rest are features
  *%>  sample_names <- data[[1]]
  class_labels <- data[[2]]
  features <- data[, -c(1, 2)]
  
  # Perform OPLS-DA
  # Convert class labels to a factor
  class_labels <- as.factor(class_labels)
  
  # Determine the optimal number of predictive components using cross-validation
  # Note: ropls automatically handles orthogonal components
  
  # Create the opls model
  set.seed(123) # For reproducibility
  c_predi <- 1
  c_ortho <- 1
  opls_model <- opls(features, class_labels, predI = c_predi, orthoI = c_ortho, crossvalI = 3, permI = 100, scaleC = "none")
  
  ############# Extract VIP scores
  
  vip_scores <- getVipVn(opls_model)
  
  # Create a data frame for VIP scores and variable names
  vip_scores_df <- data.frame(
    Variable = colnames(features),
    VIP = vip_scores
  )
  
  # Sort the data frame by VIP scores in descending order
  vip_scores_df <- vip_scores_df[order(vip_scores_df$VIP, decreasing = TRUE), ]
  
  # Save the VIP scores to an Excel file
  write.xlsx(vip_scores_df, sprintf("O2PLS_DA/%s/data/%s/OPLS-DA/%s_%s_VIP.xlsx", folder, omic, omic, i), rowNames = FALSE)
  
  # Print the first few rows of sorted VIP scores
  print(head(vip_scores_df, 10))
  
  ############# Extract the summary plot
  png(sprintf("O2PLS_DA/%s/data/%s/OPLS-DA/%s_%s.png", folder, omic, omic, i), width = 1600, height = 1200, res = 300)
  plot(opls_model, type = "summary")
  dev.off()
  
  # Extract summaryDF from the opls_model object
  summary_df <- opls_model@summaryDF
  
  # Convert the summaryDF to a data frame (if necessary)
  summary_df <- as.data.frame(summary_df)
  
  # Add a column for the dataset name
  summary_df$Dataset <- i
  
  # Append the metrics to the master data frame
  all_metrics_df <- rbind(all_metrics_df, summary_df)
  
  # Define the output file path for the Excel file
  output_excel_path <- sprintf("O2PLS_DA/%s/data/%s/OPLS-DA/%s_%s_metrics.xlsx", folder, omic, omic, i)
  
  # Write the summaryDF to an Excel file
  write.xlsx(summary_df, output_excel_path, rowNames = TRUE)
  
}

# Save the master metrics data frame to an Excel file
write.xlsx(all_metrics_df, sprintf("O2PLS_DA/%s/data/%s/OPLS-DA/all_metrics_%s_%s.xlsx", folder, omic, c_predi, c_ortho), rowNames = TRUE)

print("end of script")



################################################################  Single OPLS-DA ################################################################
  
file_path <- sprintf("O2PLS_DA/%s/data/metabolome/metabolome_%s.xlsx", folder, metabolome_dataset)
data <- read_excel(file_path)
  
# Prepare the data
#The first column is sample names, the second column is class labels, the rest are features
sample_names <- data[[1]]
class_labels <- data[[2]]
features <- data[, -c(1, 2)]
  
# Perform OPLS-DA
# Convert class labels to a factor
class_labels <- as.factor(class_labels)
  
# Determine the optimal number of predictive components using cross-validation
# Note: ropls automatically handles orthogonal components
  
# Create the opls model
set.seed(123) # For reproducibility
opls_model <- opls(features, class_labels, log10L = FALSE, predI = 1, orthoI = 1, crossvalI = 5, permI = 100, scaleC = "pareto")


############# Extract VIP scores
vip_scores <- getVipVn(opls_model)
  
# Create a data frame for VIP scores and variable names
vip_scores_df <- data.frame(
  Variable = colnames(features),
  VIP = vip_scores
)
  
# Sort the data frame by VIP scores in descending order
vip_scores_df <- vip_scores_df[order(vip_scores_df$VIP, decreasing = TRUE), ]
  
# Save the VIP scores to a CSV file
write.xlsx(vip_scores_df, sprintf("O2PLS_DA/%s/data/metabolome/OPLS-DA/metabolome_%s.xlsx", folder, metabolome_dataset), rowNames = FALSE)
  
# Print the first few rows of sorted VIP scores
print(head(vip_scores_df, 10))


############# Extract the summary plot
png(sprintf("O2PLS_DA/%s/data/metabolome/OPLS-DA/metabolome_%s_plot.png", folder, metabolome_dataset), width = 1600, height = 1200, res = 300)
plot(opls_model, type = "summary")
dev.off()


# Extract summaryDF from the opls_model object
summary_df <- opls_model@summaryDF

# Convert the summaryDF to a data frame (if necessary)
summary_df <- as.data.frame(summary_df)

# Define the output file path for the Excel file
output_excel_path <- sprintf("O2PLS_DA/%s/data/metabolome/OPLS-DA/metabolome_%s_metrics.xlsx", folder, metabolome_dataset)

# Write the summaryDF to an Excel file
write.xlsx(summary_df, output_excel_path, rowNames = TRUE)
