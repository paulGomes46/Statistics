
# Center both datasets in order to perform O2PLS
# Both datasets are equalized to have a same total variance
# Note : datasets with more variables will have the total variance shared between more variables, and that might have a negative impact on OPLS-DA


import pandas as pd

# Load datasets again
folder = "Li_Oi-DAP_DAM"
df1 = pd.read_excel(f'O2PLS_DA/{folder}/data/proteome/proteome_KNN_CV_PS.xlsx')
df2 = pd.read_excel(f'O2PLS_DA/{folder}/data/metabolome/metabolome_log2_PS.xlsx')

# Separate sample names, class labels, and features
sample_names_df1 = df1.iloc[:, 0]
class_labels_df1 = df1.iloc[:, 1]
features_df1 = df1.iloc[:, 2:]

sample_names_df2 = df2.iloc[:, 0]
class_labels_df2 = df2.iloc[:, 1]
features_df2 = df2.iloc[:, 2:]

# Ensure all features are numeric
features_df1 = features_df1.apply(pd.to_numeric, errors='coerce')
features_df2 = features_df2.apply(pd.to_numeric, errors='coerce')

# Drop rows with NaN values in features
features_df1 = features_df1.dropna()
features_df2 = features_df2.dropna()

# Calculate sum of squares for each dataset
ss_df1 = (features_df1 ** 2).sum().sum()
ss_df2 = (features_df2 ** 2).sum().sum()

print(ss_df1, ss_df2)

# Determine scaling factors
scaling_factor_df1 = (ss_df2 / ss_df1) ** 0.5
scaling_factor_df2 = (ss_df1 / ss_df2) ** 0.5

print(scaling_factor_df1, scaling_factor_df2)

# Scale the datasets
scaled_df1 = features_df1 * scaling_factor_df1
scaled_df2 = features_df2  # features_df2 is left unchanged

# Verify sum of squares
ss_scaled_df1 = (scaled_df1 ** 2).sum().sum()
ss_scaled_df2 = (scaled_df2 ** 2).sum().sum()

print(ss_scaled_df1, ss_scaled_df2)

# Combine scaled features with sample names and class labels
df1_scaled = pd.concat([sample_names_df1, class_labels_df1, scaled_df1], axis=1)
df2_scaled = pd.concat([sample_names_df2, class_labels_df2, scaled_df2], axis=1)

# Save scaled datasets to new Excel files
df1_scaled.to_excel(f'O2PLS_DA/{folder}/data/proteome/proteome_centered.xlsx', index=False)
df2_scaled.to_excel(f'O2PLS_DA/{folder}/data/metabolome/metabolome_centered.xlsx', index=False)
