import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.sparse import issparse
import matplotlib.pyplot as plt
import seaborn as sns

adata = ad.read_h5ad("C:/Users/Tycho/Desktop/SchoolTU/year3/q4_RP/data/0fce5dd5-bcec-4288-90b3-19a16b45ad16.h5ad", backed='r')
# adata = ad.read_h5ad("test_data.h5ad", backed='r')

# Cell types
# celltype = 'CD14-positive monocyte'
# celltype = 'naive thymus-derived CD4-positive, alpha-beta T cell'
# celltype = 'CD16-positive, CD56-dim natural killer cell, human'
# celltype = 'central memory CD4-positive, alpha-beta T cell'
# celltype = 'naive thymus-derived CD8-positive, alpha-beta T cell'
# celltype = 'CD8-positive, alpha-beta cytotoxic T cell'
# celltype = 'CD14-low, CD16-positive monocyte'
# celltype = 'CD8-positive, alpha-beta memory T cell'
# celltype = 'naive B cell'
# celltype = 'memory B cell'
# celltype = 'gamma-delta T cell'
# celltype = 'effector memory CD4-positive, alpha-beta T cell'
# celltype = 'mucosal invariant T cell'
# celltype = 'CD4-positive, alpha-beta T cell'
# celltype = 'T cell'
# celltype = 'monocyte'
# celltype = 'CD1c-positive myeloid dendritic cell'
# celltype = 'CD4-positive, alpha-beta cytotoxic T cell'
# celltype = 'regulatory T cell'
# celltype = 'platelet'
# celltype = 'natural killer cell'
# celltype = 'B cell'
# celltype = 'CD16-negative, CD56-bright natural killer cell, human'
# celltype = 'mature B cell'
# celltype = 'plasmacytoid dendritic cell'
# celltype = 'CD8-positive, alpha-beta T cell'
# celltype = 'plasma cell'
# celltype = 'CD141-positive myeloid dendritic cell'
# celltype = 'double negative T regulatory cell'
# celltype = 'conventional dendritic cell'
celltype = 'innate lymphoid cell'
# celltype = 'dendritic cell'
# celltype = 'erythrocyte'

# celltypes = [
#     'CD14-positive monocyte',
#     'naive thymus-derived CD4-positive, alpha-beta T cell',
#     'CD16-positive, CD56-dim natural killer cell, human',
#     'central memory CD4-positive, alpha-beta T cell',
#     'naive thymus-derived CD8-positive, alpha-beta T cell',
#     'CD8-positive, alpha-beta cytotoxic T cell',
#     'CD14-low, CD16-positive monocyte',
#     'CD8-positive, alpha-beta memory T cell',
#     'naive B cell',
#     'memory B cell',
#     'gamma-delta T cell',
#     'effector memory CD4-positive, alpha-beta T cell',
#     'mucosal invariant T cell',
#     'CD4-positive, alpha-beta T cell',
#     'T cell',
#     'monocyte',
#     'CD1c-positive myeloid dendritic cell',
#     'CD4-positive, alpha-beta cytotoxic T cell',
#     'regulatory T cell',
#     'platelet',
#     'natural killer cell',
#     'B cell',
#     'CD16-negative, CD56-bright natural killer cell, human',
#     'mature B cell',
#     'plasmacytoid dendritic cell',
#     'CD8-positive, alpha-beta T cell',
#     'plasma cell',
#     'CD141-positive myeloid dendritic cell',
#     'double negative T regulatory cell',
#     'conventional dendritic cell',
#     'innate lymphoid cell',
#     'dendritic cell',
    # 'erythrocyte'
# ]

# Create a subset of a cell type
subset = adata[adata.obs['cell_type'] == celltype].to_memory()
print(f"{celltype} shape: {subset.shape}")
print(f"{celltype} cell type: {subset.obs['cell_type'].unique()}")

# Calculate the total number of cells per age
age_counts_cells = subset.obs['development_stage'].value_counts()
# Extract numeric part from 'development_stage'
age_counts_cells.index = age_counts_cells.index.str.extract('(\d+)').astype(int).squeeze()  # Extract numeric part
age_counts = age_counts_cells.sort_index()  # Sort by age

# Remove all the genes with less than the number of cells of expression
non_zero_genes = np.array(subset.raw.X.sum(axis=0)).flatten() > subset.shape[0]
subset = subset[:, non_zero_genes]
print(f"Shape after removing low expression genes: {subset.shape}")

# Ensure 'donor_id' is in subset.columns
assert 'donor_id' in subset.obs.columns, "Column 'donor_id' not found!"

# Get unique donors and count cells per donor
donors = subset.obs['donor_id'].unique()
print(f"Found {len(donors)} donors.")

# Compute mean expression per donor
mean_expression = pd.DataFrame(
    data=np.zeros((len(donors), subset.n_vars)),
    index=donors,
    columns=subset.var_names
)

for donor in donors:
    donor_mask = subset.obs['donor_id'] == donor
    donor_cells = subset[donor_mask]
    mean_expression.loc[donor] = donor_cells.X.mean(axis=0)  # Average across cells


# Aggregate metadata (take first occurrence per donor)
donor_metadata = (
    subset.obs
    .groupby('donor_id')
    .first()
    .loc[donors]  # Preserve donor order
)

# Add the number of cells that donor had
donor_metadata['n_cells'] = subset.obs['donor_id'].value_counts().loc[donors].values

# Build new AnnData
donor_adata = sc.AnnData(
    X=mean_expression,  # or mean_expression.values for dense
    obs=donor_metadata,
    var=subset.var.copy(),
    uns=subset.uns.copy()
)

# Verify
print(f"New shape: {donor_adata.shape}")

# Calculate the total number of donors per age
age_counts_donors = donor_adata.obs['development_stage'].value_counts()
# Extract numeric part from 'development_stage'
age_counts_donors.index = age_counts_donors.index.str.extract('(\d+)').astype(int).squeeze()  # Extract numeric part
age_counts_donors = age_counts_donors.sort_index()  # Sort by age


# Plot the average number of cells per age
average_number_cells = age_counts_cells / age_counts_donors

# Plot the average number of cells per age
plt.figure(figsize=(10, 5))
plt.bar(average_number_cells.index, average_number_cells.values)
plt.xlabel('Age')
plt.ylabel('Average number of cells')
title = 'Average number of cells per age of ' + celltype
plt.title(title)
# Save the figure to folder figures
folder = "figures/average_cells_per_age/"
path = folder + title.replace(' ', '_').lower() + '.png'
plt.savefig(path)
plt.close()


# Remove all donors with less than 10% of the median number of cells
median_cells = donor_adata.obs['n_cells'].median()
print(f"Median number of cells: {median_cells}")
threshold = median_cells * 0.1
threshold = int(threshold)
print(f"Threshold: {threshold}")
donor_adata = donor_adata[donor_adata.obs['n_cells'] > threshold]
print(f"Shape after removing donors with less than threshold number of cells: {donor_adata.shape}")

# Show a distribution of the age of the donors
string_age = donor_adata.obs['development_stage'].astype(str)
string_age = string_age.str.extract('(\d+)').astype(int).squeeze()  # Extract numeric part
sns.histplot(string_age, bins=57)
plt.xlabel('Age')
plt.ylabel('Count')
title = 'Age Distribution of ' + celltype
plt.title(title)
# plt.show()
# Save the figure to folder figures
folder = "figures/distributions/"
path = folder + title.replace(' ', '_').lower() + '.png'
plt.savefig(path)

plt.close()



# Based on the age distribution create two subsets. One for the young and one for the old donors Young donors are those with age <= 36
# Old donors are those with age >= 47
# Based on the age distribution create three subsets. One for the young, one for the old and one for the middle aged donors
# Young donors are those with age <= 30
# Middle aged donors are those with age >= 40 and <= 50
# Old donors are those with age >= 60
young_donors = donor_adata[donor_adata.obs['development_stage'].astype(str).str.extract('(\d+)').astype(int).squeeze() <= 30]
middle_donors = donor_adata[(donor_adata.obs['development_stage'].astype(str).str.extract('(\d+)').astype(int).squeeze() >= 40) & (donor_adata.obs['development_stage'].astype(str).str.extract('(\d+)').astype(int).squeeze() <= 50)]
old_donors = donor_adata[donor_adata.obs['development_stage'].astype(str).str.extract('(\d+)').astype(int).squeeze() >= 60]
print(f"Young donors shape: {young_donors.shape}")
print(f"Middle donors shape: {middle_donors.shape}")
print(f"Old donors shape: {old_donors.shape}")


# Remove all the genes with zero expression in every individual donor
non_zero_genes_young = np.array(young_donors.X.sum(axis=0)).flatten() > 0
non_zero_genes_middle = np.array(middle_donors.X.sum(axis=0)).flatten() > 0
non_zero_genes_old = np.array(old_donors.X.sum(axis=0)).flatten() > 0
non_zero_genes = non_zero_genes_young & non_zero_genes_old & non_zero_genes_middle
young_donors = young_donors[:, non_zero_genes]
middle_donors = middle_donors[:, non_zero_genes]
old_donors = old_donors[:, non_zero_genes]
print(f"Young donors shape after removing zero expression genes: {young_donors.shape}")
print(f"Middle donors shape after removing zero expression genes: {middle_donors.shape}")
print(f"Old donors shape after removing zero expression genes: {old_donors.shape}")

# Save the subsets
folder = "subsets/"
young_path = folder + "{}_young_donors.h5ad".format(celltype)
middle_path = folder + "{}_middle_donors.h5ad".format(celltype)
old_path = folder + "{}_old_donors.h5ad".format(celltype)
young_donors.write_h5ad(young_path)
middle_donors.write_h5ad(middle_path)
old_donors.write_h5ad(old_path)