import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.sparse import issparse
import matplotlib.pyplot as plt
import seaborn as sns

adata = ad.read_h5ad("C:/Users/Tycho/Desktop/SchoolTU/year3/q4_RP/data/0fce5dd5-bcec-4288-90b3-19a16b45ad16.h5ad", backed='r')
celltype = 'erythrocyte'

# Create a subset of a cell type
subset = adata[adata.obs['cell_type'] == celltype].to_memory()
print(f"{celltype} shape: {subset.shape}")
print(f"{celltype} cell type: {subset.obs['cell_type'].unique()}")

# Remove all the genes with zero expression in every cell
non_zero_genes = np.array(subset.X.sum(axis=0)).flatten() > 0
subset = subset[:, non_zero_genes]
print(f"Shape after removing zero expression genes: {subset.shape}")

# Ensure 'donor_id' is in adata.obs
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

# Build new AnnData
donor_adata = sc.AnnData(
    X=mean_expression,  # or mean_expression.values for dense
    obs=donor_metadata,
    var=subset.var.copy(),
    uns=subset.uns.copy()
)

# Verify
print(f"New shape: {donor_adata.shape}")  # Should be (72, 14301)

# Show a distrubution of the age of the donors
string_age = donor_adata.obs['development_stage'].astype(str)
string_age = string_age.str.extract('(\d+)').astype(int).squeeze()  # Extract numeric part
sns.histplot(string_age, bins=54)
plt.xlabel('Age')
plt.ylabel('Count')
title = 'Age Distribution of ' + celltype
plt.title(title)
# plt.show()
# Save the figure to folder figures
path = 'figures/' + title.replace(' ', '_').lower() + '.png'
plt.savefig(path)


# Based on the age distribution create two subsets. One for the young and one for the old donors Young donors are those with age < 30
# Old donors are those with age >= 50
young_donors = donor_adata[donor_adata.obs['development_stage'].astype(str).str.extract('(\d+)').astype(int).squeeze() < 30]
old_donors = donor_adata[donor_adata.obs['development_stage'].astype(str).str.extract('(\d+)').astype(int).squeeze() >= 50]
print(f"Young donors shape: {young_donors.shape}")
print(f"Old donors shape: {old_donors.shape}")


# Save the subsets to folder subsets
young_donors.write_h5ad('subsets/young_donors.h5ad')
old_donors.write_h5ad('subsets/old_donors.h5ad')
