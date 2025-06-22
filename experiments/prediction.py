import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error
import scanpy as sc
import anndata as ad
from scipy import sparse
import os

def train_cell_type_specific_clocks(adata_path, output_dir, cell_type, genes_to_keep=None):
    """
    Final corrected version that properly creates:
    1. Model files in application-compatible format
    2. Correctly shaped imputation data
    3. Performance metrics
    """
    # Load data and subset by cell type
    adata_all = ad.read_h5ad(adata_path, backed='r')
    adata = adata_all[adata_all.obs['cell_type'] == cell_type].to_memory()
    
    # Filter genes if specified
    if genes_to_keep is not None:
        valid_genes = [gene for gene in genes_to_keep if gene in adata.var_names]
        print(f"Filtering genes - keeping {len(valid_genes)}/{len(genes_to_keep)} specified genes")
        adata = adata[:, valid_genes]
    
    # Preprocessing
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Prepare matrices
    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    y = np.array([int(stage.split('-')[0]) for stage in adata.obs['development_stage']])
    groups = adata.obs['donor_id'].values
    
    # Initialize storage for models and metrics
    all_models = []
    maes = []
    corrs = []
    
    # 5-fold CV
    cv = GroupKFold(n_splits=5)
    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y, groups)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        model = ElasticNet(alpha=1, l1_ratio=0.5, random_state=42)
        model.fit(X_train, y_train)
        
        # Create series with coefficients and intercept
        coef_series = pd.Series(model.coef_, index=adata.var_names)
        coef_series['intercept'] = model.intercept_
        all_models.append(coef_series)
        
        # Evaluate
        preds = model.predict(X_test)
        maes.append(mean_absolute_error(y_test, preds))
        corrs.append(np.corrcoef(y_test, preds)[0,1])
        print(f"Fold {fold+1}: MAE = {maes[-1]:.2f}, r = {corrs[-1]:.3f}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Save the 5-fold model (genes as columns, 5 coefficient rows)
    model_df = pd.DataFrame(all_models)
    model_filename = f"{cell_type.replace('/', '_')}_models5.csv"
    model_df.to_csv(f"{output_dir}/{model_filename}", index=False)
    
    # 2. Save imputation data
    impute_data = pd.DataFrame(
        data=[X.mean(axis=0)],
        columns=adata.var_names
    )
    impute_filename = f"Impute_avg_{cell_type.replace('/', '_')}.csv"
    impute_data.to_csv(f"{output_dir}/{impute_filename}", index=False)
    
    # 3. Save performance metrics
    metrics = pd.DataFrame({
        'cell_type': [cell_type],
        'mean_mae': [np.mean(maes)],
        'mean_r': [np.mean(corrs)],
        'n_features': [np.sum(np.any(all_models != 0, axis=0))],
        'genes_used': [str(genes_to_keep if genes_to_keep else 'all')]
    })
    metrics.to_csv(f"{output_dir}/clock_performance.csv", index=False)

    print(f"\nTraining complete! Files saved in {output_dir}:")
    print(f"- {model_filename} (5-fold models)")
    print(f"- {impute_filename} (imputation data)")
    print("- clock_performance.csv (metrics)")

if __name__ == "__main__":
    # Read the important genes from a file every line contains a gene
    cell_type = "effector memory CD4-positive, alpha-beta T cell"
    important_genes_file = f"results/{cell_type}_filtered_hubs.txt"
    if os.path.exists(important_genes_file):
        with open(important_genes_file, 'r') as f:
            important_genes = [line.strip() for line in f if line.strip()]
    
    train_cell_type_specific_clocks(
        adata_path="../h5ad/0fce5dd5-bcec-4288-90b3-19a16b45ad16.h5ad",
        output_dir="trained_clocks",
        cell_type=cell_type,
        genes_to_keep=important_genes
    )