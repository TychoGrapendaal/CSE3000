# CSE3000
All the code I used during the research project

# Experimental Setup

All results shown in the paper can be reproduced by running the Python code in this repository.

Most results are generated using **Jupyter notebooks** with:  
- **Python 3.10.11**
- **Anndata 0.11.4**

---

## Experiment Order

1. **Experiment 0: Data and Preprocessing**  
   - First, download the data and place it in the `h5ad` folder.  
   - Run `create_subset_9_best.py` to generate the required subsets.  

2. **Experiment 1: Distribution of the Donors**  
   - Run the notebook `experiments/distribution_cell_types_overlap.ipynb` to generate Figure 1 and 2 of my paper.

3. **Experiment 2: Correlation Analysis Across Age Groups**  
   - Run the notebook `experiments/correlation_analysis_experiment.ipynb` to generate:  
     - The correlation matrix  
     - The table of biggest absolute differences (for a specific cell type)  
   - Modify the `cell_type` variable in the notebook for each of the following cell types and run it:
     - CD8-positive, alpha-beta T cell
     - CD8-positive, alpha-beta memory T cell
     - CD4-positive, alpha-beta T cell
     - central memory CD4-positive, alpha-beta T cell
     - effector memory CD4-positive, alpha-beta T cell
     - gamma-delta T cell
     - regulatory T cell
     - double negative T regulatory cell
     - innate lymphoid cell

4. **Experiment 3: Gene Networks**  
   - Run the notebook `experiments/network_node_distribution.ipynb` to generate Figure 3.

5. **Experiment 4: Comparing Cell Types**  
   - Run the notebook `experiments/compare_celltypes.ipynb` to generate Figure 4. 

---

## Notes  
- Ensure all dependencies (e.g., Python 3.10.11, Anndata 0.11.4) are installed before running the experiments.  
- Run the experiments in the specified order for consistent results. 