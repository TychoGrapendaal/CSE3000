# CSE3000
All the code I used during the research project

# Experimental Setup

All results shown in the paper can be reproduced by running the Python code in this repository.

All results are generated using **Jupyter notebooks** or **Python** files with:
- **Python 3.10.11**
- **Anndata 0.11.4**
- **Numpy 1.24.3**
- **Pandas 1.5.3**
- **Scipy 1.10.1**
- **Scanpy 1.11.1**
- **Matplotlib 3.7.1**
- **Seaborn 0.13.2**
- **Sklearn 1.3.0**
- **Powerlaw 1.5**

---

## Experiment Order

1. **Experiment 0: Data and Preprocessing**  
   - First, download the data and place it in the [h5ad](h5ad) folder.
   - Run [create_subset_9_best.py](create_subset_9_best.py) to generate the required subsets.

2. **Experiment 1: Distribution of the Donors**  
   - Run the notebook [experiments/distribution_cell_types_overlap.ipynb](experiments/distribution_cell_types_overlap.ipynb)

3. **Experiment 2: Correlation Analysis Across Age Groups**  
   - Run the notebook [experiments/correlation_analysis_experiment.ipynb](experiments/correlation_analysis_experiment.ipynb) to generate:  
     - The correlation matrix
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
   - Run the notebook [experiments/network_node_distribution.ipynb](experiments/network_node_distribution.ipynb).

5. **Experiment 4: Comparing Cell Types**  
   - Run the notebook [experiments/compare_celltypes.ipynb](experiments/compare_celltypes.ipynb).

6. **Experiment 5: Prediction**  
   - Run the notebook [experiments/hubs_single_cell.ipynb](experiments/hubs_single_cell.ipynb) to retrieve the hub genes for the specified cell type.
   - Run the notebook [experiments/prediction.py](experiments/prediction.py) to generate the linear regression models for each cell type
   - Run the notebook [experiments/apply_clocks_general.py](experiments/apply_clocks_general.py) to apply the clocks to the data. (This file is created by Enikő Zakar-Polyák, Attila Csordas, Róbert Pálovics, and Csaba Kerepesi. Profiling the transcriptomic age of single-cells in humans. Commun. Biol., 7(1):1397, October 2024. https://github.com/polyake/scaging)

---

## Notes  
- Ensure all dependencies (e.g., Python 3.10.11, Anndata 0.11.4) are installed before running the experiments.  
- Run the experiments in the specified order for consistent results. 