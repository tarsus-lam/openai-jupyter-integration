# scRNA-Seq Analysis Using OpenAI

## Purpose
This repository demonstrates the feasibility of using the OpenAI API and GPT models to perform bioinformatics analyses, with a focus on single-cell RNA sequencing (scRNA-seq) analysis. Utilizing simple prompts, users can execute major steps of scRNA-seq analysis and obtain insights with minimal to no background in bioinformatics. This project also illustrates how to connect to the OpenAI API using the `openai` Python package.

## Package Structure
The `scAnalysis` package includes the following modules:

scAnalysis

├── init.py

├── read_and_merge_data.py

├── preprocess_and_normalize.py

├── identify_top_genes.py

├── cluster_top_genes.py

├── rank_and_visualize_markers.py

├── annotate_clusters.py

├── validate_canonical_markers.py

├── analyze_differential_expression.py


## Analysis Steps
The Jupyter Notebook in this repository performs the following steps for scRNA-seq analysis, with all Python code obtained from the OpenAI API outputs:
1. OpenAI API Integration Setup
2. Load required packages
3. Load scRNA-Seq data and Merge the datasets
4. Remove cells with too few genes or too many genes, and genes detected in too few cells
5. Normalize the gene expression measurements to account for differences in sequencing depth
6. Log-transform the data for downstream analysis (data processing)
7. Select genes that show high variation across cells (most informative for clustering)
8. Scale the data to have zero mean and unit variance
9. Perform Principal Component Analysis (PCA)
10. Run clustering algorithms to identify distinct groups of cells
11. Identifying marker genes to enhance visualization of distinct cellular populations
12. Identify and confirm the cell types associated with each cluster (validation of cellular identities)
13. Create mappings from clusters to identified cell types
14. Validate the cell types by plotting expression levels of known canonical markers
15. Identify genes that are differentially expressed between different cell populations or conditions

## Setup
### Creating a Conda Environment
To set up your environment for running the analysis, follow these steps:

1. Create a new Conda environment:
   ```bash
   conda create -n scrna_analysis python=3.11.5
   ```
2. Activate the environment:
   ```bash
   conda activate scrna_analysis
   ```
3. Install the required packages:
   ```bash
   conda install jupyter=7.0.6 leidenalg=0.10.1 openai=1.6.1 re=2.2.1 scanpy=1.9.6 numpy=1.26.3 pandas=2.1.4 matplotlib=3.8.2 seaborn=0.13.1 scipy=1.11.4 scikit-learn=1.3.2 scikit-misc=0.3.1
   ```
   The required packages and versions are also listed in `requirements.txt`

## Configuration
Store your OpenAI API key in a .env file at the base directory of the repository in the format:

```
OPENAI_API_KEY="[INSERT KEY HERE]"
```

## Usage
You can run the analysis either by executing the Jupyter Notebook openai_jupyter_integration.ipynb or by running the Python code output by GPT through main.py.

```
Perform scRNA-seq analysis.

positional arguments:
  treatment     Path to the treated sample for analysis.
  control       Path to the control (untreated) sample.

options:
  -h, --help    show this help message and exit
  --save_table  Save tables produced in the analysis.
  --save_plot   Save plots produced in the analysis.
```

## Data Description
- The data were obtained from 10x Genomics and produce by CellRanger. It comprises A549 lung carcinoma cells that expressed dCas9-KRAB and were transduced with a pool containing 93 total sgRNAs. Selected cells for each condition were individually frozen, then thawed and counted for analysis.
- The dataset can be obtained by searching `5k A549, Lung Carcinoma Cells, No Treatment Transduced with a CRISPR Pool` on the 10x Genomics website. The treatment sample `Gene Expression - Feature / cell matrix (raw)` is under the 'Inputs/Library' tab, and the control sample `Gene Expression - Feature / cell matrix (per-sample)` is under the 'No Treatment' tab.
- Included in this repo are the corresponding HDF5 files, stored in the `data/lung_control` and `data/lung_treatment` folders.

## Disclaimers
- The Jupyter Notebook openai_jupyter_integration.ipynb is not guaranteed to be reproducible as intended since GPT outputs may vary, providing different analysis methods each time.
- To reproduce the analysis in this repository, please use main.py.
- Note that the analysis is not comprehensive and has not undergone rigorous validation of the identified cell types. This project mainly serves as a demonstration of using OpenAI's API for bioinformatics analyses.

## Author
Tarsus Lam
