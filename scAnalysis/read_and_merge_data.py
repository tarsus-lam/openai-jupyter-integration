'''
This module utilizes Scanpy to read two sets of single-cell RNA-seq data representing samples from treated and control conditions.
The samples are stored in the 'lung_treatment' and 'lung_control' variables, and they are then merged into a single dataset named 'merged'.
The samples are annotated with a new column 'status', indicating whether each cell belongs to the treatment or control group.
'''

import scanpy as sc

def read_and_merge_samples(treatment, control):
    # Read in the samples
    lung_treatment = sc.read_10x_mtx(treatment)
    lung_control = sc.read_10x_mtx(control)
    
    # Merge the datasets
    merged = lung_treatment.concatenate(lung_control, join='outer', index_unique=None)
    
    # Annotate the samples
    merged.obs['status'] = ['treatment'] * lung_treatment.shape[0] + ['control'] * lung_control.shape[0]

    return merged
