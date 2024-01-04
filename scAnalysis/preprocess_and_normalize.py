'''
This module performs basic preprocessing on the merged scRNA-seq dataset.
It filters out low-quality cells and genes, considering a minimum gene count of 200 per cell and a minimum presence in 3 cells per gene.
QC metrics are calculated, and the data is normalized to a target sum of 1e4.
Finally, logarithmic scaling is applied to the dataset.
'''

import scanpy as sc

def preprocess_data(merged):
    # Filter out low-quality cells and genes
    sc.pp.filter_cells(merged, min_genes=200)
    sc.pp.filter_genes(merged, min_cells=3)

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(merged)

    # Normalize the data
    sc.pp.normalize_total(merged, target_sum=1e4)

    # Logarithmic scaling
    sc.pp.log1p(merged)

    return merged
