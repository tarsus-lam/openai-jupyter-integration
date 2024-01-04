'''
This module defines the number of top differentially expressed genes to keep (top_n_de_genes).
It then calculates differential expression using the t-test_overestim_var method.
The results are stored in the de_results dictionary,
where each cell type has a DataFrame containing the top differentially expressed genes, their scores, and adjusted p-values.
Finally, these results are consolidated into a single DataFrame, 'de_results_df,'
and saved to 'tables/differential_expression_results.csv'.
'''

import numpy as np
import pandas as pd
import scanpy as sc

def differential_expression(merged, save = False):
    # Define the number of top differentially expressed genes to keep
    top_n_de_genes = 10

    # Calculate differential expression
    sc.tl.rank_genes_groups(merged, groupby='Cell Type', method='t-test_overestim_var')

    # Initialize an empty dictionary to store the results
    de_results = {}

    # Loop over each cell type and get the top differentially expressed genes
    for cell_type in merged.obs['Cell Type'].unique():
        cell_type_marker_genes = pd.DataFrame(merged.uns['rank_genes_groups']['names'][cell_type])[:top_n_de_genes]
        cell_type_marker_scores = pd.DataFrame(merged.uns['rank_genes_groups']['scores'][cell_type])[:top_n_de_genes]
        cell_type_marker_pvals = pd.DataFrame(merged.uns['rank_genes_groups']['pvals_adj'][cell_type])[:top_n_de_genes]
        
        # Combine the marker genes, scores, and p-values into a single DataFrame
        cell_type_de_results = pd.concat([cell_type_marker_genes, cell_type_marker_scores, cell_type_marker_pvals], axis=1)
        cell_type_de_results.columns = ['Gene', 'Score', 'Adjusted P-value']
        
        # Save the results for the cell type in the dictionary
        de_results[cell_type] = cell_type_de_results

    # Create a new DataFrame to combine the results for all cell types
    de_results_df = pd.concat(de_results, names=['Cell Type']).reset_index(level='Cell Type')
    de_results_df.head(10)
    
    # Save the DataFrame to the 'tables' folder
    if save:
        de_results_df.to_csv('tables/differential_expression_results.csv', index=False)

    return merged
    
