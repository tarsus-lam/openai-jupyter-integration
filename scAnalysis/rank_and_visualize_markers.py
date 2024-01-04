'''
This modules sets the number of marker genes to display per cluster.
Marker genes are then ranked using the Wilcoxon test based on the clusters obtained from the Leiden algorithm.
The top N marker genes for each cluster are visualized and saved as '_marker_genes.png'.
The results are converted into a DataFrame and saved to 'tables/marker_genes_per_cluster.csv'.
Finally, the DataFrame is displayed.
'''

import numpy as np
import pandas as pd
import scanpy as sc

def rank_gene_groups(merged, save_plot = False, save_table = False):

    # Set the number of marker genes to display per cluster
    top_n_marker_genes = 10
    
    # Rank marker genes and plot the top N marker genes for each cluster
    sc.tl.rank_genes_groups(adata=merged, groupby='leiden', method='wilcoxon')
    if save_plot:
        sc.pl.rank_genes_groups(adata=merged, n_genes=top_n_marker_genes, save='_marker_genes.png')
    else:
        sc.pl.rank_genes_groups(adata=merged, n_genes=top_n_marker_genes)

    # Convert marker genes results to DataFrame
    results_df = pd.DataFrame(merged.uns['rank_genes_groups']['names'])
    results_df.head(10)

    # Save the DataFrame to the 'tables' folder
    if save_table:
        results_df.to_csv('tables/marker_genes_per_cluster.csv', index=False)

    return merged
