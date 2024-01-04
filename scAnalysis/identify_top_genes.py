'''
This module identifies highly variable genes using the Seurat flavor of the sc.pp.highly_variable_genes function,
focusing on the top 2000 genes. Subsequently, two visualizations are generated.
The first one visualizes the identified highly variable genes, saved as '_highly_variable_genes.png'.
The second plot displays the highest expressed genes (top 20), saved as '_highest_expr_genes.png'.
'''

import scanpy as sc

def get_highly_variable_genes(merged, save = False):
    # Identify highly variable genes
    sc.pp.highly_variable_genes(merged, flavor='seurat', n_top_genes=2000)
    
    # Visualize highly variable genes
    sc.pl.highly_variable_genes(merged, save='_highly_variable_genes.png')
    
    # Display and save figure
    if save:
        sc.pl.highest_expr_genes(merged, n_top=20, save='_highest_expr_genes.png')
    else:
        sc.pl.highest_expr_genes(merged, n_top=20)
    
    return merged
