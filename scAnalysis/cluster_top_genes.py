'''
This module prepares the data for clustering by scaling, performing principal component analysis (PCA),
creating a neighborhood graph, and calculating UMAP coordinates.
Clustering is then executed using the Leiden algorithm with a resolution of 0.1.
The resulting cluster relationships are visualized on a UMAP plot, and the figure is saved as '_clustering.png'.
'''

import scanpy as sc

def perform_clustering(merged, save = False):
    # Prepare the data for clustering
    sc.pp.scale(merged) # Scale the data
    sc.pp.pca(merged) # Perform dimensionality reduction
    sc.pp.neighbors(merged)
    sc.tl.umap(merged)

    # Perform clustering
    sc.tl.leiden(merged, resolution=0.1)

    # Visualize cluster relationships
    if save:
        sc.pl.umap(merged, color='leiden', save='_clustering.png')
    else:
        sc.pl.umap(merged, color='leiden')

    return merged
