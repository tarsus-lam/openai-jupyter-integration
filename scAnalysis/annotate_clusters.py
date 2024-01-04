'''
This module establishes a mapping between the clusters and cell types based on the previous analysis.
A DataFrame, 'cluster_df,' is created to consolidate information about clusters, top marker genes, and corresponding cell types.
This DataFrame is saved to 'tables/clusters_marker_genes_cell_types.csv'.
Additionally, a new column, 'Cell Type,' is added to the original data, identifying the cell type associated with each cluster.
'''


import numpy as np
import pandas as pd
import scanpy as sc

def map_marker_genes(merged, save = False):
    
    # Define the cluster to cell type mapping
    cluster_cell_type_mapping = {
        '0': 'Mitochondrial-associated cells',
        '1': 'Long non-coding RNA expressing cells',
        '2': 'Transcriptionally active cells',
        '3': 'Metabolic active cells',
        '4': 'Secretory cells',
        '5': 'Cell cycle and proliferation cells'
    }

    # Create DataFrame of clusters, top marker genes, and cell types
    cluster_marker_genes = {
        'Cluster': ['0', '1', '2', '3', '4', '5'],
        'Marker Genes': [
            'BANF1, NUTF2, MRPL12, COA3, CCNI',
            'MALAT1, MT-ND1, AC127164.1, SHANK2-AS1, MRTFA-AS1',
            'ITGA9-AS1, GRINA, TCEA2, NAE1, SMS',
            'EXOSC7, ACOT9, AC025754.2, ATP5ME, AFG3L2',
            'PDIA6, ASPH, FOXJ1, OCIAD1, MTRNR2L12',
            'CDCA5, PIGK, ITSN1, AC015871.1, MAPK9'
        ],
        'Cell Type': [cluster_cell_type_mapping[str(i)] for i in range(6)]
    }

    # Create the DataFrame
    cluster_df = pd.DataFrame(cluster_marker_genes)
    cluster_df

    # Save the DataFrame to the 'tables' folder
    if save:
        cluster_df.to_csv('tables/clusters_marker_genes_cell_types.csv', index=False)


    # Add a new column with the identified cell types to the original data
    merged.obs['Cell Type'] = [cluster_cell_type_mapping[str(i)] for i in merged.obs['leiden']]

    return merged
