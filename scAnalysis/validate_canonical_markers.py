'''
This module defines a dictionary (canonical_markers) containing known canonical markers for various cell types.
It then identifies the intersecting genes between the scRNA-seq dataset and these canonical markers.
The dataset is filtered to include only these intersecting genes.
Finally, a dot plot is generated to visualize the expression levels of the known canonical markers across different cell types,
and the plot is saved as '_canonical_markers.png'.
'''

import numpy as np
import pandas as pd
import scanpy as sc

def validate_top_markers(merged, save = False):
# Define a dictionary of known canonical markers for each cell type
    canonical_markers = {
        'Mitochondrial-associated cells': ['BANF1', 'NUTF2', 'MRPL12', 'COA3', 'CCNI'],
        'Long non-coding RNA expressing cells': ['MALAT1', 'MT-ND1', 'AC127164.1', 'SHANK2-AS1', 'MRTFA-AS1'],
        'Transcriptionally active cells': ['ITGA9-AS1', 'GRINA', 'TCEA2', 'NAE1', 'SMS'],
        'Metabolic active cells': ['EXOSC7', 'ACOT9', 'AC025754.2', 'ATP5ME', 'AFG3L2'],
        'Secretory cells': ['PDIA6', 'ASPH', 'FOXJ1', 'OCIAD1', 'MTRNR2L12'],
        'Cell cycle and proliferation cells': ['CDCA5', 'PIGK', 'ITSN1', 'AC015871.1', 'MAPK9']
    }

    # Get the intersecting genes between the dataset and canonical markers
    intersecting_genes = set(merged.var_names) & set([gene for genes in canonical_markers.values() for gene in genes])

    # Filter the dataset with the intersecting genes
    filtered_merged = merged[:, list(intersecting_genes)]

    # Plot the expression levels of known canonical markers
    if save:
        sc.pl.dotplot(filtered_merged, groupby='Cell Type', var_names=list(intersecting_genes), color_map='Blues', save='_canonical_markers.png')
    else:
        sc.pl.dotplot(filtered_merged, groupby='Cell Type', var_names=list(intersecting_genes), color_map='Blues')
    
    return merged
