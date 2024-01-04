"""
Perform scRNA-seq analysis using OpenAI outputs. Obtain top markers and cell line clusters and perform differential expression analysis.
"""

from scAnalysis.read_and_merge_data import read_and_merge_samples
from scAnalysis.preprocess_and_normalize import preprocess_data
from scAnalysis.identify_top_genes import get_highly_variable_genes
from scAnalysis.cluster_top_genes import perform_clustering
from scAnalysis.rank_and_visualize_markers import rank_gene_groups
from scAnalysis.annotate_clusters import map_marker_genes
from scAnalysis.validate_canonical_markers import validate_top_markers
from scAnalysis.analyze_differential_expression import differential_expression

import argparse
import scanpy as sc

def parse_arguments():
    parser = argparse.ArgumentParser(description='Perform scRNA-seq analysis.')
    parser.add_argument('treatment', type=str, help='Path to the treated sample for analysis.')
    parser.add_argument('control', type=str, help='Path to the control (untreated) sample.')
    parser.add_argument('--save_table', action='store_true', help='Save tables produced in the analysis.')
    parser.add_argument('--save_plot', action='store_true', help='Save plots produced in the analysis.')
	
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    try:
        adata = read_and_merge_samples(args.treatment, args.control)
        adata = preprocess_data(adata)
        adata = get_highly_variable_genes(adata, args.save_plot)
        adata = perform_clustering(adata, args.save_plot)
        adata = rank_gene_groups(adata, args.save_plot, args.save_table)
        adata = map_marker_genes(adata, args.save_table)
        adata = validate_top_markers(adata, args.save_plot)
        adata = differential_expression(adata, args.save_table)
	
        # Save the fully processed data
        adata.write('data/processed_data.h5ad')        
	
    except Exception as e:
        print(f"An error occurred: {e}")
        # Handle or log the exception as needed

if __name__ == '__main__':
    main()
