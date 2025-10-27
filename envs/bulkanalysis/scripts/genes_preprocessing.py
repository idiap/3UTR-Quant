# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# -----------------------------------------------------------------------------

# Script to filter and normalize gene expression across samples from a kallisto quantification 
# on the whole transcriptome
import argparse
import logging
import os
import sys
import yaml
import numpy as np
from bulkanalysis.plot import plot_gene_expression_distribution
from bulkanalysis.filtering import gmm_threshold, filter_samples_based_on_thresholds, remove_outliers_samples
from kallistomanager.loader import merge_kallisto_counts, from_whole_transcriptome_quantification_to_coding_genes
from bulkanalysis.normalization import quantile_normalize
from bulkanalysis.dge import create_design_matrix_for_DGE



def format_genes_from_kallisto_whole_transcriptome(df_counts_path,
                                             gene_names_path,
                                             gene_info_path, data_list=None):
    merge_kallisto_counts(data_path=df_counts_path, data_list=data_list, save_dir=df_counts_path)
    df_genes = from_whole_transcriptome_quantification_to_coding_genes(df_counts_path=os.path.join(df_counts_path, 
                                                                                                   "counts_matrix.csv"),
                                                                                            gene_names_path=gene_names_path,
                                                                                            gene_info_path=gene_info_path)
    return df_genes

def readjust_treatments(actual_samples, treatments):
    for treatment, samples in treatments.items():
        for sample in samples:
            if sample not in actual_samples:
                samples.remove(sample)
        treatments[treatment] = samples
    return treatments


def main():
    FLAGS = argparse.ArgumentParser(description='Length analysis for bulk RNAseq data.')
    FLAGS.add_argument('--config_file', help='Location of config file')
    FLAGS.add_argument(
        "-df_counts_path",
        "--df_counts_path",
        help="<Required> Path to kallisto counts",
        required=True,
    )    
    FLAGS.add_argument(
        "-gene_names",
        "--gene_names_path",
        help="<Required> Path to gene names",
        required=True,
    )
    FLAGS.add_argument(
        "-gene_info",
        "--gene_info_path",
        help="<Required> Path to gene info",
        required=True,
    )
    FLAGS.add_argument(
        "-o",
        "--path_to_results",
        help="<Required> Path to results",
        required=True,
    )

    args = FLAGS.parse_args()
    with open(args.config_file, 'r') as stream:
        config = yaml.safe_load(stream)

    if not os.path.exists(args.path_to_results):
        os.makedirs(args.path_to_results)

    old_stdout = sys.stdout

    log_file = open(os.path.join(args.path_to_results, "message.log"),"w")

    sys.stdout = log_file


    if not os.path.exists(os.path.join(args.path_to_results, 'Figures')):
        os.makedirs(os.path.join(args.path_to_results, 'Figures'))

    path_to_figures = os.path.join(os.path.join(args.path_to_results, 'Figures'))

    if 'data_list' in config.keys():
        data_list = config['data_list']
    else: 
        data_list = None

    if config['data_origin'] == 'kallisto_whole_transcriptome':
        df_genes = format_genes_from_kallisto_whole_transcriptome(df_counts_path=args.df_counts_path,
                                                                           gene_names_path=args.gene_names_path,
                                                                           gene_info_path=args.gene_info_path, 
                                                                           data_list=data_list)
    else:
        raise NotImplementedError("This script does not fit with the provided data origin")
    
    # Here add possibility to aggregate SRR ids from similar samples
    if 'aggregation' in config.keys():
        for key, values in config['aggregation'].items():
            df_genes[key] = df_genes[values].sum(axis=1)
            df_genes.drop(columns=values, inplace=True)

    
    print(f"Before filtering, we have {len(df_genes)} genes.")

    df_genes = df_genes[df_genes.sum(axis=1) != 0]
    plot_gene_expression_distribution(df_genes=df_genes, 
                                      title="Gene expression before filtering", 
                                      save_file=os.path.join(path_to_figures,
                                                             f"gene_exp_before_filtering.{config['figures_extension']}"))
    

    df_genes_log = np.log2(df_genes+1)
    ## Compute the thresholds for reliably expressed genes for each sample
    thresholds = {}
    for col in df_genes_log.columns:

        threshold = gmm_threshold(data= np.array(df_genes_log[df_genes_log[col] != 0][col]), 
                                        quantile_threshold=0.95,
                                        background_quantile=True, 
                                        density=True, show_plot=False, data_name=col, 
                                        saving_folder=path_to_figures, 
                                        extension=config['figures_extension'],
                                        return_int=False)
        thresholds[col] = threshold
    
    df_genes = remove_outliers_samples(df_genes=df_genes,
                                       thresholds=thresholds,
                                       save_path=path_to_figures, 
                                       extension=config['figures_extension'])
    
    treatments = readjust_treatments(actual_samples=df_genes.columns, treatments=config['treatments'])
    print(f"All treatments: {config['treatments']}")
    print(f"All thresholds: {thresholds}")
    print(f"Actual samples: {df_genes.columns}")
    print(f"Actual treatments: {treatments}")

    # Remove the items in thresholds that are not in the actual samples
    thresholds = {key: value for key, value in thresholds.items() if key in df_genes.columns}
    print(f"Actual thresholds: {thresholds}")

    df_genes_filtered = filter_samples_based_on_thresholds(df_genes=df_genes, 
                                                           treatments=treatments, 
                                                           thresholds=thresholds, 
                                                           pct_in_treatment=config['pct_in_treatment'], 
                                                           save_file=os.path.join(path_to_figures, 
                                                                                  f"gene_exp_after_filtering.{config['figures_extension']}"))
   

    

    create_design_matrix_for_DGE(df_counts_genes=df_genes_filtered, treatments=treatments, path_to_results=args.path_to_results)
    df_filtered_normalized = quantile_normalize(df_genes_filtered)
    plot_gene_expression_distribution(df_genes=df_filtered_normalized,
                                      title="Gene expression after filtering and quantile normalization",
                                      save_file=os.path.join(path_to_figures,
                                                             f"gene_exp_after_filtering_and_normalization.{config['figures_extension']}"))
    
    print(f"After filtering and normalization, we have {len(df_filtered_normalized)} remaining genes.")

    df_filtered_normalized.to_csv(os.path.join(args.path_to_results, 'df_genes_filtered_normalized.csv'))


    sys.stdout = old_stdout

    log_file.close()

if __name__ == "__main__":
    main()
