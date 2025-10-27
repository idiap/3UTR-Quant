# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Zhi Ming Xu <zhiming.xu@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# Description: Identifies expressed genes based on kallisto output. 
# -----------------------------------------------------------------------------

#Filter for expressed genes using kallisto output and GMM
rule gene_quant:
    input:
        config["results_dir"]+'df_counts.csv'
    output:
        gene_out="../results/genes_preprocessing/df_genes_filtered_normalized.csv"
    params:
        df_counts_path=config["bam_dir"]+"kallisto/output/",
        config_file="../config/genes_preprocessing_template.yaml",
        out_dir="../results/genes_preprocessing/",
        gene_names=config["gene_names_path"],
        gene_info=config["gene_info_path"]
    conda:
        "bulkanalysis"
    shell:
        """
            python3 ../envs/bulkanalysis/scripts/genes_preprocessing.py --config_file {params.config_file} --df_counts_path {params.df_counts_path} --gene_names_path {params.gene_names} --gene_info_path {params.gene_info} --path_to_results {params.out_dir}
        """
