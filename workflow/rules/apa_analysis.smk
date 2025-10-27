# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Zhi Ming Xu <zhiming.xu@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# Description: Calulcation of relative 3'UTR isoform usage.
# -----------------------------------------------------------------------------

#Calculates PUD and RUD based on normalized counts
rule apa_quant:
    input:
        df_counts=config["results_dir"]+'/df_counts.csv',
        annot_focus=config["annot_focus"],
        annot=config["annot"],
        gene_out="../results/genes_preprocessing/df_genes_filtered_normalized.csv"
    output:
        df_out="../results/3UTR_Analysis/df_counts_filtered_corrected_normalized.csv",
        df_out_pairs="../results/3UTR_Analysis/table_of_pairs_simplified.csv",
    params:
        out_path="../results/3UTR_Analysis/",
        config_file="../config/bulkapa_template.yaml"
    conda:
        "apatools"
    shell:
        """
            python3 ../envs/apatools/scripts/bulkapa.py --config_file {params.config_file} --df_counts_3UTR_path {input.df_counts} --df_filtered_genes_path {input.gene_out} --gtf_file_focus {input.annot_focus} --gtf_file_no_focus {input.annot} --path_to_results {params.out_path}
        """