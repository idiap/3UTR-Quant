# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# -----------------------------------------------------------------------------

import argparse
import os
import sys

import matplotlib
import numpy as np
import pandas as pd
import yaml
from bulkanalysis.filtering import filter_samples_based_on_thresholds
from bulkanalysis.filtering import gmm_threshold
from bulkanalysis.loader import bulk_aggregate_featureCounts_outputs
from bulkanalysis.normalization import quantile_normalize
from bulkanalysis.plot import plot_gene_expression_distribution

from apatools.apanalysis import create_table_of_pairs
from apatools.apanalysis import get_significant_shifts_given_fisher_and_scores
from apatools.apanalysis import plot_score_boxplot
from apatools.filtering import add_length_info
from apatools.filtering import bulk_correction
from apatools.filtering import check_proximal_abundances_in_raw_isoforms_counts
from apatools.filtering import correct_counts_for_small_proximals
from apatools.filtering import get_multi_utrs_isoforms_from_reliably_expressed_genes
from apatools.filtering import remove_length_outliers
from apatools.filtering import stats_raw

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


def readjust_treatments(actual_samples, treatments):
    for treatment, samples in treatments.items():
        for sample in samples:
            if sample not in actual_samples:
                samples.remove(sample)
        treatments[treatment] = samples
    return treatments


def load_from_featureCounts(path_to_counts, data_list: str = None, sample_names=None):
    if os.path.exists(path_to_counts):
        df_counts = pd.read_csv(path_to_counts, index_col=0)
    elif data_list is None:
        raise ValueError(
            "Enter either a valid path to the matrix of counts, \
                         or a list of .txt files output from featureCounts to be aggregated."
        )
    else:
        df_counts = bulk_aggregate_featureCounts_outputs(
            files=data_list, names=sample_names, save_path=path_to_counts
        )

    return df_counts


def main():
    FLAGS = argparse.ArgumentParser(description="Length analysis for bulk RNAseq data.")
    FLAGS.add_argument("--config_file", help="Location of config file")
    FLAGS.add_argument(
        "-df_3UTR",
        "--df_counts_3UTR_path",
        help="<Required> Path to df counts (3UTR)",
        required=True,
    )    
    FLAGS.add_argument(
        "-df_genes",
        "--df_filtered_genes_path",
        help="<Required> Path to df counts (genes)",
        required=True,
    )
    FLAGS.add_argument(
        "-gtf_focus",
        "--gtf_file_focus",
        help="<Required> Path to gtf file (focus)",
        required=True,
    )
    FLAGS.add_argument(
        "-gtf_file_no_focus",
        "--gtf_file_no_focus",
        help="<Required> Path to gtf file (no focus)",
        required=True,
    )
    FLAGS.add_argument(
        "-path_to_results",
        "--path_to_results",
        help="<Required> Path to results",
        required=True,
    )

    args = FLAGS.parse_args()
    with open(args.config_file, "r") as stream:
        config = yaml.safe_load(stream)

    if not os.path.exists(args.path_to_results):
        os.makedirs(args.path_to_results)

    if not os.path.exists(os.path.join(args.path_to_results, "Figures")):
        os.makedirs(os.path.join(args.path_to_results, "Figures"))

    path_to_figures = os.path.join(os.path.join(args.path_to_results, "Figures"))

    old_stdout = sys.stdout
    log_file = open(os.path.join(args.path_to_results, "message.log"), "w")
    sys.stdout = log_file

    # Load
    if "data_list" in config.keys():
        data_list = config["data_list"]
    else:
        data_list = None
    if "data_list" in config.keys():
        data_list = config["data_list"]
        if "sample_names" not in config:
            sample_names = [
                file.split("_")[-1].split(".")[0] for file in config["data_list"]
            ]
        else:
            sample_names = config["sample_names"]
    else:
        data_list = None
        sample_names = None

    df_counts = load_from_featureCounts(
        path_to_counts=args.df_counts_3UTR_path,
        data_list=data_list,
        sample_names=sample_names,
    )
    df_filtered_genes = pd.read_csv(args.df_filtered_genes_path, index_col=0)
    print(f"df_counts columns {df_counts.columns}")
    print(f"df_genes_columns {df_filtered_genes.columns}")
    df_counts = df_counts[df_filtered_genes.columns]

    treatments = readjust_treatments(
        actual_samples=list(df_counts.columns), treatments=config["treatments"]
    )

    print(f"Before filtering, we have {len(df_counts)} isoforms.")
    reps = df_counts.columns

    # Correct counts for isoforms whose annotation is smaller than 300nt length
    df_counts = correct_counts_for_small_proximals(
        df_counts=df_counts,
        reps=reps,
        gtf_file_quantification=args.gtf_file_focus,
        focus_nt=300,
        minimal_required_length=100,
    )

    # Remove the genes with outlier lengths (probable annotation error)
    df_counts = remove_length_outliers(df_counts[reps])

    stats_raw(
        df_counts,
        save_file=os.path.join(
            path_to_figures, f"reads_per_sample.{config['figures_extension']}"
        ),
    )

    # Keep the reliably expressed multi-UTR genes only
    df_counts_iso = get_multi_utrs_isoforms_from_reliably_expressed_genes(
        df_counts, df_filtered_genes
    )

    print(
        f"Before filtering, we have {len(df_counts_iso)} isoforms and {len(df_counts_iso.groupby('gene').sum())} genes."
    )

    df_counts_iso = add_length_info(
        df_counts_iso, path_gtf_no_length_focus=args.gtf_file_no_focus
    )

    _, df_counts_iso = check_proximal_abundances_in_raw_isoforms_counts(
        df_counts_iso=df_counts_iso,
        path_to_results=args.path_to_results,
        reps=reps,
    )

    # Correct the counts according to bulk profile
    df_counts_iso_corr = bulk_correction(
        df_counts_iso,
        columns_to_correct=list(reps),
        replace_col=True,
        gene_id="gene",
        length_id="length",
        clip=True,
    )
    print(df_counts_iso_corr.columns)

    df_counts_iso_log = np.log2(
        df_counts_iso_corr.drop(columns=["gene", "length", "strand"], axis=1).set_index(
            "transcript_id"
        )
        + 1
    )
    # Compute the thresholds for reliably expressed ISOFORMS for each sample
    thresholds = {}
    for col in df_counts_iso_log.columns:
        threshold = gmm_threshold(
            data=np.array(df_counts_iso_log[df_counts_iso_log[col] != 0][col]),
            quantile_threshold=config["quantile_threshold"],
            background_quantile=config["background_quantile"],
            density=True,
            show_plot=False,
            data_name=col,
            saving_folder=path_to_figures,
            extension=config["figures_extension"],
            return_int=False,
        )
        thresholds[col] = threshold
    print(thresholds)
    if "FOSL2.6" in df_counts_iso_corr.index:
        print("We have FOSL2.6 before filtering")
    else:
        print("We do not have FOSL2.6 before filtering")

    df_counts_iso_corr_filtered = filter_samples_based_on_thresholds(
        df_genes=df_counts_iso_corr,
        treatments=treatments,
        thresholds=thresholds,
        pct_in_treatment=config["pct_in_treatment"],
        save_file=os.path.join(
            path_to_figures,
            f"isoform_exp_after_filtering.{config['figures_extension']}",
        ),
    )

    if "FOSL2.6" in df_counts_iso_corr_filtered.index:
        print("We have FOSL2.6 after filtering")
    else:
        print("We do not have FOSL2.6 after filtering")

    print(f"{df_counts_iso_corr_filtered.shape[0]} isoforms left after filtering.")

    df_counts_iso_corr_filtered.set_index("transcript_id", inplace=True)
    print(df_counts_iso_corr_filtered.head())

    if "quantile_normalization" in config.keys():
        quantile_norm = config["quantile_normalization"]
    else:
        quantile_norm = True

    if quantile_norm:
        df_counts_iso_corr_normalized = quantile_normalize(
            df_counts_iso_corr_filtered[reps]
        )
        plot_gene_expression_distribution(
            df_counts_iso_corr_normalized,
            title="Isoform expression after filtering and normalization",
            save_file=os.path.join(
                path_to_figures,
                f"isoform_expression_after_filtering_quantile.{config['figures_extension']}",
            ),
        )
        df_counts_iso_corr_normalized["gene"] = [
            x[:-2] for x in df_counts_iso_corr_normalized.index
        ]

        df_counts_iso_corr_normalized = add_length_info(
            df_counts_iso_corr_normalized,
            path_gtf_no_length_focus=args.gtf_file_no_focus,
        )

    else:
        df_counts_iso_corr_normalized = df_counts_iso_corr_filtered.copy()

    table_of_pairs = create_table_of_pairs(
        df_counts_iso_corr_normalized.reset_index(),
        treatments=treatments,
        replicates=False,
    )

    print(f"We analyze {len(table_of_pairs)} pairs of isoforms.")

    for score in ["PUD", "RUD", "DUD"]:
        plot_score_boxplot(
            table_of_pairs=table_of_pairs,
            treatments=treatments,
            score=score,
            statistic="wilcoxon",
            save_file=os.path.join(
                path_to_figures,
                f"{score}_boxplot_distributions.{config['figures_extension']}",
            ),
        )

    table_of_pairs.to_csv(
        os.path.join(args.path_to_results, "table_of_pairs_simplified.csv")
    )

    if "thresh_FDR" in config.keys():
        thresh_FDR = config["thresh_FDR"]
    else:
        thresh_FDR = 0.01

    if "thresh_PUD" in config.keys():
        thresh_PUD = config["thresh_PUD"]
    else:
        thresh_PUD = 0.15

    if "thresh_RUD" in config.keys():
        thresh_RUD = config["thresh_RUD"]
    else:
        thresh_RUD = 1

    if "thresh_pval" in config.keys():
        thresh_pval = config["thresh_pval"]
    else:
        thresh_pval = None

    table_of_pairs = get_significant_shifts_given_fisher_and_scores(
        table_of_pairs,
        treatments=treatments,
        path_to_results=args.path_to_results,
        scores=["PUD", "RUD"],
        extension=config["figures_extension"],
        thresh_FDR=thresh_FDR,
        thresh_PUD=thresh_PUD,
        thresh_RUD=thresh_RUD,
        thresh_pval=thresh_pval,
    )

    table_of_pairs.to_csv(os.path.join(args.path_to_results, "table_of_pairs.csv"))

    print(
        f"{df_counts_iso_corr_filtered.shape[0]} isoforms will be saved in the table."
    )
    df_counts_iso_corr_normalized.to_csv(
        os.path.join(
            args.path_to_results, "df_counts_filtered_corrected_normalized.csv"
        )
    )

    sys.stdout = old_stdout

    log_file.close()


if __name__ == "__main__":
    main()
