# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# -----------------------------------------------------------------------------

import glob
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def human_format(num):
    num = float("{:.3g}".format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return "{}{}".format(
        "{:f}".format(num).rstrip("0").rstrip("."), ["", "K", "M", "B", "T"][magnitude]
    )


def merge_kallisto_counts(save_dir: str, data_list: list = None, data_path: str = None):
    """Concatenate kallisto quant outputs to create a counts matrix.

    Args:
        data_path: Path to the data folder where each folder is named with
        the sample identifier and contains a abundance.tsv file with kallisto's
        quantification results
        save_dir: Path to folder where the concatenated counts matrix will
        be saved
    """
    if data_list is not None:
        kal_dirs = data_list
    elif data_path is not None:
        data_path = os.path.join(data_path, "*/*.tsv")
        kal_dirs = glob.glob(data_path)
    else:
        raise ValueError(
            "Please provide either a valid list of data or a path to data."
        )
    sample_id = [path.split("/")[-2] for path in kal_dirs]
    df = pd.read_csv(kal_dirs[0], sep="\t", index_col=0)
    df_length = pd.DataFrame()
    df_eff_length = pd.DataFrame()
    df_est_counts = pd.DataFrame()
    df_tpm = pd.DataFrame()

    for dir, sample in zip(kal_dirs, sample_id):
        if not dir.endswith(".csv"):
            df = pd.read_csv(dir, sep="\t", index_col=0)
            df_length = pd.concat(
                [df_length, df.length.rename(sample)], axis=1, join="outer"
            )
            df_eff_length = pd.concat(
                [df_eff_length, df.eff_length.rename(sample)], axis=1, join="outer"
            )
            df_est_counts = pd.concat(
                [df_est_counts, df.est_counts.rename(sample)], axis=1, join="outer"
            )
            df_tpm = pd.concat([df_tpm, df.tpm.rename(sample)], axis=1, join="outer")

    df_length.to_csv(f"{save_dir}/gene_length_matrix.csv")
    df_eff_length.to_csv(f"{save_dir}/gene_eff_length_matrix.csv")
    df_est_counts.to_csv(f"{save_dir}/counts_matrix.csv")
    df_tpm.to_csv(f"{save_dir}/tpm_matrix.csv")


def load_raw(data_path: str, metadata_path: str, raw: bool = True):
    """Load counts matrix in csv format to AnnData.

    Args:
        data_path: Path to data folder where the counts_matrix.csv or tmp_matrix.csv is stored
        metadata_path: Path to folder where the metadata SraRunTable.txt is stored
        raw: Whether to use the raw counts or the counts normalized by length (TPM).

    Returns:
        AnnData: Transcripts counts matrix
    """

    metadata = pd.read_csv(f"{metadata_path}/SraRunTable.txt")
    metadata = metadata.set_index("Run")

    if raw:
        adata = sc.read_csv(
            f"{data_path}/counts_matrix.csv", delimiter=",", first_column_names=True
        )
    else:
        adata = sc.read_csv(
            f"{data_path}/tpm_matrix.csv", delimiter=",", first_column_names=True
        )

    adata = adata.transpose()
    adata.obs = pd.concat([adata.obs, metadata[["Cell_Line", "Treatment"]]], axis=1)

    return adata


def from_whole_transcriptome_quantification_to_coding_genes(
    df_counts_path: str, gene_names_path: str, gene_info_path: str
):
    """From a quantification on the whole trancsriptome, returns the dataframe of counts
    aggregated by coding genes.

    Args:
        df_counts_counts_path: path to the kallisto counts matrix
        gene_names_path: path to the file mapping the gene_ids contained in the gtf file to the gene_names.
        The entry should be #geneId and geneName, as when downloaded from GENCODE genes UCSC tables.
        gtf_path: GTF file for the genome annotation used for quantification with kallisto; from ENCODE.
    """
    # Load
    df_counts = pd.read_csv(df_counts_path, index_col=0)
    gene_names = pd.read_csv(gene_names_path, sep="\t")
    gene_info = pd.read_csv(gene_info_path, sep="\t")[["Symbol", "type_of_gene"]]

    # Add gene column
    df_counts["gene"] = (
        df_counts.reset_index()["index"].apply(lambda x: x.split("|")[1]).values
    )

    # Groupby genes
    df_genes = df_counts.groupby("gene").sum()
    df_genes = (
        gene_names.drop_duplicates(keep="first")
        .merge(df_genes, left_on="#geneId", right_index=True)
        .groupby("geneName")
        .sum()
    )

    # Coding genes
    coding_genes = gene_info[gene_info["type_of_gene"] == "protein-coding"][
        "Symbol"
    ].values
    df_genes = df_genes[df_genes.index.isin(coding_genes)]
    return df_genes


def add_gene_info(adata: ad.AnnData, kallisto_path: str, biomart_path: str):
    """Since the Kallisto only quantifies transcripts against the transcriptome
    and not the genome using Ensembl transcripts notations, we had to download
    from Ensembl other IDs and external information.
    Also the Kallisto index we used for quantification is downlodable from
    github with a transcriptome to gene annotation, we just would like to check
    if both notations are equivalent.

    Args:
        adata: Transcripts counts matrix
        kallisto_path: Path to the data folder where kallisto index is stored
        biomart_path: Path to tab delimited file with biomart info

    Returns:
        AnnData: Transcripts counts matrix
    """

    biomart_tab = pd.read_csv(biomart_path, sep="\t", header=0)
    gene_id_tab = pd.read_csv(
        f"{kallisto_path}/transcripts_to_genes.txt",
        sep="\t",
        names=["Transcript_ID", "Gene_ID", "Gene_name"],
    )

    kallisto = gene_id_tab.set_index("Transcript_ID").rename_axis(None, axis=0)
    biomart = biomart_tab.set_index("Transcript stable ID version").rename_axis(
        None, axis=0
    )
    check_same = pd.merge(kallisto, biomart, left_index=True, right_index=True)
    check_same = check_same.rename(
        columns={
            "Gene_ID": "Kallisto_gene_id",
            "Gene_name": "Kallisto_gene_name",
            "Gene stable ID": "bioMart_gene_id",
            "Gene stable ID version": "bioMart_gene_id_ver",
            "Transcript stable ID": "bioMart_transcript_id",
            "Gene name": "bioMart_gene_name",
            "Transcript name": "bioMart_transcript_name",
            "Gene type": "bioMart_gene_type",
            "Transcript type": "bioMart_transcript_type",
            "Chromosome/scaffold name": "bioMart_chr",
        }
    )

    to_add = check_same[
        [
            "Kallisto_gene_id",
            "Kallisto_gene_name",
            "bioMart_gene_id",
            "bioMart_gene_id_ver",
            "bioMart_transcript_id",
            "bioMart_gene_name",
            "bioMart_transcript_name",
            "bioMart_gene_type",
            "bioMart_transcript_type",
            "bioMart_chr",
        ]
    ]

    to_add = to_add.drop_duplicates(keep="first")
    to_add = to_add[
        ~to_add.index.duplicated(keep="first")
    ]  # Ensure that no duplicates remain

    adata.var = adata.var.join(to_add)
    adata.obs.columns = adata.obs.columns.str.replace(" ", "_")
    adata.var.columns = adata.var.columns.str.replace(" ", "_")

    return adata


def agg_per_gene(adata: ad.AnnData):
    """Since the analysis in the original paper is gene-centric, we need to
    aggregate counts by summing them up by gene. To not alter the transcript
    let's create another AnnData object with gene raw gene counts.

    Args:
        adata: Transcripts counts matrix

    Returns:
        AnnData: Gene counts matrix
    """
    adata = adata.copy()

    gene_index = adata.var["Kallisto_gene_id"]
    gene_index = gene_index.fillna(adata.var["bioMart_gene_id"])
    gene_index = gene_index.apply(
        lambda x: x.split(".")[0] if isinstance(x, str) else x
    )

    adata.X[
        adata.X == 0.0
    ] = np.nan  # Remove transcripts that do not seem to belong to any gene.

    gene_temp = pd.DataFrame(
        adata.X, index=adata.obs.index, columns=gene_index.squeeze()
    )
    gene_temp = gene_temp.loc[:, gene_temp.columns.notna()]
    gene_temp = gene_temp.groupby(gene_temp.columns, axis=1, sort=False).sum()

    adata_gene = sc.AnnData(gene_temp)
    adata_gene.obs = adata.obs.copy()

    gene_info = adata.var.drop(
        columns=[
            "bioMart_gene_id_ver",
            "bioMart_transcript_id",
            "bioMart_transcript_name",
            "bioMart_transcript_type",
        ]
    )
    gene_info.index = gene_index.squeeze()
    gene_info = gene_info.drop(columns=["Kallisto_gene_id"]).drop_duplicates()
    gene_info = gene_info.drop(
        index=gene_info[gene_info.index.duplicated(keep=False)].isnull().any(1).index
    )
    gene_info = gene_info[~gene_info.index.duplicated(keep="first")]

    adata_gene.var = adata_gene.var.join(gene_info)
    adata_gene.var["bioMart_chr"] = adata_gene.var["bioMart_chr"].apply(str)
    adata_gene.var["bioMart_chr"] = adata_gene.var["bioMart_chr"].replace("nan", np.nan)
    adata_gene.var.Kallisto_gene_name = adata_gene.var.Kallisto_gene_name.replace(
        "nan", np.nan
    ).fillna(adata_gene.var.index.to_series())

    return adata_gene
