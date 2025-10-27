import csv
import itertools

import pandas as pd
from gtfparse import read_gtf


def reorder_gtf(df_gtf):
    return df_gtf[
        [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        ]
    ]


def focus_annotation_focus(
    annotation, focus_type: str = "end", nt: int = 300, save_path: str = None
):
    def end_focus_coordinates(strand, start, end, nt):
        if (end - start) <= nt:
            start_focus = start
            end_focus = end
        else:
            if strand == "+":
                start_focus = end - nt
                end_focus = end
            else:
                start_focus = start
                end_focus = start + nt
        return start_focus, end_focus

    annotation_focus = annotation.copy()
    if focus_type == "end":
        annotation_focus["start_focus"] = annotation_focus.apply(
            lambda row: end_focus_coordinates(
                row["strand"], row["start"], row["end"], nt
            )[0],
            axis=1,
        )
        annotation_focus["end_focus"] = annotation_focus.apply(
            lambda row: end_focus_coordinates(
                row["strand"], row["start"], row["end"], nt
            )[1],
            axis=1,
        )
    else:
        annotation_focus["end_focus"] = annotation_focus.apply(
            lambda row: end_focus_coordinates(
                row["strand"], row["start"], row["end"], nt
            )[0],
            axis=1,
        )
        annotation_focus["start_focus"] = annotation_focus.apply(
            lambda row: end_focus_coordinates(
                row["strand"], row["start"], row["end"], nt
            )[1],
            axis=1,
        )
    annotation_focus.drop(["start", "end"], axis=1, inplace=True)
    annotation_focus.rename(
        columns={"start_focus": "start", "end_focus": "end"}, inplace=True
    )
    annotation_focus = format_dataframe_as_gtf(annotation_focus)

    if save_path is not None:
        annotation_focus.to_csv(
            save_path,
            index=False,
            sep="\t",
            quoting=csv.QUOTE_NONE,
            escapechar="\\",
            header=False,
        )
    return annotation_focus


def adjust_positions_prox(strand, start_annot, end_annot, predicted_proximal_APA):
    if strand == "+":
        start = start_annot
        end = predicted_proximal_APA
    else:
        start = predicted_proximal_APA
        end = end_annot
    return start, end


def reannotate_transcripts(
    df,
    gene_id: str = "gene_id",
    strand_col: str = "strand",
    start_col: str = "start",
    end_col: str = "end",
):
    df.reset_index(drop=True, inplace=True)
    multiple_pas_genes = list(
        df.groupby("gene_id")
        .count()[df.groupby("gene_id").count()["seqname"] > 1]
        .index
    )
    for gene in multiple_pas_genes:
        strand = list(df[df[gene_id] == gene][strand_col].values)[0]
        if strand == "-":
            for j, (i, row) in enumerate(
                df[df[gene_id] == gene]
                .sort_values(start_col, ascending=False)
                .iterrows()
            ):
                df.loc[i, "transcript_id"] = f"{df.loc[i, 'gene_id']}.{str(j+1)}"
        else:
            for j, (i, row) in enumerate(
                df[df[gene_id] == gene].sort_values(end_col, ascending=True).iterrows()
            ):
                df.loc[i, "transcript_id"] = f"{df.loc[i, 'gene_id']}.{str(j+1)}"
    return df


def get_gtf_annotation_from_DaPars_output(df: pd.DataFrame, save_path=None):
    df["gene_id"] = df["Gene"].apply(lambda x: x.split("|")[1])
    df["strand"] = df["Gene"].apply(lambda x: x.split("|")[-1])
    df["start"] = df["Loci"].apply(lambda x: int(x.split("-")[0].split(":")[1]))
    df["end"] = df["Loci"].apply(lambda x: int(x.split("-")[-1]))
    df_dapars_proxs = df.copy()
    df_dapars_proxs["start_"] = df_dapars_proxs.apply(
        lambda row: adjust_positions_prox(
            row["strand"], row["start"], row["end"], row["Predicted_Proximal_APA"]
        )[0],
        axis=1,
    )
    df_dapars_proxs["end_"] = df_dapars_proxs.apply(
        lambda row: adjust_positions_prox(
            row["strand"], row["start"], row["end"], row["Predicted_Proximal_APA"]
        )[1],
        axis=1,
    )
    df_dapars_proxs.drop(["start", "end"], axis=1, inplace=True)
    df_dapars_proxs.rename(columns={"start_": "start", "end_": "end"}, inplace=True)

    df_dapars_dist = df.copy()

    dapars_annotation = pd.concat(
        [df_dapars_proxs.reset_index(drop=True), df_dapars_dist.reset_index(drop=True)]
    )

    dapars_annotation["seqname"] = dapars_annotation["Loci"].apply(
        lambda x: x.split(":")[0]
    )
    dapars_annotation["feature"] = "exon"
    dapars_annotation["source"] = "rtracklayer"
    dapars_annotation["frame"] = "."
    dapars_annotation["score"] = "."

    dapars_annotation["transcript_id"] = dapars_annotation["gene_id"] + ".1"
    dapars_annotation = reannotate_transcripts(dapars_annotation)

    dapars_annotation["attributes"] = dapars_annotation.apply(
        lambda row: f"uniqueID \"{row['transcript_id']}\"; \
        iso {str(row['transcript_id'].split('.')[-1])}; \
        transcript_id \"{row['transcript_id']}\"; \
        gene_id \"{row['gene_id']}\"; \
        gene_name \"{row['gene_id']}\"; \
        gID \"{row['gene_id']}\"; \
        ord {str(row['transcript_id'].split('.')[-1])}",
        axis=1,
    )
    dapars_annotation = dapars_annotation[
        [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        ]
    ]
    dapars_annotation.drop_duplicates(["start", "end"], keep="first", inplace=True)

    if save_path is not None:
        dapars_annotation.to_csv(
            save_path,
            index=False,
            sep="\t",
            quoting=csv.QUOTE_NONE,
            escapechar="\\",
            header=False,
        )
    return dapars_annotation


def format_dataframe_as_gtf(df):
    def create_attribute_string(row, cols):
        for i, col in enumerate(cols):
            if i == 0:
                attribute = f'{col} "{row[col]}";'
            else:
                attribute = f'{attribute} {col} "{row[col]}";'
        return attribute

    cols = df.columns
    attributes_cols = [
        col
        for col in cols
        if col
        not in [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
        ]
    ]

    df_formatted = df.copy()
    df_formatted["attributes"] = df_formatted.apply(
        lambda row: create_attribute_string(row, attributes_cols), axis=1
    )
    df_formatted.drop(
        attributes_cols,
        axis=1,
        inplace=True,
    )

    df_formatted["score"] = "."
    df_formatted["frame"] = "."

    df_formatted = reorder_gtf(df_formatted)
    return df_formatted


def remove_overlaps_in_annotations(annotation_path, save_path: str = None):
    annotation = read_gtf(annotation_path)
    annotation_pos = annotation[annotation["strand"] == "+"]
    annotation_neg = annotation[annotation["strand"] == "-"]

    for gene in list(set(list(annotation_pos["gene_id"].values))):
        # positive strand
        df_gene = annotation_pos[annotation_pos["gene_id"] == gene].sort_values(
            "start", ascending=True
        )
        for j, (i, row) in enumerate(df_gene.iterrows()):
            if j == 0:
                first_end = row["end"]
            else:
                if row["start"] < first_end:
                    annotation_pos.drop(i, inplace=True)
                else:
                    first_end = row["end"]

    # negative strand
    for gene in list(set(list(annotation_neg["gene_id"].values))):
        df_gene = annotation_neg[annotation_neg["gene_id"] == gene].sort_values(
            "start", ascending=False
        )
        for j, (i, row) in enumerate(df_gene.iterrows()):
            if j == 0:
                first_end = row["start"]
            else:
                if row["end"] > first_end:
                    annotation_neg.drop(i, inplace=True)
                else:
                    first_end = row["start"]

    annotation_filtered = format_dataframe_as_gtf(
        pd.concat([annotation_pos, annotation_neg])
    )

    if save_path is not None:
        annotation_filtered.to_csv(
            save_path,
            index=False,
            sep="\t",
            quoting=csv.QUOTE_NONE,
            escapechar="\\",
            header=False,
        )

    return annotation_filtered


# Get the isoforms who overlap with ground truth
def check_overlap(df_, row_to_test):
    df = df_.drop(row_to_test.name)
    overlap = df.apply(
        lambda x: (x["start"] >= row_to_test["start"])
        & (x["start"] <= row_to_test["end"])
        | (x["end"] >= row_to_test["start"]) & (x["end"] <= row_to_test["end"]),
        axis=1,
    )
    if len(overlap) == 0:
        list_overlaps = []
    else:
        list_overlaps = list(df[overlap]["uniqueID"].values)
    return row_to_test["uniqueID"], list_overlaps


def keep_annotated_isoforms(myUTRs, ends_to_check, end_threshold=50, strand="+"):
    annotated_isoforms = []
    if strand == "+":
        for end in ends_to_check:
            is_close_to_end = (myUTRs["end"] >= end - end_threshold) & (
                myUTRs["end"] <= end + end_threshold
            )
            row_close_to_end = myUTRs[is_close_to_end]

            if not row_close_to_end.empty:
                annotated_isoforms.append(row_close_to_end["uniqueID"].values[0])

        return myUTRs[myUTRs["uniqueID"].isin(annotated_isoforms)].sort_values(
            "start", ascending=False
        )
    else:
        for end in ends_to_check:
            is_close_to_end = (myUTRs["start"] >= end - end_threshold) & (
                myUTRs["start"] <= end + end_threshold
            )
            row_close_to_end = myUTRs[is_close_to_end]

            if not row_close_to_end.empty:
                annotated_isoforms.append(row_close_to_end["uniqueID"].values[0])

        return myUTRs[myUTRs["uniqueID"].isin(annotated_isoforms)].sort_values(
            "start", ascending=True
        )


def get_ground_truth_one_gene(
    gene_name, refSeq, ensembl, myUTRs_focus, end_threshold=50, strand="+"
):
    if strand == "+":
        end_name = "end"
        ascending = False
    else:
        end_name = "start"
        ascending = True
    ref_seq_end = refSeq[
        (refSeq["gene_id"] == gene_name) & (refSeq["feature"] == "3UTR")
    ][end_name]
    ensembl_end = ensembl[
        (ensembl["gene_name"] == gene_name) & (ensembl["feature"] == "UTR")
    ][end_name]
    myUTRS = myUTRs_focus[myUTRs_focus["gene_id"] == gene_name]

    # Concatenate the two series
    concatenated_series = pd.concat([ref_seq_end, ensembl_end])

    # sort the values in descending order and get the unique values
    ends_to_check = concatenated_series.sort_values(ascending=ascending).unique()

    # Create ground truth -- keep annotated isoforms
    annots_all = keep_annotated_isoforms(
        myUTRS, ends_to_check, end_threshold, strand=strand
    )

    annots_ensembl = keep_annotated_isoforms(
        myUTRS,
        ensembl_end.sort_values(ascending=ascending).unique(),
        end_threshold,
        strand=strand,
    )

    # Get the isoforms who overlap with ground truth
    overlapping_isoforms = {}
    for _, row in annots_all.iterrows():
        key, vals = check_overlap(myUTRS, row)
        overlapping_isoforms[key] = vals
    # to_remove = list(itertools.chain.from_iterable(to_remove))
    # Do not throw away the annotated isoforms even if they overlap
    overlapping_not_annots = [
        to_rem
        for to_rem in list(
            itertools.chain.from_iterable(list(overlapping_isoforms.values()))
        )
        if to_rem not in annots_all["uniqueID"].values
    ]
    overlapping_annots = [
        to_rem
        for to_rem in list(
            itertools.chain.from_iterable(list(overlapping_isoforms.values()))
        )
        if to_rem in annots_all["uniqueID"].values
    ]

    # If there are overlapping annotations; keep the Ensembl one;
    # or the longest one if they belong to the same annotation
    isoforms_to_keep = []
    isoforms_to_remove = []
    for overlap_iso in overlapping_annots:
        list_to_test = [overlap_iso] + overlapping_isoforms[overlap_iso]
        for i, iso in enumerate(list_to_test):
            if i == 0:
                iso_to_keep = iso
            else:
                is_ensembl_keep = iso_to_keep in annots_ensembl["uniqueID"].values
                is_ensembl_to_test = iso in annots_ensembl["uniqueID"].values
                if strand == "+":
                    is_longer_to_test = (
                        myUTRS[myUTRS["uniqueID"] == iso]["end"].values[0]
                        > myUTRS[myUTRS["uniqueID"] == iso_to_keep]["end"].values[0]
                    )
                else:
                    is_longer_to_test = (
                        myUTRS[myUTRS["uniqueID"] == iso]["start"].values[0]
                        < myUTRS[myUTRS["uniqueID"] == iso_to_keep]["start"].values[0]
                    )
                if is_ensembl_keep:
                    if is_ensembl_to_test:
                        if is_longer_to_test:
                            iso_to_keep = iso
                else:
                    if is_ensembl_to_test:
                        iso_to_keep = iso
        isoforms_to_keep.append(iso_to_keep)
        isoforms_to_remove += [iso for iso in list_to_test if iso != iso_to_keep]
    isoforms_to_keep = list(set(isoforms_to_keep))
    isoforms_to_remove = list(set(isoforms_to_remove))

    # Remove overlapping isoforms not annotated
    ground_truth = myUTRS[~myUTRS["uniqueID"].isin(overlapping_not_annots)].sort_values(
        "start", ascending=ascending
    )
    # Remove overlapping isoforms annotated
    ground_truth = ground_truth[
        ~ground_truth["uniqueID"].isin(isoforms_to_remove)
    ].sort_values("start", ascending=ascending)

    # Now remove the remaining overlapping isoforms

    if strand == "+":
        for j, (i, row) in enumerate(ground_truth.iterrows()):
            if j == 0:
                last_start = row["start"]
            else:
                if row["end"] > last_start:
                    ground_truth.drop(i, inplace=True)
                else:
                    last_start = row["start"]
    else:
        for j, (i, row) in enumerate(ground_truth.iterrows()):
            if j == 0:
                last_start = row["end"]
            else:
                if row["start"] < last_start:
                    ground_truth.drop(i, inplace=True)
                else:
                    last_start = row["end"]

    return (
        ground_truth,
        myUTRS[~myUTRS["uniqueID"].isin(ground_truth["uniqueID"].values)],
    )


def get_ground_truth_one_gene_side_prox(
    gene_name, refSeq, ensembl, myUTRs_focus, end_threshold=50, strand="+"
):
    if strand == "+":
        end_name = "end"
        ascending = True
    else:
        end_name = "start"
        ascending = False
    ref_seq_end = refSeq[
        (refSeq["gene_id"] == gene_name) & (refSeq["feature"] == "3UTR")
    ][end_name]
    ensembl_end = ensembl[
        (ensembl["gene_name"] == gene_name) & (ensembl["feature"] == "UTR")
    ][end_name]
    myUTRS = myUTRs_focus[myUTRs_focus["gene_id"] == gene_name]

    # Concatenate the two series
    concatenated_series = pd.concat([ref_seq_end, ensembl_end])

    # sort the values in descending order and get the unique values
    ends_to_check = concatenated_series.sort_values(ascending=ascending).unique()

    # Create ground truth -- keep annotated isoforms
    annots_all = keep_annotated_isoforms(
        myUTRS, ends_to_check, end_threshold, strand=strand
    )

    annots_ensembl = keep_annotated_isoforms(
        myUTRS,
        ensembl_end.sort_values(ascending=ascending).unique(),
        end_threshold,
        strand=strand,
    )

    # Get the isoforms who overlap with ground truth
    overlapping_isoforms = {}
    for _, row in annots_all.iterrows():
        key, vals = check_overlap(myUTRS, row)
        overlapping_isoforms[key] = vals
    # to_remove = list(itertools.chain.from_iterable(to_remove))
    # Do not throw away the annotated isoforms even if they overlap
    overlapping_not_annots = [
        to_rem
        for to_rem in list(
            itertools.chain.from_iterable(list(overlapping_isoforms.values()))
        )
        if to_rem not in annots_all["uniqueID"].values
    ]
    overlapping_annots = [
        to_rem
        for to_rem in list(
            itertools.chain.from_iterable(list(overlapping_isoforms.values()))
        )
        if to_rem in annots_all["uniqueID"].values
    ]

    # If there are overlapping annotations; keep the Ensembl one;
    # or the longest one if they belong to the same annotation
    isoforms_to_keep = []
    isoforms_to_remove = []
    for overlap_iso in overlapping_annots:
        list_to_test = [overlap_iso] + overlapping_isoforms[overlap_iso]
        for i, iso in enumerate(list_to_test):
            if i == 0:
                iso_to_keep = iso
            else:
                is_ensembl_keep = iso_to_keep in annots_ensembl["uniqueID"].values
                is_ensembl_to_test = iso in annots_ensembl["uniqueID"].values
                if strand == "+":
                    is_longer_to_test = (
                        myUTRS[myUTRS["uniqueID"] == iso]["end"].values[0]
                        > myUTRS[myUTRS["uniqueID"] == iso_to_keep]["end"].values[0]
                    )
                else:
                    is_longer_to_test = (
                        myUTRS[myUTRS["uniqueID"] == iso]["start"].values[0]
                        < myUTRS[myUTRS["uniqueID"] == iso_to_keep]["start"].values[0]
                    )
                if is_ensembl_keep:
                    if is_ensembl_to_test:
                        if is_longer_to_test:
                            iso_to_keep = iso
                else:
                    if is_ensembl_to_test:
                        iso_to_keep = iso
        isoforms_to_keep.append(iso_to_keep)
        isoforms_to_remove += [iso for iso in list_to_test if iso != iso_to_keep]
    isoforms_to_keep = list(set(isoforms_to_keep))
    isoforms_to_remove = list(set(isoforms_to_remove))

    # Remove overlapping isoforms not annotated
    ground_truth = myUTRS[~myUTRS["uniqueID"].isin(overlapping_not_annots)].sort_values(
        "start", ascending=ascending
    )
    # Remove overlapping isoforms annotated
    ground_truth = ground_truth[
        ~ground_truth["uniqueID"].isin(isoforms_to_remove)
    ].sort_values("start", ascending=ascending)

    # Now remove the remaining overlapping isoforms

    if strand == "+":
        for j, (i, row) in enumerate(ground_truth.iterrows()):
            if j == 0:
                first_end = row["end"]
            else:
                if row["start"] < first_end:
                    ground_truth.drop(i, inplace=True)
                else:
                    first_end = row["end"]
    else:
        for j, (i, row) in enumerate(ground_truth.iterrows()):
            if j == 0:
                first_end = row["start"]
            else:
                if row["end"] > first_end:
                    ground_truth.drop(i, inplace=True)
                else:
                    first_end = row["start"]

    return (
        ground_truth,
        myUTRS[~myUTRS["uniqueID"].isin(ground_truth["uniqueID"].values)],
    )


def keep_annotated_and_remove_overlaps(
    annotation_path,
    refSeq_path,
    ensembl_path,
    remove_overlaps_starting_from_dist=True,
    save_path: str = None,
):
    # Read the annotations
    annotation = read_gtf(annotation_path).to_pandas()
    refSeq = read_gtf(refSeq_path).to_pandas()
    ensembl = read_gtf(ensembl_path).to_pandas()

    # Separate into positive and negative strand
    annotation_pos = annotation[annotation["strand"] == "+"].copy()
    annotation_neg = annotation[annotation["strand"] == "-"].copy()

    if remove_overlaps_starting_from_dist:
        # positive strand
        for gene in list(set(list(annotation_pos["gene_id"].values))):
            _, inverse_ground_truth = get_ground_truth_one_gene(
                gene, refSeq, ensembl, annotation_pos, end_threshold=30, strand="+"
            )
            annotation_pos.drop(inverse_ground_truth.index, inplace=True)

        # negative strand
        for gene in list(set(list(annotation_neg["gene_id"].values))):
            _, inverse_ground_truth = get_ground_truth_one_gene(
                gene, refSeq, ensembl, annotation_neg, end_threshold=30, strand="-"
            )
            annotation_neg.drop(inverse_ground_truth.index, inplace=True)

    else:
        # positive strand
        for gene in list(set(list(annotation_pos["gene_id"].values))):
            _, inverse_ground_truth = get_ground_truth_one_gene_side_prox(
                gene, refSeq, ensembl, annotation_pos, end_threshold=30, strand="+"
            )
            annotation_pos.drop(inverse_ground_truth.index, inplace=True)

        # negative strand
        for gene in list(set(list(annotation_neg["gene_id"].values))):
            _, inverse_ground_truth = get_ground_truth_one_gene_side_prox(
                gene, refSeq, ensembl, annotation_neg, end_threshold=30, strand="-"
            )
            annotation_neg.drop(inverse_ground_truth.index, inplace=True)

    annotation_filtered = format_dataframe_as_gtf(
        pd.concat([annotation_pos, annotation_neg])
    )

    if save_path is not None:
        annotation_filtered.to_csv(
            save_path,
            index=False,
            sep="\t",
            quoting=csv.QUOTE_NONE,
            escapechar="\\",
            header=False,
        )

    return annotation_filtered
