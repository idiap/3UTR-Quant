# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Zhi Ming Xu <zhiming.xu@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# Description: Remove overlapping 3'UTR annotation based on order of priority
# -----------------------------------------------------------------------------

rule process_annot:
    input:
        annot_focus=config["annot_focus_raw"],
    output:
        annot_filt=config["annot_focus"]
    params:
        ensembl=config["ensembl"],
        ref_seq=config["refseq"]
    conda:
        "bulkanalysis"
    shell:
        """
            python ../envs/genomemanager/scripts/create_focused_ground_truth_no_overlaps.py -a {input.annot_focus} -r {params.ref_seq} -e {params.ensembl} -o {output.annot_filt}
        """
