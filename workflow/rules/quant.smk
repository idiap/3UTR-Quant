# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Zhi Ming Xu <zhiming.xu@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# Description: Use featurecounts to quantify reads within 3'UTRs
# -----------------------------------------------------------------------------

#Run featurecounts
rule quant_1:
    input:
        bam=config["bam_dir"]+'{SRR_ID}.Aligned.sortedByCoord.out.bam',
        annot_filt=config["annot_focus"]
    output:
        counts_out=config["results_dir"]+'hg38_UTRs/myUTRs_hg38_{SRR_ID}.txt'
    params:
        log=config["results_dir"]+"hg38_UTRs/{SRR_ID}.log"
    conda:
        "subread"
    shell:
        """
            featureCounts -a {input.annot_filt} -o {output.counts_out} -g transcript_id -f -O -M -s 2 -p -T 12 --verbose {input.bam} > {params.log}
        """

#Process featurecounts output
rule quant_2:
    input:
        expand(config["results_dir"]+'hg38_UTRs/myUTRs_hg38_{SRR_ID}.txt',SRR_ID=SRR_ID)    
    output:
        out=config["results_dir"]+'df_counts.csv'
    params:
        all_files=' '.join(expand(config["results_dir"]+'hg38_UTRs/myUTRs_hg38_{SRR_ID}.txt',SRR_ID=SRR_ID)),
        all_names=expand("{SRR_ID}",SRR_ID=SRR_ID)
    conda:
        "bulkanalysis"
    shell:
        """
            python3 ../envs/bulkanalysis/scripts/aggregate_featureCounts_output.py -f {params.all_files} -n {params.all_names} -s {output.out}
        """
