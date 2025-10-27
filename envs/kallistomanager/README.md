<!--
SPDX-FileCopyrightText: 2023 Idiap Research Institute <contact@idiap.ch>

SPDX-FileContributor: <Author last name> <Author first name> firstname.lastname@idiap.ch

SPDX-License-Identifier: LicenseRef-ProjectNameID
-->

# Project

This is a package to help with the analysis of kallisto outputs.
One particularly useful application is to load and merge outputs of different samples.

From `path_to_kallisto_output`; the output path of `kallisto quant` (`-o` argument), you can merge all the samples outputs and read it as a pandas DataFrame as follows:

``` python
        merge_kallisto_counts(os.path.join(`path_to_kallisto_output`, "*"),  `path_to_kallisto_output`)

```

## Installation of this package in another project

To install this library in another Python project, execute simply:

```bash
pip install git+ssh://git@gitlab.idiap.ch:genomics/genomics/toolboxes/kallistomanager.git
```

## State of the project

This package is under development. Some functions have been developed by Christelle Schneuwly for analysis of kallisto outputs from SMART-seq data.
Note: the poetry environment has to be solved. In particular scanpy is not installed due to some conflicts.

# How to use the package

## Coverage

This is to visualize in IGV the coverage output from kallisto quant, when we used a transcriptome annotation at the transcript level and not at the chromosome level.

**Step 1:** Include the arrgument `--pseudobam` when running kallisto quant as follow:

```bash
kallisto quant --index path_to_kallisto_index -o path_to_output --pseudobam file1.fastq file2.fastq
```

**Step 2:** Sort the bam file using `samtools` as follow:

```bash
samtools sort pseudoalignments.bam -o pseudoalignments_sorted.bam
```

**Step 3:** Use `igvtools` to create a first coverage file, but based on transcripts instead of chromosomes. Then convert it into a bedfile:

```bash
igvtools count pseudoalignments_sorted.bam transcripts.tdf mygenome.chrom.sizes
igvtools tdftobedgraph transcripts.tdf transcripts.bed
```

**Step 4:** Sort the resulting bed file using `bedtools`:

```bash
bedtools sort transcripts.bed -g genome.transcripts.sizes > transcripts_sorted.bed
```

**Step 5:** Use the script `scripts/coverage.py` to convert the coverage file from transcripts to chromosomes:

```bash
python scripts/coverage.py --annotation_path path_to_gtf --bed_transcripts_path transcripts_sorted.bed --output chromosomes_coverage.bed
```

**Step 6:** Sort the resulting bedfile:

```bash
bedtools sort chromosomes_coverage.bed -g genome.chrom.sizes > chromosomes_coverage_sorted.bedgraph
```
