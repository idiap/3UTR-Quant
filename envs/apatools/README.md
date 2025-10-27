# Project

This project aims at creating a package with multiple tools to analyze Alternative Polyadenylation (APA) from RNA sequencing data (bulk and single-cell).

## State

This project is under development.


## Installation

To install this library in another Python project, execute simply:

```bash
pip install git+https://gitlab.idiap.ch/genomics/genomics/toolboxes/apatools.git
```
or
```bash
git clone https://gitlab.idiap.ch/genomics/genomics/toolboxes/apatools.git
cd apatools
pip install .
```

## Scripts
### `bulkapa`
`bulkapa.py` aims at performing 3' UTR analysis of bulk RNA sequencing data. To run it, udpate the file `config/bulkapa_template.yaml` with the right paths. You can then run the script as follows:
```bash
python3 scripts/bulkapa.py --config_file config/bulkapa_template.yaml
```

The config file should contain:
- `df_counts_3UTR_path`: Path to the matrix of 3' UTR isoforms (output of quantification).
- `gtf_file_focus: GTF` file used for quantification of the 3' UTR isoforms -- usually it corresponds to a focus on the last nucleotides of the original 3' UTR isoform annotation.
- `gtf_file_no_focus`: Original GTF file describing the 3' UTR isoforms (no focus on the last nucleotides).
- `df_filtered_genes_path`: Path to the reliably expressed genes dataframe. Should contain gene names in the first column. This matrix is used to filter the reliably expressed 3' UTR isoforms.
- `extension_figures`: Extension you want to save your figures with, e.g "pdf", "png",...
- `treatments`: dictionnary with name of the treatments in keys and list of corresponding sample names in keys.
- `path_to_results`: Directory where to save the results.
