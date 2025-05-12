# QIIME2 Microbiome Analysis Workflow

This repository contains a complete QIIME2 workflow for analyzing microbiome data from paired-end sequencing.

## Prerequisites

- QIIME2 installed (recommended version: 2023.5 or later)
- Paired-end fastq.gz files in the `data/` directory
- Properly formatted metadata file (`metadata.tsv`)

## Setup Instructions

1. **Activate QIIME2 environment**:
   ```bash
   conda activate qiime2-2023.5  # or your QIIME2 environment name
   ```

2. **Prepare your data**:
   - Place all your fastq.gz files in the `data/` directory
   - Files should be named with the pattern: `SAMPLEID_1.fastq.gz` and `SAMPLEID_2.fastq.gz`
   - Update the `metadata.tsv` file with your sample information

3. **Run the workflow**:
   ```bash
   bash qiime2_workflow.sh
   ```

## Workflow Steps

1. **Data Import**: Imports paired-end fastq.gz files using a manifest file
2. **Quality Control**: Uses DADA2 for denoising, chimera removal, and ASV identification
3. **Phylogenetic Analysis**: Builds a phylogenetic tree for diversity analyses
4. **Taxonomic Classification**: Assigns taxonomy using the Silva 138 classifier
5. **Diversity Analysis**: Calculates alpha and beta diversity metrics
6. **Statistical Tests**: Performs statistical tests on diversity metrics based on metadata categories

## Output Files

All output files are stored in the `qiime2_output/` directory:

- `qiime2_output/table.qza`: Feature table (ASV counts)
- `qiime2_output/rep-seqs.qza`: Representative sequences
- `qiime2_output/taxonomy.qza`: Taxonomic classifications
- `qiime2_output/rooted-tree.qza`: Phylogenetic tree
- `qiime2_output/diversity/`: Directory containing diversity metrics

Visualizations are stored in `qiime2_output/visualizations/`:

- `demux-paired-end.qzv`: Quality plots for raw sequences
- `table.qzv`: Feature table summary
- `rep-seqs.qzv`: Representative sequences
- `taxonomy.qzv`: Taxonomy summary
- `taxa-bar-plots.qzv`: Taxonomic bar plots
- `alpha-rarefaction.qzv`: Rarefaction curves
- Various alpha and beta diversity visualizations

## Interpreting Results

1. **Rarefaction Curves**: Check if curves plateau to ensure sufficient sequencing depth
2. **Alpha Diversity**: Compare diversity metrics across groups using statistical tests
3. **Beta Diversity**: Examine community composition differences between groups using PERMANOVA tests

## Troubleshooting

- If the script fails at any step, check the error message for details
- Ensure QIIME2 is properly installed and activated
- Verify that your metadata file is properly formatted
- Check that your fastq.gz files follow the expected naming pattern

## Customization

You can modify the following parameters in the script:

- DADA2 trimming parameters (trim-left-f, trim-left-r, trunc-len-f, trunc-len-r)
- Sampling depth for diversity analyses
- Taxonomic classifier

