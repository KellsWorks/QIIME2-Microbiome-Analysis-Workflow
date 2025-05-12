#!/bin/bash

# QIIME2 Microbiome Analysis Workflow
# This script processes fastq.gz files and performs a complete QIIME2 analysis

# Create necessary directories
mkdir -p qiime2_output
mkdir -p qiime2_output/visualizations
mkdir -p qiime2_output/temp

# Exit on error
set -e

# Check if QIIME2 is activated
if ! command -v qiime &> /dev/null; then
    echo "QIIME2 is not activated. Please activate your QIIME2 environment first."
    echo "Example: conda activate qiime2-2023.5"
    exit 1
fi

echo "Step 1: Preparing manifest file for importing data..."
# Create a manifest file for importing the data
echo "sample-id,absolute-filepath,direction" > qiime2_output/temp/manifest.csv

# Find all forward reads and create manifest entries
for FWD in data/*_1.fastq.gz; do
    SAMPLE=$(basename $FWD _1.fastq.gz)
    REV=${FWD/_1.fastq.gz/_2.fastq.gz}

    if [ -f "$REV" ]; then
        echo "$SAMPLE,$(readlink -f $FWD),forward" >> qiime2_output/temp/manifest.csv
        echo "$SAMPLE,$(readlink -f $REV),reverse" >> qiime2_output/temp/manifest.csv
    fi
done

echo "Manifest file created at qiime2_output/temp/manifest.csv"

echo "Step 2: Importing data into QIIME2..."
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path qiime2_output/temp/manifest.csv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path qiime2_output/demux-paired-end.qza

qiime demux summarize \
  --i-data qiime2_output/demux-paired-end.qza \
  --o-visualization qiime2_output/visualizations/demux-paired-end.qzv

echo "Imported data visualization created. Check demux-paired-end.qzv to determine trimming parameters."

# Check if the import was successful
if [ ! -f "qiime2_output/demux-paired-end.qza" ]; then
    echo "Error: Import failed. The demux-paired-end.qza file was not created."
    exit 1
fi

echo "Step 3: Running DADA2 for quality control and denoising..."
# Check if DADA2 plugin is available
if qiime --help | grep -q "dada2"; then
    qiime dada2 denoise-paired \
      --i-demultiplexed-seqs qiime2_output/demux-paired-end.qza \
      --p-trim-left-f 0 \
      --p-trim-left-r 0 \
      --p-trunc-len-f 250 \
      --p-trunc-len-r 250 \
      --o-table qiime2_output/table.qza \
      --o-representative-sequences qiime2_output/rep-seqs.qza \
      --o-denoising-stats qiime2_output/denoising-stats.qza \
      --p-n-threads 0
else
    echo "Error: DADA2 plugin is not available. Please install it with:"
    echo "conda install -c bioconda -c conda-forge -c qiime2 -c defaults qiime2-dada2-plugin"
    exit 1
fi

# Check if DADA2 was successful
if [ ! -f "qiime2_output/denoising-stats.qza" ] || [ ! -f "qiime2_output/table.qza" ] || [ ! -f "qiime2_output/rep-seqs.qza" ]; then
    echo "Error: DADA2 processing failed. Required output files are missing."
    exit 1
fi

qiime metadata tabulate \
  --m-input-file qiime2_output/denoising-stats.qza \
  --o-visualization qiime2_output/visualizations/denoising-stats.qzv

qiime feature-table summarize \
  --i-table qiime2_output/table.qza \
  --o-visualization qiime2_output/visualizations/table.qzv \
  --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data qiime2_output/rep-seqs.qza \
  --o-visualization qiime2_output/visualizations/rep-seqs.qzv

echo "DADA2 processing complete. Feature table and representative sequences generated."

echo "Step 4: Building phylogenetic tree..."
# Check if phylogeny plugin is available
if qiime --help | grep -q "phylogeny"; then
    qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences qiime2_output/rep-seqs.qza \
      --o-alignment qiime2_output/aligned-rep-seqs.qza \
      --o-masked-alignment qiime2_output/masked-aligned-rep-seqs.qza \
      --o-tree qiime2_output/unrooted-tree.qza \
      --o-rooted-tree qiime2_output/rooted-tree.qza

    # Check if phylogenetic tree construction was successful
    if [ ! -f "qiime2_output/rooted-tree.qza" ]; then
        echo "Error: Phylogenetic tree construction failed."
        exit 1
    fi

    echo "Phylogenetic tree construction complete."
else
    echo "Error: Phylogeny plugin is not available. Please install it with:"
    echo "conda install -c bioconda -c conda-forge -c qiime2 -c defaults qiime2-phylogeny-plugin"
    exit 1
fi

echo "Step 5: Assigning taxonomy..."
# Check if feature-classifier plugin is available
if qiime --help | grep -q "feature-classifier"; then
    # Download the classifier if it doesn't exist
    if [ ! -f "classifier.qza" ]; then
        echo "Downloading Silva 138 99% classifier..."
        wget -O "classifier.qza" "https://data.qiime2.org/2023.5/common/silva-138-99-nb-classifier.qza"
    fi

    # Check if classifier download was successful
    if [ ! -f "classifier.qza" ]; then
        echo "Error: Failed to download the classifier. Please download it manually and place it in the current directory."
        echo "You can download it from: https://data.qiime2.org/2023.5/common/silva-138-99-nb-classifier.qza"
        exit 1
    fi

    qiime feature-classifier classify-sklearn \
      --i-classifier classifier.qza \
      --i-reads qiime2_output/rep-seqs.qza \
      --o-classification qiime2_output/taxonomy.qza

    # Check if taxonomy assignment was successful
    if [ ! -f "qiime2_output/taxonomy.qza" ]; then
        echo "Error: Taxonomy assignment failed."
        exit 1
    fi

    qiime metadata tabulate \
      --m-input-file qiime2_output/taxonomy.qza \
      --o-visualization qiime2_output/visualizations/taxonomy.qzv

    qiime taxa barplot \
      --i-table qiime2_output/table.qza \
      --i-taxonomy qiime2_output/taxonomy.qza \
      --m-metadata-file metadata.tsv \
      --o-visualization qiime2_output/visualizations/taxa-bar-plots.qzv

    echo "Taxonomic classification complete."
else
    echo "Error: Feature-classifier plugin is not available. Please install it with:"
    echo "conda install -c bioconda -c conda-forge -c qiime2 -c defaults qiime2-feature-classifier"
    exit 1
fi

echo "Step 6: Performing diversity analyses..."

# Check if diversity plugin is available
if qiime --help | grep -q "diversity"; then
    # Determine sampling depth for rarefaction
    # We'll use a more dynamic approach to determine sampling depth
    echo "Determining appropriate sampling depth..."

    # Default sampling depth if we can't determine it dynamically
    SAMPLING_DEPTH=10000

    # Try to get the minimum frequency from the feature table
    if [ -f "qiime2_output/visualizations/table.qzv" ]; then
        echo "Please check qiime2_output/visualizations/table.qzv to determine an appropriate sampling depth."
        echo "The sampling depth should be less than or equal to the minimum frequency per sample."
        echo "Using default sampling depth of $SAMPLING_DEPTH for now."
        echo "You can modify this value in the script if needed."
    fi

    echo "Using sampling depth: $SAMPLING_DEPTH"

    qiime diversity core-metrics-phylogenetic \
      --i-phylogeny qiime2_output/rooted-tree.qza \
      --i-table qiime2_output/table.qza \
      --p-sampling-depth $SAMPLING_DEPTH \
      --m-metadata-file metadata.tsv \
      --output-dir qiime2_output/diversity

    # Check if diversity analysis was successful
    if [ ! -d "qiime2_output/diversity" ] || [ ! -f "qiime2_output/diversity/faith_pd_vector.qza" ]; then
        echo "Error: Diversity analysis failed."
        exit 1
    fi

    echo "Core diversity metrics calculated."

    echo "Step 7: Generating rarefaction curves..."
    qiime diversity alpha-rarefaction \
      --i-table qiime2_output/table.qza \
      --i-phylogeny qiime2_output/rooted-tree.qza \
      --p-max-depth $SAMPLING_DEPTH \
      --m-metadata-file metadata.tsv \
      --o-visualization qiime2_output/visualizations/alpha-rarefaction.qzv

    # Check if rarefaction curve generation was successful
    if [ ! -f "qiime2_output/visualizations/alpha-rarefaction.qzv" ]; then
        echo "Error: Rarefaction curve generation failed."
        exit 1
    fi

    echo "Rarefaction curves generated."
else
    echo "Error: Diversity plugin is not available. Please install it with:"
    echo "conda install -c bioconda -c conda-forge -c qiime2 -c defaults qiime2-diversity-plugin"
    exit 1
fi

echo "Step 8: Performing alpha diversity statistical tests..."

# Check if required files exist
if [ ! -f "qiime2_output/diversity/faith_pd_vector.qza" ] || [ ! -f "qiime2_output/diversity/shannon_vector.qza" ]; then
    echo "Error: Required diversity files are missing. Cannot perform alpha diversity tests."
    exit 1
fi

# Extract metadata columns from metadata.tsv
echo "Checking metadata columns..."
METADATA_COLUMNS=$(head -n 1 metadata.tsv | tr '\t' '\n' | grep -v "#" | grep -v "SampleID" | grep -v "description")
echo "Found metadata columns: $METADATA_COLUMNS"

# Check if trial_point column exists
if grep -q "trial_point" metadata.tsv; then
    echo "Performing alpha diversity tests for trial_point..."
    qiime diversity alpha-group-significance \
      --i-alpha-diversity qiime2_output/diversity/faith_pd_vector.qza \
      --m-metadata-file metadata.tsv \
      --o-visualization qiime2_output/visualizations/faith-pd-group-significance-trial-point.qzv \
      --m-metadata-column trial_point

    qiime diversity alpha-group-significance \
      --i-alpha-diversity qiime2_output/diversity/shannon_vector.qza \
      --m-metadata-file metadata.tsv \
      --o-visualization qiime2_output/visualizations/shannon-group-significance-trial-point.qzv \
      --m-metadata-column trial_point
else
    echo "Warning: 'trial_point' column not found in metadata. Skipping related tests."
fi

# Check if sex column exists
if grep -q "sex" metadata.tsv; then
    echo "Performing alpha diversity tests for sex..."
    qiime diversity alpha-group-significance \
      --i-alpha-diversity qiime2_output/diversity/faith_pd_vector.qza \
      --m-metadata-file metadata.tsv \
      --o-visualization qiime2_output/visualizations/faith-pd-group-significance-sex.qzv \
      --m-metadata-column sex

    qiime diversity alpha-group-significance \
      --i-alpha-diversity qiime2_output/diversity/shannon_vector.qza \
      --m-metadata-file metadata.tsv \
      --o-visualization qiime2_output/visualizations/shannon-group-significance-sex.qzv \
      --m-metadata-column sex
else
    echo "Warning: 'sex' column not found in metadata. Skipping related tests."
fi

echo "Alpha diversity statistical tests completed."

echo "Step 9: Performing beta diversity statistical tests..."

# Check if required files exist
if [ ! -f "qiime2_output/diversity/unweighted_unifrac_distance_matrix.qza" ] || [ ! -f "qiime2_output/diversity/weighted_unifrac_distance_matrix.qza" ]; then
    echo "Error: Required distance matrix files are missing. Cannot perform beta diversity tests."
    exit 1
fi

# Check if trial_point column exists
if grep -q "trial_point" metadata.tsv; then
    echo "Performing beta diversity tests for trial_point..."
    qiime diversity beta-group-significance \
      --i-distance-matrix qiime2_output/diversity/unweighted_unifrac_distance_matrix.qza \
      --m-metadata-file metadata.tsv \
      --m-metadata-column trial_point \
      --o-visualization qiime2_output/visualizations/unweighted-unifrac-trial-point-significance.qzv \
      --p-pairwise

    qiime diversity beta-group-significance \
      --i-distance-matrix qiime2_output/diversity/weighted_unifrac_distance_matrix.qza \
      --m-metadata-file metadata.tsv \
      --m-metadata-column trial_point \
      --o-visualization qiime2_output/visualizations/weighted-unifrac-trial-point-significance.qzv \
      --p-pairwise
else
    echo "Warning: 'trial_point' column not found in metadata. Skipping related tests."
fi

# Check if sex column exists
if grep -q "sex" metadata.tsv; then
    echo "Performing beta diversity tests for sex..."
    qiime diversity beta-group-significance \
      --i-distance-matrix qiime2_output/diversity/unweighted_unifrac_distance_matrix.qza \
      --m-metadata-file metadata.tsv \
      --m-metadata-column sex \
      --o-visualization qiime2_output/visualizations/unweighted-unifrac-sex-significance.qzv \
      --p-pairwise

    qiime diversity beta-group-significance \
      --i-distance-matrix qiime2_output/diversity/weighted_unifrac_distance_matrix.qza \
      --m-metadata-file metadata.tsv \
      --m-metadata-column sex \
      --o-visualization qiime2_output/visualizations/weighted-unifrac-sex-significance.qzv \
      --p-pairwise
else
    echo "Warning: 'sex' column not found in metadata. Skipping related tests."
fi

echo "Beta diversity statistical tests completed."

echo "QIIME2 analysis workflow completed successfully!"
echo "All results are stored in the qiime2_output directory."
echo "Visualizations can be viewed using 'qiime tools view' command or at https://view.qiime2.org/"
echo ""
echo "Key steps and values used in this analysis:"
echo "1. Data import: Used manifest file format for paired-end sequences"
echo "2. Quality control: DADA2 with trim-left-f=0, trim-left-r=0, trunc-len-f=250, trunc-len-r=250"
echo "3. Taxonomy: Silva 138 99% classifier"
echo "4. Diversity analysis: Sampling depth = $SAMPLING_DEPTH"
echo ""
echo "To interpret the results:"
echo "1. Rarefaction curves (alpha-rarefaction.qzv): Check if curves plateau to ensure sufficient sequencing depth"
echo "2. Alpha diversity: Compare diversity metrics across groups (trial_point, sex) using statistical tests"
echo "3. Beta diversity: Examine community composition differences between groups using PERMANOVA tests"
