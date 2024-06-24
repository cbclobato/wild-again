#!/bin/bash

# Import raw sequencing data into QIIME2

# Activate QIIME2 environment
conda activate qiime2-2022.2

# Set paths
raw_data_dir="data/raw"
output_dir="data/qiime2"

# Import data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $raw_data_dir \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path $output_dir/raw_data.qza
