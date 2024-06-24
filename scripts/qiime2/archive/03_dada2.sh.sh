#!/bin/bash

# Set variables
INPUT_TABLE=data/qiime2/table-dada2.qza
INPUT_REFS=data/qiime2/silva-138-99-seqs.qza
OUTPUT_CLASS=data/qiime2/taxonomy.qza

# Load QIIME2 module
module load qiime2

# Perform taxonomic assignment
qiime feature-classifier classify-sklearn \
  --i-classifier $INPUT_REFS \
  --i-reads $INPUT_TABLE \
  --o-classification $OUTPUT_CLASS

# Export taxonomic assignments as TSV file
qiime tools export \
  --input-path $OUTPUT_CLASS \
  --output-path data/taxa

# Convert TSV file to BIOM format
biom convert \
  -i data/taxa/taxonomy.tsv \
  -o data/taxa/taxonomy.biom \
  --table-type "Taxonomy"