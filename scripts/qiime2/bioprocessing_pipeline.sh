#!/usr/bin/bash

# Setup

## Activate QIIME2 environment
conda activate qiime2-2023.2

## Set paths
export POOL=CB3 #1, 2, 3
export OUTPUT_DIR="outputs/qiime2"
export DATASET="data/raw/${POOL}"
export BC="data/metadata/${POOL}"
export METADATA="data/metadata"
export REFDB="data/db"

## Create output folders from rootcd ..
mkdir -p "${OUTPUT_DIR}/1-demux/CB1"
mkdir -p "${OUTPUT_DIR}/1-demux/CB2"
mkdir -p "${OUTPUT_DIR}/1-demux/CB3"
mkdir -p "${OUTPUT_DIR}/2-combine/CB1"
mkdir -p "${OUTPUT_DIR}/2-combine/CB2"
mkdir -p "${OUTPUT_DIR}/2-combine/CB3"
mkdir -p "${OUTPUT_DIR}/3-import/CB1"
mkdir -p "${OUTPUT_DIR}/3-import/CB2"
mkdir -p "${OUTPUT_DIR}/3-import/CB3"
mkdir -p "${OUTPUT_DIR}/4-dada2"
mkdir -p "${OUTPUT_DIR}/5-taxonomy"
mkdir -p "${OUTPUT_DIR}/6-filter1"
mkdir -p "${OUTPUT_DIR}/7-phylogeny"
mkdir -p "${OUTPUT_DIR}/8-filter2"
mkdir -p "${OUTPUT_DIR}/9-export"

###########################################################################################################
# 1-Demultiplexing

## Round 1
cutadapt \
  --cores 10 \
  --pair-adapters \
  -g "file:${BC}/bcpr-fw.fasta" \
  -G "file:${BC}/bcpr-rv.fasta" \
  -o "${OUTPUT_DIR}/1-demux/${POOL}/round1-{name}-1.fastq.gz" \
  -p "${OUTPUT_DIR}/1-demux/${POOL}/round1-{name}-2.fastq.gz" \
  "${DATASET}/forward.fastq.gz" "${DATASET}/reverse.fastq.gz" \
  --no-indels \
  --errors 0.05 #&> cutadapt-1.log

## Round 2
cutadapt \
  --cores 10 \
  --pair-adapters \
  -g "file:${BC}/bcpr-fw.fasta" \
  -G "file:${BC}/bcpr-rv.fasta" \
  -o "${OUTPUT_DIR}/1-demux/${POOL}/round2-{name}-1.fastq.gz" \
  -p "${OUTPUT_DIR}/1-demux/${POOL}/round2-{name}-2.fastq.gz" \
  "${OUTPUT_DIR}/1-demux/${POOL}/round1-unknown-2.fastq.gz" "${OUTPUT_DIR}/1-demux/${POOL}/round1-unknown-1.fastq.gz" \
  --no-indels \
  --errors 0.05 #&> cutadapt-2.log

###########################################################################################################
# 2-Combine

chmod +x "${BC}/cat.sh"

bash "${BC}/cat.sh"

###########################################################################################################
# 3-Import

## Check manifest
touch "${BC}/manifest.csv"

qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-path "${BC}/manifest.csv" \
  --input-format PairedEndFastqManifestPhred33 \
  --output-path "${OUTPUT_DIR}/3-import/${POOL}/demux.qza"

qiime demux summarize \
  --i-data "${OUTPUT_DIR}/3-import/${POOL}/demux.qza" \
  --o-visualization "${OUTPUT_DIR}/3-import/${POOL}/demux.qzv"

qiime tools view "${OUTPUT_DIR}/3-import/${POOL}/demux.qzv"

###########################################################################################################
# 4-Dada2

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${OUTPUT_DIR}/3-import/${POOL}/demux.qza" \
  --p-trim-left-f 19 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 139 \
  --p-trunc-len-r 212 \
  --p-n-threads 10 \
  --o-table "${OUTPUT_DIR}/4-dada2/${POOL}/table.qza" \
  --o-representative-sequences "${OUTPUT_DIR}/4-dada2/${POOL}/rep-seqs.qza" \
  --o-denoising-stats "${OUTPUT_DIR}/4-dada2/${POOL}/denoising-stats.qza" \
  --verbose

## Look for substancial drops. Percentage of non-chimeric input should be > 85%
qiime metadata tabulate \
  --m-input-file "${OUTPUT_DIR}/4-dada2/${POOL}/denoising-stats.qza" \
  --o-visualization "${OUTPUT_DIR}/4-dada2/${POOL}/denoising-stats.qzv"

qiime tools view "${OUTPUT_DIR}/4-dada2/${POOL}/denoising-stats.qzv"

## Merge the different outputted feature tables and representative sequences
qiime feature-table merge \
  --i-tables "${OUTPUT_DIR}/4-dada2/CB1/table.qza" \
  --i-tables "${OUTPUT_DIR}/4-dada2/CB2/table.qza" \
  --i-tables "${OUTPUT_DIR}/4-dada2/CB3/table.qza" \
  --o-merged-table "${OUTPUT_DIR}/4-dada2/merged_table.qza"

qiime feature-table merge-seqs \
 --i-data "${OUTPUT_DIR}/4-dada2/CB1/rep-seqs.qza" \
 --i-data "${OUTPUT_DIR}/4-dada2/CB2/rep-seqs.qza" \
 --i-data "${OUTPUT_DIR}/4-dada2/CB3/rep-seqs.qza" \
 --o-merged-data "${OUTPUT_DIR}/4-dada2/merged_rep-seqs.qza"

## Verify number of samples 
qiime feature-table summarize \
  --i-table "${OUTPUT_DIR}/4-dada2/merged_table.qza" \
  --o-visualization "${OUTPUT_DIR}/4-dada2/merged_table.qzv" \
  --m-sample-metadata-file "${METADATA}/metadata.tsv" # .tsv required

qiime tools view "${OUTPUT_DIR}/4-dada2/merged_table.qzv"

## Check sequence length and number of ASVs. Sequences should start with “TAC….” and finish with "..GG".
qiime feature-table tabulate-seqs \
  --i-data "${OUTPUT_DIR}/4-dada2/merged_rep-seqs.qza" \
  --o-visualization "${OUTPUT_DIR}/4-dada2/merged_rep-seqs.qzv"

qiime tools view "${OUTPUT_DIR}/4-dada2/merged_rep-seqs.qzv"

###########################################################################################################
# 5-Taxonomy

## 16S with vsearch
qiime feature-classifier classify-consensus-vsearch \
  --i-query "${OUTPUT_DIR}/4-dada2/merged_rep-seqs.qza" \
  --i-reference-reads "${REFDB}/silva-138-99-seqs.qza" \
  --i-reference-taxonomy "${REFDB}/silva-138-99-tax.qza" \
  --p-threads 10 \
  --o-classification "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-vsearch.qza" \
  --o-search-results "${OUTPUT_DIR}/5-taxonomy/16S-tophits-vsearch.qza"

## Import updated SILVA DB with Peribacillus
  qiime tools import \
  --input-path "${REFDB}/dna-sequences.fasta" \
  --output-path "${REFDB}/dna-sequences.qza" \
  --type 'FeatureData[Sequence]'

    qiime tools import \
  --input-path "${REFDB}/taxonomy.tsv" \
  --output-path "${REFDB}/taxonomy.qza" \
  --type 'FeatureData[Taxonomy]'

## 16S with blast
qiime feature-classifier classify-consensus-blast \
  --i-query "${OUTPUT_DIR}/4-dada2/merged_rep-seqs.qza" \
  --i-reference-reads "${REFDB}/dna-sequences.qza" \
  --i-reference-taxonomy "${REFDB}/taxonomy.qza" \
  --o-classification "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-blast.qza" \
  --o-search-results "${OUTPUT_DIR}/5-taxonomy/16S-tophits-blast.qza"


qiime metadata tabulate \
    --m-input-file "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-vsearch.qza" \
    --o-visualization "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-vsearch.qzv"

qiime tools view "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-vsearch.qzv"

qiime metadata tabulate \
    --m-input-file "${OUTPUT_DIR}/5-taxonomy/16S-tophits-vsearch.qza" \
    --o-visualization "${OUTPUT_DIR}/5-taxonomy/16S-tophits-vsearch.qzv"

qiime tools view "${OUTPUT_DIR}/5-taxonomy/16S-tophits-vsearch.qzv"

###########################################################################################################
# 6-Filter

## Keep only bacteria. 
qiime taxa filter-table \
  --i-table "${OUTPUT_DIR}/4-dada2/merged_table.qza" \
  --i-taxonomy "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-vsearch.qza" \
  --p-include bacteria \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table "${OUTPUT_DIR}/6-filter/bac_table.qza"

qiime feature-table summarize \
  --i-table "${OUTPUT_DIR}/6-filter/bac_table.qza" \
  --o-visualization "${OUTPUT_DIR}/6-filter/bac_table.qzv" \
  --m-sample-metadata-file "${METADATA}/metadata.tsv" # .tsv required

qiime tools view "${OUTPUT_DIR}/6-filter/bac_table.qzv" 

## Filter samples with low-read number (<1000 reads) 
qiime feature-table filter-samples \
  --i-table "${OUTPUT_DIR}/6-filter/bac_table.qza" \
  --p-min-frequency 1000 \
  --o-filtered-table "${OUTPUT_DIR}/6-filter/bac_1000_table.qza"

qiime feature-table summarize \
  --i-table "${OUTPUT_DIR}/6-filter/bac_1000_table.qza" \
  --o-visualization "${OUTPUT_DIR}/6-filter/bac_1000_table.qzv" \
  --m-sample-metadata-file "${METADATA}/metadata.tsv" # .tsv required

qiime tools view "${OUTPUT_DIR}/6-filter/bac_1000_table.qzv" # pick --p-max-depth = 5607

## Filter rep-seqs based on the table
qiime feature-table filter-seqs \
  --i-data "${OUTPUT_DIR}/4-dada2/merged_rep-seqs.qza" \
  --i-table "${OUTPUT_DIR}/6-filter/bac_1000_table.qza" \
  --o-filtered-data "${OUTPUT_DIR}/6-filter/bac_1000_rep-seqs.qza"

qiime feature-table tabulate-seqs \
  --i-data "${OUTPUT_DIR}/6-filter/bac_1000_rep-seqs.qza" \
  --o-visualization "${OUTPUT_DIR}/6-filter/bac_1000_rep-seqs.qzv"

qiime tools view "${OUTPUT_DIR}/6-filter/bac_1000_rep-seqs.qzv"

###########################################################################################################
# 7-Phylogeny

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${OUTPUT_DIR}/6-filter/bac_1000_rep-seqs.qza" \
  --output-dir "${OUTPUT_DIR}/7-phylogeny"

###########################################################################################################
# 8-Export

qiime tools export \
  --input-path "${OUTPUT_DIR}/6-filter/bac_1000_table.qza" \
  --output-path "${OUTPUT_DIR}/8-export/"

qiime tools export \
  --input-path "${OUTPUT_DIR}/6-filter/bac_1000_rep-seqs.qza" \
  --output-path "${OUTPUT_DIR}/8-export/"

qiime tools export \
  --input-path "${OUTPUT_DIR}/7-phylogeny/tree.qza" \
  --output-path "${OUTPUT_DIR}/8-export/tree"

qiime tools export \
  --input-path "${OUTPUT_DIR}/7-phylogeny/rooted_tree.qza" \
  --output-path "${OUTPUT_DIR}/8-export/rooted-tree"

## Convert
biom convert \
   -i "${OUTPUT_DIR}/8-export/feature-table.biom" \
   -o "${OUTPUT_DIR}/8-export/feature-table.txt" \
   --to-tsv
