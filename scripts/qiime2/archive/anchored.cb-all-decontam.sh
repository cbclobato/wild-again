#!/usr/bin/bash

# Setup

## Activate QIIME2 environment
conda info --envs
conda activate qiime2-2022.2
pwd
cd /home/jmachadofreitas/Documents/Carolina/22-amp-seq-seedlibrary/metabarcoding/

## Set variables
export JOBNAME="anchored"
export POOL=CB-all  #CB1 #CB2 #CB3 #CB-all 
export DATASET="datasets/${POOL}"
export REFERENCES="datasets/"
export METADATA="metadata/${POOL}"
export OUTPUT_DIR="results/${JOBNAME}/${POOL}"

## Set paths
mkdir -p "${OUTPUT_DIR}/1-demux"
mkdir -p "${OUTPUT_DIR}/2-combined"
mkdir -p "${OUTPUT_DIR}/3-import"
mkdir -p "${OUTPUT_DIR}/4-dada2"
mkdir -p "${OUTPUT_DIR}/5-taxonomy"
mkdir -p "${OUTPUT_DIR}/6-filter1"
mkdir -p "${OUTPUT_DIR}/7-decontam"
mkdir -p "${OUTPUT_DIR}/8-phylogeny"
mkdir -p "${OUTPUT_DIR}/8-phylogeny-1000"
mkdir -p "${OUTPUT_DIR}/9-rarefaction"
mkdir -p "${OUTPUT_DIR}/10-filter2"
mkdir -p "${OUTPUT_DIR}/12-alpha-group"
mkdir -p "${OUTPUT_DIR}/12-beta-group"
mkdir -p "${OUTPUT_DIR}/rstudio"

chmod +x

###########################################################################################################

# 1-Demultiplexing

## Round 1
cutadapt \
  --cores 8 \
  --pair-adapters \
  -g "file:${METADATA}/^bcpr-fw.fasta" \
  -G "file:${METADATA}/^bcpr-rv.fasta" \
  -o "${OUTPUT_DIR}/1-demux/round1-{name}-1.fastq.gz" \
  -p "${OUTPUT_DIR}/1-demux/round1-{name}-2.fastq.gz" \
  "${DATASET}/forward.fastq.gz" "${DATASET}/reverse.fastq.gz" \
  --no-indels \
  --errors 0.05 &> cutadapt-1.log

## Round 2
cutadapt \
  --cores 8 \
  --pair-adapters \
  -g "file:${METADATA}/^bcpr-fw.fasta" \
  -G "file:${METADATA}/^bcpr-rv.fasta" \
  -o "${OUTPUT_DIR}/1-demux/round2-{name}-1.fastq.gz" \
  -p "${OUTPUT_DIR}/1-demux/round2-{name}-2.fastq.gz" \
  "${OUTPUT_DIR}/1-demux/round1-unknown-2.fastq.gz" "${OUTPUT_DIR}/1-demux/round1-unknown-1.fastq.gz" \
  --no-indels \
  --errors 0.05 &> cutadapt-2.log

###########################################################################################################

# 2-Combine
chmod +x "${METADATA}/cat.sh" \
bash "${METADATA}/cat.sh"

## Check if manifest exists
touch "manifests/${JOBNAME}.csv"

###########################################################################################################

# 3-Import

qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-path "manifests/${POOL}/${JOBNAME}.csv" \
  --input-format PairedEndFastqManifestPhred33 \
  --output-path "${OUTPUT_DIR}/3-import/demux.qza"

## Quality Check
qiime demux summarize \
  --i-data "${OUTPUT_DIR}/3-import/demux.qza" \
  --o-visualization "${OUTPUT_DIR}/3-import/demux.qzv"
qiime tools view "${OUTPUT_DIR}/3-import/demux.qzv"

###########################################################################################################

# 4-Filter and Denoise
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${OUTPUT_DIR}/3-import/demux.qza" \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 142 \
  --p-trunc-len-r 213 \
  --p-n-threads 6 \
  --o-table "${OUTPUT_DIR}/4-dada2/table.qza" \
  --o-representative-sequences "${OUTPUT_DIR}/4-dada2/rep-seqs.qza" \
  --o-denoising-stats "${OUTPUT_DIR}/4-dada2/denoising-stats.qza" \
  --verbose

## Look for substancial drops. Percentage of input non-chimeric > 85%.
qiime metadata tabulate \
  --m-input-file "${OUTPUT_DIR}/4-dada2/denoising-stats.qza" \
  --o-visualization "${OUTPUT_DIR}/4-dada2/denoising-stats.qzv"
qiime tools view "${OUTPUT_DIR}/4-dada2/denoising-stats.qzv"

## Verify number of samples, and sequencing depth.
qiime feature-table summarize \
  --i-table "${OUTPUT_DIR}/4-dada2/table.qza" \
  --o-visualization "${OUTPUT_DIR}/4-dada2/table.qzv" \
  --m-sample-metadata-file "metadata/${POOL}/metadata.tsv" # mapping-file.tsv optional
qiime tools view "${OUTPUT_DIR}/4-dada2/table.qzv"

## Check Sequence length, and should start with “TAC….” and finish with "..GG", and number of ASVs.
qiime feature-table tabulate-seqs \
  --i-data "${OUTPUT_DIR}/4-dada2/rep-seqs.qza" \
  --o-visualization "${OUTPUT_DIR}/4-dada2/rep-seqs.qzv"
qiime tools view "${OUTPUT_DIR}/4-dada2/rep-seqs.qzv"

###########################################################################################################

# 5-Taxonomy

## 16S
qiime feature-classifier classify-consensus-vsearch \
  --i-query "${OUTPUT_DIR}/4-dada2/rep-seqs.qza" \
  --i-reference-reads "${REFERENCES}/silva-138-99-seqs.qza" \
  --i-reference-taxonomy "${REFERENCES}/silva-138-99-tax.qza" \
  --p-threads 7 \
  --o-classification "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-vsearch.qza"

qiime metadata tabulate \
    --m-input-file "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-vsearch.qza" \
    --o-visualization "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-vsearch.qzv"
qiime tools view "${OUTPUT_DIR}/5-taxonomy/taxonomy-vsearch.qzv"

###########################################################################################################

# 6-Filtering 1

### Filter host-related ASVs from Table
qiime taxa filter-table \
  --i-table "${OUTPUT_DIR}/4-dada2/table.qza" \
  --i-taxonomy "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-vsearch.qza" \
  --p-exclude "mitochondria,chloroplast" \
  --o-filtered-table "${OUTPUT_DIR}/6-filter1/bac-table.qza"

qiime feature-table summarize \
  --i-table "${OUTPUT_DIR}/6-filter1/bac-table.qza" \
  --o-visualization "${OUTPUT_DIR}/6-filter1/bac-table.qzv" \
  --m-sample-metadata-file "metadata/${POOL}/metadata.tsv" # mapping-file.tsv optional
qiime tools view "${OUTPUT_DIR}/6-filter1/bac-table.qzv" 

# **TODO**
### Filter samples with low-read number (<1000 reads) from Table
qiime feature-table filter-samples \
  --i-table "${OUTPUT_DIR}/6-filter1/bac-table.qza" \
  --p-min-frequency 1000 \
  --o-filtered-table "${OUTPUT_DIR}/6-filter1/bac-1000-table.qza"

qiime feature-table summarize \
  --i-table "${OUTPUT_DIR}/6-filter1/bac-1000-table.qza" \
  --o-visualization "${OUTPUT_DIR}/6-filter1/bac-1000-table.qzv" \
  --m-sample-metadata-file "metadata/${POOL}/metadata.tsv" # mapping-file.tsv optional
qiime tools view "${OUTPUT_DIR}/6-filter1/bac-1000-table.qzv"  # pick --p-max-depth

### Filter host-related Sequences
qiime feature-table filter-seqs \
  --i-data "${OUTPUT_DIR}/4-dada2/rep-seqs.qza" \
  --i-table "${OUTPUT_DIR}/6-filter1/bac-1000-table.qza" \
  --o-filtered-data "${OUTPUT_DIR}/6-filter1/bac-1000-rep-seqs.qza"

qiime feature-table tabulate-seqs \
  --i-data "${OUTPUT_DIR}/6-filter1/bac-1000-rep-seqs.qza" \
  --o-visualization "${OUTPUT_DIR}/6-filter1/bac-1000-rep-seqs.qzv" 
qiime tools view "${OUTPUT_DIR}/6-filter1/bac-1000-rep-seqs.qzv" # get sequences if blast is needed

###########################################################################################################

# **TODO**
#### Go to R to run decontam ####

# Export or download from .qzv files

 qiime tools export \
  --input-path "${OUTPUT_DIR}/6-filter1/bac-table.qza" \
  --output-path "${OUTPUT_DIR}/rstudio/phyloseq-obj/exp-bac-table.biom"

 qiime tools export \
  --input-path "${OUTPUT_DIR}/6-filter1/bac-rep-seqs.qza" \
  --output-path "${OUTPUT_DIR}/rstudio/phyloseq-obj/exp-bac-rep-seqs.fasta"

# Convert

## .csv to .tsv
cat "${OUTPUT_DIR}/rstudio/taxonomy_decontam_r.csv" | sed 's/,/\t/g' > "${OUTPUT_DIR}/rstudio/taxonomy_decontam_r.tsv" 

## .biom to .tsv (in r sample names should have _ intead of - and first column should be named OTUID)
biom convert \
  -i "${OUTPUT_DIR}/rstudio/phyloseq-obj/exp-bac-table.biom" \
  -o "${OUTPUT_DIR}/rstudio/phyloseq-obj/exp-bac-table.tsv" \
  --to-tsv 

## .tsv to .biom
biom convert \
  -i "${OUTPUT_DIR}/rstudio/asv_decontam_r.tsv" \
  -o "${OUTPUT_DIR}/rstudio/asv_decontam_r.biom" \
  --table-type="OTU table" \
  --to-hdf5 

# Import

 qiime tools import \
  --type FeatureTable[Frequency] \
  --input-format "BIOMV210Format" \
  --input-path "${OUTPUT_DIR}/rstudio/asv_decontam_r.biom" \
  --output-path "${OUTPUT_DIR}/7-decontam/asv-table-decontam.qza"

 qiime tools import \
  --type FeatureData[Taxonomy] \
  --input-format "HeaderlessTSVTaxonomyFormat" \
  --input-path "${OUTPUT_DIR}/rstudio/taxonomy_decontam_r.tsv" \
  --output-path "${OUTPUT_DIR}/7-decontam/taxonomy-decontam.qza"

## Visualizations
qiime metadata tabulate \
  --m-input-file "${OUTPUT_DIR}/7-decontam/taxonomy-decontam.qza" \
  --o-visualization "${OUTPUT_DIR}/7-decontam/taxonomy-decontam.qzv"
qiime tools view "${OUTPUT_DIR}/7-decontam/taxonomy-decontam.qzv"

qiime feature-table summarize \
  --i-table "${OUTPUT_DIR}/7-decontam/asv-table-decontam.qza" \
  --m-sample-metadata-file "metadata/${POOL}/metadata.tsv" \
  --o-visualization "${OUTPUT_DIR}/7-decontam/feature-table-decontam.qzv"
qiime tools view "${OUTPUT_DIR}/7-decontam/feature-table-decontam.qzv"

###########################################################################################################

#### Back to Qiime2 ####

# 8-Phylogeny on non-rarefied and rarefied dataset

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${OUTPUT_DIR}/6-filter1/bac-rep-seqs.qza" \
  --output-dir "${OUTPUT_DIR}/8-phylogeny"

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${OUTPUT_DIR}/6-filter1/bac-1000-rep-seqs.qza" \
  --output-dir "${OUTPUT_DIR}/8-phylogeny-1000"

# **TODO**
## Export trees
qiime tools export \
  --input-path "${OUTPUT_DIR}/8-phylogeny-1000/tree.qza" \
  --output-path "${OUTPUT_DIR}/rstudio/tree"

qiime tools export \
  --input-path "${OUTPUT_DIR}/8-phylogeny-1000/rooted_tree.qza" \
  --output-path "${OUTPUT_DIR}/rstudio/rooted-tree"

# **TODO**
## Prune tree (keep just core taxa)
source /macqiime/configs/bash_profile.txt
filter_tree.py 
   -i Subset1-16S-V4-rooted-tree.nwk 
   -t 16SV4_Core_Microbiome.txt 
   -o Core-taxa-16S-V4-rooted-tree.tre

###########################################################################################################

# 9-Rarefaction

## Check rarefaction curves - change pmax-depth depending on the dataset
qiime diversity alpha-rarefaction \
  --i-table "${OUTPUT_DIR}/6-filter1/bac-1000-table.qza" \
  --i-phylogeny "${OUTPUT_DIR}/8-phylogeny-1000/rooted_tree.qza" \
  --p-min-depth 1000 \
  --p-max-depth 30000 \
  --m-metadata-file "metadata/${POOL}/metadata.tsv" \
  --o-visualization "${OUTPUT_DIR}/9-rarefaction/alpha-rarefaction-1000.qzv"

# **TODO**
qiime diversity alpha-rarefaction \
  --i-table "${OUTPUT_DIR}/6-filter1/bac-1000-table.qza" \
  --p-max-depth 8348 \
  --p-steps 20 \
  --o-visualization "${OUTPUT_DIR}/9-rarefaction/alpha-rarefaction-1000-persample.qzv"

# **TODO**
## Choose the rarefaction level and rarefy the table
qiime feature-table rarefy \
  --i-table "${OUTPUT_DIR}/6-filter1/bac-1000-table.qza" \
  --p-sampling-depth 8696 \
  --o-rarefied-table "${OUTPUT_DIR}/9-rarefaction/bac-1000-table-rare8696.qza"

# **TODO**
## Filter the rep-seqs to keep only the ASVs present in the rarefied ASV table 
 qiime feature-table filter-seqs \
  --i-data "${OUTPUT_DIR}/6-filter1/bac-1000-rep-seqs.qza" \
  --i-table "${OUTPUT_DIR}/9-rarefaction/bac-1000-table-rare8696.qza" \
  --o-filtered-data "${OUTPUT_DIR}/9-rarefaction/bac-1000-rep-seqs-rare8696.qza"

# **TODO**
## Check summaries
qiime feature-table summarize \
  --i-table "${OUTPUT_DIR}/9-rarefaction/bac-1000-table-rare8696.qza" \
  --o-visualization "${OUTPUT_DIR}/9-rarefaction/bac-1000-table-rare8696.qzv" \
  --m-sample-metadata-file "metadata/${POOL}/metadata.tsv"

###########################################################################################################

# Export unrarefied and rarefied datasets

## ASV tables 
qiime tools export \
  --input-path "${OUTPUT_DIR}/6-filter1/bac-1000-table.qza" \
  --output-path "${OUTPUT_DIR}/rstudio/bac-1000-table.biom"

qiime tools export \
  --input-path "${OUTPUT_DIR}/9-rarefaction/bac-1000-table-rare8696.qza" \
  --output-path "${OUTPUT_DIR}/rstudio//bac-1000-table-rare8696.biom"

## Rep-seqs 
qiime tools export \
  --input-path "${OUTPUT_DIR}/6-filter1/bac-1000-rep-seqs.qza" \
  --output-path "${OUTPUT_DIR}/rstudio/bac-1000-rep-seqs"

  qiime tools export \
  --input-path "${OUTPUT_DIR}/9-rarefaction/bac-1000-rep-seqs-rare8696.qza" \
  --output-path "${OUTPUT_DIR}/rstudio/bac-1000-rep-seqs-rare8696"

## Convert format
biom convert \
  -i "${OUTPUT_DIR}/rstudio/bac-1000-table.biom" \
  -o "${OUTPUT_DIR}/rstudio/bac-1000-table.txt" \
  --to-tsv

biom convert \
  -i "${OUTPUT_DIR}/rstudio/bac-1000-table-rare8696.biom" \
  -o "${OUTPUT_DIR}/rstudio/bac-1000-table-rare8696.txt" \
  --to-tsv

## Remove white space from taxonomy file and re-import
D_0__Bacteria;D_1__Bacteroidetes;D_2__Ignavibacteria;D_3__OPB56;D_4__uncultured bacterium 
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Carnobacteriaceae;D_5__Alloiococcus;D_6__bacterium
D_0__Bacteria;D_1__Chloroflexi;D_2__Anaerolineae;D_3__Anaerolineales;D_4__Anaerolineaceae;D_5__uncultured;D_6__uncultured bacterium 
D_0__Bacteria;D_1__Bacteroidetes;D_2__Ignavibacteria;D_3__OPB56;D_4__uncultured bacterium 
D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales;D_4__Carnobacteriaceae;D_5__Alloiococcus;D_6__bacterium 

# Import
 qiime tools import \
  --type FeatureData[Taxonomy] \
  --input-format "HeaderlessTSVTaxonomyFormat" \
  --input-path "${OUTPUT_DIR}/rstudio/bac-1000-table.txt" \
  --output-path "${OUTPUT_DIR}/rstudio/bac-1000-table-clean.qza" \

###########################################################################################################

# 10-Filter clean database

## Remove Unassigned, Archaea and Eukaryotic ASVs in the table -- Filter clean database

### Unrarefied data
qiime taxa filter-table \
  --i-table "${OUTPUT_DIR}/6-filter1/bac-1000-table.qza" \
  --i-taxonomy "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-vsearch.qza" \
  --p-exclude d__Eukaryota,Archaea,Unassigned  \
  --o-filtered-table "${OUTPUT_DIR}/10-filter2/bac-1000-table-filt1.qza"

### Rarefied data
qiime taxa filter-table \
  --i-table "${OUTPUT_DIR}/9-rarefaction/bac-1000-table-rare8696.qza" \
  --i-taxonomy "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-vsearch.qza" \
  --p-exclude d__Eukaryota,Archaea,Unassigned \
  --o-filtered-table "${OUTPUT_DIR}/10-filter2/rarefied-data/bac-1000-table-rare8696-filt1.qza"

qiime feature-table summarize \
  --i-table "${OUTPUT_DIR}/10-filter2/rarefied-data/bac-1000-table-rare8696-filt1.qza" \
  --o-visualization "${OUTPUT_DIR}/10-filter2/rarefied-data/bac-1000-table-rare8696-filt1.qzv" \
  --m-sample-metadata-file "metadata/${POOL}/metadata.tsv"
qiime tools view "${OUTPUT_DIR}/10-filter2/rarefied-data/bac-1000-table-rare8696-filt3.qzv" 


## Filter low rare ASVs

### Unrarefied data
qiime feature-table filter-features \
  --i-table "${OUTPUT_DIR}/10-filter2/bac-1000-table-filt1.qza" \
  --p-min-frequency 20 \
  --p-min-samples 2 \
  --o-filtered-table "${OUTPUT_DIR}/10-filter2/bac-1000-table-filt2.qza"

## Filter rep-seqs based on table

### Unrarefied data
qiime feature-table filter-seqs \
  --i-data "${OUTPUT_DIR}/6-filter1/bac-1000-rep-seqs.qza" \
  --i-table "${OUTPUT_DIR}/10-filter2/bac-1000-table-filt2.qza" \
  --o-filtered-data "${OUTPUT_DIR}/10-filter2/bac-1000-rep-seqs-filt2.qza"

### Rarefied data
qiime feature-table filter-seqs \
  --i-data "${OUTPUT_DIR}/9-rarefaction/bac-1000-rep-seqs-rare8696.qza" \
  --i-table "${OUTPUT_DIR}/10-filter2/rarefied-data/bac-1000-table-rare8696-filt1.qza" \
  --o-filtered-data "${OUTPUT_DIR}/10-filter2/rarefied-data/bac-1000-rep-seqs-rare8696-filt1.qza"

## Remove ASVs that are too short for 16S (< 200 bp)

### Unrarefied data
qiime feature-table filter-seqs \
  --i-data "${OUTPUT_DIR}/10-filter2/bac-1000-rep-seqs-filt2.qza" \
  --m-metadata-file "${OUTPUT_DIR}/10-filter2/bac-1000-rep-seqs-filt2.qza" \
  --p-where "length(Sequence) > 200" \
  --o-filtered-data "${OUTPUT_DIR}/10-filter2/bac-1000-rep-seqs-filt3.qza"

qiime feature-table tabulate-seqs \
  --i-data "${OUTPUT_DIR}/10-filter2/bac-1000-rep-seqs-filt3.qza" \
  --o-visualization "${OUTPUT_DIR}/10-filter2/bac-1000-rep-seqs-filt3.qzv"
qiime tools view "${OUTPUT_DIR}/10-filter2/bac-1000-rep-seqs-filt3.qzv"

## Filter table based on rep seq file

### Unrarefied data
qiime feature-table filter-features \
  --i-table "${OUTPUT_DIR}/10-filter2/bac-1000-table-filt2.qza" \
  --m-metadata-file "${OUTPUT_DIR}/10-filter2/bac-1000-rep-seqs-filt3.qza" \
  --o-filtered-table "${OUTPUT_DIR}/10-filter2/bac-1000-table-filt3.qza"

qiime feature-table summarize \
  --i-table "${OUTPUT_DIR}/10-filter2/bac-1000-table-filt3.qza" \
  --o-visualization "${OUTPUT_DIR}/10-filter2/bac-1000-table-filt3.qzv" \
  --m-sample-metadata-file "metadata/${POOL}/metadata.tsv"
qiime tools view "${OUTPUT_DIR}/10-filter2/bac-1000-table-filt3.qzv" 

## Filter phylogeny 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${OUTPUT_DIR}/10-filter2/bac-1000-rep-seqs-filt3.qza" \
  --output-dir "${OUTPUT_DIR}/8-phylogeny-1000-p"

qiime tools export \
  --input-path "${OUTPUT_DIR}/8-phylogeny-1000-p/tree.qza" \
  --output-path "${OUTPUT_DIR}/rstudio/phyloseq-obj-1000-p/tree"

qiime tools export \
  --input-path "${OUTPUT_DIR}/8-phylogeny-1000-p/rooted_tree.qza" \
  --output-path "${OUTPUT_DIR}/rstudio/phyloseq-obj-1000-p/rooted-tree"

###########################################################################################################

# Check taxonomy bar graph
qiime taxa barplot \
  --i-table "${OUTPUT_DIR}/10-filter2/bac-1000-table-filt3.qza" \
  --i-taxonomy "${OUTPUT_DIR}/5-taxonomy/16S-taxonomy-vsearch.qza" \
  --m-metadata-file "metadata/${POOL}/metadata.tsv" \
  --o-visualization "${OUTPUT_DIR}/10-filter2/barplot-1000-p.qzv"
qiime tools view "${OUTPUT_DIR}/10-filter2/barplot-1000-p.qzv"

###########################################################################################################

# 11-Core Diversity

qiime diversity core-metrics-phylogenetic \
  --i-table "${OUTPUT_DIR}/10-filter2/bac-1000-table-filt3.qza" \
  --i-phylogeny "${OUTPUT_DIR}/8-phylogeny-1000-p/rooted_tree.qza" \
  --p-sampling-depth 8696 \
  --m-metadata-file "metadata/${POOL}/metadata.tsv" \
  --output-dir "${OUTPUT_DIR}/11-core-diversity-metrics-1000-p"

# Group Significance
## Alpha - Visually and statistically compare groups of alpha diversity values.
### Richness
qiime diversity alpha-group-significance \
  --i-alpha-diversity "${OUTPUT_DIR}/11-core-diversity-metrics-1000-p/shannon_vector.qza" \
  --m-metadata-file "metadata/${POOL}/metadata.tsv" \
  --o-visualization "${OUTPUT_DIR}/12-alpha-group/shannon-1000-p.qzv"

## Evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity "${OUTPUT_DIR}/11-core-diversity-metrics-1000-p/evenness_vector.qza" \
  --m-metadata-file "metadata/${POOL}/metadata.tsv" \
  --o-visualization "${OUTPUT_DIR}/12-alpha-group/evenness-1000-p.qzv"

## Beta - Determine whether groups of samples are significantly different from one another using a permutation-based statistical test.
### Change --m-metadata-column value to compare different features
### Permanova use option --p-pairwise  or --p-no-pairwise
qiime diversity beta-group-significance \
  --i-distance-matrix "${OUTPUT_DIR}/11-core-diversity-metrics-1000-p/bray_curtis_distance_matrix.qza" \
  --m-metadata-file "metadata/${POOL}/metadata.tsv" \
  --m-metadata-column cultivar \
  --p-method "permanova" \
  --p-pairwise \
  --o-visualization "${OUTPUT_DIR}/12-beta-group/comparing-cultivar-1000-p.qzv"

### Run .qzv plugin for each alpha diversity result
files="${OUTPUT_DIR}/11-core-diversity-metrics-1000-p"/*vector.qza
for result in $files; do 
outname=${result/_vector.qza/_group_significance.qzv}
qiime diversity alpha-group-significance \
--i-alpha-diversity $result \
--m-metadata-file "metadata/${POOL}/metadata.tsv" \
--o-visualization $outname
done

### Diversity results generated by study - export results in 1 file
#### Extraction diversity result by study from zip files (ex qza)
#### Transform the .qza in .zip files
#### might be necessary to run: chmod +x diversity-index-extraction-script.sh 
#### copy the script in "extraction" folder and execute the script - it will create a .txt file with all div results
#### Run the following command to execute the script - don't forget to change the name of the 1st file in the script: ./diversity-index-extraction-script.sh 

###########################################################################################################

# Export final

qiime tools export \
  --input-path "${OUTPUT_DIR}/10-filter2/bac-1000-table-filt3.qza" \
  --output-path "${OUTPUT_DIR}/rstudio/phyloseq-obj-1000-p"

 qiime tools export \
   --input-path "${OUTPUT_DIR}/10-filter2/bac-1000-rep-seqs-filt3.qza" \
   --output-path "${OUTPUT_DIR}/rstudio/phyloseq-obj-1000-p"

biom convert \
   -i "${OUTPUT_DIR}/rstudio/phyloseq-obj-1000-p/feature-table.biom" \
   -o "${OUTPUT_DIR}/rstudio/phyloseq-obj-1000-p/feature-table.txt" \
   --to-tsv

###########################################################################################################

# 12-Software for predicting functional abundances based only on marker gene sequences

source activate Picrust2
picrust2_pipeline.py \
   -s Subset3_16S-V4-MiSeq_rep-seq-FINAL-rarefied.fasta \
   -i Subset3-16S-V4-MiSeq_table-FINAL-rarefied.biom \
   -o picrust2_out_16S_V4_MiSeq_Subset3 \
   -p 7

