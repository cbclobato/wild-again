#!/bin/bash

conda activate genome_env

# Check whats installed within an environment and install necessary tools
conda list # | grep flye
conda install -c defaults -c conda-forge -c bioconda pbtk

# Setup
export SAMPLE=pacbio
export OUTPUT_DIR="outputs/${SAMPLE}"
export DATASET="data/raw/${SAMPLE}"

# FLYE
## Convert to .FASTA
bam2fasta -o "${DATASET}/Peribacillus" \
"${DATASET}/m64291e_240314_122818.subreads.bam"

## Check number of subreads with report (CLR reads = 370273)
grep -c ">" "${DATASET}/Peribacillus.fasta"

## Run
### Estimated genome size: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_024169475.1/
flye \
--pacbio-raw "${DATASET}/Peribacillus.fasta" \
--out-dir "${OUTPUT_DIR}/flye" \
--genome-size 5.5m \
--threads 32 \
--debug

# CHECKM (informs about the quality of the genome assembly in terms of completeness and contamination)
checkm lineage_wf \
-x fasta \
--reduced_tree "${OUTPUT_DIR}/flye" \
"${OUTPUT_DIR}/flye/checkm" \
-t 20

# gtdbt (for taxonomy)
gtdbtk identify \
--genome_dir "${OUTPUT_DIR}/flye" \
--out_dir "${OUTPUT_DIR}/flye/gtdbk" \
-x fasta \
--cpus 10

gtdbtk align \
--identify_dir "${OUTPUT_DIR}/flye/gtdbk" \
--out_dir "${OUTPUT_DIR}/flye/gtdbk_align" \
--cpus 10

gtdbtk classify \
--genome_dir "${OUTPUT_DIR}/flye" \
--align_dir "${OUTPUT_DIR}/flye/gtdbk_align" \
--out_dir "${OUTPUT_DIR}/flye/gtdbk_classify" \
-x fasta \
--cpus 10

# DRAM
DRAM.py annotate \
-i "${OUTPUT_DIR}/flye/*.fasta" \
-o "${OUTPUT_DIR}/flye/DRAM"

DRAM.py distill \
-i "${OUTPUT_DIR}/flye/DRAM/annotations.tsv" \
-o "${OUTPUT_DIR}/flye/DRAM/genome_summaries" \
--trna_path "${OUTPUT_DIR}/flye/DRAM/trnas.tsv" \
--rrna_path "${OUTPUT_DIR}/flye/DRAM/rrnas.tsv"

# PLABASE: https://plabase.cs.uni-tuebingen.de/pb/plabase.php > PGPT-Pred > Upload genes.faa from DRAM > balastp+hmmer
