#!/bin/bash

# $(OUTPUT_DIR)/%.fastq.gz
POOL=all
cutadapta \
  --errors 2 \
  --no-indels \
  --pair-adapters \
  -g "file:metadata/${POOL}-bc-fw.fasta" \
  -G "file:metadata/${POOL}-bc-rv.fasta" \
  -o "${OUTPUT_DIR}/demux/{name}-1.fastq.gz" \
  -p "${OUTPUT_DIR}/demux/{name}-2.fastq.gz" \
  "raw/${POOL}/forward.fastq.gz" "raw/${POOL}/reverse.fastq.gz" \
  --discard-untrimmed \
  --cores 6 \
  --json cutadapt.json
