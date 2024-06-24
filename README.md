## Wild again: Recovery of a beneficial Cannabis seed endophyte from low domestication genotypes
Authors: Carolina Lobato,  João Machado de Freitas, Daniel Habich, Gabriele Berg, Tomislav Cernava

Institute of Environmental Biotechnology (UBT), Graz University of Technology, Graz (Austria)
contact: cbotocourinhalobato@tugraz.at

### Project structure
```text
project/
├── data/
│       └── metadata/
├── scripts/
│   ├── qiime2/
│   │   └── bioprocessing_pipeline.sh
│   ├── r/
│   │   ├── Setup.Rmd
│   │   ├── Figure1.Rmd
│   │   ├── Figure2.Rmd
│   │   ├── Figure3.Rmd
│   │   ├── Figure4.Rmd
│   │   └── Figure5.Rmd
│   └── utils/
│       ├── csv2fasta.sh
│       ├── csv2tsv.sh
│       ├── qiime2r.sh
│       ├── install.R
│       └── plot_composition_v2.R
├── README.md
└── LICENSE
```

### Details   
- data/
  - metadata/  This subdirectory contains:
      - fasta barcode files (bcpr-fw.fasta and bcpr-rv.fasta)
      - the concatenating files (cat.sh),
      - the manifest files (manifest.csv),
      - the sample information for the analysis (metadata.csv)

- scripts/
  - qiime2/  This subdirectory contains the pipeline used for:
      - demultiplexing with CUTADAPT v4.2,
      - importing into QIIME2 v2023.5 and the further bioinformatic processing steps using the DADA2 pipeline and the VSEARCH algorithm using the SILVA v138 reference database, which generated the feature table, taxonomy file, representative sequences and phylogenetic tree
      - exporting from QIIME2
    
  - r/  This subdirectory contains the scripts used in R to create the phyloseq objects, preprocess the data and prepare the figures in the manuscript.
    
  -  utils/  This subdirectory contains utility scripts that are used by other scripts in the project, such as:
      - csv2fasta.sh for converting .csv to .fasta format,
      - csv2tsv.sh for converting .csv to .tsv format
      - qiime2r.sh for converting .biom to .tsv format and back,
      - install.R for installing the necessary packages in R.
      - plot_composition_v2 modified microbiome::plot_composition function for running when the microbiome version is above 1.6.
   
### Further content
The 16S rRNA gene amplicon raw FASTQ files were deposited at the European Nucleotide Archive (ENA; https://www.ebi.ac.uk/ena) under the accession number PRJEB64469.

The assembled genome of Bacillus frigotolerans, with the associated annotations, was deposited in the National Center for Biotechnology Information (NCBI; https://www.ncbi.nlm.nih.gov/) under accession number PRJNA1113337.

### References
