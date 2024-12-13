## Wild again: Recovery of a beneficial Cannabis seed endophyte from low domestication genotypes
Authors: Carolina Lobato, João Machado de Freitas, Daniel Habich, Isabella Kögl, Gabriele Berg & Tomislav Cernava

[Institute of Environmental Biotechnology (UBT)](https://www.tugraz.at/institute/ubt/home/) — TU Graz


### Project structure
```text
project/
├── data/
│   └── metadata/
│       ├── CB1/
│       │   ├── bcpr-fw.fasta
│       │   ├── bcpr-rv.fasta
│       │   ├── cat.sh
│       │   └── manifest.csv
│       ├── CB2/
│       │   ├── bcpr-fw.fasta
│       │   ├── bcpr-rv.fasta
│       │   ├── cat.sh
│       │   └── manifest.csv
│       ├── CB3/
│       │   ├── bcpr-fw.fasta
│       │   ├── bcpr-rv.fasta
│       │   ├── cat.sh
│       │   └── manifest.csv
│       ├── metadata.csv
│       ├── pouches-ind.tsv
│       ├── pouches-exp.tsv
│       ├── field23.tsv
│       └── PLaBase.tsv
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
│       ├── plot_composition_v2.R
│       ├── umap/
│       │   ├── data.R
│       │   ├── project.R
│       │   ├── cluster-analysis.R
│       │   └── run-all.R
│       └── biomarkers/
│           ├── data.R
│           ├── feature-importance.R
│           ├── sv-importance-patch.R
│           ├── train-eval.R 
│           └── run-all.R
├── outputs/
│   ├── qiime2/
│   └── r/
├── README.md
└── LICENSE
```

### Details   
- data/
  - metadata/  This subdirectory contains:
      - fasta barcode files for each pool (bcpr-fw.fasta and bcpr-rv.fasta).
      - the concatenating files for each pool (cat.sh).
      - the manifest files for each pool (manifest.csv).
      - the sample information for the metabarcoding analysis (metadata.csv).
      - the metadata and measurements at individual-level in lab trials (pouches-ind.tsv)
      - the metadata and measurments related to each replicate experiment in lab trials (pouches-exp.tsv)
      - the metadata and measurements of the field trials (field23.tsv)
      - the output gene class table from PLaBase (PLaBase.tsv)

- scripts/
  - qiime2/  This subdirectory contains the pipeline (bioprocessing_pipeline.sh) used for:
      - demultiplexing with CUTADAPT v4.2.
      - importing into QIIME2 v2023.5 and the further bioinformatic processing steps using the DADA2 pipeline and the VSEARCH algorithm using the SILVA v138 reference database, which generated the feature table, taxonomy file, representative sequences and phylogenetic tree.
      - exporting from QIIME2.
    
  - r/  This subdirectory contains the scripts used in R to create the phyloseq objects, preprocess the data and prepare the figures in the manuscript.
    
  -  utils/  This subdirectory contains utility scripts that are used by other scripts in the project, such as:
      - csv2fasta.sh for converting .csv to .fasta format.
      - csv2tsv.sh for converting .csv to .tsv format.
      - qiime2r.sh for converting .biom to .tsv format and back,
      - install.R for installing the necessary packages in R.
      - plot_composition_v2 modified microbiome::plot_composition function for running when the microbiome version is above 1.6.
      - umap/ contains the scrips used for beta diversity representation with UMAP shown in Figure2.
      - biomarkers/ contains the scripts used for biomarker assessment shown in Figure3.
   
- outputs/ contains qiime2 and r saved outputs.
   
### Further content
The 16S rRNA gene amplicon raw FASTQ files were deposited in [ENA](https://www.ebi.ac.uk/ena) under the accession number PRJEB64469.

The assembled genome of *Bacillus frigotolerans*, with the associated annotations, was deposited in [NCBI](https://www.ncbi.nlm.nih.gov/) under accession number PRJNA1113337.

### References
Lobato, C., de Freitas, J.M., Habich, D. et al. Wild again: recovery of a beneficial Cannabis seed endophyte from low domestication genotypes. Microbiome 12, 239 (2024). https://doi.org/10.1186/s40168-024-01951-5
