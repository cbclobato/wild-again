## Linux
# sudo apt install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
#  libharfbuzz-dev libfribidi-dev fzf cmake libudunits2-dev libgdal-dev

# igraph (phyloseq)
# sudo apt-get install gfortran libblas-dev liblapack-dev

# ggpattern
# sudo apt install libudunits2-dev libgdal-dev

## R
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install(c(
   "phyloseq",
   "microbiome",
   "metagenomeSeq",
   "biomformat",
   "treeio",
   "Biostrings",
   "decontam"
))

# Palettes:
# https://github.com/EmilHvitfeldt/r-color-palettes

install.packages(c(
  "tidyverse",
  "devtools",
  "remotes",
  "curl",
  "vegan",
  "metacoder",
  "genefilter",
  "DESeq2",
  "RColorBrewer",
  "ggsignif",
  "pheatmap",
  "gdtools",
  "ggrepel",
  "metacoder",
  "stats",
  "gplots",
  "VennDiagram",
  "adespatial",
  "reshape2",
  "multcompView",
  "spaa",
  "tidymodels",
  "tidytext",
  "treemap",
  "asremlPlus",
  "recipes",
  "png",
  "umap",
  "plotly",
  "network",
  "ggnetwork",
  "GGally",
  "intergraph",
  "RMThreshold".
  "ggtext",
  "UpSetR",
  repos = "http://R-Forge.R-project.org"
))

remotes::install_github(c(
   "vmikk/metagMisc",
   "cpauvert/psadd",
   "umerijaz/microbiomeSeq",
   "gauravsk/ranacapa",
   "coolbutuseless/ggpattern",
   "KarstensLab/microshades",
   "mattmar/rasterdiv"
))

devtools::install_github(c(
   "microsud/microbiomeutilities",
   "kassambara/ggpubr",
   "mikemc/speedyseq",
   "pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",
   "trinker/qdapTools",
   "helixcn/spaa",
   "GuillemSalazar/EcolUtils",
   "zdk123/SpiecEasi",
   "david-barnett/microViz",
   "briatte/ggnet")
)

# From summaries. Needed?
# install.packages('Rmisc', dependencies = TRUE)
# install.packages('rsvg')
# remotes::install_github('coolbutuseless/ggsvg')