# UMAP
# conda create -n umap -c conda-forge umap-learn

# Imports ======================================================================
rm(list = ls())

library(RColorBrewer)


source("scripts/r/umap/data.R")
source("scripts/r/umap/project.R")
source("scripts/r/umap/cluster-analysis.R")

# Config =======================================================================

# General
seed <- 1234
set.seed(seed)
input_path <- "data/umap"
output_path <- "outputs/umap"
dir.create(input_path, FALSE, TRUE)
dir.create(output_path, FALSE, TRUE)

# Preprocessing
true_false <- c(TRUE)
# true_false <- c(TRUE, FALSE)

# Data
name <- "domestication"
normalized <- TRUE
domestication_order <-
  c("Landrace", "Selected line",  "Cross hybrid", "Inbred line", "NA")
dom2_order <- c("High", "Low", "NA")

# Eval projection
evaluate <- TRUE

# UMAP
cfg = umap.defaults
cfg$n_components = 2
cfg$metric = "braycurtis"
cfg$n_neighbors = 45
cfg$min_dist = 0.8
cfg$densmap = TRUE
cfg$n_epochs = 500
cfg$init = "pca"
cfg$random_seed = 5432
# cfg$n_jobs=1

# Figure
label_name <- "domestication"
# label_name <- "dom2"


# Preprocessing ================================================================
print("Preprocessing...")
# true_false <- c(TRUE, FALSE)

## Import metadata file
metadata <- read.table(
  "outputs/Tax4Fun/output_Ref99NR/Fun-Met-adj.csv",
  sep = ";",
  header = TRUE,
  fill = TRUE,
  row.names = "KO",
  na.strings = "")

## Import feature table
asv_tab <- read.table("outputs/Tax4Fun/output_Ref99NR/Fun-OTU-adj.txt",
                      header = TRUE,
                      row.names = "KO")

## Import taxonomy file
tax_tab <- read.table(
  "outputs/Tax4Fun/output_Ref99NR/Fun-Tax.tsv",
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  row.names = "KO")

## Convert to phyloseq object
# ps <- merge_phyloseq(
#   metadata,
#   tax_tab,
#   asv_tab)
# ps # 8949 taxa and 451 samples
# save(ps,
#      file = "outputs/Tax4Fun/output_Ref99NR/ps-tax.RData")
#load("outputs/Tax4Fun/output_Ref99NR/ps-tax.RData")

dataset <- list(data = as.data.frame(t(asv_tab)), metadata = metadata)


# Project ASV counts ===========================================================
print("Run UMAP...")
print(paste("Dataset:", name))
df_umap <- fit_transform(dataset, cfg)


# Post-processing ==============================================================
label_sym <- sym(label_name)
df_postproc <-
  df_umap %>%
  mutate(
    domestication = fct_relevel(replace_na(domestication, "NA"), domestication_order),
    dom2 = fct_relevel(replace_na(dom2, "NA"), dom2_order),
    chemotype2 = as.factor(replace_na(chemotype2, "NA"))
  ) %>%
  rename(
    Domestication = !!label_sym,
    Chemotype = chemotype2,
    Genotype = genotype
  )

# save(df_postproc, file = "outputs/Tax4Fun/output_Ref99NR/UMAP.RData")

# Cluster analysis =============================================================
source("scripts/r/umap/cluster-analysis.R")
embeds <- df_postproc[,c("UMAP 1", "UMAP 2")]
labels <- df_postproc[,"Domestication"]
p <- silhouette_plot(embeds, labels)
print(p)


# Figures ======================================================================
p1 <-
  ggplot(df_postproc,
         aes(
           x = `UMAP 1`,
           y = `UMAP 2`,
           color = Domestication,
           shape = Chemotype
         )) +
  geom_point() +
  scale_color_discrete(type = palette_domestication) +
  scale_shape_manual(values = c(19, 1, 4)) + 
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme_classic()
print(p1)

#plotly::ggplotly(p1)

palette_genotype <- colorRampPalette(brewer.pal(46, "Paired"))(46)
p2 <- ggplot(df_postproc, aes(x = `UMAP 1`, y = `UMAP 2`, color = Genotype)) +
  geom_point() +
  scale_color_discrete(type = palette_genotype) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme_classic()
print(p2)

#plotly::ggplotly(p2)

# Query ========================================================================
# domestication_val = "Cross hybrid"
# chemotype_val = "Low"
# df_postproc %>%
#   filter(`UMAP 2` < 8 & `UMAP 1` < 10 &
#            Domestication == "Cross hybrid" &
#            Chemotype == "Low")


# End ==========================================================================
print("Complete!")
