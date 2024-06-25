# UMAP
# conda create -n umap -c conda-forge umap-learn

# Imports ======================================================================
rm(list = ls())

library(phyloseq)

source("scripts/utils/umap/data.R")
source("scripts/utils/umap/project.R")

# Config =======================================================================

# General
seed <- 1234
set.seed(seed)
input_path <- "outputs/r/fig_2/umap"
output_path <- "outputs/r/fig_2/umap"
dir.create(input_path, FALSE, TRUE)
dir.create(output_path, FALSE, TRUE)

# Preprocessing
true_false <- c(TRUE)

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
cfg$min_dist = .5
cfg$densmap = TRUE
cfg$n_epochs = 500
cfg$init = "pca"
cfg$random_seed = 5432
# cfg$n_jobs=1

# Figure
label_name <- "domestication"

# Preprocessing ================================================================
print("Preprocessing...")
# true_false <- c(TRUE, FALSE)


load("outputs/r/setup/bh.RData")
dataset <- preprocess_umap_asv(bh,
                               name = name,
                               normalize = normalized)


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

save(df_postproc, file = "outputs/r/fig_2/UMAP.RData")

# End ==========================================================================
print("Complete!")
