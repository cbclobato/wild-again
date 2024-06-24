# sudo apt install libharfbuzz-dev libfribidi-dev
# System dependencies:
# sudo apt install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

# install.packages(c("tidyverse", "tidymodels"))

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18")

# System dependencies:
# sudo apt install libcurl4-openssl-dev
# sudo apt install libxml2-dev
# sudo apt-get install gfortran libblas-dev liblapack-dev
# sudo apt-get -y install libglpk-dev

# BiocManager::install("phyloseq")
# BiocManager::install("metagenomeSeq")

# install.packages("devtools")
# devtools::install_github("vmikk/metagMisc")

library(phyloseq)
library(tidyverse)
library(broom)

# SHAP

# Save feature importance datasets

## Definitions

### Get Dataset

get_dataset <- function(ps, dataset = "domestication") {
  datasets <- c(
    "domestication"
  )
  if (dataset == "domestication") {
    # domestication ALL
    dataset_name <- "domestication"
    label_name <- "domestication"
    categories <- c(
      "Landrace",
      "Selection",
      "Cross hybrid",
      "Inbred line"
    )
    dataset_name <- dataset
    label_name <- dataset
    ps_subset <- ps # TODO: remove NA in biomarkers preprocessing
  } else {
    stop(
      paste(
        "\nUnknown dataset:", dataset, "\n",
        "List of datasets:\n\t", paste(datasets, collapse = "\n\t")
      )
    )
  }
  list(
    ps=ps_subset,
    name=dataset,
    label_name=label_name,
    categories=categories
  )
}


# Find Significant
find_significant <- function(ps, target, normalize = TRUE) {
  target <- enquo(target)
  if (normalize) {
    ps <- transform_sample_counts(ps, function(x) x/sum(x))
  }
  if (taxa_are_rows(ps)) {
    ps <- t(ps)
  }
  metadata <- as(sample_data(ps), "data.frame")
  counts_table <- as.data.frame(otu_table(ps))
  counts_table$target <- metadata %>% pull(!!target)
  significant_table <-
    counts_table %>%
    gather(ASV, val, -target) %>%
    group_by(ASV) %>%
    do(tidy(kruskal.test(val ~ target, .))) %>%
    filter(p.value < 0.05)
  return(significant_table)
}

# Preprocess
preprocess <- function(
    ps,
    dataset,
    normalize = FALSE,
    filter_features = FALSE,
    output_path = "data/biomarkers",
    verbose = TRUE
) {
  
  dir.create(output_path, FALSE, TRUE)
  
  if (verbose) {
    print_str <- paste0(
      "Dataset: ", dataset,
      ", Normalize: ", normalize,
      ", Filter: ", filter_features
    )
    print(print_str)
  }
  
  dataset <- get_dataset(ps, dataset)
  ps <- dataset$ps
  dataset_name <- dataset$name
  label_name <- dataset$label_name
  categories <- dataset$categories
  
  
  if (verbose == 2) {
    print_str <- paste0(
      "Dataset: ", dataset_name,
      ", Target: ", label_name,
      ", Categories: ", paste0(categories, collapse = ", ")
    )
    print(print_str)
  }
  
  if (normalize) {
    library(metagMisc)
    ps <- metagMisc::phyloseq_transform_css(ps, norm = TRUE, log = TRUE)
  }
  
  
  if (taxa_are_rows(ps)) {
    ps <- t(ps)
  }
  
  asv_table <-
    ps %>%
    otu_table %>%
    as.data.frame %>%
    rownames_to_column("sample_id")
  
  if (filter_features) {
    significant_table <- find_significant(ps, label_name)
    significant <- pull(significant_table, ASV)
    asv_table <- select(asv_table, sample_id, all_of(significant))
  }
  
  # metadata <- as(sample_data(ps), "data.frame")
  metadata <-
    ps %>%
    sample_data %>%
    as.matrix %>%
    data.frame %>%
    rownames_to_column("sample_id")
  
  
  # Join
  all_data <- inner_join(asv_table, metadata, by="sample_id")
  
  
  # Select Label
  label_sym <- sym(label_name)
  dataset <-
    all_data %>%
    select(genotype, starts_with("ASV", ignore.case = FALSE), !!label_name) %>%
    filter(!is.na(!!label_sym)) %>%
    mutate(
      "{{ label_sym }}" := factor(!!label_sym, levels=categories),
      genotype := factor(genotype)
    )
  
  if (verbose) {
    # Summary
    paste("Number of observations:", dim(dataset)[1])
    paste("Number of features:", dim(dataset)[2])
    paste("Number of genotypes (cultivars):", length(unique(dataset$genotype))) #
  }
  
  # Save
  prefix <- "all-"
  if (filter_features) {
    prefix  <- ""
  }
  
  suffix  <- ""
  if (normalize) {
    suffix  <- "-norm"
  }
  
  
  file_path <- file.path(
    output_path,
    paste0(prefix, dataset_name, suffix, ".RData")
  )
  print(file_path)
  save(dataset_name, dataset, label_name, categories, file=file_path)
  
  # Return
  list(
    name = dataset_name,
    dataset = dataset,
    label_name = label_name,
    categories = categories
  )
}


# UMAP

# Project ASVs such that samples are observations 
preprocess_umap_asv <- function(
    ps,
    dataset,
    normalize = TRUE,
    output_path = "data/umap"
) {
  
  dir.create(output_path, FALSE, TRUE)
  
  dataset <- get_dataset(ps, dataset)
  ps <- dataset$ps
  dataset_name <- dataset$name
  
  if (normalize) {
    library(metagMisc)
    ps <- metagMisc::phyloseq_transform_css(ps, norm = TRUE, log = TRUE)
  }
  
  if (taxa_are_rows(ps)) {
    ps <- t(ps)
  }
  
  asv_table <-
    ps %>%
    otu_table %>%
    as.data.frame %>%
    rownames_to_column("sample_id")
  
  metadata <-
    ps %>%
    sample_data %>%
    as.matrix %>%
    data.frame %>%
    rownames_to_column("sample_id")
  
  # Save
  prefix <- ""
  suffix <- ""
  if (normalize) {
    suffix <- "-norm"
  }
  
  file_path <- file.path(
    output_path,
    paste0(prefix, dataset_name, suffix, ".csv")
  )
  print(file_path)
  write.csv(asv_table, file_path, row.names = FALSE)
  
  file_path <- file.path(
    output_path,
    paste0(dataset_name, "-metadata.csv")
  )
  print(file_path)
  write.csv(metadata, file_path, row.names = FALSE)
  
  
  # Return
  list(asv_table = asv_table, metadata = metadata)
}


# Project Samples such that ASVs are observations
preprocess_umap_samples <- function(
    ps,
    dataset,
    normalize = TRUE,
    output_path = "data/umap"
) {
  
  dir.create(output_path, FALSE, TRUE)
  
  dataset <- get_dataset(ps, dataset)
  ps <- dataset$ps
  dataset_name <- dataset$name
  
  if (normalize) {
    library(metagMisc)
    ps <- metagMisc::phyloseq_transform_css(ps, norm = TRUE, log = TRUE)
  }
  
  if (taxa_are_rows(ps)) {
    ps <- t(ps)
  }
  
  asv_table <-
    ps %>%
    otu_table %>%
    as.data.frame %>%
    rownames_to_column("sample_id")
  
  sample_data <-
    ps %>%
    otu_table %>%
    as.data.frame %>%
    t %>%
    as.data.frame %>%
    rownames_to_column("ASV")
  
  sample_table <-
    sample_data %>%
    select(ASV, starts_with("C", ignore.case = FALSE))

  # Save
  prefix <- ""
  suffix <- ""
  if (normalize) {
    suffix <- "-norm"
  }
  
  file_path <- file.path(
    output_path,
    paste0("t_", prefix, dataset_name, suffix, ".csv")
  )
  print(file_path)
  write.csv(sample_table, file_path, row.names = FALSE)

  # Return
  list(sample_table = sample_table)
}


