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


get_dataset <- function(ps, name = "domestication") {
  # Get Dataset
  names <- c("domestication", "domestication-2")
  if (name == "domestication") {
    # domestication MULTICLASS
    dataset_name <- "domestication"
    label_name <- "domestication"
    categories <- c("Landrace",
                    "Selection",
                    "Cross hybrid",
                    "Inbred line")
    dataset_name <- name
    ps_subset <- ps
  } else if (name == "domestication-2") {
    # domestication BINARY
    dataset_name <- "domestication-2"
    label_name <- "dom2"
    categories <- c("High", "Low")
    dataset_name <- name
    ps_subset <- ps
  } else {
    stop(paste(
      "\nUnknown dataset:",
      name,
      "\n",
      "List of datasets:\n\t",
      paste(names, collapse = "\n\t")
    ))
  }
  list(
    ps = ps_subset,
    name = name,
    label_name = label_name,
    categories = categories
  )
}


# Find Significant
find_significant <- function(ps, target, normalize = TRUE) {
  target <- enquo(target)
  if (normalize) {
    ps <- transform_sample_counts(ps, function(x)
      x / sum(x))
  }
  if (taxa_are_rows(ps)) {
    ps <- t(ps)
  }
  metadata <- as(sample_data(ps), "data.frame")
  counts_table <- as.data.frame(otu_table(ps))
  counts_table$target <- metadata %>% pull(!!target)
  significant_table <-
    counts_table %>%
    gather(ASV, val,-target) %>%
    group_by(ASV) %>%
    do(tidy(kruskal.test(val ~ target, .))) %>%
    filter(p.value < 0.05)
  return(significant_table)
}


# Preprocess
preprocess <- function(ps,
                       name,
                       normalize = FALSE,
                       filter_features = FALSE,
                       output_path = "data/biomarkers",
                       verbose = TRUE,
                       save_output = TRUE) {
  normalized <- normalize
  filtered <- filter_features
  
  if (verbose) {
    cat("\nDataset: ",
        name,
        ", Normalize: ",
        normalize,
        ", Filter: ",
        filter_features)
  }
  
  dataset <- get_dataset(ps, name)
  ps <- dataset$ps
  name <- dataset$name
  label_name <- dataset$label_name
  categories <- dataset$categories
  
  if (verbose) {
    cat(
      "\nDataset: ",
      name,
      "\nTarget: ",
      label_name,
      "\nCategories: ",
      paste0(categories, collapse = ", ")
    )
  }
  
  if (normalize) {
    library(metagMisc)
    ps <-
      metagMisc::phyloseq_transform_css(ps, norm = TRUE, log = TRUE)
  }
  
  
  if (taxa_are_rows(ps)) {
    ps <- t(ps)
  }
  
  # Get count table
  asv_table <-
    ps %>%
    otu_table %>%
    as.data.frame %>%
    rownames_to_column("sample_id")
  
  
  if (filter_features) {
    significant_table <- find_significant(ps, label_name)
    significant <- pull(significant_table, ASV)
    asv_table <- dplyr::select(asv_table, sample_id, all_of(significant))
  }
  
  # Get metadata
  metadata <-
    ps %>%
    sample_data %>%
    as.matrix %>%
    data.frame %>%
    rownames_to_column("sample_id")
  
  # Join
  all_data <- inner_join(asv_table, metadata, by = "sample_id")
  
  # Select Label
  label_sym <- sym(label_name)
  data <-
    all_data %>%
    dplyr::select(genotype,
           starts_with("ASV", ignore.case = FALSE),
           !!label_name) %>%
    filter(!is.na(!!label_sym)) %>%
    mutate(
      "{{ label_sym }}" := factor(!!label_sym, levels = categories),
      genotype := factor(genotype)
    )
  
  if (verbose) {
    # Summary
    paste("\nNumber of observations:", dim(data)[1])
    paste("Number of features:", dim(data)[2])
    paste("Number of genotypes (cultivars):", length(unique(data$genotype))) #
  }
  
  output <- list(
    name = name,
    data = data,
    label_name = label_name,
    categories = categories,
    normalized = normalized,
    filtered = filtered
  )
  
  if (save_output) {
    # Save
    dir.create(output_path, FALSE, TRUE)
    prefix <- "all-"
    if (filter_features) {
      prefix  <- ""
    }
    
    suffix  <- ""
    if (normalize) {
      suffix  <- "-norm"
    }
    
    file_path <- file.path(output_path,
                           paste0(prefix, name, suffix, ".RData"))
    save(name,
         data,
         label_name,
         categories,
         normalized,
         filtered,
         file = file_path)
  }
  
  # Return
  output
}
