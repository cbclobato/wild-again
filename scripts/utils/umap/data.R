library(phyloseq)
library(tidyverse)


# Project ASVs
preprocess_umap_asv <- function(ps,
                                name,
                                normalize = TRUE,
                                output_path = "data/umap",
                                save = TRUE) {
  dir.create(output_path, FALSE, TRUE)
  
  if (normalize) {
    library(metagMisc)
    ps <-
      metagMisc::phyloseq_transform_css(ps, norm = TRUE, log = TRUE)
  }
  
  if (taxa_are_rows(ps)) {
    ps <- t(ps)
  }
  
  asv_table <-
    ps %>%
    otu_table %>%
    as.data.frame # %>%
    # rownames_to_column("sample_id")
  
  metadata <-
    ps %>%
    sample_data %>%
    as.matrix %>%
    data.frame #%>%
    # rownames_to_column("sample_id")
  
  # Save
  if (save) {
    prefix <- ""
    suffix <- ""
    if (normalize) {
      suffix <- "-norm"
    }
    
    file_path <- file.path(output_path,
                           paste0(prefix, name, suffix, ".csv"))
    write.csv(asv_table, file_path, row.names = TRUE)
    
    file_path <- file.path(output_path,
                           paste0(name, "-metadata.csv"))
    write.csv(metadata, file_path, row.names = TRUE)
  }
  
  # Return
  list(data = asv_table, metadata = metadata)
}


# ASVs observations
preprocess_umap_samples <- function(ps,
                                    dataset,
                                    normalize = TRUE,
                                    output_path = "data/umap") {
  dir.create(output_path, FALSE, TRUE)
  
  dataset <- get_dataset(ps, dataset)
  ps <- dataset$ps
  dataset_name <- dataset$name
  
  if (normalize) {
    library(metagMisc)
    ps <-
      metagMisc::phyloseq_transform_css(ps, norm = TRUE, log = TRUE)
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
  
  file_path <- file.path(output_path,
                         paste0("t_", prefix, dataset_name, suffix, ".csv"))
  write.csv(sample_table, file_path, row.names = FALSE)
  
  # Return
  list(sample_table = sample_table)
}
