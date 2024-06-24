library(phyloseq)
library(tidyverse)
library(broom)


get_dataset <- function(ps, dataset = "chemotype") {
  datasets <- c(
      "cultivar",
      "chemotype",
      "hemp",
      "flower",
      "marijuana-chemotype",
      "marijuana-domestication"
  )
  if (dataset == "cultivar") {
    # All cultivars
    dataset_name <- dataset
    label_name <- dataset
    categories <- sort(unique(sample_data(ps)$cultivar))
    ps_subset <- ps
  } else if (dataset == "chemotype") {
     # chem (Subset B, n=8)
    dataset_name <- "chemotype"
    label_name <- "chemotype2"
    categories <- c("High", "Low")
    ps_subset <-
      subset_samples(
        ps,
        !is.na(chemotype)
        & inflorescence == "Dioecious"
        & domestication == "Cross hybrid"
        
      )  
  } else if (dataset == "hemp") {
    # dom-hemp (Subset C, n=26)
    dataset_name <- "hemp"
    label_name <- "domestication"
    # categories <- c("Hybrid", "Landrace", "Wild/Feral")
    categories <- c("Inbred line", "Landrace", "Cross hybrid", "Selection")
    ps_subset <-
      subset_samples(
        ps,
        !is.na(domestication)
        & chemotype2 == "Low"
      )
  } else if (dataset == "flower") {
    # Flower (Subset A, n=13)
    dataset_name <- "flower"
    label_name <- "inflorescence"
    # categories <- c("Monoecious", "Dioecious")
    categories <- c("Monoecious", "Subdioecious", "Dioecious")
    ps_subset <-
      subset_samples(
        ps,
        !is.na(inflorescence)
        & chemotype == "Low"
        & domestication %in% c("Cross hybrid", "Inbred line")
      )
  } else if (dataset == "marijuana-chemotype") {
    dataset_name <- "marijuana-chemotype"
    label_name <- "chemotype"  # CBD-rich, CBD/THC-rich, THC-rich
    categories <- c("THC-rich", "CBD/THC-rich", "CBD-rich")
    ps_subset <-
      subset_samples(
        ps,
        !is.na(domestication) 
        & chemotype2 == "High"
      )
  } else if (dataset == "marijuana-domestication") {
    # dom-chem-mj (Subset D, n=9)
    dataset_name <- "marijuana-domestication"
    label_name <- "domestication"
    # categories <- c("Hybrid", "Landrace", "Unstable hybrid", "Hybrid S1")
    categories <- c(
      # "Unstable hybrid",
      "Segregating hybrid",
      "Cross hybrid S1",
      "Landrace",
      "Cross hybrid",
      "Selection"
    )
    ps_subset <-
      subset_samples(
        ps,
        !is.na(domestication)
        & chemotype2 == "High"
      )
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


preprocess <- function(
    ps,
    dataset,
    normalize = FALSE,
    filter_features = FALSE,
    output_path = "outputs/r/310_bac_coocur",
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
    # ps <- transform_sample_counts(ps, function(x) x/sum(x))
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
    select(cultivar, starts_with("ASV", ignore.case = FALSE), !!label_name) %>%
    filter(!is.na(!!label_sym)) %>%
    mutate(
      "{{ label_sym }}" := factor(!!label_sym, levels=categories),
      cultivar := factor(cultivar)
    )
  
  if (verbose) {
    # Summary
    paste("Number of observations:", dim(dataset)[1])
    paste("Number of features:", dim(dataset)[2])
    paste("Number of cultivars:", length(unique(dataset$cultivar)))    
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
  save(dataset_name, dataset, label_name, categories, file=file_path)
  
  # Return
  list(
    name = dataset_name,
    dataset = dataset,
    label_name = label_name,
    categories = categories
  )
}


## Samples are observations
preprocess_umap_asv <- function(
    ps,
    dataset,
    normalize = TRUE,
    output_path = "outputs/r/310_bac_coocur"
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
  write.csv(asv_table, file_path, row.names = FALSE)
  
  file_path <- file.path(
    output_path,
    paste0(dataset_name, "-metadata.csv")
  )
  write.csv(metadata, file_path, row.names = FALSE)
  

  # Return
  list(asv_table = asv_table, metadata = metadata)
}


### Project Samples
## ASVs are observations
preprocess_umap_samples <- function(
    ps,
    dataset,
    normalize = TRUE,
    output_path = "outputs/r/310_bac_coocur"
  ) {
  
  dir.create(output_path, FALSE, TRUE)
  
  dataset <- get_dataset(ps, dataset)
  ps <- dataset$ps
  dataset_name <- dataset$name
  
  
  if (normalize) {
    library(metagMisc)
    ps <- metagMisc::phyloseq_transform_css(ps, norm = TRUE, log = TRUE) 
    # ps <- transform_sample_counts(ps, function(x) x/sum(x))
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
  write.csv(sample_table, file_path, row.names = FALSE)

  # Return
  list(sample_table = sample_table)
}


dir.create("data/r/310_bac_coocur", FALSE, TRUE)

dataset <- "cultivar"
normalize <- TRUE
filter_features <-TRUE

load("outputs/r/204_setup/bh.RData")

# Saves the table
output <- preprocess_umap_samples(
  bh,
  dataset = dataset,
  normalize = normalize,
)
output
output[["sample_table"]]


#install.packages("factoextra")
#install.packages("ggfortify")
library(factoextra)
library(corrplot)
library(ggfortify)

df <- output[["sample_table"]]
summary(df)
head(df)
row.names(df) <- df$ASV
numerical_data <- df[,2:452]
head(numerical_data)
res.pca <- prcomp(numerical_data, scale = TRUE)
print(res.pca)
summary(res.pca)
eig.val <- get_eigenvalue(res.pca); eig.val
fviz_eig(res.pca, col.var = "blue")
var <- get_pca_var(res.pca); var
head(var$cos2)
corrplot(var$cos2, is.corr = FALSE)
fviz_cos2(res.pca, choice = "var", axes = 1:2)

fviz_pca_var(res.pca,
             col.var = "cos2", # Color by the quality of representation
             label = "var",
             gradient.cols = c("darkorchid4", "gold", "darkorange"),
             repel = TRUE)

# Contributions of variables to PC1
a <- fviz_contrib(res.pca, choice = "var", axes = 1)
# Contributions of variables to PC2
b <- fviz_contrib(res.pca, choice = "var", axes = 2)
gridExtra::grid.arrange(a, b, ncol = 2, 
                        top = 'Contribution of the variables to the first two PCs')

ind <- get_pca_ind(res.pca); ind
name = "ASV 3055: Bacillus_sp."
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             label = "none",
             pointsize = 3,
             gradient.cols = c("darkorchid4", "gold", "darkorange"),
             # select.ind = list(contrib = 20),
             repel = TRUE) +
  theme_classic()

fviz_pca_biplot(res.pca, 
                label = "var", 
                habillage=df$Species,
                addEllipses=TRUE, ellipse.level=0.95)

autoplot(res.pca, 
         loadings = TRUE,
         loadings.colour = 'darkorchid4', 
         loadings.label = TRUE, 
         loadings.label.size = 3)

kmeans <- eclust(df, k = 4)
autoplot(res.pca, data = kmeans, colour = "cluster")

