# Feature Importance
#install.packages("discrim")
#install.packages("xgboost")
#install.packages("shapviz")

# Imports ======================================================================
rm(list = ls())

# library(foreach)
# library(doParallel)

source("scripts/r/biomarkers/data.R")
source("scripts/r/biomarkers/train-eval.R")
source("scripts/r/biomarkers/feature-importance.R")

# cluster <- makeCluster(2)
# registerDoParallel(cluster)


# Config =======================================================================

# General
seed <- 12345
input_path <- "data/biomarkers"
output_path <- "outputs/biomarkers"

# Preprocessing
has_preprocess <- TRUE
true_false <- c(TRUE)
# true_false <- c(TRUE, FALSE)

# Data
dataset_names <- c("domestication")
# dataset_names <- c("domestication", "domestication-2")
filtered <- TRUE
normalized <- TRUE

# Models
# model_names <- c("xgboost")
model_names <- c("lda", "xgboost")
eval_model <- FALSE

# Feature importance
top_n <- 10
fill <- "#fca50a"
label_color <- "#769C64"# , "#4A623E", "#94C47D"

dir.create(input_path, FALSE, TRUE)
dir.create(output_path, FALSE, TRUE)
set.seed(seed)

# Preprocessing ================================================================
print("Preprocessing...")
# true_false <- c(TRUE, FALSE)


if (has_preprocess) {
  # foreach(dataset_name = dataset_names) %dopar% {
  for (dataset_name in dataset_names) {
    for (normalize in true_false) {
      for (filter_features in true_false) {
        load("outputs/r/supplementary/bh.RData")
        dataset <- preprocess(bh,
                              dataset_name,
                              normalize = normalize,
                              filter_features = filter_features)
        print("ASV 3055: Bacillus_sp." %in% names(dataset$data))
      }
    }
  }
}


# Train and Eval ===============================================================
print("Training and evaluation...")
results <- list()
# results <- foreach(dataset_name = dataset_names) %dopar% {
for (dataset_name in dataset_names) {
  print(paste("Dataset:", dataset_name))
  
  # list: name, dataset, categories, label_name, normalized, filtered
  dataset <- load_preprocessed_data(dataset_name,
                                    filtered = filtered,
                                    normalized = normalized)
  
  for (model_name in model_names) {
    print(paste("Model:", model_name))
    
    # Train and eval
    results[[dataset_name]][[model_name]] <-
      train_eval(dataset, model_name, eval = eval_model, seed = seed)
    results[[dataset_name]][["categories"]] <- dataset$categories
  }
}

# Feature Importance ===========================================================
print("Feature Importance and Biomarker discovery...")
results.imp <- list()
for (dataset_name in dataset_names) {
  source("scripts/r/biomarkers/feature-importance.R")
  result <- results[[dataset_name]]
  lda_fit_engine  <- result[["lda"]]$fit_engine
  xgb_fit_engine  <- result[["xgboost"]]$fit_engine
  data_mat <- result[["xgboost"]]$data_mat
  categories <- results[[dataset_name]][["categories"]]
  results.imp[[dataset_name]] <-
    feature.importance(lda_fit_engine,
                       xgb_fit_engine,
                       data_mat,
                       categories,
                       top_n = top_n)
}


# Figures ======================================================================
print("Plots...")
mycolors <- c("#ce232a","#918bc0","#ff6682")
plots <- list()
for (dataset_name in dataset_names) {
  df_plot <- results.imp[[dataset_name]] %>% 
    filter(Kind != "Gain") 
save(df_plot, file = "outputs/r/fig_3/biomarkers.RData")   

  p <-
    ggplot(df_plot,
           aes(
             x = Importance,
             y = ASV_ordered,
             color = Kind,
             fill = Kind
           )) +
    # geom_bar(stat="identity", position="dodge") +
    geom_point(position = position_dodge(width = 0.8))
  if (dataset_name == "domestication") {
    p <- p +
      facet_wrap(~ Domestication, 
                 scales = "free", 
                 nrow = 1)
  }
  plots[[dataset_name]] <- p +
    scale_y_reordered() +
    scale_color_manual(values = mycolors) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          strip.background = element_blank(),
          strip.text = element_text(size = 10,
                                    face = "bold"),
          legend.position = "bottom",
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm")) +
    guides(fill = guide_legend(nrow = 1)) + 
    labs(x = "", 
         y = "")
  # ggsave
}
print(plots[["domestication"]])


# End ==========================================================================
# stopCluster(cluster)
print("Complete!")
