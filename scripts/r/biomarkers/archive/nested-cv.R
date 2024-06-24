library(phyloseq)
library(tidyverse)
library(tidymodels)


evaluation <- function(dataset, seed = 1234) {
  # Model  and evaluation
  
  name <- dataset[["name"]]
  data <- dataset[["data"]]
  label_name <- dataset[["label_name"]]
  categories <- dataset[["categories"]]
  filtered <- dataset[["filtered"]]
  normalized <- dataset[["normalized"]]
  
  label_sym <- sym(label_name)
  if (length(categories) == 2) {
    .pred_categories <- paste0(".pred_", categories[1])
  } else {
    .pred_categories <- paste0(".pred_", categories)
  }
  
  # Create Resamples
  print("Setup...")
  set.seed(seed)
  
  # Pre-processing Recipe
  data_recipe <-
    recipe(data) %>%
    update_role(starts_with("ASV"), new_role = "predictor") %>%
    update_role(!!label_sym, new_role = "outcome") %>%
    step_rm(genotype) %>%
    # Remove any zero variance predictors
    step_zv(all_predictors()) %>%
    # Center numeric predictors
    step_center(all_numeric_predictors()) %>%
    prep()
  
  
  # Specify Metrics
  model_metrics <-
    metric_set(roc_auc, accuracy, kap, f_meas, sens, spec, recall, precision)
  
  # Specify Model
  model_spec <-
    boost_tree(
      trees = 5000,
      min_n = 3,
      tree_depth = 8,
      learn_rate = 0.02,
      loss_reduction = 1e-8,
      sample_size = 0.8,
    ) %>%
    set_engine("xgboost") %>%
    set_mode("classification")
  
  # Workflow
  wf <-
    workflow() %>%
    add_recipe(data_recipe) %>%
    add_model(model_spec)
  
  eval_plot <- NULL
  eval_results <- NULL
  if (eval_model) {
    print("Evaluate with Resamples...")
    set.seed(seed)
    
    # Nested CV
    n_folds <- 10
    resamples <- nested_cv(
      train_dat,
      outside = vfold_cv(
        v = n_folds,
        repeats = 1,
        strata = !!label_sym,
        pool = 0.1
      ),
      inside = bootstraps(times = 25)
    )
    
    
    # Evaluation
    resample_results <-
      wf %>%
      fit_resamples(
        resamples,
        control = control_resamples(save_pred = TRUE, save_workflow = TRUE),
        metrics = model_metrics
      )
    
    # Save results RData
    prefix_name_suffix <- get_basename(name, filtered, normalized)
    file_path <- file.path(output_path,
                           paste0(prefix_name_suffix, "-resample-results.RData"))
    save(resample_results, file = file_path)
    
    # Summarize Metrics
    resample_results_summary <-
      collect_metrics(resample_results, summarize = TRUE)
    print(resample_results_summary)
    
    # ROC curve
    roc_resampled <-
      resample_results %>%
      unnest(.predictions) %>%
      roc_curve(truth = !!label_sym,!!!.pred_categories)
    p <- autoplot(roc_resampled)
    eval_plot <- p
    
    # Save figure
    basename <- get_basename(name, filtered, normalized)
    file_path <- file.path(output_path,
                           paste0(basename, "-resample-results.png"))
    ggsave(file_path, plot = p, bg = "white")
  }
  
  # Fit model
  print("Train Model...")
  
  set.seed(seed)
  model_fit <- wf %>% fit(data = data)
  
  new_data <- bake(prep(data_recipe),
                   has_role("predictor"),
                   new_data = data)
  
  data_mat <- bake(
    prep(data_recipe),
    has_role("predictor"),
    new_data = data,
    composition = "matrix"
  )
  
  fit_engine <- extract_fit_engine(model_fit)
  
  list(
    data_mat = data_mat,
    fit_engine = fit_engine,
    eval_results = eval_results,
    eval_plot = eval_plot
  )
}


get_basename <- function(name, filtered, normalized) {
  prefix <- "all-"
  if (filtered) {
    prefix <- ""
  }
  
  suffix <- ""
  if (normalized) {
    suffix <- "-norm"
  }
  
  prefix_name_suffix <- paste0(prefix, name, suffix)
  prefix_name_suffix
}


load_preprocessed_data <- function(name, filtered, normalized) {
  # Load preprocessed data
  # Returns:
  #   list: name, data, categories, label_name, normalized, filtered
  basename <- get_basename(name, filtered, normalized)
  load(file.path(input_path, paste0(basename, ".RData")),
       envir = tmp.env <- new.env())
  dataset <- as.list(tmp.env)
  dataset
}


