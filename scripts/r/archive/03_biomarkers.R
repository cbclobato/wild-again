# library(doParallel)
library(phyloseq)
library(tidyverse)
library(tidymodels)
library(tidytext)
library(colorspace)
library(ggtext)
library(glue)
library(cowplot)
library(ggbeeswarm)
library(ggtext)
library(glue)
library(patchwork)

source("scripts/r/archive/03_sv-importance-patch.R")

# cluster <- makeCluster(4)
# registerDoParallel(cluster)

filtered <- TRUE
normalized <- TRUE
input_path <- "data/biomarkers"
output_path <- "outputs/biomarkers"
dir.create(output_path, FALSE, TRUE)


dataset_names <- c("domestication")


eval_model <- F

meanAbsPlots <- list()
collectedMetrics <- list()
for (dataset_name in dataset_names) {
# foreach(dataset_name=datasets) %dopar% {
  
  print(paste("Dataset:", dataset_name))
  prefix <- "all-"
  if (filtered) {
    prefix <- "" 
  }
  
  suffix <- ""
  if (normalized) {
    suffix <- "-norm" 
  }
  
  prefix_dataset_suffix <- paste0(prefix, dataset_name, suffix)
  input_file_path <- 
    file.path(input_path, paste0(prefix_dataset_suffix, ".RData"))
  load(input_file_path)
  
  
  label_sym <- sym(label_name)
  if (length(categories) == 2) {
    .pred_categories <- paste0(".pred_", categories[1])
  } else {
    .pred_categories <- paste0(".pred_", categories)  
  }
  
  ## Create Resamples
  print("Setup...")
  set.seed(1234)
  n_folds <- 10
  resamples <-
    vfold_cv(
      dataset,
      v = n_folds,
      repeats = 5,
      strata = !!label_sym,
      pool = 0.1
    )
  
  ## Pre-processing Recipe
  dataset_recipe <-
    recipe(dataset) %>%
    update_role(starts_with("ASV"), new_role = "predictor") %>%
    update_role(!!label_sym, new_role="outcome") %>%
    step_rm(genotype) %>%
    # Remove any zero variance predictors
    step_zv(all_predictors()) %>%
    # Center numeric predictors
    step_center(all_numeric_predictors()) %>%
    prep()
  
  
  ## Specify Metrics
  model_metrics <- metric_set(
    roc_auc, accuracy, kap, f_meas, sens, spec, recall, precision
  )
  
  ## Specify Model
  model_spec <- 
    boost_tree(
      trees = 1000,
      min_n = 8,
      tree_depth = 8,
      learn_rate = 0.02,
      loss_reduction = 1e-8,
      sample_size = 0.8,
    ) %>%
    set_engine("xgboost") %>%
    set_mode("classification")
  
  ## Workflow
  wf <-
    workflow() %>%
    add_recipe(dataset_recipe) %>%
    add_model(model_spec)
  
  if (eval_model) {
    ## Evaluate with Resamples
    print("Evaluate with Resamples...")
    
    resample_results <-
      wf %>%
      fit_resamples(
        resamples,
        control = control_resamples(save_pred = TRUE, save_workflow=TRUE),
        metrics = model_metrics
      )
    
    file_path <- file.path(
      output_path,
      paste0(prefix_dataset_suffix, "-resample-results.RData")
    )
    save(resample_results, file=file_path)
    
    ## Summarize Metrics
    collectedMetrics[[dataset_name]] <- collect_metrics(resample_results, summarize = TRUE)
    
    # ROC curve
    roc_resampled <-
      resample_results %>%
      unnest(.predictions) %>%
      roc_curve(
        truth = !!label_sym, !!!.pred_categories
      )
    p <- autoplot(roc_resampled)
    file_path <- file.path(
      output_path,
      paste0(prefix_dataset_suffix, "-resample-results.png")
    )
    ggsave(file_path, plot = p, bg = "white")
  }
  
  
  ## Train Model
  print("Train Model...")
  
  # Fit Model on all data
  shap_dataset <- dataset
  model_fit <- wf %>% fit(data = shap_dataset)
  
  eval_data <- shap_dataset
  new_data <- bake(
    prep(dataset_recipe),
    has_role("predictor"),
    new_data=eval_data,
  )
  
  data_mat <- bake(
    prep(dataset_recipe),
    has_role("predictor"),
    new_data=eval_data,
    composition="matrix"
  )
  
  fit_engine <- extract_fit_engine(model_fit)
  
  ## Plot Feature Importance per Class
  print("Feature Importance and Biomarker Selection...")
  
  # load("data/spec-gen-1000-cultivar.RData")
  # load("data/core/gen-core-rel.RData")
  # load("outputs/r/305_core_flex/gen-core.RData")
  # load("outputs/r/305_core_flex/gen-core-rel.RData")
  load("outputs/r/archive/305_core_flex/gen-core-rel.RData")
  core.list <- row.names(gen.core.rel)
  # ecotype <- res3  %>% rownames_to_column("ASV")
  
  ## Parameters
  dataset <- eval_data
  # obj <- shap
  top_n <- 20
  bar_width <- 2/3
  fill <- "#fca50a"
  
  ## Create explanations for all classes
  label_levels <- levels(dataset[[label_name]])
  n_levels <- length(label_levels)
  
  ## >>> Only Selection:
  idx <- 1
  matrix_list <- list()
  df_list <- list()
  obj <- shapviz::shapviz(fit_engine, data_mat)
  beeswarm <- sv_importance_v2(obj, kind = "beeswarm", max_display=20)
  
  file_path <- file.path(output_path, paste0(dataset_name, "-selection-beeswarm.png"))
  ggsave(
    file_path, plot=beeswarm, bg = "white", width =  2400, height = 1000, units = "px"
  )
  
  S <- shapviz::get_shap_values(obj)
  df_list[[1]] <- base::transform(
    as.data.frame.table(S, responseName = "SHAP"),
    feature = factor(Var2),
    label = label_levels[idx]
  )
  s.df <- bind_rows(df_list)
  levels(s.df$label) <- c("Selection")
  
  ## Plot
  label_color <- "#769C64"  # "#4A623E", #94C47D"
  group.imp <- 
    s.df %>%
    group_by(label, feature) %>%
    summarise(
      SHAP_n=n(),
      SHAP_mean=mean(SHAP),
      SHAP_sd=sd(SHAP),
      SHAP_se=sd(SHAP) / sqrt(SHAP_n),
      absSHAP_mean=mean(abs(SHAP)),
      absSHAP_sd=sd(abs(SHAP)),
      absSHAP_se=sd(abs(SHAP)) / sqrt(SHAP_n),
      .groups="drop"
    )  %>%
    mutate(
      is.core = feature %in% core.list,
      color = if_else(is.core, label_color, "black"),
      feature_core = glue("<span style='color:{color}'>{feature}</span>")
    )
  # slice_max(absSHAP_mean, n=10, with_ties=FALSE) %>%
  # mutate(
  # ASV = fct_reorder(ASV, absSHAP_mean),
  # ASV_reorder = reorder_within(ASV, absSHAP_mean, Class)
  # )
  
  top_asvs <-
    group.imp %>% 
    group_by(feature) %>%
    summarize(absSHAP_mean=max(absSHAP_mean)) %>%
    slice_max(absSHAP_mean, n=top_n, with_ties=FALSE) %>%
    pull(feature)
  
  group.imp <-
    group.imp %>%
    filter(feature %in% top_asvs)  %>%
    mutate(feature_core = fct_reorder(feature_core, absSHAP_mean, .fun=max))
  
  # s.df <-
  #   s.df %>%
  #   filter(feature %in% top_asvs)  %>%
  #   mutate(feature = fct_reorder(feature, abs(SHAP), .fun=max))
  
  
  ## Colors
  # palette_name <- "Dark 3"
  palette_name <- "Dynamic"
  # n_colors <- length(unique(with(group.imp, ASV)))
  # many_colors <- qualitative_hcl(n_colors, palette = palette_name)
  n_colors <- length(unique(with(group.imp, label)))
  many_colors <- qualitative_hcl(n_colors, palette = palette_name)
  # many_colors <- divergingx_hcl(n_colors, palette = palette_name)
  
  p1 <- 
    ggplot(
      group.imp,
      aes(
        x = absSHAP_mean, 
        y = feature_core,
        color=label,
      )
    ) +
    geom_pointrange(
      aes(
        # xmin = SHAP_mean + qt(0.025, SHAP_n - 1) * SHAP_se,
        # xmax = SHAP_mean + qt(0.975, SHAP_n - 1) * SHAP_se,
        xmin = absSHAP_mean - 1.96 * absSHAP_se,
        xmax = absSHAP_mean + 1.96 * absSHAP_se,
      ),
      size=0.3,
      position=position_dodge(width=0.5),
      # show.legend = FALSE
    ) +
    scale_colour_manual(name="", values = many_colors) +
    # facet_grid(rows = vars(Class), scales="free") +
    labs(x = "mean(|SHAP value|)", y = element_blank()) +
    theme_minimal() +
    theme(
      axis.text.y = element_markdown(),
    ) +
    scale_y_reordered()
  
  if (dataset_name %in% c("flower", "hemp", "marijuana-domestication")) {
    p1 <- p1 +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.direction = "horizontal"
      ) +
      guides(color = guide_legend(title = "", nrow = 1))
  }
  
  ## Plots
  meanAbsPlots[[dataset_name]] <- p1
  
  # p <- plot_grid(
  #   p1 + guides(color = FALSE),
  #   p2 +
  #     theme(
  #       axis.text.y = element_blank(),
  #       axis.ticks.y = element_blank(),
  #       axis.title.y = element_blank()
  #     ),
  #   labels = "AUTO",
  #   label_size=12
  # ); p
  
  ## Save figures
  print("Saving figures...")
  
  p <- meanAbsPlots[["hemp"]]
  ggsave("outputs/biomarkers/hemp-selection.png", plot=p, bg = "white", width =  2400, height = 1000, units = "px")
  
  p <- meanAbsPlots[["domestication"]]
  ggsave("outputs/biomarkers/domestication-selection.png", plot=p, bg = "white", width =  2400, height = 1000, units = "px")
  
  # file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-abs.png"))
  # ggsave(file_path, plot = p1, bg = "white", width =  2400, height = 1000, units = "px")
  # 
  # file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-abs.svg"))
  # ggsave(file_path, plot = p1, bg = "white", width =  2400, height = 1000, units = "px")
  # 
  # file_path <- file.path(output_path, paste0(prefix_dataset_suffix, ".png"))
  # ggsave(file_path, plot = p2, bg = "white", width =  2400, height = 1000, units = "px")
  # 
  # file_path <- file.path(output_path, paste0(prefix_dataset_suffix, ".svg"))
  # ggsave(file_path, plot = p2, bg = "white", width =  2400, height = 1000, units = "px")
}
stopCluster(cluster)
print("Complete!")

## Chemotype + Marijuana Chemotype

# p_grid <- (
#   meanAbsPlots[["chemotype"]] + meanAbsPlots[["marijuana-chemotype"]]
# ) +
#   plot_annotation(tag_levels = 'A')
# p_grid
# prefix_dataset_suffix <- paste0(prefix, "chemotypes", suffix)
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-abs.png"))
# ggsave(
#   file_path,
#   plot = p_grid,
#   bg = "white",
#   width = 30,  # 30
#   height = 12, # 12
#   units = "cm",
#   scale = 1,
#   dpi = 600
# )
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-abs.svg"))
# ggsave(
#   file_path,
#   plot = p_grid,
#   bg = "white",
#   width = 30,
#   height = 12,
#   units = "cm",
#   scale = 1,
#   dpi = 600
# )
# 
# 
# ## Flower + Domestication
# 
# p_grid <- (
#   (meanAbsPlots[["hemp"]] + meanAbsPlots[["marijuana-domestication"]]) + meanAbsPlots[["flower"]]
# ) + plot_annotation(tag_levels = 'A')
# p_grid
# prefix_dataset_suffix <- paste0(prefix, "flower-domestication", suffix)
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-abs.png"))
# ggsave(
#   file_path,
#   plot = p_grid,
#   bg = "white",
#   width = 48,  # 30
#   height = 12, # 12
#   units = "cm",
#   scale = 1,
#   dpi = 600
# )
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-abs.svg"))
# ggsave(
#   file_path,
#   plot = p_grid,
#   bg = "white",
#   width = 48,
#   height = 12,
#   units = "cm",
#   scale = 1,
#   dpi = 600
# )
# 
# 
# ## Flower + Domestication V2
# 
# p_grid <- (
#   meanAbsPlots[["flower"]] + (meanAbsPlots[["hemp"]] / meanAbsPlots[["marijuana-domestication"]])
# ) + plot_annotation(tag_levels = 'A')
# p_grid
# prefix_dataset_suffix <- paste0(prefix, "flower-domestication", suffix)
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-abs-v2.png"))
# ggsave(
#   file_path,
#   plot = p_grid,
#   bg = "white",
#   width = 40,  # 30
#   height = 14, # 12
#   units = "cm",
#   scale = 1,
#   dpi = 600
# )
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-abs-v2.svg"))
# ggsave(
#   file_path,
#   plot = p_grid,
#   bg = "white",
#   width = 40,
#   height = 14,
#   units = "cm",
#   scale = 1,
#   dpi = 600
# )

## Old
# pair <- c("hemp", "marijuana-domestication")
# p1 <- 
#   plot_grid(
#     plotlist=meanAbsPlots[pair],
#     labels = "AUTO",
#     label_size=12
#   )

# p2 <-
#   plot_grid(
#     plotlist=meanPlots[pair],
#     labels = "AUTO",
#     label_size=12
#   )
# prefix_dataset_suffix <- paste0(prefix, "chemotypes", suffix)
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, ".png"))
# ggsave(
#   file_path,
#   plot = p1,
#   bg = "white",
#   width = 30,
#   height = 10,
#   units = "cm",
#   scale = 1,
#   dpi = 600
# )
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, ".svg"))
# ggsave(
#   file_path,
#   plot = p1,
#   bg = "white",
#   width = 30,
#   height = 10,
#   units = "cm",
#   scale = 1,
#   dpi = 600
# )


## ROC plots

# filtered <- TRUE
# normalized <- TRUE
# output_path <- "outputs/biomarkers"
# 
# 
# dataset_names <- c(
#   "chemotype",
#   "hemp",
#   "flower",
#   "marijuana-chemotype",
#   "marijuana-domestication"
# )
# 
# roc_plots <- list()
# for (dataset_name in dataset_names) {
#   print(paste("Dataset:", dataset_name))
#   prefix <- "all-"
#   if (filtered) {
#     prefix <- ""
#   }
# 
#   suffix <- ""
#   if (normalized) {
#     suffix <- "-norm"
#   }
# 
#   prefix_dataset_suffix <- paste0(prefix, dataset_name, suffix)
#   input_file_path <-
#     file.path(input_path, paste0(prefix_dataset_suffix, ".RData"))
#   load(input_file_path)
# 
# 
#   label_sym <- sym(label_name)
#   if (length(categories) == 2) {
#     .pred_categories <- paste0(".pred_", categories[1])
#   } else {
#     .pred_categories <- paste0(".pred_", categories)
#   }
# 
#   print(paste("Dataset:", dataset_name))
#   prefix <- "all-"
#   if (filtered) {
#     prefix <- ""
#   }
# 
#   suffix <- ""
#   if (normalized) {
#     suffix <- "-norm"
#   }
# 
#   prefix_dataset_suffix <- paste0(prefix, dataset_name, suffix)
# 
#   file_path <- file.path(
#     output_path,
#     paste0(prefix_dataset_suffix, "-resample-results.RData")
#   )
#   load(file_path)
# 
#   print(dataset_name)
#   print(collect_metrics(resample_results, summarize = TRUE)[6,])
# 
#   roc_resampled <-
#     resample_results %>%
#     unnest(.predictions) %>%
#     roc_curve(
#       truth = !!label_sym, !!!.pred_categories
#     )
#   p <- autoplot(roc_resampled)
# 
#   roc_plots[[dataset_name]] <- p
# }
# 
# p <-
#   plot_grid(
#   plotlist=roc_plots,
#   label_size = 12,
#   labels="AUTO"
# )
# 
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-roc.png"))
# ggsave(file_path, plot = p, bg = "white", width =  2400, height = 1000, units = "px")
# 
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-roc.svg"))
# ggsave(file_path, plot = p, bg = "white", width =  2400, height = 1000, units = "px")


obj <- shapviz::shapviz(fit_engine, data_mat, which_class=5)

sv_importance_v2(obj, kind = "beeswarm") + theme_classic()
