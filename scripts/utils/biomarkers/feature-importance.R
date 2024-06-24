library(tidyverse)
library(tidytext)


source("scripts/r/biomarkers/sv-importance-patch.R")

# xgboost::xgb.importance(model=result[["xgboost"]]$fit_engine) # Gain
feature.importance <-
  function(lda_fit_engine,
           xgb_fit_engine,
           data_mat,
           categories,
           top_n = 15) {
    n_categories <- length(categories)
    multiclass <- n_categories > 2
    
    # TODO: Load core per class
    # load("outputs/r/archive/305_core_flex/gen-core-rel.RData")
    # core.list <- row.names(gen.core.rel)
    
    lda_means <-
      as.data.frame(t(lda_fit_engine$means), check.names = FALSE)
    
    
    xgb_gain <- xgboost::xgb.importance(model = xgb_fit_engine) %>%
      dplyr::select(Feature, Gain) %>%
      rename(ASV = Feature) %>%
      mutate(
        `Selection` = Gain,
        `Landrace` = Gain,
        `Cross hybrid` = Gain,
        `Inbred line` = Gain
      ) %>%
      dplyr::select(-Gain)
    
    sv <- shapviz::shapviz(xgb_fit_engine, data_mat)
    
    abs_shap <-
      as.data.frame(shapviz::sv_importance(sv, kind = "no")) %>%
      rownames_to_column("ASV")
    
    f_shap <- as.data.frame(sv_importance_v2(sv, kind = "no")) %>%
      rownames_to_column("ASV")
    
    if (multiclass) {
      replace_vec <- character()
      for (i in  seq_len(n_categories)) {
        replace_vec[[categories[i]]] <- paste("Class", i, sep = "_")
      }
      
      # From https://github.com/xia-lab/MicrobiomeAnalystR/blob/master/R/general_anal.R
      # lines 738:740
      class_no <- length(categories)
      lda_means$max <- apply(lda_means[, 1:class_no], 1, max)
      
      lda_means$min <- apply(lda_means[, 1:class_no], 1, min)
      
      lda_means <- lda_means %>%
        rownames_to_column("ASV") %>%
        mutate(ASV = str_replace_all(ASV, "`", ""))
      
      lda_means$LDAScore <-
        signif(log10(1 + abs(lda_means$max - lda_means$min) / 2), digits = 3)
      
      lda_means <- lda_means %>% 
        dplyr::select(ASV, LDAScore) %>%
        mutate(
          `Selection` = LDAScore,
          `Landrace` = LDAScore,
          `Cross hybrid` = LDAScore,
          `Inbred line` = LDAScore
        ) %>%
        dplyr::select(-LDAScore)
        
      
      lda_means <- lda_means %>%
        pivot_longer(
          # -starts_with("ASV"),
          cols = all_of(categories),
          names_to = "Domestication",
          values_to = "Importance"
        ) %>%
        mutate(Kind = "LDA Score")
      
      
      f_shap <- f_shap %>%
        rename(all_of(replace_vec)) %>%
        pivot_longer(
          cols = all_of(categories),
          names_to = "Domestication",
          values_to = "Importance"
        ) %>%
        mutate(Kind = "Score")
      
      abs_shap <- abs_shap %>%
        rename(all_of(replace_vec)) %>%
        pivot_longer(
          cols = all_of(categories),
          names_to = "Domestication",
          values_to = "Importance"
        ) %>%
        mutate(Kind = "Mean(|SHAP|)")
      
      
      xgb_gain <- xgb_gain %>%
        pivot_longer(
          cols = all_of(categories),
          names_to = "Domestication",
          values_to = "Importance"
        ) %>%
        mutate(Kind = "Gain")
      
      df <- rbind(lda_means, abs_shap, f_shap, xgb_gain)
      scores <- unique(df$Kind)
      
      df_plot <- df %>%
        mutate(Domestication = factor(Domestication, levels = categories)) %>%
        pivot_wider(
          id_cols = c("ASV", "Domestication"),
          names_from = "Kind",
          values_from = "Importance"
        ) %>%
        group_by(Domestication) %>%
        slice_max(order_by = Score,
                  n = top_n,
                  with_ties = FALSE) %>%
        mutate(ASV_ordered = reorder_within(ASV, Score, Domestication)) %>%
        pivot_longer(cols = all_of(scores),
                     names_to = "Kind",
                     values_to = "Importance")
      
    } else {
      lda_means <- lda_means %>%
        rownames_to_column("ASV") %>%
        mutate(ASV = str_replace_all(ASV, "`", ""))
      # TODO: ...
      
      df_plot <-
        inner_join(f_shap , abs_shap, by = "ASV") %>%
        rename(Score = imp, `Mean(|SHAP|)` = imp_abs) %>%
        slice_max(order_by = Score,
                  n = top_n,
                  with_ties = FALSE) %>%
        mutate(ASV_ordered = fct_reorder(ASV, Score)) %>%
        pivot_longer(
          cols = c("Score", "Mean(|SHAP|)"),
          names_to = "Kind",
          values_to = "Importance"
        )
    }
    
    df_plot
  }
