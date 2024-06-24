## Chemotype

> Chemotype across dioecious hybrid

```{r}
library(phyloseq)
load("data/phyloseq/bh.RData")
dataset_name <- "chemotype"
label_name <- "chemotype2"  # High, Low
categories <- c("High", "Low")
bh_subset <-
  subset_samples(
    bh,
    !is.na(chemotype)
    & inflorescence == "Dioecious"
    & domestication == "Hybrid"
  )
dim(sample_data(bh_subset))
```


## Hemp Domestication

> Domestication across dioecious hemp


```{r}
library(phyloseq)
load("data/phyloseq/bh.RData")
dataset_name <- "hemp"
label_name <- "domestication"  # Wild/Feral, Landrace, Hybrid
categories <- c("Hybrid", "Landrace", "Wild/Feral")
bh_subset <-
  subset_samples(
    bh,
    !is.na(domestication)
    & inflorescence %in% c("Dioecious", "Subdioecious")
    & chemotype2 == "Low"
  )
```


## Flower

> Flower type across industrial hemp varieties

```{r}
library(phyloseq)
load("data/phyloseq/bh.RData")
dataset_name <- "flower"
label_name <- "inflorescence"  # Monoecious, Dioecious
categories <- c("Monoecious", "Dioecious")
bh_subset <-
  subset_samples(
    bh,
    !is.na(inflorescence)
    & application == "Industrial hemp"
  )
```


## Marijuana Chemotype

```{r}
library(phyloseq)
load("data/phyloseq/bh.RData")
dataset_name <- "marijuana-chemotype"
# label_name <- "domestication"  # Landrace, Unstable Hybrid, Hybrid, Hybrid S1
# categories <- c("Landrace", "Unstable Hybrid", "Hybrid", "Hybrid S1")
label_name <- "chemotype"  # CBD-rich, CBD/THC-rich, THC-rich
categories <- c("THC-rich", "CBD/THC-rich", "CBD-rich")
bh_subset <-
  subset_samples(
    bh,
    !is.na(domestication) 
    & chemotype2 == "High"
  )
```


## Marijuana Domestication

```{r}
library(phyloseq)
load("data/phyloseq/bh.RData")
dataset_name <- "marijuana-domestication"
label_name <- "domestication"
categories <- c("Hybrid", "Landrace", "Unstable Hybrid", "Hybrid S1")
bh_subset <-
  subset_samples(
    bh,
    !is.na(domestication) 
    & chemotype2 == "High"
  )
```


```{r}
load("data/spec-gen-1000-cultivar.RData")  # res3

## Ecotype
ecotype <- 
  res3 %>%
  rownames_to_column("ASV")

## ASVs
asv_rarefied <- ecotype %>% arrange(ASV) %>% pull(ASV) 
significant <- 
  ecotype %>% 
  filter(sign == "SPECIALIST" | sign == "GENERALIST") %>%
  arrange(ASV) %>%
  pull(ASV)

if (only_significant) {
  prefix <- ""
  ecotype <- filter(ecotype, sign == "SPECIALIST" | sign == "GENERALIST")
}

```




abundance_table <- 
  bh_subset %>%
  otu_table %>%
  data.frame %>%
  t %>%
  as.data.frame %>%
  rownames_to_column("sample_id")


# if (with_feature_selection) {
#   significant <- pull(significant_table, ASV)
#   abundance_table <- select(abundance_table, id, all_of(significant))
# }


# df <- data.frame(filter_features = NA, normalize = NA, results = NA)
# new_row <- c(
#   filter_features = filter_features,
#   normalize = normalize,
#   result = result
# )
# df <- rbind(df, new_row)


set.seed(1234)

splits <- 
  group_initial_split(
    dataset,
    group = cultivar,
    prop = 3/4,
    strata = !!label_sym
  )


# dataset <- select(dataset, -!!sym(label_name))
# dataset_name <- "cultivar"
# label_name <- "cultivar"
# categories <- levels(dataset$cultivar)
# pos.outcome <- categories[1]


## Check

for (dataset_name in datasets) {
  print(dataset_name)
  load(paste0("data/feature-importance/", dataset_name, ".RData"))
  label_sym <- sym(label_name)
  .pred_categories <- paste0(".pred_", categories)
  set.seed(1234)
  #print(table(pull(dataset, !!label_sym), pull(dataset, cultivar)))
  splits <- 
    initial_split(
      dataset,
      prop = 3/4,
      strata = !!label_sym
    )
  train_data <- training(splits)
  test_data <- testing(splits)
  
  n_folds <- 10
  resamples <-
    vfold_cv(
      train_data,
      v = n_folds,
      repeats = 5,
      strata = !!label_sym,
      pool = 0.1
    )
  # print(table(pull(test_data, !!label_sym)))
}


set.seed(1234)

# splits <- 
#   group_initial_split(
#   dataset,
#   group = cultivar,
#   prop = 4/5,
#   strata = !!label_sym
# )
# train_data <- training(splits)
# test_data <- testing(splits)



n_times <- 25
if (dataset_name == "marijuana-chemotype") {
  resamples <-
    group_bootstraps(
      dataset,
      group = cultivar,
      times = n_times,
    )
} else {
  resamples <-
    group_bootstraps(
      dataset,
      group = cultivar,
      times = n_times,
      strata = !!label_sym  # marijuana chemotype
    )
}

## HPARAMS models

### Chemotype

```{r}
# Gradient Boosted Trees
model_spec <- 
  boost_tree(
    trees = 122,
    min_n = 3,
    tree_depth = 8,
    learn_rate = 0.0658,
    loss_reduction = 1.6e-1,
    sample_size = 0.34,
  ) %>%
  set_engine("xgboost") %>%
  set_mode("classification")
```


### Hemp

```{r}
# Gradient Boosted Trees
model_spec <- 
  boost_tree(
    trees = 31,
    min_n = 5,
    tree_depth = 1,
    learn_rate = 0.08,
    loss_reduction = 4e-10,
    sample_size = 0.8,
  ) %>%
  set_engine("xgboost") %>%
  set_mode("classification")
```


### Flower

```{r}
# Gradient Boosted Trees
model_spec <- 
  boost_tree(
    trees = 1703,
    min_n = 16,
    tree_depth = 12,
    learn_rate = 0.0028,
    loss_reduction = 1.14e-10,
    sample_size = 0.92,
  ) %>%
  set_engine("xgboost") %>%
  set_mode("classification")
```


### Marijuana Chemotype

```{r}
# Gradient Boosted Trees
model_spec <- 
  boost_tree(
    trees = 974,
    min_n = 8,
    tree_depth = 10,
    learn_rate = 0.035,
    loss_reduction = 2e-7,
    sample_size = 0.65,
  ) %>%
  set_engine("xgboost") %>%
  set_mode("classification")
```


### Marijuana Domestication

```{r}
# Gradient Boosted Trees
model_spec <- 
  boost_tree(
    trees = 807,
    min_n = 2,
    tree_depth = 5,
    learn_rate = 0.008,
    loss_reduction = 2e-10,
    sample_size = 0.29,
  ) %>%
  set_engine("xgboost") %>%
  set_mode("classification")
```


```{r}
dir.create("data/umap", FALSE, TRUE)

load("data/phyloseq/bh.RData")


filename <- paste0("data/umap/", prefix, dataset_name, suffix, ".csv")
write.csv(asv_table, filename, row.names = FALSE)

filename <- paste0("data/umap/", dataset_name, "-metadata.csv")
write.csv(metadata, filename, row.names = FALSE)




# Save
filename <- paste0("data/umap/", prefix, dataset_name, normalization, "-t.csv")
write.csv(abundance_transposed, filename, row.names = FALSE)

filename <- paste0("data/umap/", prefix, "ecotype.csv")
write.csv(filtered_ecotype, filename, row.names = FALSE)
```



## FlashWeave only One-Hot encodes a metadata variable
##  if it has 3 or more categories
n_groups <- length(cultivars)
if (n_groups == 2) {
  # Create a recipe with dummy variables
  recipe_obj <- recipes::recipe( ~ ., data = metadata_table) %>%
    step_dummy(all_nominal(), -all_outcomes(), one_hot = FALSE)
  
  # Prepare the data using the recipe
  metadata_table <- 
    prep(recipe_obj, training = metadata_table) %>%
    bake(new_data = metadata_table)
}




#  shap

matrix_list <- list()
for (idx in 1:2) {
  obj <- shapviz(fit_engine, data_mat, which_class=idx)
  matrix_list[[idx]] <- 
    as.data.frame(obj[["S"]]) %>%
    mutate(Class = label_levels[idx])
}
s.df <- bind_rows(matrix_list)
top_n <- 10
#obj <- shapviz(fit_engine, data_mat, which_class=idx)
#s.df <-
#  (baseline + as.data.frame(obj[["S"]])) %>%
#  mutate(Class = label_levels[1])

# s.df <- 
#   s.df %>%
#   mutate(
#     Class = fct_relabel(Class, ~ paste0(label_levels, collapse = " / "))
#   )



if (length(label_levels) == 2) {
  group.imp <- 
    s.df %>%
    pivot_longer(
      starts_with("ASV", ignore.case=FALSE),
      names_to = "ASV",
      values_to = "SHAP"
    ) %>
    group_by(Class, ASV) %>%
    summarise(
      absSHAP_n=n(),
      SHAP_mean=mean(SHAP),
      absSHAP_mean=mean(abs(SHAP)),
      absSHAP_sd=sd(abs(SHAP)),
      absSHAP_se=sd(abs(SHAP)) / sqrt(absSHAP_n),
      # .groups="drop_last"
      .groups="keep"
    ) %>%
    slice_max(absSHAP_mean, n=top_n, with_ties=FALSE) %>%
    ungroup() %>%
    mutate(
      ASV_reorder = reorder_within(ASV, absSHAP_mean, Class)
    )
} else {
  
}
  
  
  # group.imp <-
  #   group.imp %>%
  #   inner_join(ecotype, by = "ASV") %>%
  #   mutate(
  #     ecotype_color = case_when(
  #       sign == "GENERALIST" ~ "#ffa600",
  #       #sign == "SPECIALIST" ~ "lightblue3",
  #       sign == "SPECIALIST" ~ "#00A9FF",
  #       .default = "gray40"
  #     ),
  #     ASV = glue("<b style='color:{ecotype_color}'>{ASV}</b>"),
  #     ASV_reorder = reorder_within(ASV, absSHAP_mean, Class)
  #   )