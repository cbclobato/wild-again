## OLD 


```{r}
set.seed(1234)

resample <- group_initial_split(dataset, group = cultivar)
unique(training(resample)$cultivar)
unique(testing(resample)$cultivar)


rsample2caret(training(resample))

groups <- d$cultivar
d <- select(d, -cultivar)
folds <- groupKFold(groups, k = 2)

# Show groups
lapply(folds, function(x, y) table(y[x]), y = groups)

# Create train and test datasets
train_data <- d %>% slice(folds$Fold1)
test_data <- d %>% slice(folds$Fold2)

# <- preProcess(train_data, method='knnImpute')

# Has overlap?
# any(train_data$cultivar %in% test_data$cultivar)

tr_ctrl <- trainControl(
  method="cv",
  summaryFunction=multiClassSummary()
)

model <- train(
  as.formula(paste(label_name, "~ .")),  # E.g. domestication ~ . ## All columns 
  data = d,
  method = "rf",
  ntree = 1,
  na.action = na.omit
)
predictions <- predict(model, test_data)
model
```

Random Forest 

382 samples
3361 predictors
5 classes: 'Hybrid', 'Hybrid S1', 'Landrace', 'Unstable hybrid', 'Wild/Feral' 

No pre-processing
Resampling: Bootstrapped (25 reps) 
Summary of sample sizes: 308, 308, 308, 308, 308, 308, ... 
Resampling results across tuning parameters:
  
  mtry  Accuracy   Kappa     
2  0.6015188  0.03658599
81  0.6214202  0.33054970
3361  0.7091208  0.48864744

Accuracy was used to select the optimal model using the largest value.
The final value used for the model was mtry = 3361.


```{r}
tr_ctrl <- trainControl(
  method = "cv",
  number = 1,
  repeats = 1,
  index = indices$index,
  indexOut = indices$indexOut,
  savePredictions = "final",
  classProbs = TRUE,
  summaryFunction = multiClassSummary()
)

<- groupKFold(group, k = n_folds)
group_kfold <- createFolds(
  group=train_data$cultivar,
  k = n_folds,
  list = TRUE,
  returnTrain = TRUE
)

tr_ctrl <- trainControl(
  method = "cv",
  number = cv_folds,
  index = group_kfold,
  savePredictions = "final",
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

model <- train(
  y ~ .,
  data = train_data,
  method = "rf",  # e.g., "rf" for random forests, "xgbTree" for XGBoost
  trControl = tr_ctrl,
  metric = "ROC",  # Use AUC-ROC as the performance metric
  preProcess = c("center", "scale") # Pre-process features by centering and scaling
)  

# Make predictions and evaluate the model
predictions <- predict(model, test_data, type = "prob")
roc_obj <- roc(test_data$y, predictions[, 2])
auc <- auc(roc_obj)
print(auc)
```

## Confusion Matrix

```{r}
confusionMatrix(
  reference = test_data[[label_name]],
  data = predictions,
  mode='everything',
  positive='MM'
)
# indices <- rsample2caret(folds)
```



## Feature Importance

```{r}
featurePlot(
  x = trainData[, 1:18],
  y = trainData[[label_name]],
  plot = "box",
  strip = strip.custom(par.strip.text=list(cex=.7)),
  scales = list(
    x = list(relation="free"),
    y = list(relation="free")
  )
)
```


```{r}
featurePlot(
  x = trainData[, 1:18],
  y = trainData[[label_name]],
  plot = "density",
  strip = strip.custom(par.strip.text=list(cex=.7)),
  scales = list(
    x = list(relation="free"),
    y = list(relation="free")
  )
)
```

## with base R
# abundance_table <- as.data.frame(t(as.matrix(otu_table(bh_raref))))
# abundance_table$id <- row.names(abundance_table)
# metadata <- as.data.frame(as.matrix(sample_data(bh_raref)))
# metadata$id <- row.names(metadata)
# ## metadata <- metadata[names(metadata) %in% c("id", label_name)]
# d <- merge(abundance_table, metadata, by="id")
# ... with(res3, sign == "SPECIALIST" | sign == "GENERALIST")
## Cross-validation

```{r}
n_folds <- length(unique(train_data$cultivar)); n_folds
folds <- group_vfold_cv(
  train_data,
  group=cultivar,
  v=n_folds,
  repeats=1,
  pool=0.1
)

model <- randomForest( ~ ., data = iris ,ntree = 100, proximity = TRUE)

tune_spec <- rand_forest(
  formula=form,  
  data = train_data, 
  mtry = tune(), 
  trees = tune(), 
  mode = "classification"
)

wf <- workflow() %>%
  add_model(tune_spec) %>%
  add_formula(f)
