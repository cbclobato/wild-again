# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("lefser")

library(phyloseq)
library(lefser)
library(tidytext)

source("scripts/r/biomarkers/data.R")

name <- "domestication"

load("outputs/r/supplementary/bh.RData")

bh.ra <- transform_sample_counts(bh, function(x) x / sum(x))
ra <- unclass(otu_table(bh.ra))
colData <- as(sample_data(bh.ra), "data.frame")


colData$dom2 <- as.factor(colData$dom2)
colData$dom2 <- relevel(colData$dom2, ref="Low")

se <- SummarizedExperiment(
  assays = list(ra=ra), colData = colData
)
# se <- relativeAb(se)

result <- lefser(se, groupCol = "dom2", lda.threshold = 0.001)
# print(result)
# lefserPlot(result)


df_plot <- as.data.frame(t(result$means), check.names = FALSE) %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-starts_with("ASV"), names_to = "Domestication", values_to = "LDA Score") %>%
  group_by(Domestication) %>%
  slice_max(abs(`LDA Score`), n = 20) %>%
  # mutate(ASV_ordered = fct_reorder(ASV, `LDA Score`))
  mutate(ASV_ordered = reorder_within(ASV, `LDA Score`, Domestication))

ggplot(df_plot, aes(x=`LDA Score`, y=ASV_ordered)) +
  geom_bar(stat="identity") + 
  scale_y_reordered() +
  facet_wrap(~ Domestication, scales = "free")


# LDA ==========================================================================
# install.packages("discrim") # for tidymodels parsnip

source("scripts/r/biomarkers/data.R")
load("outputs/r/supplementary/bh.RData")

name <- "domestication"
dataset <- preprocess(bh,
                      name,
                      normalize = TRUE,
                      filter_features = TRUE)
print(dataset)
df <- dataset$data
training <- subset(df, select = -genotype)
result.lda <- MASS::lda(domestication~., training)
plot(result.lda)

library(pROC)

predictions <- predict(result.lda, training, type="prob", method="debiased")
labels <- predictions$class
probs <- predictions$posterior
result.roc <- multiclass.roc(labels, probs)
result.roc


plot(ciobj, type="shape", col="#1c61b6AA") # plot as a blue shape
plot(ci(rocobj, of="thresholds", thresholds="best")) # add one threshold
rocobj <- plot.roc(labels, probs,
                   main="Confidence intervals", percent=TRUE,
                   ci=TRUE, # compute AUC (of AUC by default)
                   print.auc=TRUE) # print the AUC (will contain the CI)
ciobj <- ci.se(rocobj, # CI of sensitivity
               specificities=seq(0, 100, 5)) # over a select set of specificities
