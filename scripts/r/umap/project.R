library(umap)
library(reticulate)
use_condaenv("umap")


fit_transform <- function(dataset, cfg) {
  obj <- umap(dataset$data, config = cfg, method = "umap-learn")
  
  umap.data <-
    data.frame(obj$layout) %>% rownames_to_column("sample_id")
  
  metadata <- dataset$metadata %>% rownames_to_column("sample_id")
  
  df <- inner_join(umap.data, metadata, by = "sample_id") %>%
    select(X1, X2, genotype, chemotype2, domestication, dom2) %>%
    rename(`UMAP 1` = X1,
           `UMAP 2` = X2, )
  df
}

palette_domestication <-
  c("#0EA5E9", "#ff6682", "#6159a4", "#feae01", "#a6a6a6")


palette_genotype <-
  c(
    "#A6CEE3",
    "#85B8D7",
    "#64A3CC",
    "#438EC0",
    "#2279B5",
    "#3F8EAA",
    "#63A8A0",
    "#87C196",
    "#ABDA8B",
    "#98D277",
    "#79C360",
    "#5AB349",
    "#3BA432",
    "#569E3F",
    "#879D5A",
    "#B89B74",
    "#E99A8F",
    "#F78685",
    "#F16667",
    "#EB4748",
    "#E52829",
    "#E62F27",
    "#EC583B",
    "#F3804F",
    "#F9A963",
    "#FDB762",
    "#FDA847",
    "#FE982C",
    "#FE8811",
    "#FA8313",
    "#ED9047",
    "#E09C7B",
    "#D3A8AF",
    "#C3AAD2",
    "#AC8DC3",
    "#9471B4",
    "#7D54A5",
    "#704599",
    "#957599",
    "#B9A499",
    "#DDD399",
    "#FDFB96",
    "#EAD27A",
    "#D7AA5F",
    "#C48143",
    "#B15928"
  )