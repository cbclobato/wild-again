# Silhoute plot ggplot 
# https://rpkgs.datanovia.com/factoextra/reference/fviz_silhouette.html
library(cluster)

read.csv("outputs/umap/")
embeds <- 

s <- silhouette(labels, dist(embeds, method="euclidean"), full=TRUE)

