library(cluster)
library(ggplot2)

palette_domestication <-
  c("#0EA5E9", "#ff6682", "#6159a4", "#feae01", "#a6a6a6")


silhouette_plot <- function(embeds, labels) {
  print(levels(labels))
  labels_int <- as.integer(labels)
  sil.obj <-
    silhouette(labels_int, dist(embeds, method = "euclidean"), full = TRUE)

  
  df <- as.data.frame(sil.obj[, 1:3], stringsAsFactors = TRUE) %>%
    rownames_to_column("row_names")
  labels_df <- as.data.frame(labels) %>% 
    rownames_to_column("row_names")
  df <- merge(df, labels_df, by = "row_names")
  df <- df[order(df$cluster, -df$sil_width), ]

  # df <- df[bar_order, ]
  # df$row_names <- factor(df$row_names, levels=bar_order)
  p <-
    ggplot(df,
           aes(
             x = fct_reorder(row_names, order(df$cluster, -df$sil_width)),
             y = sil_width,
             color = labels,
             fill = labels
           )) +
    geom_bar(stat = "identity") +
    geom_hline(
      yintercept = mean(df$sil_width),
      linetype = "dashed",
      color = "#ce232a"
    ) +
    scale_fill_discrete(type = palette_domestication) +
    scale_color_discrete(type = palette_domestication) +
    ylim(c(NA, 1)) +
    labs(
      y = "Si", #Silhouette width 
      x = "",
      color = "Domestication",
      fill = "Domestication",
      # title = paste0(
      #   "Clusters silhouette plot ",
      #   "\nAverage silhouette width: ",
        round(mean(df$sil_width),
              2)
      ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.text.y = (element_text(angle = 0,
                                  hjust = 0.5,
                                  size = 12)),
      axis.title.y = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.position = "none") # c(0.05, 0.9)
  
  list(df = df, p = p)
}
