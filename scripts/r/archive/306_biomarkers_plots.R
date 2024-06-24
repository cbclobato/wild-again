## Plot 1
label_size <- 10
prefix_dataset_suffix <- paste0(prefix, "two", suffix)
# datasets <- c("chemotype", "marijuana-chemotype",
#               "hemp", "marijuana-domestication")
p1 <-
  plot_grid(
    meanAbsPlots[["chemotype"]],
    meanAbsPlots[["marijuana-chemotype"]],
    labels = "AUTO",
    ncol = 2,
    label_size=label_size
  )
file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-abs.png"))
ggsave(
  file_path,
  plot = p1,
  bg = "white",
  width = 30,
  height = 10,
  units = "cm",
  scale = 1,
  dpi = 600
)
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-abs.svg"))
# ggsave(
#   file_path,
#   plot = p1,
#   bg = "white",
#   width = 30,
#   height = 15,
#   units = "cm",
#   scale = 1,
#   dpi = 600
# )

# p2 <-
#   plot_grid(
#     plotlist=meanPlots[datasets],
#     labels = "AUTO",
#     label_size=label_size
#   )
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, ".png"))
# ggsave(
#   file_path,
#   plot = p2,
#   bg = "white",
#   width = 30,
#   height = 15,
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
#   height = 15,
#   units = "cm",
#   scale = 1,
#   dpi = 600
# )

## Plot 2: Flower
prefix_dataset_suffix <- paste0(prefix, "three", suffix)
# p1 <-
#   plot_grid(
#     plotlist=meanAbsPlots[datasets],
#     labels = "AUTO",
#     label_size=12
#   )

left_col <-
  plot_grid(
    meanAbsPlots[["hemp"]] + theme(axis.title.x = element_blank()),
    meanAbsPlots[["marijuana-domestication"]],
    labels = c("A", "B"),
    label_size=label_size,
    ncol=1
  )
p2 <- plot_grid(
  left_col,
  meanAbsPlots[["flower"]],  
  labels = c('', 'C'),
  ncol=2,
  label_size=label_size,
  rel_widths = c(1.2, 1)
)

file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-abs.png"))
ggsave(
  file_path,
  plot = p2,
  bg = "white",
  width = 30,
  height = 10,
  units = "cm",
  scale = 1,
  dpi = 600
)
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, "-abs.svg"))
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

# p2 <- meanPlots[datasets]
# file_path <- file.path(output_path, paste0(prefix_dataset_suffix, ".png"))
# ggsave(
#   file_path,
#   plot = p1,
#   bg = "white",
#   width = 20,
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
#   width = 20,
#   height = 10,
#   units = "cm",
#   scale = 1,
#   dpi = 600
# )