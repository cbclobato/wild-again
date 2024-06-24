library(phyloseq)

plot_composition_v2 <- function (x, sample.sort = NULL, otu.sort = NULL, x.label = "sample", 
          plot.type = "barplot", verbose = FALSE, average_by = NULL, 
          group_by = NULL, ...) 
{
  Sample <- Abundance <- Taxon <- Group <- Tax <- horiz <- value <- scales <- ID <- meta <- OTU <- taxic <- otu.df <- taxmat <- new.tax <- NULL
  if (!is.null(x@phy_tree)) {
    x@phy_tree <- NULL
  }
  xorig <- x
  if (verbose) {
    message("Pick the abundance matrix taxa x samples")
  }
  abu <- abundances(x)
  if (verbose) {
    message("Average the samples by group")
  }
  group <- NULL
  if (!is.null(average_by)) {
    dff <- as.data.frame(t(abu))
    dff$group <- sample_data(x)[[average_by]]
    if (is.numeric(dff$group) || is.character(dff$group)) {
      dff$group <- factor(dff$group, levels = sort(unique(dff$group)))
    }
    dff <- dff %>% filter(!is.na(group))
    dff$group <- droplevels(dff$group)
    av <- aggregate(. ~ group, data = dff, mean)
    rownames(av) <- as.character(av$group)
    av$group <- NULL
    abu <- t(av)
  }
  if (verbose) {
    message("Sort samples")
  }
  if (is.null(sample.sort) || sample.sort == "none" || !is.null(average_by)) {
    sample.sort <- colnames(abu)
  }
  else if (length(sample.sort) == 1 && sample.sort %in% taxa(xorig)) {
    tax <- sample.sort
    sample.sort <- rev(sample_names(x)[order(abundances(x)[tax, 
    ])])
  }
  else if (length(sample.sort) == 1 && sample.sort %in% names(sample_data(x)) && 
           is.null(average_by)) {
    sample.sort <- rownames(sample_data(x))[order(sample_data(x)[[sample.sort]])]
  }
  else if (all(sample.sort %in% sample_names(x)) & is.null(average_by)) {
    sample.sort <- sample.sort
  }
  else if (length(sample.sort) == 1 && sample.sort == "neatmap") {
    sample.sort <- neatsort(x, method = "NMDS", distance = "bray", 
                            target = "sites", first = NULL)
  }
  else if (is.vector(sample.sort) && length(sample.sort) > 
           1) {
    sample.sort <- sample_names(x)[sample.sort]
  }
  else if (!sample.sort %in% names(sample_data(x))) {
    warning(paste("The sample.sort argument", sample.sort, 
                  "is not included in sample_data(x). \n            Using original sample ordering."))
    sample.sort <- sample_names(x)
  }
  if (is.null(otu.sort)) {
    otu.sort <- taxa(x)
  }
  else if (length(otu.sort) == 1 && otu.sort == "none")  {
    otu.sort <- taxa(x)
  }
  else if (length(otu.sort) == 1 && otu.sort == "abundance2") {
    otu.sort <- rev(c(rev(names(sort(rowSums(abu)))[seq(1, 
                                                        nrow(abu), 2)]), names(sort(rowSums(abu)))[seq(2, 
                                                                                                       nrow(abu), 2)]))
  }
  else if (length(otu.sort) == 1 && otu.sort == "abundance") {
    otu.sort <- rev(names(sort(rowSums(abu))))
  }
  else if (length(otu.sort) == 1 && otu.sort %in% colnames(tax_table(x))) {
    otu.sort <- rownames(sample_data(x))[order(tax_table(x)[[otu.sort]])]
  }
  else if (all(otu.sort %in% taxa(x))) {
    otu.sort <- otu.sort
  }
  else if (length(otu.sort) == 1 && otu.sort == "neatmap") {
    otu.sort <- neatsort(x, method = "NMDS", distance = "bray", 
                         target = "species", first = NULL)
  }
  dfm <- psmelt(otu_table(abu, taxa_are_rows = TRUE))
  names(dfm) <- c("Tax", "Sample", "Abundance")
  dfm$Sample <- factor(dfm$Sample, levels = sample.sort)
  dfm$Tax <- factor(dfm$Tax, levels = otu.sort)
  if (!is.null(group_by)) {
    if (!is.null(average_by)) {
      dfm$Group <- meta(x)[[group_by]][match(as.character(dfm$Sample), 
                                             meta(x)[[average_by]])]
    }
    else {
      dfm$Group <- meta(x)[[group_by]][match(as.character(dfm$Sample), 
                                             sample_names(x))]
    }
  }
  if (x.label %in% colnames(sample_data(x)) & is.null(average_by)) {
    meta <- sample_data(x)
    dfm$xlabel <- as.vector(unlist(meta[as.character(dfm$Sample), 
                                        x.label]))
    if (is.factor(meta[, x.label])) {
      lev <- levels(meta[, x.label])
    }
    else {
      lev <- unique(as.character(unname(unlist(meta[, 
                                                    x.label]))))
    }
    dfm$xlabel <- factor(dfm$xlabel, levels = lev)
  }
  else {
    dfm$xlabel <- dfm$Sample
  }
  if (verbose) {
    message("Construct the plots")
  }
  if (plot.type == "barplot") {
    p <- make_barplot1(dfm, group_by)
  }
  else if (plot.type == "heatmap") {
    p <- make_heatmap1(x, otu.sort, sample.sort, verbose)
  }
  else if (plot.type == "lineplot") {
    p <- make_lineplot1(dfm)
  }
  else {
    stop("plot.type argument not recognized")
  }
  p
}


make_heatmap1 <- function (x, otu.sort, sample.sort, verbose=FALSE) {
  
  if (verbose) { message("Constructing the heatmap.") }
  
  # Taxa x samples otu matrix
  otu <- abundances(x)
  
  # Remove NAs after the transform
  otu <- otu[rowMeans(is.na(otu)) < 1, colMeans(is.na(otu)) < 1]
  
  
  
  
  otu.sort <- otu.sort[otu.sort %in% rownames(otu)]
  sample.sort <- sample.sort[sample.sort %in% colnames(otu)]
  tmp <- melt(otu[otu.sort, sample.sort])
  
  p <- heat(tmp, colnames(tmp)[[1]], colnames(tmp)[[2]],
            colnames(tmp)[[3]])
  
  p
}

make_barplot1 <- function (dfm, group_by) {
  
  Tax <- Sample <- Abundance <- NULL
  
  # Provide barplot
  dfm <- dfm %>% arrange(Tax)  # Show Taxs always in the same order
  dfm$Tax <- factor(dfm$Tax, levels = unique(dfm$Tax))
  
  p <- ggplot(dfm, aes(x=Sample, y=Abundance, fill=Tax)) +
    geom_bar(position="stack", stat="identity") +
    scale_x_discrete(labels=dfm$xlabel, breaks=dfm$Sample)
  
  # Name appropriately
  p <- p + labs(y = "Abundance")
  
  # Rotate horizontal axis labels, and adjust
  p <- p + theme(axis.text.x=element_text(angle=90, vjust=0.5,
                                          hjust=0))
  p <- p + guides(fill=guide_legend(reverse=FALSE))
  
  if (!is.null(group_by)) {
    p <- p + facet_grid(.~Group, drop = TRUE,
                        space = "free", scales = "free") 
    
  }
  p
}

make_lineplot1 <- function (dfm) {
  
  Tax <- Sample <- Abundance <- NULL
  
  dfm <- dfm %>% arrange(Tax)  # Show Taxs always in the same order
  p <- ggplot(dfm, aes(x=Sample, y=Abundance, color=Tax, group = Tax)) +
    geom_point() +
    geom_line() +
    scale_x_discrete(labels=dfm$xlabel, breaks=dfm$Sample) +
    labs(y = "Abundance") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0)) +
    guides(fill=guide_legend(reverse=FALSE))
  
  p
  
}
