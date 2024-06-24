#==============================================================================
# Helper functions to produce the beeswarm plots
# Patch for https://github.com/ModelOriented/shapviz/blob/main/R/sv_importance.R
#==============================================================================

# Example
# ggplot(iris, aes(Species, Sepal.Width)) +
#   geom_point(position = position_bee(), aes(color = Species))

# Beeswarm position
position_bee <- function(width = NULL, adjust = NULL) {
  ggplot2::ggproto(NULL, PositionBee, width = width, adjust = adjust)
}

PositionBee <- ggplot2::ggproto(
  "PositionBee",
  ggplot2::Position,
  required_aes = c("x", "y"),
  
  setup_params = function(self, data) {
    list(
      width = if (!is.null(self$width))
        self$width
      else
        ggplot2::resolution(data$y, zero = FALSE) * 0.4,
      adjust = if (!is.null(self$adjust))
        self$adjust
      else
        0.5
    )
  },
  
  compute_panel = function(self, data, params, scales) {
    data <- ggplot2::flip_data(data, params$flipped_aes)
    y_jit <-
      ave2(data$x,
           g = data$y,
           FUN = shifter,
           adjust = params$adjust)
    data <- ggplot2::transform_position(
      data,
      trans_y = function(y)
        y + y_jit * params$width
    )
    ggplot2::flip_data(data, params$flipped_aes)
  }
)

# Shift values according to their density in the unit interval by quasi-random numbers
shifter <- function(y, ...) {
  if (length(y) == 1L) {
    return(0)
  }
  dens <- stats::density(y, ...)
  dens_y <- dens[["y"]] / max(dens[["y"]])
  shift <-
    halton_sequence(length(y))[rank(y, ties.method = "first")] - 0.5
  2 * shift * stats::approx(dens[["x"]], dens_y, xout = y)[["y"]]
}

# "stats::ave" for grouping variable "g" and additional arguments ...
ave2 <- function(x, g = NULL, FUN = mean, ...) {
  if (is.null(g)) {
    x[] <- FUN(x, ...)
  } else {
    split(x, g) <- lapply(split(x, g), FUN, ...)
  }
  x
}

# First n values of the 1-dimensional Halton sequence (= van der Corput sequence)
# https://en.wikipedia.org/wiki/Halton_sequence
halton_sequence <- function(n, b = 2) {
  vapply(seq_len(n), halton, FUN.VALUE = 0.0)
}

# i-th element of above sequence
halton <- function(i, b = 2) {
  f <- 1
  r <- 0
  while (i > 0) {
    f <- f / b
    r <- r + f * (i %% b)
    i <- trunc(i / b)
  }
  r
}

# Adds non-null titles "nms" to list of ggplots
add_titles <- function(a_list, nms = NULL) {
  if (is.null(nms)) {
    return(a_list)
  }
  mapply(function(p, nm)
    p + ggplot2::ggtitle(nm), a_list, nms, SIMPLIFY = FALSE)
}

#' SHAP Importance Plots
#'
#' This function provides two types of SHAP importance plots: a bar plot
#' and a beeswarm plot (sometimes called "SHAP summary plot").
#' The two types of plots can also be combined.
#'
#' The bar plot shows SHAP feature importances, calculated as the average absolute SHAP
#' value per feature. The beeswarm plot displays SHAP values per feature, using min-max
#' scaled feature values on the color axis. Non-numeric features are transformed
#' to numeric by calling [data.matrix()] first. For both types of plots, the features
#' are sorted in decreasing order of importance.
#'
#' @param object An object of class "(m)shapviz".
#' @param kind Should a "bar" plot (the default), a "beeswarm" plot, or "both" be shown?
#'   Set to "no" in order to suppress plotting. In that case, the sorted
#'   SHAP feature importances of all variables are returned.
#' @param max_display Maximum number of features (with highest importance) to plot.
#'   Set to `Inf` to show all features. Has no effect if `kind = "no"`.
#' @param fill Color used to fill the bars (only used if bars are shown).
#' @param bar_width Relative width of the bars (only used if bars are shown).
#' @param bar_type For "mshapviz" objects with `kind = "bar"`: How should bars be
#'   represented? The default is "dodge" for dodged bars. Other options are "stack",
#'   "wrap", or "separate" (via "patchwork"). Note that "separate" is currently
#'   the only option that supports `show_numbers = TRUE`.
#' @param bee_width Relative width of the beeswarms.
#' @param bee_adjust Relative bandwidth adjustment factor used in
#'   estimating the density of the beeswarms.
#' @param viridis_args List of viridis color scale arguments. The default points to the
#'   global option `shapviz.viridis_args`, which corresponds to
#'   `list(begin = 0.25, end = 0.85, option = "inferno")`. These values are passed to
#'   [ggplot2::scale_color_viridis_c()]. For example, to switch to standard viridis,
#'   either change the default with `options(shapviz.viridis_args = list())` or set
#'   `viridis_args = list()`.
#' @param color_bar_title Title of color bar of the beeswarm plot. Set to `NULL`
#'   to hide the color bar altogether.
#' @param show_numbers Should SHAP feature importances be printed? Default is `FALSE`.
#' @param format_fun Function used to format SHAP feature importances
#'   (only if `show_numbers = TRUE`). To change to scientific notation, use
#'   `function(x) = prettyNum(x, scientific = TRUE)`.
#' @param number_size Text size of the numbers (if `show_numbers = TRUE`).
#' @param ... Arguments passed to [ggplot2::geom_bar()] (if `kind = "bar"`) or to
#'   [ggplot2::geom_point()] otherwise. For instance, passing `alpha = 0.2` will produce
#'   semi-transparent beeswarms, and setting `size = 3` will produce larger dots.
#' @returns
#'   A "ggplot" (or "patchwork") object representing an importance plot, or - if
#'   `kind = "no"` - a named numeric vector of sorted SHAP feature importances
#'   (or a matrix in case of an object of class "mshapviz").
#' @examples
#' X_train <- data.matrix(iris[, -1])
#' dtrain <- xgboost::xgb.DMatrix(X_train, label = iris[, 1], nthread = 1)
#' fit <- xgboost::xgb.train(data = dtrain, nrounds = 10, nthread = 1)
#' x <- shapviz(fit, X_pred = X_train)
#' sv_importance_v2(x)
#' sv_importance_v2(x, kind = "no")
#' sv_importance_v2(x, kind = "beeswarm", show_numbers = TRUE)
#' @seealso \code{\link{sv_interaction}}
#' @export
sv_importance_v2 <- function(object, ...) {
  UseMethod("sv_importance_v2")
}

#' @describeIn sv_importance_v2
#'   Default method.
#' @export
sv_importance_v2.default <- function(object, ...) {
  stop("No default method available.")
}

#' @describeIn sv_importance_v2
#'   SHAP importance plot for an object of class "shapviz".
#' @export
sv_importance_v2.shapviz <-
  function(object,
           kind = c("bar", "beeswarm", "both", "no"),
           max_display = 15L,
           fill = "#fca50a",
           bar_width = 2 / 3,
           bee_width = 0.4,
           bee_adjust = 0.5,
           viridis_args = getOption("shapviz.viridis_args"),
           color_bar_title = "Feature value",
           show_numbers = FALSE,
           format_fun = format_max,
           number_size = 3.2,
           ...) {
    print("shapviz")
    stopifnot("format_fun must be a function" = is.function(format_fun))
    kind <- match.arg(kind)
    
    # imp <- .get_imp(shapviz::get_shap_values(object))
    imp <- .get_imp(shapviz::get_shap_values(object),
                    shapviz::get_feature_values(object))
    
    if (kind == "no") {
      return(imp)
    }
    
    # Deal with too many features
    if (ncol(object) > max_display) {
      imp <- imp[seq_len(max_display)]
    }
    ord <- names(imp)
    object <- object[, ord]  # not required for kind = "bar"
    
    # ggplot will need to work with data.frame
    imp_df <- data.frame(feature = factor(ord, rev(ord)), value = imp)
    is_bar <- kind == "bar"
    if (is_bar) {
      p <- ggplot2::ggplot(imp_df, ggplot2::aes(x = value, y = feature)) +
        ggplot2::geom_bar(fill = fill,
                          width = bar_width,
                          stat = "identity",
                          ...) +
        ggplot2::labs(x = "mean(|SHAP value|)", y = ggplot2::element_blank())
    } else {
      # Prepare data.frame for beeswarm plot
      S <- shapviz::get_shap_values(object)
      
      X <- .scale_X(shapviz::get_feature_values(object))
      df <- transform(
        as.data.frame.table(S, responseName = "value"),
        feature = factor(Var2, levels = rev(ord)),
        color = as.data.frame.table(X)$Freq
      )
      
      p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = feature))
      if (kind == "both") {
        p <- p +
          ggplot2::geom_bar(
            data = imp_df,
            fill = fill,
            width = bar_width,
            stat = "identity"
          )
      }
      p <- p +
        ggplot2::geom_vline(xintercept = 0, color = "darkgray") +
        ggplot2::geom_point(ggplot2::aes(color = color),
                            position = position_bee(width = bee_width, adjust = bee_adjust),
                            ...) +
        .get_color_scale(
          viridis_args = viridis_args,
          bar = !is.null(color_bar_title),
          ncol = length(unique(df$color))   # Special case of constant feature values
        ) +
        ggplot2::labs(x = "SHAP value",
                      y = ggplot2::element_blank(),
                      color = color_bar_title) +
        ggplot2::theme(legend.box.spacing = grid::unit(0, "pt"))
    }
    if (show_numbers) {
      p <- p +
        ggplot2::geom_text(
          data = imp_df,
          ggplot2::aes(
            x = if (is_bar)
              value + max(value) / 60
            else
              min(df$value) - diff(range(df$value)) / 20,
            label = format_fun(value)
          ),
          hjust = !is_bar,
          size = number_size
        ) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.05 + c(0.12 *
                                                                                  !is_bar, 0.09 * is_bar)))
    }
    p
  }

#' @describeIn sv_importance_v2
#'   SHAP importance plot for an object of class "mshapviz".
#' @export
sv_importance_v2.mshapviz <-
  function(object,
           kind = c("bar", "beeswarm", "both", "no"),
           max_display = 15L,
           fill = "#fca50a",
           bar_width = 2 / 3,
           bar_type = c("dodge", "stack", "facets", "separate"),
           bee_width = 0.4,
           bee_adjust = 0.5,
           viridis_args = getOption("shapviz.viridis_args"),
           color_bar_title = "Feature value",
           show_numbers = FALSE,
           format_fun = format_max,
           number_size = 3.2,
           ...) {
    kind <- match.arg(kind)
    bar_type <- match.arg(bar_type)
    
    # All other cases are done via {patchwork}
    if (kind %in% c("bar", "no") && bar_type != "separate") {
      x <- shapviz::get_shap_values(object)
      w <- shapviz::get_feature_values(object)
      imp <- .get_imp(x, w)
      
      if (kind == "no") {
        return(imp)
      }
      if (nrow(imp) > max_display) {
        imp <- imp[seq_len(max_display), , drop = FALSE]
      }
      ord <- rownames(imp)
      imp_df <- data.frame(feature = factor(ord, rev(ord)), utils::stack(as.data.frame(imp)))
      
      if (bar_type %in% c("dodge", "stack")) {
        imp_df <- transform(imp_df, ind = factor(ind, rev(levels(ind))))
        if (is.null(viridis_args)) {
          viridis_args <- list()
        }
        p <-
          ggplot2::ggplot(imp_df, ggplot2::aes(x = values, y = feature)) +
          ggplot2::geom_bar(
            ggplot2::aes(fill = ind),
            width = bar_width,
            stat = "identity",
            position = bar_type,
            ...
          ) +
          ggplot2::labs(fill = ggplot2::element_blank()) +
          do.call(ggplot2::scale_fill_viridis_d, viridis_args) +
          ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))
      } else {
        # facets
        p <-
          ggplot2::ggplot(imp_df, ggplot2::aes(x = values, y = feature)) +
          ggplot2::geom_bar(fill = fill,
                            width = bar_width,
                            stat = "identity",
                            ...) +
          ggplot2::facet_wrap("ind")
      }
      p <- p +
        ggplot2::xlab("mean(|SHAP value|)") +
        ggplot2::ylab(ggplot2::element_blank())
      return(p)
    }
    
    # Now, patchwork
    plot_list <- lapply(
      object,
      FUN = sv_importance_v2,
      # Argument list (simplify via match.call() or some rlang magic?)
      kind = kind,
      max_display = max_display,
      fill = fill,
      bar_width = bar_width,
      bee_width = bee_width,
      bee_adjust = bee_adjust,
      viridis_args = viridis_args,
      color_bar_title = color_bar_title,
      show_numbers = show_numbers,
      format_fun = format_fun,
      number_size = number_size,
      ...
    )
    if (kind == "no") {
      return(plot_list)
    }
    plot_list <-
      add_titles(plot_list, nms = names(object))  # see sv_waterfall()
    patchwork::wrap_plots(plot_list)
  }

# Helper functions
.min_max_scale <- function(z, na.rm = TRUE) {
  r <- range(z, na.rm = na.rm)
  d <- diff(r)
  if (is.na(d) || d == 0) {
    z[!is.na(z)] <- 0.5
    return(z)
  }
  (z - r[1L]) / (r[2L] - r[1L])
}

.standard_scale <- function(z, i, na.rm = TRUE) {
  z_mat <- as.matrix(z)
  (z[, i] - mean(z_mat, na.rm = na.rm)) / sd(z_mat, na.rm = na.rm)
}

# Original
# .get_imp <- function(z) {
#   if (is.matrix(z)) {
#     return(sort(colMeans(abs(z)), decreasing = TRUE))
#   }
#   # list/mshapviz
#   imp <- sapply(z, function(x) colMeans(abs(x)))
#   imp[order(-rowSums(imp)), ]
# }


.get_imp <- function(x, w) {
  # Patch for https://github.com/ModelOriented/shapviz/blob/main/R/sv_importance.R
  
  # Weighted average
  colWMean <- function (x, w) {
    col_names <- colnames(x)
    weighted_averages <- c()
    for (i in 1:ncol(x)) {
      xc <- x
      wc <- w
      # wc[x[,i] < 0,] <- NA
      ws <- .standard_scale(wc, i)
      # wc <- ws
      
      # wc[w[,i] < 0 | x[,i] < 0] <- NA # feature values
      # xc[w[,i] < 0 | x[,i] < 0, i] <- NA # shap values
      # wc[ws < 0] <- NA # feature values
      xc[ws < 0, i] <- NA # shap values
      
      wss <- .min_max_scale(ws)
      # xs <- .min_max_scale(xc[,i])
      xs <- xc[, i]
      # ws <- wc
      weighted_average <-
        sum(xs * wss, na.rm = TRUE) / sum(wss, na.rm = TRUE)
      # weighted_average <- mean(xc[,i], na.rm=TRUE)
      
      weighted_averages <- c(weighted_averages, weighted_average)
    }
    names(weighted_averages) <- col_names
    weighted_averages
  }
  
  if (is.matrix(x)) {
    return(sort(colWMean(x, w), decreasing = TRUE))
  }
  
  # list/mshapviz and if bar
  imp <- mapply(function(x, w)
    colWMean(x, w), x, w)
  imp[order(-rowSums(imp)), ]
}


.scale_X <- function(X) {
  X_scaled <- apply(data.matrix(X), 2L, FUN = .min_max_scale)
  if (nrow(X) == 1L)
    t(X_scaled)
  else
    X_scaled
}

# ncol < 2 treats the special case of constant feature values (e.g., if n = 1)
.get_color_scale <- function(viridis_args,
                             bar = TRUE,
                             ncol = 2L) {
  if (bar) {
    viridis_args_plus <-
      list(
        breaks = if (ncol >= 2L)
          0:1
        else
          0.5,
        labels = if (ncol >= 2L)
          c("Low", "High")
        else
          "Avg",
        guide = ggplot2::guide_colorbar(
          barwidth = 0.4,
          barheight = 8,
          title.theme = ggplot2::element_text(
            angle = 90,
            hjust = 0.5,
            vjust = 0
          ),
          title.position = "left"
        )
      )
  } else {
    viridis_args_plus <- list(guide = "none")
  }
  return(do.call(
    ggplot2::scale_color_viridis_c,
    c(viridis_args, viridis_args_plus)
  ))
}


#' Number Formatter
#'
#' Formats a numeric vector in a way that its largest absolute value determines
#' the number of digits after the decimal separator. This function is helpful in
#' perfectly aligning numbers on plots. Does not use scientific formatting.
#'
#' @param x A numeric vector to be formatted.
#' @param digits Number of significant digits of the largest absolute value.
#' @param ... Further arguments passed to [format()], e.g., `big.mark = "'"`.
#' @returns A character vector of formatted numbers.
#' @export
#' @examples
#' x <- c(100, 1, 0.1)
#' format_max(x)
#'
#' y <- c(100, 1.01)
#' format_max(y)
#' format_max(y, digits = 5)
format_max <- function(x, digits = 4L, ...) {
  mx <- trunc(log10(max(abs(x), na.rm = TRUE))) + 1L
  x_rounded <- round(x, pmax(0L, digits - mx))
  format(x_rounded, scientific = FALSE, trim = TRUE, ...)
}
