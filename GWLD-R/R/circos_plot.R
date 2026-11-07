.circos_numeric <- function(x, name) {
  value <- suppressWarnings(as.numeric(as.character(x)))
  if (any(!is.finite(value))) stop(name, " must contain finite numeric values")
  value
}

# Use the circle's upper-right corner rather than a separate legend column.
.circos_legend_layout <- function(labels = c(0, 0.25, 0.5, 0.75, 1), fontsize = 12, family = "serif") {
  x_scale <- diff(graphics::grconvertX(c(0, 1), "user", "ndc"))
  y_scale <- diff(graphics::grconvertY(c(0, 1), "user", "ndc"))
  text <- textGrob(labels, gp = gpar(fontsize = fontsize, fontfamily = family))
  label_width <- convertWidth(grobWidth(text) + unit(2, "mm"), "npc", valueOnly = TRUE) / x_scale
  label_height <- convertHeight(grobHeight(text), "npc", valueOnly = TRUE) / y_scale
  # Keep the complete legend outside the chromosome-label envelope.
  left <- 1.12 - 0.06 - label_width
  bottom <- sqrt(max(0, 1.04^2 - max(0, left)^2))
  bar_height <- max(0.05, min(0.45, 0.92 - bottom - label_height))
  list(x = graphics::grconvertX(1.12, "user", "ndc"),
       y = graphics::grconvertY(0.92, "user", "ndc"),
       width = unit(0.06 * x_scale, "npc"),
       height = unit(bar_height * y_scale, "npc"))
}

.circos_color_function <- function(limits, colors) {
  if (length(colors) == 1L) colors <- rep(colors, 2L)
  circlize::colorRamp2(seq(limits[1], limits[2], length.out = length(colors)), colors, space = "RGB")
}

#' circos snp
#'
#' @param data chromosome/position matrix or data.frame, or custom chromosome ranges
#' @details Numeric coordinate strings are converted internally; chromosome names
#' are preserved, including nonnumeric names.
#' @param ylim Range of data on y-axis
#' @param cell.padding only affect current track
#' @param track.height Height of the track. It is the percentage to the radius of the unit circles
#' @param cex.lab a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the defaul
#' @param bg.border Color for the border of the plotting regions
#' @param family the name of a font family for drawing text
#' @param legend.space Space reserved to the right of the circle, in circle-radius units; use zero when no legend is needed.
#' @param ... other graphical parameters
#'
#' @import circlize
#' @importFrom graphics par rect text
#' @rdname snpcircos
#' @examples NULL
#' @export

circos.ideogram <- function(data, ylim = c(0, 1), cell.padding = c(0.02, 0, 0.02, 0),
                            track.height = 0.1, cex.lab = 0.6, bg.border = NA, family = NULL, legend.space = 0, ...) {
  if (!is.numeric(legend.space) || length(legend.space) != 1L ||
      !is.finite(legend.space) || legend.space < 0) stop("legend.space must be a finite nonnegative number")
  ### Use CHROM, POS, and ID metadata; only the chromosome and position columns are required.
  if (inherits(data, "custom")) {
    # data structure: chr, start, end
    df <- as.data.frame(data)
    df[, 2] <- .circos_numeric(df[, 2], "Chromosome starts")
    df[, 3] <- .circos_numeric(df[, 3], "Chromosome ends")
  } else {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (ncol(data) < 2L || !nrow(data)) stop("Supply chromosome and position columns")
    data[, 2] <- .circos_numeric(data[, 2], "SNP positions")
    df <- data.frame()
    for (i in unique(data[, 1])) {
      tmp <- data[data[, 1] == i, ]
      df <- rbind(df, c(i, 0, max(tmp[, 2])))
    }
    dimnames(df) <- list(paste0("chr", df[, 1]), c("chr", "start", "end"))
  }
  df[, 1] <- as.character(df[, 1])
  if (anyNA(df[, 1]) || any(!nzchar(df[, 1]))) stop("Chromosome names must not be missing")
  df[, 2] <- .circos_numeric(df[, 2], "Chromosome starts")
  df[, 3] <- .circos_numeric(df[, 3], "Chromosome ends")
  # Set the font.
  family <- ifelse(is.null(family), "serif", family)
  # Reset the circos plotting state.
  circos.clear()
  # Shift the circle left without reducing its radius.
  circos.par(cell.padding = cell.padding, canvas.xlim = c(-0.92, 1.08 + legend.space),
             canvas.ylim = c(-1, 1))
  ## Use factor levels to preserve the chromosome order.
  chromosome_order <- if (all(grepl("^[0-9]+$", df[, 1])))
    df[order(as.numeric(df[, 1])), 1] else sort(df[, 1])
  circos.initialize(factors = factor(df[, 1], levels = chromosome_order), xlim = df[, 2:3])
  ## Draw the axes and labels.
  circos.track(
    ylim = ylim, track.height = track.height, bg.border = bg.border,
    panel.fun = function(x, y, ...) {
      sector.index <- get.current.sector.index()
      # Draw the axis.
      circos.axis(h = "bottom", labels = F)
      # Draw the chromosome label.
      xlim <- get.cell.meta.data("xlim")
      circos.text(mean(xlim), 0.6,
        labels = sector.index, facing = "downward",
        cex = cex.lab, niceFacing = TRUE, family = family
      )
    }
  )
}

#' link snp
#'
#' @param data  A data frame with 7 columns;Chr_1, Pos_1, ID_1, Chr_2, Pos_2, ID_2 value
#' @details Character matrices returned by circos/calc_circos are accepted directly.
#' Coordinates and LD values are converted internally; empty pair tables are accepted.
#' @param ylim  Range of data on y-axis
#' @param color.theme The color theme in the legend
#' @param track.height Height of the track
#' @param legend.show logical whether to show the legend
#' @param legend.position NULL places the legend at the upper right of the circle; an explicit c(x, y) uses normalized device coordinates.
#' @param legend.fontsize the font size of legend
#' @param legend.limits A numeric vector of length two providing limits of the scale
#' @param legend.breaks A numeric vector of positions
#' @param legend.labels A vector giving labels (must be same length as legend.breaks)
#' @param legend.n.breaks An integer guiding the number of major breaks
#' @param bg.border Color for the border of the plotting regions
#' @param bg.col Background color for the plotting regions
#' @param family The font family of legend
#' @param lwd The line width
#' @param lty The line type
#' @param ... other graphical parameters
#'
#' @import viridisLite
#' @importFrom grDevices colorRampPalette grey.colors rainbow col2rgb rgb
#' @export
#'
#' @examples NULL
# The seven columns are CHROM_1, POS_1, ID_1, CHROM_2, POS_2, ID_2, and Value.
# Color links according to their LD values.
circos.linksnp <- function(data, ylim = c(0, 1), color.theme = "viridis", track.height = 0.05, legend.show = TRUE,
                           legend.position = NULL, legend.fontsize = 12, legend.limits = NULL,
                           legend.breaks = NULL, legend.labels = NULL, legend.n.breaks = NULL, family = NULL,
                           bg.border = "black", bg.col = NULL, lwd = par("lwd"), lty = par("lty"), ...) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (ncol(data) != 7L) stop("Circos pairs must have seven columns")
  if (!nrow(data)) {
    circos.clear()
    return(invisible(NULL))
  }
  data[, 1] <- as.character(data[, 1])
  data[, 4] <- as.character(data[, 4])
  data[, 2] <- .circos_numeric(data[, 2], "First SNP positions")
  data[, 5] <- .circos_numeric(data[, 5], "Second SNP positions")
  data[, 7] <- .circos_numeric(data[, 7], "LD values")
  # Set the number of legend colors.
  n.breaks <- if (is.null(legend.n.breaks)) 100L else .whole_number(legend.n.breaks, "legend.n.breaks")
  # Use a palette from viridisLite.
  if (color.theme %in% c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako")) {
    link.colors <- do.call(color.theme, list(n.breaks, alpha = 1, begin = 0, end = 1, direction = -1))
  }
  if (is.null(legend.limits)) {
    # Find the minimum and maximum LD values.
    limits <- c(floor(min(data[, 7]) * 10) / 10, ceiling(max(data[, 7]) * 10) / 10)
  } else {
    limits <- legend.limits
  }
  if (limits[1] == limits[2]) limits <- limits + c(-0.05, 0.05)
  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits)) || diff(limits) <= 0)
    stop("legend.limits must contain two increasing finite numbers")
  colorkey <- .circos_color_function(limits, link.colors)(data[, 7])
  for (i in seq_len(nrow(data))) {
    circos.link(data[i, 1], data[i, 2], data[i, 4], data[i, 5],
      col = as.character(colorkey[i]),
      lwd = lwd,
      lty = lty
    )
  }
  ## Select colors for the chromosome blocks.
  if (is.null(bg.col)) {
    bg.col <- rainbow(length(get.all.sector.index()))
  } else {
    bg.col <- bg.col
  }
  # Draw the chromosome blocks.
  circos.track(
    ylim = ylim, track.height = track.height,
    bg.col = bg.col, bg.border = bg.border,
    bg.lty = par("lty"), bg.lwd = par("lwd")
  )

  # Set the font.
  family <- ifelse(is.null(family), "serif", family)
  # Optionally draw the legend.
  if (legend.show) {
    if (is.null(legend.breaks)) {
      legend.breaks <- seq(limits[1], limits[2], length.out = 5)
    } else {
      legend.breaks <- legend.breaks
    }
    if (is.null(legend.labels)) {
      labels <- legend.breaks
    } else {
      labels <- legend.labels
    }
    layout <- .circos_legend_layout(labels, legend.fontsize, family)
    if (!is.null(legend.position)) {
      if (!is.numeric(legend.position) || length(legend.position) != 2L ||
          any(!is.finite(legend.position))) stop("legend.position must be NULL or two finite numbers")
      layout$x <- legend.position[1]
      layout$y <- legend.position[2]
    }
    legend.just <- if (is.null(legend.position)) c("right", "top") else c("center", "center")
    legend <- circos.legendlink(
      x = layout$x, y = layout$y, width = layout$width, height = layout$height,
      just = legend.just, legend.limits = limits,
      legend.breaks = legend.breaks, legend.n.breaks = n.breaks, legend.col = link.colors,
      legend.labels = labels, legend.family = family, legend.fontsize = legend.fontsize, ...
    )
    grid.draw(legend)
  }
  circos.clear() ## Reset the circos plotting state.
}

#' circos legend
#'
#' @param x numeric
#' @param y numeric
#' @param legend.limits A numeric vector of length two providing limits of the scale
#' @param legend.n.breaks An integer guiding the number of major breaks.
#' @param legend.breaks A numeric vector of positions
#' @param legend.labels A vector giving labels (must be same length as legend.breaks)
#' @param legend.col A specification for the default plotting color
#' @param legend.fontsize numeric; the font size of legend
#' @param legend.family The font family of legend
#' @param just Legend justification, passed to grid::viewport.
#' @param ... Additional arguments passed to color_legend.
#'
#' @param width Width of the legend color bar, as a grid unit.
#' @param height Height of the legend color bar, as a grid unit.
#' @details Uses the built-in color_legend component for the color bar, ticks and label layout.
#' @return A grid grob containing the legend.
#' @export
#'

circos.legendlink <- function(x = 0.88, y = 0.5, legend.limits = c(0, 1), legend.n.breaks = 100,
                              legend.breaks = pretty(legend.limits, n = 4), legend.labels = legend.breaks,
                              legend.col = viridisLite::viridis(legend.n.breaks, direction = -1),
                              legend.fontsize = 12, legend.family = "serif",
                              width = unit(0.035, "snpc"), height = unit(0.25, "snpc"), just = c("center", "center"), ...) {
  if (!is.numeric(legend.limits) || length(legend.limits) != 2L ||
      any(!is.finite(legend.limits)) || diff(legend.limits) <= 0)
    stop("legend.limits must contain two increasing finite numbers")
  legend.n.breaks <- .whole_number(legend.n.breaks, "legend.n.breaks", 1)
  if (!is.numeric(legend.breaks) || length(legend.labels) != length(legend.breaks))
    stop("legend.breaks must be numeric and match legend.labels")
  if (length(legend.col) != legend.n.breaks) stop("Supply one color per legend interval")
  keep <- is.finite(legend.breaks) & legend.breaks >= legend.limits[1] & legend.breaks <= legend.limits[2]
  grob <- color_legend(
    col_fun = .circos_color_function(legend.limits, legend.col), limits = legend.limits,
    at = legend.breaks[keep], labels = legend.labels[keep],
    grid_width = width, legend_height = height, n = max(2L, legend.n.breaks),
    labels_gp = gpar(fontsize = legend.fontsize, fontfamily = legend.family), ...
  )
  grob$vp <- viewport(x = x, y = y, width = grob$vp$width, height = grob$vp$height,
                       just = just, name = "legend_vp")
  grob
}
