# Spread nearby labels while keeping their order and their tick positions intact.
.legend_label_positions <- function(target, size, axis_length, gap) {
  if (length(target) < 2L) return(target)
  separation <- (utils::head(size, -1L) + utils::tail(size, -1L)) / 2 + gap
  placed <- target
  for (i in 2:length(placed)) placed[i] <- max(placed[i], placed[i - 1L] + separation[i - 1L])
  placed[length(placed)] <- min(placed[length(placed)], axis_length)
  for (i in rev(seq_len(length(placed) - 1L))) placed[i] <- min(placed[i], placed[i + 1L] - separation[i])
  placed
}

#' Continuous color legend
#'
#' A standalone grid legend with physical sizing, accurate ticks and automatic
#' spacing of nearby labels. No ComplexHeatmap dependency is required.
#' @param col_fun Function mapping numeric values to R colors.
#' @param limits Two increasing finite scale limits.
#' @param at Numeric tick positions. Out-of-range ticks are omitted.
#' @param labels Labels corresponding to at; expressions are supported.
#' @param title Optional legend title.
#' @param direction Vertical or horizontal color bar.
#' @param legend_height Requested vertical bar height, as a grid unit.
#' @param legend_width Requested horizontal bar width, as a grid unit.
#' @param grid_width Vertical bar thickness, as a grid unit.
#' @param grid_height Horizontal bar thickness, as a grid unit.
#' @param labels_gp Graphical parameters for labels.
#' @param title_gp Graphical parameters for the title.
#' @param border Border and tick color; NA uses white ticks without a border.
#' @param n Number of color intervals.
#' @param gap Space between the bar and labels, as a grid unit.
#' @param label_gap Minimum spacing between labels, as a grid unit.
#' @return A grid grob. Draw it with grid::grid.draw(); its viewport can be
#' positioned with grid::editGrob(). The bar grows when needed to fit labels.
#' @export
color_legend <- function(col_fun, limits = c(0, 1), at = pretty(limits, n = 4), labels = at,
                         title = NULL, direction = c("vertical", "horizontal"),
                         legend_height = unit(40, "mm"), legend_width = unit(40, "mm"),
                         grid_width = unit(4, "mm"), grid_height = unit(4, "mm"),
                         labels_gp = gpar(fontsize = 12, fontfamily = "serif"),
                         title_gp = gpar(fontsize = 12, fontface = "bold"),
                         border = NA, n = 256L, gap = unit(2, "mm"), label_gap = unit(1, "mm")) {
  direction <- match.arg(direction)
  if (!is.function(col_fun)) stop("col_fun must map numeric values to colors")
  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits)) || diff(limits) <= 0)
    stop("limits must contain two increasing finite numbers")
  if (!is.numeric(at) || length(at) != length(labels)) stop("at must be numeric and match labels")
  n <- .whole_number(n, "n", 2L)
  keep <- is.finite(at) & at >= limits[1] & at <= limits[2]
  at <- at[keep]
  labels <- labels[keep]
  order <- order(at)
  at <- at[order]
  labels <- labels[order]
  if (anyDuplicated(at)) stop("Legend tick positions must be unique")
  mm <- function(x) unit(x, "mm")
  width_mm <- function(x) convertWidth(x, "mm", valueOnly = TRUE)
  height_mm <- function(x) convertHeight(x, "mm", valueOnly = TRUE)
  text_width <- vapply(seq_along(labels), function(i) width_mm(grobWidth(textGrob(labels[i], gp = labels_gp))), numeric(1))
  text_height <- vapply(seq_along(labels), function(i) height_mm(grobHeight(textGrob(labels[i], gp = labels_gp))), numeric(1))
  vertical <- direction == "vertical"
  size <- if (vertical) text_height else text_width
  spacing <- if (vertical) height_mm(label_gap) else width_mm(label_gap)
  gap_mm <- if (vertical) width_mm(gap) else height_mm(gap)
  requested <- if (vertical) height_mm(legend_height) else width_mm(legend_width)
  thickness <- if (vertical) width_mm(grid_width) else height_mm(grid_height)
  if (any(!is.finite(c(spacing, gap_mm, requested, thickness))) ||
      spacing < 0 || gap_mm < 0 || requested <= 0 || thickness <= 0)
    stop("Legend sizes must be positive and spacing must be nonnegative")
  minimum <- if (length(size) > 1L) sum((utils::head(size, -1L) + utils::tail(size, -1L)) / 2 + spacing) else 0
  bar_length <- max(requested, minimum)
  ticks <- (at - limits[1]) / diff(limits) * bar_length
  positions <- .legend_label_positions(ticks, size, bar_length, spacing)
  tw <- max(c(0, text_width))
  th <- max(c(0, text_height))
  content_width <- if (vertical) thickness + gap_mm + tw else bar_length + tw
  content_height <- if (vertical) bar_length + th else thickness + gap_mm + th
  has_title <- !is.null(title) && length(title) > 0L
  title_grob <- if (has_title) textGrob(title, gp = title_gp, name = "legend_title") else nullGrob()
  title_width <- if (has_title) width_mm(grobWidth(title_grob)) else 0
  title_height <- if (has_title) height_mm(grobHeight(title_grob)) else 0
  full_width <- max(content_width, title_width)
  full_height <- content_height + if (has_title) gap_mm + title_height else 0
  offset <- (full_width - content_width) / 2
  colors <- col_fun(seq(limits[1], limits[2], length.out = n))
  if (length(colors) != n || anyNA(colors)) stop("col_fun must return one nonmissing color per value")
  grDevices::col2rgb(colors)
  centers <- (seq_len(n) - 0.5) / n * bar_length
  tick_color <- if (is.na(border)) "white" else border
  if (vertical) {
    pad <- th / 2
    bar <- rectGrob(x = mm(offset), y = mm(pad + centers), width = mm(thickness),
                    height = mm(bar_length / n), just = c("left", "center"),
                    gp = gpar(fill = colors, col = NA), name = "legend_bar")
    ticks_grob <- if (length(ticks)) segmentsGrob(x0 = mm(offset + rep(c(0, thickness), each = length(ticks))),
      x1 = mm(offset + rep(c(thickness * 0.15, thickness * 0.85), each = length(ticks))),
      y0 = mm(pad + rep(ticks, 2L)), y1 = mm(pad + rep(ticks, 2L)),
      gp = gpar(col = tick_color, lwd = 0.7), name = "legend_ticks") else nullGrob(name = "legend_ticks")
    labels_grob <- if (length(labels)) textGrob(labels, x = mm(offset + thickness + gap_mm), y = mm(pad + positions),
                             just = "left", gp = labels_gp, name = "legend_labels") else nullGrob(name = "legend_labels")
    leaders <- if (length(ticks)) segmentsGrob(x0 = mm(offset + thickness), x1 = mm(offset + thickness + gap_mm * 0.75),
      y0 = mm(pad + ticks), y1 = mm(pad + positions),
      gp = gpar(col = ifelse(abs(positions - ticks) > 1e-7, "grey40", NA)), name = "legend_leaders") else nullGrob(name = "legend_leaders")
    box <- rectGrob(x = mm(offset), y = mm(pad), width = mm(thickness), height = mm(bar_length),
                     just = c("left", "bottom"), gp = gpar(fill = NA, col = border), name = "legend_border")
  } else {
    pad <- tw / 2
    bottom <- th + gap_mm
    bar <- rectGrob(x = mm(offset + pad + centers), y = mm(bottom), width = mm(bar_length / n),
                    height = mm(thickness), just = c("center", "bottom"),
                    gp = gpar(fill = colors, col = NA), name = "legend_bar")
    ticks_grob <- if (length(ticks)) segmentsGrob(x0 = mm(offset + pad + rep(ticks, 2L)), x1 = mm(offset + pad + rep(ticks, 2L)),
      y0 = mm(bottom + rep(c(0, thickness), each = length(ticks))),
      y1 = mm(bottom + rep(c(thickness * 0.15, thickness * 0.85), each = length(ticks))),
      gp = gpar(col = tick_color, lwd = 0.7), name = "legend_ticks") else nullGrob(name = "legend_ticks")
    labels_grob <- if (length(labels)) textGrob(labels, x = mm(offset + pad + positions), y = mm(th / 2),
                             gp = labels_gp, name = "legend_labels") else nullGrob(name = "legend_labels")
    leaders <- if (length(ticks)) segmentsGrob(x0 = mm(offset + pad + ticks), x1 = mm(offset + pad + positions),
      y0 = mm(bottom), y1 = mm(th + gap_mm * 0.25),
      gp = gpar(col = ifelse(abs(positions - ticks) > 1e-7, "grey40", NA)), name = "legend_leaders") else nullGrob(name = "legend_leaders")
    box <- rectGrob(x = mm(offset + pad), y = mm(bottom), width = mm(bar_length), height = mm(thickness),
                     just = c("left", "bottom"), gp = gpar(fill = NA, col = border), name = "legend_border")
  }
  if (has_title) {
    title_grob$x <- mm(full_width / 2)
    title_grob$y <- mm(content_height + gap_mm + title_height / 2)
  }
  gTree(children = gList(bar, ticks_grob, leaders, box, labels_grob, title_grob),
        vp = viewport(width = mm(full_width), height = mm(full_height), name = "legend_vp"), name = "legend")
}
