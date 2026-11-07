#' linkage disequilibrium heatmap
#'
#' @param data genotype matrix/data.frame, or a reader result with Genotype and Info
#' @param method measure linkage disequilibrium method
#' @param SnpPosition data.frame or matrix
#' @param SnpName vector
#' @param color a preset name or a nonempty vector of valid R colors
#' @param showLDvalues logical
#' @param title heatmap's title
#' @param label.size label's font size of heatmap
#' @param family Font for heatmap
#' @param LDmatrix optional precomputed LD matrix; named variants must match the supplied IDs and order
#' @param cores thread number
#' @param ... NULL
#'
#' @import grid
#' @return HeatMap
#' @export
#' @examples NULL
HeatMap <- function(data = NULL, method = "r^2", SnpPosition = NULL, SnpName = NULL, cores = 1,
                    color = NULL, showLDvalues = FALSE, title = "Pairwise LD", label.size = 10, family = NULL,
                    LDmatrix = NULL, ...) {
  method <- .method_name(method)
  cached <- !is.null(LDmatrix)
  if (is.list(data) && !is.data.frame(data)) {
    if (is.null(data$Genotype) || is.null(data$Info)) stop("Reader result must contain Genotype and Info")
    if (is.null(SnpPosition)) SnpPosition <- data$Info[, 1:2, drop = FALSE]
    if (is.null(SnpName)) SnpName <- data$Info[, 3L]
    data <- data$Genotype
  }
  if (method == "r^2") {
    breaks <- c(0, 1)
    labels <- seq(0, 1, 0.2)
    legend.title <- expression(paste(r^2, " Color Key"))
  } else if (method == "D'") {
    breaks <- c(0, 1)
    labels <- seq(0, 1, 0.2)
    legend.title <- "D' Color Key"
  } else if (method == "D") {
    breaks <- c(-1, 1)
    labels <- seq(-1, 1, 0.5)
    legend.title <- "D Color Key"
  } else if (method == "RMI") {
    breaks <- c(-0.1, 1.6)
    labels <- seq(0, 1.5, 0.3)
    legend.title <- "RMI Color Key"
  } else if (method == "MI") {
    breaks <- c(0, 1.6)
    labels <- seq(0, 1.5, 0.3)
    legend.title <- "MI Color Key"
  } else {
    stop("Invalid LD measurement")
  }
  if (is.null(LDmatrix)) {
    if (is.null(data)) stop("Supply genotypes or a precomputed LDmatrix")
    LDmatrix <- GWLD(codegeno(data), method = method, cores = cores)
  }
  if (!is.matrix(LDmatrix) || !is.numeric(LDmatrix) || !nrow(LDmatrix) || nrow(LDmatrix) != ncol(LDmatrix))
    stop("LDmatrix must be a nonempty square numeric matrix")
  if (!is.null(data) && ncol(data) != nrow(LDmatrix)) stop("LDmatrix and genotype variant counts differ")
  if (is.null(SnpName)) SnpName <- colnames(LDmatrix)
  if (!is.null(SnpName) && length(SnpName) != nrow(LDmatrix)) stop("SnpName length must match LDmatrix")
  if (!is.null(SnpName)) {
    SnpName <- as.character(SnpName)
    if (anyNA(SnpName)) stop("SnpName must not contain missing variant IDs")
    if (cached && anyDuplicated(SnpName)) stop("Cached LDmatrix matching requires unique variant IDs")
    for (ids in dimnames(LDmatrix)) {
      if (!is.null(ids) && !identical(as.character(ids), SnpName))
        stop("LDmatrix variant IDs/order do not match SnpName or reader metadata; align them explicitly")
    }
    if (is.null(rownames(LDmatrix))) rownames(LDmatrix) <- SnpName
    if (is.null(colnames(LDmatrix))) colnames(LDmatrix) <- SnpName
  }
  # The SnpPosition CHORM, POS;
  # Distances should not be calculated between different chromosomes
  if (!is.null(SnpPosition)) {
    if (is.null(dim(SnpPosition))) {
      if (length(SnpPosition) != nrow(LDmatrix)) stop("SnpPosition length must match LDmatrix")
    } else {
      if (nrow(SnpPosition) != nrow(LDmatrix) || ncol(SnpPosition) < 2L)
        stop("SnpPosition must have one chromosome/position row per variant")
      Chr <- SnpPosition[, 1]
      SnpPosition <- if (length(unique(Chr)) > 1) NULL else SnpPosition[, 2]
    }
    if (!is.null(SnpPosition)) {
      SnpPosition <- suppressWarnings(as.numeric(as.character(SnpPosition)))
      if (any(!is.finite(SnpPosition))) stop("SnpPosition must contain finite numeric positions")
    }
  }

  # Create the plot viewports.
  mainvp <- viewport(
    width = unit(0.8, "snpc"), height = unit(0.8, "snpc"),
    name = "mainvp"
  )
  flipVP <- viewport(
    width = unit(0.85, "snpc"), height = unit(0.85, "snpc"),
    y = 0.6, angle = -45, name = "flipVP"
  )
  # Select the color palette.
  if (!is.null(color) && (!is.character(color) || !length(color) || anyNA(color)))
    stop("color must be a preset name or a nonempty color vector")
  if (is.null(color)) {
    color <- grey.colors(100)
  } else if (length(color) == 1L && color == "blueTored") {
    colors <- colorRampPalette(c("red", "purple", "blue"), space = "rgb")
    color <- colors(100)
  } else if (length(color) == 1L && color == "YellowTored") {
    colors <- colorRampPalette(c("yellow", "orange", "red"), space = "rgb")
    color <- colors(100)
  } else {
    color <- color
  }
  if (!is.character(color) || !length(color) || anyNA(color)) stop("color must be a preset name or a nonempty color vector")
  tryCatch(col2rgb(color), error = function(e) stop("Invalid heatmap color: ", conditionMessage(e)))
  # Select the color palette.
  mybreak <- seq(breaks[1], breaks[2], length.out = length(color) + 1)
  colcut <- as.character(cut(LDmatrix, mybreak, labels = as.character(color), include.lowest = T, right = T))
  # Create the heatmap cells.
  rectcol <- ifelse(nrow(LDmatrix) <= 60, "white", NA)
  ImageRect <- MakeImageRect(dim(LDmatrix)[1], dim(LDmatrix)[2], colcut,
    name = "heatmapRect", byrow = FALSE, rectlwd = 0.5,
    rectcol = rectcol
  )
  # Set the plot font.
  family <- ifelse(is.null(family), "serif", family)
  # Optionally display LD values in the cells.
  if (showLDvalues) {
    ImageText <- MakeImageText(dim(LDmatrix)[1], dim(LDmatrix)[2],
      round(LDmatrix, digits = 2),
      name = "heatmaptext",
      flip = TRUE, family = family
    )
  } else {
    ImageText <- NULL
  }
  # Place the cells and LD values in the rotated viewport.
  heatMap <- gTree(children = gList(ImageRect, ImageText), vp = flipVP, name = "heatMap")
  # Set the main plot title.
  ## Place the main title above the heatmap.
  main_title <- textGrob(title, 0.5, 1.05,
    gp = gpar(cex = 1.2, fontfamily = family),
    name = "main_title"
  )
  ## Place the region label above the heatmap.
  if (is.null(SnpPosition)) {
    mapLabel <- NULL
  } else {
    mapLabel <- paste("Physical Length: ", round(((max(SnpPosition) - min(SnpPosition)) / 1000), 1), " Kb", sep = "")
  }
  label_title <- textGrob(mapLabel, 0.5, 0.95,
    gp = gpar(cex = 1, fontfamily = family),
    just = "centre", name = "label_title"
  )
  title <- gTree(children = gList(main_title, label_title), name = "title")
  # Create the SNP labels.
  label <- HeatmapLabels(LDmatrix, SnpPosition, SnpName, vp = flipVP, label.size = label.size, family = family)
  # Create the legend.
  legend <- HeatmapLegend(0.8, 0,
    legend.colors = color, legend.title, breaks = breaks, n.breaks = length(color),
    labels = labels, legend.family = family
  )
  LDheatmap <- gTree(
    children = gList(heatMap, legend, label, title),
    vp = mainvp, name = "LDheatmap"
  )
  grid.newpage()
  grid.draw(LDheatmap)

  invisible(structure(list(
    LDmatrix = LDmatrix, LDheatMap = LDheatmap,
    SnpName = SnpName, SnpPosition = SnpPosition,
    Color = color
  ), class = "HeatMap"))
}

## Drawing methods for HeatMap objects.
#' @importFrom grid grid.draw
#' @export
grid.draw.HeatMap <- function(x, recording = TRUE) {
  grid.newpage()
  grid.draw(x$LDheatMap)
}

#' @export
print.HeatMap <- function(x, ...) {
  grid.draw(x)
}
