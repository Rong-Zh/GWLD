#' @import grid
MakeImageRect <- function(nrow, ncol, cols, name, byrow = TRUE, rectlwd = NULL, rectcol = NA) {
  xx <- (1:ncol) / ncol
  yy <- (1:nrow) / nrow
  # Create coordinate pairs, repeating columns for row-wise filling or rows for column-wise filling.
  if (byrow) {
    right <- rep(xx, nrow)
    top <- rep(yy, each = ncol)
  } else {
    right <- rep(xx, each = nrow)
    top <- rep(yy, ncol)
  }
  # Use NA for the rectangle border color.
  # rectcol <- NA
  rectGrob(
    x = right, y = top,
    width = 1 / ncol, height = 1 / nrow,
    just = c("right", "top"),
    gp = gpar(lwd = rectlwd, col = rectcol, fill = cols),
    name = name
  )
}

MakeImageText <- function(nrow, ncol, cols, name, flip = FALSE, family = family) {
  cols <- as.character(cols)
  cols[is.na(cols)] <- ""
  xx <- (1:ncol) / ncol
  yy <- (1:nrow) / nrow

  # Fill cells in the order used to construct the image.
  if (flip) {
    right <- rep(xx, each = nrow)
    top <- rep(yy, ncol)
  } else {
    right <- rep(xx, nrow)
    top <- rep(yy, each = ncol)
  }
  textGrob(cols,
    x = right - 0.5 / nrow, y = top - 0.5 / ncol, rot = 45,
    gp = gpar(cex = 0.45, fontfamily = family),
    name = name
  )
}

HeatmapLegend <- function(x, y, legend.colors, legend.title, breaks, n.breaks,
                          labels, legend.family) {
  # Normalize to the range [0, 1].
  labels_pos <- seq(0, 1, length.out = n.breaks)
  legend_breaks <- seq(breaks[1], breaks[2], length.out = n.breaks + 1)
  pos <- labels_pos[cut(labels, legend_breaks, labels = FALSE, right = TRUE, include.lowest = T)]
  right <- rep((1:n.breaks) / n.breaks, each = 1)
  top <- rep(1, n.breaks)
  legend_title <- textGrob(legend.title,
    x = 0.5, y = 1.3, name = "legend_title", just = "bottom",
    gp = gpar(cex = 0.85, fontfamily = legend.family)
  )
  legend_rect <- rectGrob(
    x = right, y = top, width = 1 / n.breaks, height = 1,
    just = c("right", "top"), name = "legend_rect", gp = gpar(fill = legend.colors, col = NA)
  )
  legend_box <- linesGrob(x = c(0, 1, 1, 0, 0), y = c(0, 0, 1, 1, 0), name = "legend_box", gp = gpar(col = "black", lwd = 0.4))
  legend_ticks <- segmentsGrob(
    x0 = pos,
    y0 = rep(0, length(pos)),
    x1 = pos,
    y1 = rep(-0.15, length(pos)),
    gp = gpar(col = "black", lwd = 0.7),
    name = "legend_ticks"
  )
  legend_text <- textGrob(labels,
    x = pos, y = rep(-0.2, length(pos)), just = "top",
    gp = gpar(cex = 0.75, fontfamily = legend.family),
    name = "legend_text"
  )
  legend_vp <- viewport(x = x, y = y, width = 0.4, height = 0.035, name = "legend_vp")
  legend <- gTree(children = gList(legend_title, legend_rect, legend_ticks, legend_box, legend_text), vp = legend_vp, name = "legend")
  legend
}

HeatmapLabels <- function(LDmatrix, pos, name, vp, label.size = label.size, family = family) {
  nsnps <- ncol(LDmatrix)
  step <- 1 / (nsnps - 1)
  snp <- ((1:nsnps - 1) + 0.5) / nsnps
  geneMapLocation <- ifelse(is.null(pos), -0.004, -0.15)
  seq.x <- c(
    0.5 * geneMapLocation + 1 / (nsnps * 2),
    1 + 0.5 * geneMapLocation - 1 / (nsnps * 2)
  )
  seq.y <- c(
    -0.5 * geneMapLocation + 1 / (nsnps * 2),
    1 - 0.5 * geneMapLocation - 1 / (nsnps * 2)
  )
  ## Set the SNP positions and spacing.
  if (is.null(pos)) {
    pos <- 1:nsnps ## Use numeric positions.
    total.dist <- max(pos) - min(pos)
    min.dist <- min(pos)
    fraction <- if (total.dist == 0) rep(0.5, nsnps) else (pos - min.dist) / total.dist
    regionx <- seq.x[1] + fraction * (seq.x[2] - seq.x[1])
    regiony <- seq.y[1] + fraction * (seq.y[2] - seq.y[1])
    diagonal <- NULL
    segments <- NULL
  } else {
    ### Draw the physical-position axis.
    diagonal <- linesGrob(seq.x, seq.y, gp = gpar(lty = 1), name = "diagonal", vp = vp)
    ### Connect matrix positions to physical positions.
    total.dist <- max(pos) - min(pos)
    min.dist <- min(pos)
    fraction <- if (total.dist == 0) rep(0.5, nsnps) else (pos - min.dist) / total.dist
    regionx <- seq.x[1] + fraction * (seq.x[2] - seq.x[1])
    regiony <- seq.y[1] + fraction * (seq.y[2] - seq.y[1])
    segments <- segmentsGrob(snp, snp, regionx, regiony, name = "segments", vp = vp)
  }
  ## Set the SNP labels.
  if (!is.null(segments) & is.null(name)) {
    ### Set marker positions and symbols, adjusting their size for the SNP count.
    symbols_pch <- pointsGrob(snp[1:nsnps], snp[1:nsnps],
      pch = "*",
      gp = gpar(cex = -0.02 * nsnps + 2, bg = "blue", col = "blue"),
      name = "symbols_pch", vp = vp
    )
    SNPnames <- NULL
  } else if (!is.null(name)) {
    ## Convert the label gap to viewport coordinates.
    name_gap <- convertWidth(grobWidth(textGrob("r")), "npc", valueOnly = TRUE) / sqrt(2)
    ind <- seq_len(nsnps)
    # Set marker positions and symbols.
    symbols_pch <- pointsGrob(snp[ind], snp[ind],
      pch = "*",
      gp = gpar(cex = -0.02 * nsnps + 2, bg = "blue", col = "blue"),
      name = "symbols_pch", vp = vp
    )
    # Keep labels far enough apart to avoid overlap and update the visible indices.
    min_dist <- convertHeight(grobHeight(textGrob("rs", gp = gpar(fontsize = label.size))), "npc", valueOnly = TRUE) * 1.1
    if (length(ind) >= 2) {
      start <- ind[1]
      res <- c(start)
      for (i in 2:length(ind)) {
        if (abs(regionx[ind[i]] - regionx[start]) >= min_dist) {
          res <- append(res, ind[i])
          start <- ind[i]
        }
      }
      ind <- res
    }
    ## Labels may extend beyond the plot when the device is too small.
    SNPnames <- textGrob(name[ind],
      just = c("left", "bottom"), rot = 90,
      regionx[ind] - sqrt(2 + 0.5) * name_gap,
      regiony[ind] + sqrt(2 + 0.5) * name_gap,
      gp = gpar(fontsize = label.size, col = "blue", fontfamily = family),
      name = "SNPnames", vp = vp
    )
  } else {
    symbols_pch <- NULL
    SNPnames <- NULL
  }
  label <- gTree(children = gList(symbols_pch, segments, diagonal, SNPnames), name = "label")
  label
}
