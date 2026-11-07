# Run from the repository root with the locally installed package.
library(GWLD)
info <- data.frame(CHROM = c("1", "2"), POS = c(100, 100))
pairs <- data.frame("1", 20, "a", "2", 70, "b", 0.2)
output <- tempfile("gwld-legend-")
dir.create(output)
for (size in list(c(1600, 1600), c(2400, 1600), c(1600, 2400), c(1000, 1000))) {
  png(file.path(output, paste0(paste(size, collapse = "x"), ".png")),
      width = size[1], height = size[2], res = 160)
  circos.ideogram(info)
  expected <- GWLD:::.circos_legend_layout()
  # The circle fills at least 80% of the shorter device dimension.
  diameter <- diff(graphics::grconvertX(c(-0.95, 0.95), "user", "ndc")) * size[1]
  stopifnot(diameter / min(size) > 0.8)
  circos.linksnp(pairs)
  legend <- grid::grid.get("legend")
  w <- grid::convertWidth(legend$vp$width, "npc", valueOnly = TRUE)
  h <- grid::convertHeight(legend$vp$height, "npc", valueOnly = TRUE)
  stopifnot(isTRUE(all.equal(as.numeric(legend$vp$x), expected$x)),
            isTRUE(all.equal(as.numeric(legend$vp$y), expected$y)),
            identical(legend$vp$justification, c(1, 1)),
            expected$x - w > 0, expected$x < 1, expected$y - h > 0.5, expected$y < 1)
  dev.off()
  png(file.path(output, paste0("manual-", paste(size, collapse = "x"), ".png")),
      width = size[1], height = size[2], res = 160)
  circos.ideogram(info)
  circos.linksnp(pairs, legend.position = c(0.8, 0.7))
  stopifnot(as.numeric(grid::grid.get("legend")$vp$x) == 0.8,
            as.numeric(grid::grid.get("legend")$vp$y) == 0.7)
  dev.off()
}
png(file.path(output, "component.png"), width = 1200, height = 1200, res = 160)
grob <- circos.legendlink(legend.breaks = c(0, 0.13, 1), legend.labels = c("0", ".13", "1"))
y <- as.numeric(grob$children$legend_ticks$y0)[1:3]
stopifnot(abs((y[2] - y[1]) / (y[3] - y[1]) - 0.13) < 1e-12)
# Clustered labels must fit without overlapping, even if the requested bar is too short.
for (direction in c("vertical", "horizontal")) {
  labels <- c("zero", "very close", "also close", "maximum")
  grob <- color_legend(circlize::colorRamp2(c(0, 1), c("yellow", "purple")),
                       at = c(0, 0.001, 0.002, 1), labels = labels, title = "RMI",
                       direction = direction, legend_height = grid::unit(1, "mm"),
                       legend_width = grid::unit(1, "mm"))
  text <- grob$children$legend_labels
  if (direction == "vertical") {
    centers <- as.numeric(text$y)
    size <- vapply(labels, function(s) grid::convertHeight(grid::grobHeight(grid::textGrob(s, gp = text$gp)), "mm", TRUE), numeric(1))
  } else {
    centers <- as.numeric(text$x)
    size <- vapply(labels, function(s) grid::convertWidth(grid::grobWidth(grid::textGrob(s, gp = text$gp)), "mm", TRUE), numeric(1))
  }
  stopifnot(all(diff(centers) >= (head(size, -1) + tail(size, -1)) / 2 + 1 - 1e-9))
  grid::grid.draw(grob)
}
grid::grid.draw(color_legend(circlize::colorRamp2(c(0, 1), c("yellow", "purple")),
                             at = numeric(), labels = character()))
stopifnot(inherits(try(circos.legendlink(legend.limits = c(1, 0)), silent = TRUE), "try-error"))
stopifnot(inherits(try(color_legend(identity, at = c(0, 0)), silent = TRUE), "try-error"))
dev.off()
stopifnot(!"ComplexHeatmap" %in% loadedNamespaces())
unlink(output, recursive = TRUE)
cat("Compact circle layout, corner/manual positioning, exact ticks and nonoverlapping vertical/horizontal labels passed without ComplexHeatmap\n")
