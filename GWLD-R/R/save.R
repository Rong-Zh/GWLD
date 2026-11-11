#' save the plots
#'
#' @param filename output filename of plot.
#' @param plots Plot to save.
#' @param width plot size.
#' @param height plot size.
#' @param dpi Plot resolution.
#' @param units "in", "cm", "mm", or "px".
#' @param device Device to use,Can either be a device function, or one of "pdf",
#' "png", "jpg", "jpeg", "tiff".
#' @param scale Multiplicative scaling factor.
#' @param ... Other arguments passed on to the graphics device function,
#' as specified by \code{device}.
#'
#' @details Dimensions are multiplied by scale and converted to inches.
#' Raster resolution is controlled by dpi; pixel dimensions are converted using dpi.
#' PDF and custom device functions receive width and height in inches.
#' Custom functions must open a graphics device and accept filename, width and height.
#' Missing width and height default independently to 6 and 5 in the requested units.
#' @name save
#' @export
#'

savefig <- function(filename, plots, width = NULL, height = NULL, dpi = 300,
                    units = NULL, device = NULL, scale = 1, ...) {
  positive <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
      stop(name, " must be one finite positive number")
    value
  }
  if (!is.character(filename) || length(filename) != 1L || is.na(filename) || !nzchar(filename))
    stop("filename must be one nonempty path")
  units <- if (is.null(units)) "in" else match.arg(units, c("in", "cm", "mm", "px"))
  dpi <- positive(dpi, "dpi")
  scale <- positive(scale, "scale")
  width <- positive(if (is.null(width)) 6 else width, "width")
  height <- positive(if (is.null(height)) 5 else height, "height")
  per.inch <- switch(units, "in" = 1, cm = 2.54, mm = 25.4, px = dpi)
  width <- positive(width * scale / per.inch, "scaled width")
  height <- positive(height * scale / per.inch, "scaled height")
  dev <- plotdev(device, filename, dpi = dpi, units = "in")
  old.dev <- grDevices::dev.cur()
  dev(filename = path.expand(filename), width = width, height = height, ...)
  output.dev <- grDevices::dev.cur()
  on.exit({
    if (output.dev %in% grDevices::dev.list()) invisible(grDevices::dev.off(output.dev))
    if (old.dev %in% grDevices::dev.list()) invisible(grDevices::dev.set(old.dev))
  })
  grid::grid.draw(plots)
}

## Graphics devices. All sizes passed here are in inches.
plotdev <- function(device, filename = NULL, dpi = 300, units = "in") {
  if (is.function(device)) {
    # Built-in device functions retain the same units/resolution as their names.
    native <- list(pdf = grDevices::pdf, png = grDevices::png,
                   jpeg = grDevices::jpeg, tiff = grDevices::tiff)
    match <- names(native)[vapply(native, identical, logical(1), device)]
    if (!length(match)) return(device)
    device <- match[1L]
  }
  if (is.null(device)) device <- tools::file_ext(filename)
  if (!is.character(device) || length(device) != 1L || is.na(device))
    stop("device must be a graphics function or one supported device name")
  devices <- list(
    pdf = function(filename, ..., version = "1.4") grDevices::pdf(file = filename, ..., version = version),
    png = function(...) grDevices::png(..., res = dpi, units = units),
    jpg = function(...) grDevices::jpeg(..., res = dpi, units = units),
    jpeg = function(...) grDevices::jpeg(..., res = dpi, units = units),
    tiff = function(...) grDevices::tiff(..., res = dpi, units = units)
  )
  dev <- devices[[tolower(device)]]
  if (is.null(dev)) stop("Unknown graphics device")
  dev
}
