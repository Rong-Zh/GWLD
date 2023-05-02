#' save the plots
#'
#' @param filename output filename of plot.
#' @param plots Plot to save.
#' @param width plot size.
#' @param height plot size.
#' @param dpi Plot resolution.
#' @param units "in", "cm", "mm", or "px".
#' @param device Device to use,Can either be a device function. one of "pdf",
#' "png", "jpg", "jepg", "tiff".
#' @param scale Multiplicative scaling factor.
#' @param ... Other arguments passed on to the graphics device function,
#' as specified by \code{device}.
#'
#' @name save
#' @export
#'

savefig<- function(filename, plots, width=NULL, height=NULL, dpi=300,
                    units=NULL, device=NULL, scale=1, ...) {

  ##设置单位("in", "cm", "mm", "px"),英寸,厘米,毫米,像素
  if(is.null(units)) {
    units <- "in"
  } else {
    units <- units
  }
  ##旧的绘图设备
  old.dev <- grDevices::dev.cur() 
  dev <- plotdev(device, filename, dpi = dpi, units)
  ##设置绘图大小
  if (is.null(width) | is.null(height)) {
    dim <- c(6, 5)
  } else {
    dim <- c(width, height)
  }
  dev(filename = filename, width = dim[1], height = dim[2], ...)
  on.exit(utils::capture.output({
    grDevices::dev.off()
    if (old.dev > 1) grDevices::dev.set(old.dev) # 保留之前的绘图窗口
  }))
  grid.draw(plots)
}

##绘图设备
plotdev <- function(device, filename = NULL, dpi = 300, units) {
  force(filename)
  force(dpi)
  pngdev <- grDevices::png
  jpegdev <- grDevices::jpeg
  tiffdev <- grDevices::tiff
  ##绘图设备列表
  devices <- list(
    pdf =  function(filename, ..., version = "1.4") grDevices::pdf(file = filename, ..., version = version),
    png =  function(...) pngdev(..., res = dpi, units = units),
    jpg =  function(...) jpegdev(..., res = dpi, units = units),
    jpeg = function(...) jpegdev(..., res = dpi, units = units),
    tiff = function(...) tiffdev(..., res = dpi, units = units)
  )
  if (is.null(device)) {
    ##使用的绘图设备
    device <- tools::file_ext(filename)
  }
  dev <- devices[[device]]
  if (is.null(dev)) {
    stop("Unknown graphics device")
  }
  dev
}
