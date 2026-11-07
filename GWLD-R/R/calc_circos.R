#' Calculate linkage disequilibrium circos
#'
#' @param geno012 a genetype matrix coded by 0,1,2, NA(NA for missing values)
#' @param snpinfo a matrix or data.frame consists of three columns(chr,pos,id)
#' @param method linkage disequilibrium measures
#' @param cores The number of threads
#' @param threshold the minimum threshold of preserve pairwise calculations
#' @param file optional output TSV path; stream in bounded batches instead of returning a matrix
#' @details Values retain double precision. Streaming files include a header.
#' Calculations check for interrupts between batches; interrupted files retain
#' completed rows. cores must be a whole number.
#'
#' @rdname calc_circos
#' @export
#'
calc_circos <- function(geno012, snpinfo, method = "RMI", cores = 1, threshold = 0.2, file = NULL) {
  result <- .Call("_GWLD_calc_circos", geno012, as.matrix(snpinfo), .method_name(method),
    .whole_number(cores, "cores"), threshold, .pair_output_path(file), PACKAGE = "GWLD")
  if (is.null(file)) result else invisible(result)
}

#' @param m  a upper.tri matrix calculated by linkage disequilibrium measures
#' @param snpinfo a matrix or data.frame consists of three columns(chr,pos,id)
#' @param threshold the minimum threshold of preserve pairwise calculations
#'
#' @rdname calc_circos
#' @export
#'
circos <- function(m, snpinfo, threshold = 0.2, file = NULL) {
  result <- .Call("_GWLD_circos", m, as.matrix(snpinfo), threshold, .pair_output_path(file), PACKAGE = "GWLD")
  if (is.null(file)) result else invisible(result)
}
