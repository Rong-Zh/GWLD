#' Calculate linkage disequilibrium decay
#'
#' @param geno012 a genetype matrix coded by 0,1,2, NA(NA for missing values)
#' @param snpinfo a matrix or data.frame consists of three columns(chr,pos,id)
#' @param method linkage disequilibrium measures
#' @param cores The number of threads
#' @param maxdist The maximum distance calculated by pairwise SNP
#' @param file optional output TSV path; stream in bounded batches instead of returning a matrix
#' @details Values retain double precision. Streaming files include a header.
#' Calculations check for interrupts between batches; interrupted files retain
#' completed rows. cores and maxdist must be whole numbers.
#'
#' @rdname calc_decay
#' @export
#'
calc_decay <- function(geno012, snpinfo, method = "RMI", cores = 1, maxdist = 300, file = NULL) {
  result <- .Call("_GWLD_calc_decay", geno012, as.matrix(snpinfo), .method_name(method),
    .whole_number(cores, "cores"), .whole_number(maxdist, "maxdist", 0L), .pair_output_path(file), PACKAGE = "GWLD")
  if (is.null(file)) result else invisible(result)
}

#' @param m  a upper.tri matrix calculated by linkage disequilibrium measures
#' @param snpinfo a matrix or data.frame consists of three columns(chr,pos,id)
#' @param maxdist The maximum distance calculated by pairwise SNP
#'
#' @rdname calc_decay
#'
#' @export
#'
decay <- function(m, snpinfo, maxdist = 300, file = NULL) {
  result <- .Call("_GWLD_decay", m, as.matrix(snpinfo), .whole_number(maxdist, "maxdist", 0L), .pair_output_path(file), PACKAGE = "GWLD")
  if (is.null(file)) result else invisible(result)
}
