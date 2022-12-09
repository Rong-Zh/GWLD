#' Calculate linkage disequilibrium circos 
#'
#' @param geno012 a genetype matrix coded by 0,1,2, NA(NA for missing values)
#' @param snpinfo a matrix or data.frame consists of three columns(chr,pos,id)
#' @param method linkage disequilibrium measures
#' @param cores The number of threads
#' @param threshold the minimum threshold of preserve pairwise calculations
#' 
#' @rdname calc_circos
#' @export
#'
calc_circos <- function(geno012, snpinfo, method="RMI", cores=1, threshold=0.2) {
  if(method %in% c("D", "D'", "r^2", "MI", "RMI")) {
    .Call("_GWLD_calc_circos", geno012, as.matrix(snpinfo), method, thread=cores, 
          threshold=threshold, PACKAGE = "GWLD")
  } else {
    stop("Unknown method! Input one of methods(D, D', r^2, MI, RMI) ")
  }
}

#' @param m  a upper.tri matrix calculated by linkage disequilibrium measures
#' @param snpinfo a matrix or data.frame consists of three columns(chr,pos,id)
#' @param threshold the minimum threshold of preserve pairwise calculations
#' 
#' @rdname calc_circos
#' @export
#'
circos <- function(m, snpinfo, threshold=0.2) {
  .Call("_GWLD_circos", m, as.matrix(snpinfo), threshold=threshold, PACKAGE = "GWLD")
}


