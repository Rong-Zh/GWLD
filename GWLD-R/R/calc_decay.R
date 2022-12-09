#' Calculate linkage disequilibrium decay
#'
#' @param geno012 a genetype matrix coded by 0,1,2, NA(NA for missing values)
#' @param snpinfo a matrix or data.frame consists of three columns(chr,pos,id)
#' @param method linkage disequilibrium measures
#' @param cores The number of threads
#' @param maxdist The maximum distance calculated by pairwise SNP
#' 
#' @rdname calc_decay
#' @export
#'
calc_decay <- function(geno012, snpinfo, method="RMI", cores=1, maxdist=300) {
  if(method %in% c("D", "D'", "r^2", "MI", "RMI")) {
   .Call("_GWLD_calc_decay", geno012, as.matrix(snpinfo), method, thread=cores, maxdist=maxdist, PACKAGE = "GWLD")
  } else {
    stop("Unknown method! Input one of methods(D, D', r^2, MI, RMI) ")
  }
}

#' @param m  a upper.tri matrix calculated by linkage disequilibrium measures
#' @param snpinfo a matrix or data.frame consists of three columns(chr,pos,id)
#' @param maxdist The maximum distance calculated by pairwise SNP
#' 
#' @rdname calc_decay
#'
#' @export
#'
decay <- function(m, snpinfo, maxdist=300) {
  .Call("_GWLD_decay", m, as.matrix(snpinfo), maxdist=maxdist, PACKAGE = "GWLD")
}


