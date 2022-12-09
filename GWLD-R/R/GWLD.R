#' GWLD: A package for Genome wide linkage disequilibrium calculate 
#'
#' @useDynLib GWLD, .registration=TRUE
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom Rcpp evalCpp
#' 
#' @aliases data.frame,ANY
#' @aliases matrix,ANY
#' @aliases numeric,numeric
#' @param g1 vector,dataframe or matrix containing genotype coded with interger(0, 1, 2, NA) 
#' @param g2 vector(ignored if g1 is a dataframe or matrix)
#' @param method methods
#' @param cores thread number
#' @param ... optional arguments (ignored)
#' @export
#' @rdname GWLD

setGeneric("GWLD", function(g1, g2, method="RMI", cores=1, ...){standardGeneric("GWLD")})

#' @rdname GWLD
setMethod("GWLD", signature = (g1="data.frame"),
          function(g1, method="r^2", cores=1, ...) {
            if(method %in% c("D", "D'", "r^2", "MI", "RMI")) {
              .Call("_GWLD_gwld", geno012=as.matrix(g1), method, thread=cores, PACKAGE = "GWLD")
            } else {
              stop("Unknown method! Input one of methods(D, D', r^2, MI, RMI) ")
            }
          })

#' @rdname GWLD
setMethod("GWLD", signature = (g1="matrix"),
          function(g1, method="RMI", cores=1, ...) {
            if(method %in% c("D", "D'", "r^2", "MI", "RMI")) {
              .Call("_GWLD_gwld", geno012=as.matrix(g1), method, thread=cores, PACKAGE = "GWLD")
            } else {
              stop("Unknown method! Input one of methods(D, D', r^2, MI, RMI) ")
            }
          })

#' @rdname GWLD
setMethod("GWLD", signature = c(g1="numeric", g2="numeric"),
          function(g1, g2, method="RMI", ...) {
            if(method %in% c("D", "D'", "r^2", "MI", "RMI")) {
              LD(g1, g2, method)
            } else if(method=="MI") {
              MI(g1, g2, cores)
            } else if(method=="RMI") {
              RMI(g1, g2, cores)
            } else {
              stop("Unknown method! Input one of methods(D, D', r^2, MI, RMI) ")
            }
          })

