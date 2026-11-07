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
#' @param g1 numeric vector, data.frame or matrix of exact 0/1/2/NA genotype codes (-1 also denotes missing)
#' @param g2 vector(ignored if g1 is a dataframe or matrix)
#' @param method RMI, MI, r2/r^2, D or Dprime/D'; lowercase names are accepted
#' @param cores thread number
#' @param ... optional arguments (ignored)
#' @export
#' @rdname GWLD

setGeneric("GWLD", function(g1, g2, method = "RMI", cores = 1, ...) {
  standardGeneric("GWLD")
})

#' @rdname GWLD
setMethod("GWLD",
  signature = c(g1 = "data.frame"),
  function(g1, method = "RMI", cores = 1, ...) {
    GWLD(as.matrix(g1), method = method, cores = cores)
  }
)

#' @rdname GWLD
setMethod("GWLD",
  signature = c(g1 = "matrix"),
  function(g1, method = "RMI", cores = 1, ...) {
    method <- .method_name(method)
    cores <- .whole_number(cores, "cores")
    res <- .Call("_GWLD_gwld", geno012 = g1, method, thread = cores, PACKAGE = "GWLD")
    dimnames(res) <- list(colnames(g1), colnames(g1))
    res
  }
)

#' @rdname GWLD
setMethod("GWLD",
  signature = c(g1 = "numeric", g2 = "numeric"),
  function(g1, g2, method = "RMI", cores = 1, ...) {
    method <- .method_name(method)
    if (method %in% c("D", "D'", "r^2")) {
      LD(g1, g2, method)
    } else if (method == "MI") {
      MI(g1, g2)
    } else if (method == "RMI") {
      RMI(g1, g2)
    } else {
      stop("Unknown method! Input one of methods(D, D', r^2, MI, RMI) ")
    }
  }
)
