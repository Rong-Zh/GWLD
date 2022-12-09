#' GWLD: a package for genome-wide linkage disequilibrium analysis
#' 
#' 
#' @param g1 vector,dataframe or matrix containing genotype coded with interger(0, 1, 2, NA) 
#' @param g2 vector(ignored if g1 is a dataframe or matrix)
#' @param ... optional arguments (ignored)
#' @param method methods
#' @param cores thread number
#' 
#' @return a matrix, row names are the same as column names.
#' 
#' @import methods
#' @export
#' @examples
#' g0 <- c('T/A', NA, 'T|T', NA, 'T|A', NA, 'T|T', 'T/A','T/T', 'T/T',
#'         'T/A', 'A|A', 'T/T', 'T|A', 'T/A', 'T|T', NA, 'T|A', 'T/A', NA)
#' 
#' g1 <- c('C/A', 'C/A', 'C/C', 'C/A', 'C/C', 'C/A', 'C/A', 'C/A',
#'        'C/A', 'C/C', 'C/A', 'A/A', 'C/A', 'A/A', 'C/A', 'C/C',
#'        'C/A', 'C/A', 'C/A', 'A/A')
#' 
#' gt <- data.frame(g0, g1)
#' geno012 <- codegeno(gt, "/|")
#' LD(geno012, method="r^2", cores=1)
#' 
#' g2 <- recode(g0,"/|")
#' g3 <- recode(g1,"/")
#' 
#' LD(g2,g3, method="r^2")
#' 
#' g4 <- data.frame(g2, g3)
#' 
#' LD(g4, method="r^2", cores=2)
#' 
#' @rdname LD

setGeneric("LD", function(g1, g2, method="r^2", cores=1, ...){standardGeneric("LD")})

#' @rdname LD
setMethod("LD", signature = (g1="data.frame"),
          function(g1, method="r^2", cores=1) {
            res <- .Call('_GWLD_LDC_Mat', as.matrix(g1), method, cores, PACKAGE = "GWLD")
            row.names(res) <- colnames(res) <- colnames(g1)
            return(res)
          })

#' @rdname LD
setMethod("LD", signature = (g1="matrix"),
          function(g1, method="r^2", cores=1) {
            res <- .Call('_GWLD_LDC_Mat', g1, method, cores, PACKAGE = "GWLD")
            row.names(res) <- colnames(res) <- colnames(g1)
            return(res)
          })

#' @rdname LD
setMethod("LD", signature = c(g1="numeric", g2="numeric"),
          function(g1, g2, method="r^2", ...) {
            .Call('_GWLD_LDC', g1, g2, method, PACKAGE = "GWLD")
          })
