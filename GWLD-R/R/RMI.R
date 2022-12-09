#' Compute pairwise linkage disequilibrium with reduced mutual information
#'
#' @param g1 vector,dataframe or matrix containing genotype coded with interger(0, 1, 2, NA)
#' @param g2 vector(ignored if g1 is a dataframe or matrix)
#' @param ... optional arguments (ignored)
#' @param cores thread number
#' 
#' @return a matrix, row names are the same as column names.
#' 
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
#' RMI(geno012, cores=2)
#' 
#' g2 <- recode(g0,"/|")
#' g3 <- recode(g1,"/")
#' 
#' RMI(g2,g3)
#' 
#' g4 <- data.frame(g2, g3)
#' 
#' RMI(g4, cores=2)
#' 
#' @rdname RMI
#' @export
setGeneric("RMI", function(g1, g2, cores=1, ...){standardGeneric("RMI")})

#' @rdname RMI
setMethod("RMI", signature(g1="data.frame"), 
          function(g1, cores=1) {
            res <- .Call('_GWLD_RMIC_Mat', as.matrix(g1), cores, PACKAGE = "GWLD")
            row.names(res) <- colnames(res) <- colnames(g1)
            return(res)
          })

#' @rdname RMI
setMethod("RMI", signature = (g1="matrix"),
          function(g1, cores=1) {
            res <- .Call('_GWLD_RMIC_Mat', g1, cores, PACKAGE = "GWLD")
            row.names(res) <- colnames(res) <- colnames(g1)
            return(res)
          })

#' @rdname RMI
setMethod("RMI", signature(g1="numeric", g2="numeric"),
          function(g1, g2) {
            .Call('_GWLD_RMIC', g1, g2, PACKAGE = "GWLD")
          })
