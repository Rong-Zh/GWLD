#' Convert character genotype matrix to integer genotype matrix
#'
#' convert genotypes matrix separated by ("", "/", "|", ...) to integer genotypes matrix(0, 1, 2, NA, NA for missing value)
#'
#' @param data character numeric matrix of genotype
#' @param ... optional arguments (ignored)
#' @param sep separator characters for character genotypes (default "/|")
#' @details Numeric matrices already coded 0/1/2/NA retain their allele direction.
#'
#' @export
#' @name codegeno
#' @rdname codegeno
#' @examples
#' g0 <- c(
#'   "T/A", NA, "T|T", NA, "T|A", NA, "T|T", "T/A", "T/T", "T/T",
#'   "T/A", "A|A", "T/T", "T|A", "T/A", "T|T", NA, "T|A", "T/A", NA
#' )
#'
#' g1 <- c(
#'   "0/0", "1/1", "0/1", "0|0", "1|1", "1/1", "1/1", "0/1", "0/0", "0/0",
#'   "1|1", "0|1", "0|0", "0|0", "1|1", "0/0", "0|0", "0|0", "1|1", "0/0"
#' )
#'
#' gt_mat <- as.matrix(data.frame(g0 = g0, g1 = g1))
#'
#' geno012 <- codegeno(gt_mat, "/|")
setGeneric("codegeno", function(data, sep = "/|", ...) {
  standardGeneric("codegeno")
})

#' @rdname codegeno
setMethod("codegeno",
  signature = (data <- "data.frame"),
  function(data, sep = "/|", ...) {
    codegeno(as.matrix(data), sep = sep)
  }
)

#' @rdname codegeno
#'
setMethod("codegeno",
  signature = (data <- "matrix"),
  function(data, sep = "/|", ...) {
    if (is.numeric(data)) {
      if (any(!is.na(data) & !(data %in% c(0, 1, 2))))
        stop("Numeric genotypes must be 0, 1, 2 or NA")
      if (is.integer(data)) return(data)
      storage.mode(data) <- "integer"
      return(data)
    }
    res <- .Call("_GWLD_code_Mat", data, sep, PACKAGE = "GWLD")
    row.names(res) <- row.names(data)
    colnames(res) <- colnames(data)
    return(res)
  }
)
