.method_name <- function(method) {
  if (!is.character(method) || length(method) != 1L || is.na(method)) stop("Supply one LD method name")
  methods <- c(rmi = "RMI", mi = "MI", r2 = "r^2", "r^2" = "r^2", d = "D", dprime = "D'", "d'" = "D'")
  selected <- unname(methods[tolower(method)])
  if (is.na(selected)) stop("Unknown method; choose RMI, MI, r2, D or Dprime")
  selected
}
.whole_number <- function(value, name, minimum = 1L) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < minimum || value != floor(value) || value > .Machine$integer.max)
    stop(name, " must be a whole number between ", minimum, " and ", .Machine$integer.max)
  as.integer(value)
}
.pair_output_path <- function(file) {
  if (is.null(file)) return("")
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file))
    stop("file must be NULL or one nonempty output path")
  path.expand(file)
}
.plink_sample_names <- function(family, individual) {
  family <- as.character(family)
  individual <- as.character(individual)
  if (anyNA(family) || anyNA(individual) || any(!nzchar(family)) || any(!nzchar(individual)))
    stop("PLINK family and individual IDs must not be missing")
  escape <- function(x) gsub(":", "%3A", gsub("%", "%25", x, fixed = TRUE), fixed = TRUE)
  names <- paste(escape(family), escape(individual), sep = ":")
  if (anyDuplicated(names)) stop("Duplicate PLINK family/individual IDs")
  names
}
