#' plink format(.ped, .map) to  matrix
#'
#' @rdname ped2matrix
#' @param FilePrefix plink format(.ped, .map) file prefix
#' @param genotype "int" (default) for 0/1/2/NA or "char" for allele strings
#' @details Integer genotypes count non-major alleles. Missing PED alleles (0)
#' become NA; frequency ties use lexical allele order. No intermediate genotype
#' strings are constructed for integer output.
#' Positions normally use numeric storage; if any exceeds 2^53, the POS column
#' uses exact decimal text. Fractional or out-of-range positions are rejected.
#' Sample names use FID:IID; percent signs and colons within each component
#' are escaped, and duplicate family/individual pairs are rejected.
#' @param ... NULL
#' @importFrom data.table fread
#' @export
ped2matrix <- function(FilePrefix, genotype = c("int", "char"), ...) {
  genotype <- match.arg(genotype)
  FilePrefix <- path.expand(FilePrefix)
  pedFile <- paste0(FilePrefix, ".ped")
  mapFile <- paste0(FilePrefix, ".map")
  ## The first six PED columns are family ID, individual ID, paternal ID, maternal ID, sex (1=male; 2=female), and phenotype.
  ## Genotypes start at column seven, with two allele columns per variant.
  ped <- data.table::fread(file = pedFile, data.table = FALSE, colClasses = "character", na.strings = NULL, quote = "")
  ## MAP files contain four columns: CHROM, SNPID, genetic distance, and POS.
  map <- data.table::fread(file = mapFile, data.table = FALSE, colClasses = "character", na.strings = NULL, quote = "")
  if (ncol(map) != 4L) stop("MAP must contain four columns")
  SnpInfo <- data.frame(CHROM = map[, 1], POS = .Call("_GWLD_parse_positions", map[[4L]], PACKAGE = "GWLD"), ID = map[, 2])
  if (!nrow(ped) || !nrow(map) || ncol(ped) != 6L + 2L * nrow(map))
    stop("PED must contain six sample fields and two alleles per MAP variant")
  samples <- ped[, 1:6, drop = FALSE]
  gt <- matrix(if (genotype == "int") NA_integer_ else NA_character_, nrow(ped), nrow(map))
  for (j in seq_len(nrow(map))) {
    first <- as.character(ped[[2L * j + 5L]])
    second <- as.character(ped[[2L * j + 6L]])
    first[is.na(first) | first == "0"] <- NA_character_
    second[is.na(second) | second == "0"] <- NA_character_
    if (genotype == "char") {
      first[is.na(first)] <- "."
      second[is.na(second)] <- "."
      gt[, j] <- paste(first, second, sep = "/")
    } else {
      alleles <- c(first, second)
      levels <- sort(unique(alleles[!is.na(alleles)]), method = "radix")
      if (length(levels) > 2L) stop("PED contains a multiallelic variant: ", map[j, 2L])
      if (length(levels)) {
        counts <- tabulate(match(alleles, levels), nbins = length(levels))
        major <- levels[which.max(counts)]
        gt[, j] <- as.integer(first != major) + as.integer(second != major)
      }
    }
  }
  # Add sample names.
  rownames(gt) <- .plink_sample_names(samples[, 1], samples[, 2])
  colnames(gt) <- map[, 2]
  ped_map <- list(Info = SnpInfo, Genotype = gt)
  return(ped_map)
}

#' read plink format(.bed, .bim, .fam)
#'
#' @rdname read.plink
#' @details Positions are validated exactly. BIM position columns normally use
#' numeric storage, switching to decimal text if any position exceeds 2^53.
#' @param FilePrefix plink format(.bed, .bim, .fam) file's prefix
#' @param ... NULL
#' @importFrom utils read.table
#' @return list
#' @export
read.plink <- function(FilePrefix, ...) {
  prefix <- path.expand(FilePrefix)
  bedfile <- paste(prefix, ".bed", sep = "")
  famfile <- paste(prefix, ".fam", sep = "")
  bimfile <- paste(prefix, ".bim", sep = "")
  bim <- read.table(bimfile, header = FALSE, sep = "", stringsAsFactors = FALSE, na.strings = NULL, quote = "", comment.char = "",
                    colClasses = c("character", "character", NA, "character", "character", "character"))
  fam <- read.table(famfile, header = FALSE, sep = "", stringsAsFactors = FALSE, na.strings = NULL, quote = "", comment.char = "",
                    colClasses = c("character", "character", NA, NA, NA, NA))
  if (ncol(fam) != 6L || ncol(bim) != 6L || !nrow(fam) || !nrow(bim))
    stop("Nonempty FAM and BIM files must each have six columns")
  bim[[4L]] <- .Call("_GWLD_parse_positions", bim[[4L]], PACKAGE = "GWLD")
  sample_names <- .plink_sample_names(fam[, 1], fam[, 2])
  nSample <- nrow(fam)
  bed <- .Call("_GWLD_read_bed", bedfile, nSample, PACKAGE = "GWLD")
  if (ncol(bed) != nrow(bim)) stop("BED/BIM variant count mismatch")
  colnames(bed) <- bim[, 2]
  rownames(bed) <- sample_names
  plink <- list(bed = bed, bim = bim, fam = fam)
  return(plink)
}

#' plink format(.bed, .bim, .fam) to matrix
#'
#' @rdname bed2matrix
#' @details Info POS is numeric unless any coordinate exceeds 2^53; then it
#' uses exact decimal text. Fractional or out-of-range positions are rejected.
#' @param FilePrefix plink format(.bed, .bim, .fam) file's prefix
#' @param ... NULL
#' @return list
#' @export
bed2matrix <- function(FilePrefix, ...) {
  plink <- read.plink(FilePrefix)
  bed <- plink$bed
  bim <- plink$bim
  SnpInfo <- data.frame(CHROM = bim[, 1], POS = bim[, 4], ID = bim[, 2])
  res <- list(Info = SnpInfo, Genotype = bed)
  return(res)
}


#' vcf file to matrix
#'
#' @rdname vcf2matrix
#' @return list
#' @param file VCF, compressed VCF (.vcf.gz), or BCF filename
#' @param threads thread budget for BGZF decompression (default 1)
#' @param ... NULL
#' @details Uses the shared HTSlib reader. Genotype is an integer matrix with
#' samples as rows and variants as columns, coded 0/1/2/NA. Sample names are
#' preserved without replacing underscores. Supports biallelic diploid GT.
#' Info POS is numeric unless a coordinate exceeds 2^53, in which case the
#' column uses exact decimal text to preserve all integer coordinates.
#' @export
vcf2matrix <- function(file, threads = 1L, ...) {
  .Call("_GWLD_read_vcf_matrix", path.expand(file), .whole_number(threads, "threads"), PACKAGE = "GWLD")
}

#'  read vcf format file
#' @rdname read.vcf
#' @param filename VCF, compressed VCF (.vcf.gz), or BCF filename
#' @param genotype genotype display type, one of "None, "char", "int" and "allele"
#' @param ... NULL
#'
#' @return matrix
#' @export
read.vcf <- function(filename, genotype = "None", ...) {
  .Call("_GWLD_read_vcf", filename, genotype, PACKAGE = "GWLD")
}
