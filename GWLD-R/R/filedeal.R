#' plink format(.ped, .map) to  matrix
#'
#' @rdname ped2matrix
#' @param FilePrefix plink format(.ped, .map) file prefix
#' @param genotype "int" (default) for 0/1/2/NA or "char" for allele strings
#' @details Integer genotypes count non-major alleles. Missing PED alleles (0)
#' become NA; frequency ties use lexical allele order. No intermediate genotype
#' strings are constructed for integer output.
#' @param ... NULL
#' @importFrom data.table fread
#' @export
ped2matrix <- function(FilePrefix, genotype = c("int", "char"), ...) {
  genotype <- match.arg(genotype)
  FilePrefix <- path.expand(FilePrefix)
  pedFile <- paste0(FilePrefix, ".ped")
  mapFile <- paste0(FilePrefix, ".map")
  ## ped文件前六列:FamilyID, IndividualID, PaternalID, MaternalID, Sex(1=male; 2=female),Phenotype
  ## 从第七列就是碱基基因型数据, 两列为一个基因型
  ped <- data.table::fread(file = pedFile, data.table = FALSE)
  ## map文件有四列, CHROM, SNPID, 0, POS
  map <- data.table::fread(file = mapFile, data.table = FALSE)
  SnpInfo <- data.frame(CHROM = map[, 1], POS = as.numeric(map[, 4]), ID = map[, 2])
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
  # 添加行名
  rownames(gt) <- paste(samples[, 1], samples[, 2], sep = ":")
  colnames(gt) <- map[, 2]
  ped_map <- list(Info = SnpInfo, Genotype = gt)
  return(ped_map)
}

#' read plink format(.bed, .bim, .fam)
#'
#' @rdname read.plink
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
  bim <- read.table(bimfile, header = FALSE, sep = "", stringsAsFactors = FALSE)
  fam <- read.table(famfile, header = FALSE, sep = "", stringsAsFactors = FALSE)
  nSample <- nrow(fam)
  bed <- .Call("_GWLD_read_bed", bedfile, nSample, PACKAGE = "GWLD")
  colnames(bed) <- bim[, 2]
  rownames(bed) <- paste(fam[, 1], fam[, 2], sep = ":")
  plink <- list(bed = bed, bim = bim, fam = fam)
  return(plink)
}

#' plink format(.bed, .bim, .fam) to matrix
#'
#' @rdname bed2matrix
#' @param FilePrefix plink format(.bed, .bim, .fam) file's prefix
#' @param ... NULL
#' @return list
#' @export
bed2matrix <- function(FilePrefix, ...) {
  plink <- read.plink(FilePrefix)
  bed <- plink$bed
  bim <- plink$bim
  SnpInfo <- data.frame(CHROM = bim[, 1], POS = as.numeric(bim[, 4]), ID = bim[, 2])
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
#' @export
vcf2matrix <- function(file, threads = 1L, ...) {
  .Call("_GWLD_read_vcf_matrix", path.expand(file), as.integer(threads), PACKAGE = "GWLD")
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
