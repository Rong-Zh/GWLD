#' plink format(.ped, .map) to  matrix
#'
#' @rdname ped2matrix
#' @param FilePrefix plink format(.ped, .map) file prefix
#' @param ... NULL
#' @importFrom data.table fread
#' @export
ped2matrix <- function(FilePrefix, ...) {
  pedFile <- paste0(FilePrefix, ".ped")
  mapFile <- paste0(FilePrefix, ".map")
  ##ped文件前六列:FamilyID, IndividualID,PaternalID, MaternalID,Sex(1=male; 2=female),Phenotype
  ##从第七列就是碱基基因型数据, 两列为一个基因型
  ped <- data.table::fread(file = pedFile, data.table = F)
  ##map文件有四列,CHROM,SNPID,0,POS
  map <- data.table::fread(file = mapFile, data.table = F)
  SnpInfo <- data.frame("CHROM"=map[,1], "POS"=map[,4], "ID"=map[,2])
  samples <- ped[,1:6]
  ped <- ped[,7:ncol(ped)]
  gt <- mapply(function(i, j) paste(ped[,i], ped[,j], sep = "/"),
               c(seq(1,ncol(ped), 2)), c(seq(2,ncol(ped), 2)),
               SIMPLIFY = T, USE.NAMES = F)
  #添加行名
  rownames(gt) <- paste(samples[,1], samples[, 2], sep=":")
  colnames(gt) <- map[,2]
  ped_map <- list(Info=as.matrix(SnpInfo), Genotype = gt)
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
read.plink <- function(FilePrefix,...) {
  prefix = path.expand(FilePrefix)
  bedfile <- paste(prefix, ".bed", sep="")
  famfile <- paste(prefix, ".fam", sep="")
  bimfile <- paste(prefix, ".bim", sep="")
  bim <- read.table(bimfile, header=FALSE, sep="", stringsAsFactors=FALSE)
  fam <- read.table(famfile, header=FALSE, sep="", stringsAsFactors=FALSE)
  nSample <- nrow(fam)
  bed <- .Call("_GWLD_read_bed", bedfile, nSample, PACKAGE="GWLD")
  colnames(bed) <- bim[,2]
  rownames(bed) <- paste(fam[,1], fam[, 2], sep=":")
  plink <- list(bed=bed, bim=bim, fam=fam)
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
  SnpInfo <- data.frame("CHROM"=bim[,1], "POS"=bim[,4], "ID"=bim[,2])
  res <- list(Info=as.matrix(SnpInfo), Genotype=bed)
  return(res)
}


#' vcf file to matrix
#'
#' @rdname vcf2matrix
#' @return list
#' @param file vcf filename
#' @param ... NULL
#' @importFrom data.table fread
#' @export
vcf2matrix <- function(file, ...) {
  
  vcf <- data.table::fread(file = file, data.table = F)
  SnpInfo <- vcf[,1:3]
  colnames(SnpInfo) <- c("CHROM", "POS", "ID")
  gt <- vcf[,10:ncol(vcf)]
  ##转换为行为样品，列为位点
  gt <- t(gt)
  colnames(gt) <- SnpInfo$ID
  rownames(gt) <- gsub("_", ":", row.names(gt))
  return(list(Info=as.matrix(SnpInfo), Genotype=gt))
}

#'  read vcf format file
#' @rdname read.vcf
#' @param filename vcf format file's name
#' @param genotype genotype display type, one of "None, "char", "int" and "allele" 
#' @param ... NULL
#'
#' @return matrix
#' @export
read.vcf <- function(filename, genotype="None", ...) {
    .Call("_GWLD_read_vcf", filename, genotype, PACKAGE="GWLD")
}
