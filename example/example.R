# install GWLD
# Install Rtools first for Windows. Visit https://cran.r-project.org/bin/.
# install LLVM and libomp first for macOS M chip.
install.packages("devtools")
devtools::install_github("Rong-Zh/GWLD/GWLD-R")

# loading package
library(GWLD)

# loading example data with base-type genotype
data("duck")
dta <- duck$SNP
SNP <- dta$genotype
Info <- dta$info
geno012 <- codegeno(dta$genotype, sep = "/")


## VCF format
# Approach 1
SNP <- read.vcf(filename = "example-3.vcf", genotype = "int")
SNP <- data.frame(SNP)
Info <- SNP[, 1:3]
SNP <- data.frame(t(SNP[, 10:ncol(SNP)]))
geno012 <- apply(SNP, 2, as.numeric)
row.names(geno012) <- row.names(SNP)
colnames(geno012) <- Info$ID

# Approach 2
SNP <- vcf2matrix("example-3.vcf")
Info <- SNP$Info
head(SNP$Genotype)
geno012 <- codegeno(SNP$Genotype, sep = "/")

## Plink format
# read plink format(ped, map)
SNP <- ped2matrix(FilePrefix = "example-2")
Info <- SNP$Info
head(SNP$Genotype)
geno012 <- codegeno(SNP$Genotype, sep = "/")

# read Plink format data(bed, bim, fam)
# Approach 1
SNP <- read.plink("example")
Info <- SNP$bim
Info <- Info[, c(1, 4, 2)]
colnames(Info) <- c("CHROM", "POS", "ID")
geno012 <- SNP$bed

# Approach 2
SNP <- bed2matrix("example")
Info <- SNP$Info
geno012 <- SNP$Genotype


# Using the plink data above as an example
# Calculating values with different methods.
r2 <- LD(geno012, cores = 1)
mi <- MI(geno012, cores = 1)
RMI <- RMI(geno012, cores = 1)
value <- GWLD(geno012, method = "RMI", cores = 1)

# Calculating the value between two loci.
r2 <- LD(geno012[, 1], geno012[, 2], method = "r^2", cores = 1)
value <- GWLD(geno012[, 1], geno012[, 2], method = "RMI", cores = 1)

# Calculate the decay from result
RMI <- RMI(geno012, cores = 1)
rmi_decay <- decay(RMI, Info)

# or calculate circos result from data
rmi_decay <- calc_decay(geno012, Info, method = "RMI")


# Calculate the decay from result
RMI <- RMI(geno012, cores = 1)
rmi_circos <- circos(RMI, Info, threshold = 0.1)

# Calculate the distance between different chromosomes
rmi_circos <- calc_circos(geno012, Info, method = "RMI", threshold = 0.1, cores = 1)


# plot Ciros
circosdata <- duck$Circos
circos.ideogram(circosdata$chr)
circos.linksnp(circosdata$linkdata)

# plot using the example data,
# It is possible to construct all or partial chromosome length data
chr_info <- Info |>
  data.frame() |>
  dplyr::mutate_at(c("CHROM", "POS"), .funs = as.numeric)

rmi_res <- rmi_circos |>
  data.frame() |>
  dplyr::mutate_at(c("CHR_1", "POS_1", "CHR_2", "POS_2", "Value"), .funs = as.numeric)

# show all chromosomes of the duck data
circos.ideogram(circosdata$chr)
circos.linksnp(rmi_res)

# show partial chromosomes of the duck data
circos.ideogram(chr_info)
circos.linksnp(rmi_res)

# construct chromosomes data
data <- data.frame(
  chr = c(1, 2, 3, 4, 5),
  start = c(0, 0, 0, 0, 0),
  end = c(164695780, 124299693, 49759293, 10956466, 60265033)
)

row.names(data) <- paste0("chr", data[,1]) 
circos.ideogram(data)
circos.linksnp(rmi_res)

# For more details, see https://doi.org/10.1093/g3journal/jkad154 Box 1
